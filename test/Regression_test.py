#!/usr/bin/env python3
# ATAC-seq Regression Model Evaluation Pipeline
# Evaluates a trained deep learning model on chromatin accessibility regression data

# ===== 1. Environment Configuration and Argument Parsing =====
import argparse
import random
import os
import sys
import numpy as np
import gc
from pyfasta import Fasta
import time
import pandas as pd
import torch
import torch.nn as nn
import math
from torch.utils.data import Dataset, DataLoader
from sklearn.metrics import mean_squared_error, r2_score
from scipy import stats

# Parse command line arguments
parser = argparse.ArgumentParser(description='ATAC-seq Regression Model Evaluation Pipeline')
parser.add_argument('fold_idx', type=int, help='Fold index for cross-validation')
parser.add_argument('input_dir', help='Directory containing input files')
parser.add_argument('output_dir', help='Directory for output files')
parser.add_argument('-g', '--genome_fasta', required=True, 
                    help='Path to genome FASTA file (union.peaks.pad1k.fa)')
parser.add_argument('-l', '--label_file', required=True, 
                    help='Path to label file (union.peaks.labels.regress.txt)')
parser.add_argument('-e', '--expression_file', required=True, 
                    help='Path to expression matrix file (final_exp_matrix.txt)')
parser.add_argument('-c', '--crossval_file', required=True, 
                    help='Path to cross-validation file (Leave_one_out_cross_validation.txt)')
parser.add_argument('-m', '--model_path', required=True, 
                    help='Path to trained model file (regress_model.pth)')
parser.add_argument('-s', '--seq_len', type=int, default=1000, 
                    help='Sequence length (default: 1000)')
parser.add_argument('-n', '--num_gene', type=int, default=2500, 
                    help='Number of genes (default: 2500)')
parser.add_argument('-b', '--batch_size', type=int, default=12000, 
                    help='Batch size for evaluation (default: 12000)')
parser.add_argument('-r', '--ratio', type=float, default=1.0, 
                    help='Sampling ratio of regions (default: 1.0)')
args = parser.parse_args()

# Set global parameters from arguments
SEQ_LEN = args.seq_len
NUM_GENE = args.num_gene
BATCH_SIZE = args.batch_size
RATIO = args.ratio

# Create output directory
os.makedirs(args.output_dir, exist_ok=True)

# Set CUDA configuration
os.environ["CUDA_LAUNCH_BLOCKING"] = "1"

# ===== 2. Function Definitions =====
def seq_to_mat(seq):
    """Convert DNA sequence to one-hot encoded matrix"""
    d = {'a': 0, 'A': 0, 'g': 1, 'G': 1, 'c': 2, 'C': 2, 't': 3, 'T': 3, 'N': 4, 'n': 4}
    mat = np.zeros((5, SEQ_LEN))
    for i in range(len(seq)):
        mat[d[seq[i]], i] = 1
    mat = mat[:4, :]
    return mat

class ChromatinAccessibilityRegressor(nn.Module):
    """Deep learning model for predicting chromatin accessibility levels"""
    def __init__(self, num_gene, seq_len):
        super(ChromatinAccessibilityRegressor, self).__init__()
        # Convolutional layers for sequence processing
        self.seq_conv1 = nn.Conv2d(1, 160, kernel_size=(4, 15))
        self.seq_pool1 = nn.MaxPool2d(kernel_size=(1, 4))
        self.seq_conv2 = nn.Conv2d(160, 160, kernel_size=(1, 12))
        self.seq_pool2 = nn.MaxPool2d(kernel_size=(1, 4))
        self.seq_conv3 = nn.Conv2d(160, 160, kernel_size=(1, 12))
        self.seq_pool3 = nn.MaxPool2d(kernel_size=(1, 4))
        
        # Fully connected layers
        self.flatten = nn.Flatten()
        self.fc1 = nn.Linear(num_gene, 256)  # For gene expression input
        self.fc2 = nn.Linear(1760, 256)      # For sequence features
        self.fc3 = nn.Linear(512, 64)        # Combined features
        self.fc4 = nn.Linear(64, 1)          # Output layer (regression)
        
        # Batch normalization and dropout
        self.BN1 = nn.BatchNorm2d(160)
        self.BN4 = nn.BatchNorm1d(64)
        self.dropout1 = nn.Dropout(0.2)
        self.dropout2 = nn.Dropout(0.5)

    def forward(self, seq_input, gexp_input):
        """Forward pass through the model"""
        # Process sequence input
        seq_input = torch.reshape(seq_input, (-1, 1, 4, SEQ_LEN))
        
        seq_output = self.seq_conv1(seq_input)
        seq_output = self.BN1(seq_output)
        seq_output = nn.functional.relu(seq_output)
        seq_output = self.seq_pool1(seq_output)
        seq_output = self.dropout1(seq_output)
        
        seq_output = self.seq_conv2(seq_output)
        seq_output = self.BN1(seq_output)
        seq_output = nn.functional.relu(seq_output)
        seq_output = self.seq_pool2(seq_output)
        seq_output = self.dropout1(seq_output)
        
        seq_output = self.seq_conv3(seq_output)
        seq_output = self.BN1(seq_output)
        seq_output = nn.functional.relu(seq_output)
        seq_output = self.seq_pool3(seq_output)
        seq_output = self.dropout1(seq_output)
        
        seq_output = self.flatten(seq_output)
        seq_output = self.fc2(seq_output)
        seq_output = nn.functional.relu(seq_output)
        
        # Process gene expression input
        gexp_input = torch.reshape(gexp_input, (-1, NUM_GENE))
        gexp_output = self.fc1(gexp_input)
        gexp_output = nn.functional.relu(gexp_output)
        
        # Combine features
        combined = torch.cat((seq_output, gexp_output), dim=1)
        combined = self.fc3(combined)
        combined = self.BN4(combined)
        combined = nn.functional.relu(combined)
        combined = self.dropout2(combined)
        
        # Output layer
        output = self.fc4(combined)
        output = nn.functional.relu(output)
        return output

def evaluate_regression_model(model, genome, label, tf_gexp_norm, test_cell_idx, region_idx, output_dir, fold_idx):
    """Evaluate the trained regression model on test data"""
    start_time = time.time()
    print("Starting model evaluation...")
    
    # Prepare results dataframe
    results = []
    
    # Process each test cell
    for cellid in test_cell_idx:
        print(f"Processing cell: {cellid}")
        
        # Split regions into batches
        num_batches = math.ceil(len(region_idx) / BATCH_SIZE)
        batch_indices = [(i * BATCH_SIZE, min((i + 1) * BATCH_SIZE, len(region_idx))) 
                        for i in range(num_batches)]
        
        # Process each batch
        for batch_start, batch_end in batch_indices:
            batch_size = batch_end - batch_start
            batch_regions = []
            
            # Prepare input tensors
            X_seq = np.zeros((batch_size, 4, SEQ_LEN))
            X_gexp = np.zeros((batch_size, NUM_GENE))
            Y_true = np.zeros(batch_size)
            
            # Fill batch data
            for i, idx in enumerate(range(batch_start, batch_end)):
                region_info = label.index[region_idx[idx]]
                batch_regions.append(region_info)
                
                # Get sequence and convert to matrix
                sequence = str(genome[region_info])
                X_seq[i] = seq_to_mat(sequence)
                
                # Get gene expression
                X_gexp[i] = tf_gexp_norm.loc[cellid].values
                
                # Get true label
                Y_true[i] = label.iloc[region_idx[idx]].loc[str(cellid)]
            
            # Convert to tensors
            X_seq_tensor = torch.tensor(X_seq, dtype=torch.float32)
            X_gexp_tensor = torch.tensor(X_gexp, dtype=torch.float32)
            
            # Move to GPU if available
            device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
            X_seq_tensor = X_seq_tensor.to(device)
            X_gexp_tensor = X_gexp_tensor.to(device)
            
            # Run model inference
            model.eval()
            with torch.no_grad():
                Y_pred = model(X_seq_tensor, X_gexp_tensor).cpu()
            
            # Extract predictions
            Y_pred_numpy = Y_pred.squeeze().numpy()
            
            # Collect results for this batch
            for i in range(batch_size):
                results.append({
                    'region': batch_regions[i],
                    'cellid': cellid,
                    'Y_pred': Y_pred_numpy[i],
                    'Y_true': Y_true[i]
                })
    
    # Save results to CSV
    results_df = pd.DataFrame(results)
    output_file = os.path.join(output_dir, f'{fold_idx}_regress.csv')
    results_df.to_csv(output_file, index=False)
    
    # Calculate evaluation metrics
    y_true = results_df['Y_true'].values
    y_pred = results_df['Y_pred'].values
    
    mse = mean_squared_error(y_true, y_pred)
    r2 = r2_score(y_true, y_pred)
    pearson_r, _ = stats.pearsonr(y_true, y_pred)
    spearman_r, _ = stats.spearmanr(y_true, y_pred)
    
    # Print evaluation metrics
    print("\nEvaluation Metrics:")
    print(f"Mean Squared Error (MSE): {mse:.4f}")
    print(f"R-squared (R2): {r2:.4f}")
    print(f"Pearson Correlation: {pearson_r:.4f}")
    print(f"Spearman Correlation: {spearman_r:.4f}")
    
    # Save metrics to file
    metrics_file = os.path.join(output_dir, f'{fold_idx}_regress_metrics.txt')
    with open(metrics_file, 'w') as f:
        f.write(f"Mean Squared Error (MSE): {mse:.4f}\n")
        f.write(f"R-squared (R2): {r2:.4f}\n")
        f.write(f"Pearson Correlation: {pearson_r:.4f}\n")
        f.write(f"Spearman Correlation: {spearman_r:.4f}\n")
    
    # Print completion message
    end_time = time.time()
    duration = (end_time - start_time) / 60
    print(f"\nEvaluation completed in {duration:.2f} minutes")
    print(f"Results saved to: {output_file}")
    print(f"Metrics saved to: {metrics_file}")

# ===== 3. Main Execution Block =====
if __name__ == "__main__":
    # Load input data
    print("Loading input data...")
    genome = Fasta(args.genome_fasta)
    label = pd.read_csv(args.label_file, sep='\t', header=0, index_col=[0])
    tf_gexp_norm = pd.read_csv(args.expression_file, sep='\t', header=0, index_col=[0])
    
    # Get cross-validation split
    with open(args.crossval_file) as f:
        crossval_lines = f.readlines()
    cellid_info = crossval_lines[args.fold_idx]
    test_cell_idx = [int(cellid) for cellid in cellid_info.strip().split('\t')[2].split(' ')]
    
    # Set random seeds for reproducibility
    random.seed(1234)
    torch.manual_seed(1234)
    
    # Prepare region indices
    region_idx = random.sample(list(range(label.shape[0])), int(label.shape[0] * RATIO))
    region_idx.sort()
    
    # Load trained model
    print("Loading trained model...")
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = torch.load(args.model_path, map_location=device)
    model = model.to(device)
    
    # Evaluate model
    evaluate_regression_model(
        model=model,
        genome=genome,
        label=label,
        tf_gexp_norm=tf_gexp_norm,
        test_cell_idx=test_cell_idx,
        region_idx=region_idx,
        output_dir=args.output_dir,
        fold_idx=args.fold_idx
    )
    
    print("Evaluation completed successfully")