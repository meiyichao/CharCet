#!/usr/bin/env python3
# ATAC-seq Regression Model Training Pipeline
# Trains a convolutional neural network to predict chromatin accessibility levels

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
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader, Subset
from torch.utils.data.sampler import WeightedRandomSampler
from torch.optim.lr_scheduler import ExponentialLR, ReduceLROnPlateau
from sklearn import metrics
from sklearn.metrics import roc_auc_score, precision_recall_curve, auc, average_precision_score

# Parse command line arguments
parser = argparse.ArgumentParser(description='ATAC-seq Regression Model Training Pipeline')
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
parser.add_argument('-s', '--seq_len', type=int, default=1000, 
                    help='Sequence length (default: 1000)')
parser.add_argument('-n', '--num_gene', type=int, default=2500, 
                    help='Number of genes (default: 2500)')
parser.add_argument('-b', '--batch_size', type=int, default=128, 
                    help='Batch size for training (default: 128)')
parser.add_argument('-p', '--patience', type=int, default=12, 
                    help='Early stopping patience (default: 12)')
parser.add_argument('-r', '--ratio', type=float, default=1.0, 
                    help='Sampling ratio of regions (default: 1.0)')
args = parser.parse_args()

# Set global parameters from arguments
SEQ_LEN = args.seq_len
NUM_GENE = args.num_gene
BATCH_SIZE = args.batch_size
PATIENCE = args.patience
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
        return output

class ChromatinDataset(Dataset):
    """Dataset class for chromatin accessibility regression data"""
    def __init__(self, training_queue, genome, label, tf_gexp_norm):
        self.training_queue = training_queue
        self.genome = genome
        self.label = label
        self.tf_gexp_norm = tf_gexp_norm

    def __len__(self):
        return len(self.training_queue)

    def __getitem__(self, idx):
        """Get a single data sample"""
        region_idx, cellid = self.training_queue[idx]
        info = self.label.index[region_idx]
        sequence = str(self.genome[info])
        X_mat = seq_to_mat(sequence)
        X_gexp = self.tf_gexp_norm.loc[cellid].values.reshape((1, -1))
        Y = self.label.iloc[region_idx].loc[str(cellid)]
        return (X_mat, X_gexp), Y

def train_model(model, data_queue, genome, label, tf_gexp_norm, output_dir):
    """Train the chromatin accessibility regression model"""
    # Split data into training and validation sets
    train_queue = data_queue[:int(0.95 * len(data_queue))]
    valid_queue = data_queue[int(0.95 * len(data_queue)):]
    print(f'{len(train_queue)} training and {len(valid_queue)} validation examples')
    
    # Create datasets
    train_dataset = ChromatinDataset(train_queue, genome, label, tf_gexp_norm)
    valid_dataset = ChromatinDataset(valid_queue, genome, label, tf_gexp_norm)
    
    # Create data loaders
    train_dataloader = DataLoader(
        train_dataset, 
        batch_size=BATCH_SIZE, 
        num_workers=8,
        drop_last=True
    )
    
    valid_dataloader = DataLoader(
        valid_dataset, 
        batch_size=BATCH_SIZE, 
        shuffle=True,
        num_workers=8
    )
    
    # Set up training components
    criterion = nn.MSELoss().cuda()
    optimizer = optim.Adam(model.parameters(), lr=0.001, weight_decay=0.005)
    model = model.cuda()
    
    # Early stopping configuration
    best_validation_loss = float('inf')
    early_stopping_counter = 0
    
    # Training loop
    start_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    print(f"Training started at: {start_time}")
    
    for epoch in range(1, 101):  # Train for maximum 100 epochs
        print(f"---------- Epoch {epoch}/100 ----------")
        
        # Training phase
        model.train()
        train_loss = 0.0
        for batch_data, batch_labels in train_dataloader:
            seq_input, gexp_input = batch_data
            seq_input = seq_input.cuda()
            gexp_input = gexp_input.cuda()
            batch_labels = batch_labels.float().unsqueeze(1).cuda()
            
            optimizer.zero_grad()
            outputs = model(seq_input, gexp_input)
            loss = criterion(outputs, batch_labels)
            loss.backward()
            optimizer.step()
            train_loss += loss.item()
        
        train_loss /= len(train_dataloader)
        
        # Validation phase
        model.eval()
        valid_loss = 0.0
        with torch.no_grad():
            for batch_data, batch_labels in valid_dataloader:
                seq_input, gexp_input = batch_data
                seq_input = seq_input.cuda()
                gexp_input = gexp_input.cuda()
                batch_labels = batch_labels.float().unsqueeze(1).cuda()
                
                outputs = model(seq_input, gexp_input)
                loss = criterion(outputs, batch_labels)
                valid_loss += loss.item()
        
        valid_loss /= len(valid_dataloader)
        
        # Early stopping check
        if valid_loss < best_validation_loss - 0.0005:
            best_validation_loss = valid_loss
            early_stopping_counter = 0
            torch.save(model, os.path.join(output_dir, 'regress_model.pth'))
            print(f"Saved improved model with validation loss: {valid_loss:.4f}")
        else:
            early_stopping_counter += 1
            if early_stopping_counter >= PATIENCE:
                print(f"Early stopping triggered after {PATIENCE} epochs without improvement")
                break
        
        print(f'Epoch {epoch} - Train Loss: {train_loss:.4f} - Valid Loss: {valid_loss:.4f}')
        print(f"Current time: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}\n")
    
    return model

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
    train_cell_idx = [int(cellid) for cellid in cellid_info.split('\t')[1].split(' ')]
    
    # Set random seeds for reproducibility
    random.seed(1234)
    torch.manual_seed(1234)
    
    # Prepare data queue
    region_idx = random.sample(list(range(label.shape[0])), int(label.shape[0] * RATIO))
    data_queue = []
    for idx in region_idx:
        for cellid in train_cell_idx:
            data_queue.append((idx, cellid))
    random.shuffle(data_queue)
    
    # Initialize and train model
    print("Initializing model...")
    model = ChromatinAccessibilityRegressor(NUM_GENE, SEQ_LEN)
    
    print("Starting training...")
    trained_model = train_model(model, data_queue, genome, label, tf_gexp_norm, args.output_dir)
    
    print("Training completed successfully")