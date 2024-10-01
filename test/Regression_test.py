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
from torch.utils.data import Dataset, DataLoader
from torch.utils.tensorboard import SummaryWriter
from torch.optim.lr_scheduler import ExponentialLR
from sklearn import metrics
from sklearn.metrics import roc_auc_score,roc_curve,precision_recall_curve,auc,average_precision_score
from scipy import stats
import math

fold_idx = int(sys.argv[1])
input_dir = sys.argv[2]
output_dir = sys.argv[3]

######################## Global Settings #######################
os.environ["CUDA_LAUNCH_BLOCKING"] = "1"
SEQ_LEN = 1000
NUM_GENE = 2500

###################### Model Implementation #####################
# One-hot coding
def seq_to_mat(seq):
    d = {'a': 0, 'A': 0, 'g': 1, 'G': 1, 'c': 2, 'C': 2, 't': 3, 'T': 3, 'N': 4, 'n': 4}
    mat = np.zeros((5, SEQ_LEN))
    for i in range(len(seq)):
        mat[d[seq[i]], i] = 1
    mat = mat[:4, :]
    return mat

class MyModel(nn.Module):
    def __init__(self, NUM_GENE, SEQ_LEN):
        super(MyModel, self).__init__()
        self.conv1_params = {'nb_filter': 160, 'filter_size': (4, 15), 'pool_size': (1, 4)}
        self.conv2_params = {'nb_filter': 160, 'filter_size': (1, 12), 'pool_size': (1, 4)}
        self.conv3_params = {'nb_filter': 160, 'filter_size': (1, 12), 'pool_size': (1, 4)}
        
        self.seq_conv1 = nn.Conv2d(1, self.conv1_params['nb_filter'], kernel_size=self.conv1_params['filter_size'])
        self.seq_pool1 = nn.MaxPool2d(kernel_size=self.conv1_params['pool_size'])
        
        self.seq_conv2 = nn.Conv2d(self.conv1_params['nb_filter'], self.conv2_params['nb_filter'], kernel_size=self.conv2_params['filter_size'])
        self.seq_pool2 = nn.MaxPool2d(kernel_size=self.conv2_params['pool_size'])
        
        self.seq_conv3 = nn.Conv2d(self.conv2_params['nb_filter'], self.conv3_params['nb_filter'], kernel_size=self.conv3_params['filter_size'])
        self.seq_pool3 = nn.MaxPool2d(kernel_size=self.conv3_params['pool_size'])
        
        self.flatten = nn.Flatten()
        self.fc1 = nn.Linear(2500, 256)
        self.fc2 = nn.Linear(1760, 256)
        self.fc3 = nn.Linear(256 + 256, 64)
        self.fc4 = nn.Linear(64, 1)
        
        self.BN1 = nn.BatchNorm2d(160)
        self.BN4 = nn.BatchNorm1d(64)
        
        self.dropout1 = nn.Dropout(0.2)
        self.dropout2 = nn.Dropout(0.5)

    def forward(self, seq_input, gexp_input):
        seq_input = seq_input.float()  
        gexp_input = gexp_input.float()
        seq_input = torch.reshape(seq_input,(-1,1,4,SEQ_LEN))
        gexp_input = torch.reshape(gexp_input,(-1,NUM_GENE))
        gexp_input = self.fc1(gexp_input)
        gexp_input = nn.functional.relu(gexp_input)
        
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

        fc1_output = nn.functional.relu(seq_output)
        fc1_output = self.fc2(fc1_output)
        
        concat_output = torch.cat((fc1_output, gexp_input), dim=1)
        fc2_output = self.fc3(concat_output)
        fc2_output = self.BN4(fc2_output)
        fc2_output = nn.functional.relu(fc2_output)
        fc2_output = self.dropout2(fc2_output)
        output = self.fc4(fc2_output)
        output = nn.functional.relu(output) 

        return output

################################## Model Evaluation##################################
def model_evaluation(model, batch_size):
    start_time = time.time()
    results_df = pd.DataFrame(columns=["region", 'cellid', 'Y_pred', 'Y_true'])
    region = []

    num_batches = math.ceil(len(region_idx) / batch_size)
    batch_indices = [(i * batch_size, min((i + 1) * batch_size, len(region_idx))) for i in range(num_batches)]

    for cellid in test_cell_idx:
        results_batch_df = pd.DataFrame(columns=["region", 'cellid', 'Y_pred', 'Y_true'])
        
        for batch_start, batch_end in batch_indices:
            batch_region = []
            X_test_mat = np.zeros((batch_end - batch_start, 4, SEQ_LEN, 1))
            X_test_gexp = np.zeros((batch_end - batch_start, NUM_GENE))
            Y_true = np.zeros(batch_end - batch_start)

            for i, idx in enumerate(range(batch_start, batch_end)):
                info = label.index[region_idx[idx]]
                batch_region.append(info)
                chrom, start, end = info.split(':')[0], int(info.split(':')[1].split('-')[0]), int(
                    info.split(':')[1].split('-')[1])
                sequence = str(genome[info])
                X_test_mat[i, :, :, 0] = seq_to_mat(sequence)
                X_test_gexp[i, :] = tf_gexp_norm.loc[cellid].values.reshape((1, -1))
                Y_true[i] = label.iloc[region_idx[idx]].loc[str(cellid)]

            X_test_mat = torch.tensor(X_test_mat, dtype=torch.float32)
            X_test_gexp = torch.tensor(X_test_gexp, dtype=torch.float32)
            Y_true = torch.tensor(Y_true, dtype=torch.float32)

            device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
            model.eval()
            Y_pred_list = []
            with torch.no_grad():
                X_batch_mat = X_test_mat.to(device)
                X_batch_gexp = X_test_gexp.to(device)
                Y_batch = Y_true.to(device)

                Y_pred_batch = model(X_batch_mat, X_batch_gexp).unsqueeze(1).cpu()
                Y_pred_list.append(Y_pred_batch)

            Y_pred = torch.cat(Y_pred_list, dim=0)

            Y_pred_numpy = Y_pred.squeeze().cpu().numpy()
            Y_true_numpy = Y_true.squeeze().cpu().numpy()

            for i in range(batch_end - batch_start):
                results_batch_df = results_batch_df.append(
                    {'region': batch_region[i], 'cellid': cellid, 'Y_pred': Y_pred_numpy[i], 'Y_true': Y_true_numpy[i]},
                    ignore_index=True)

        results_df = pd.concat([results_df, results_batch_df], ignore_index=True)

    results_df.to_csv('%s/%d_regress.csv' % (output_dir,fold_idx), index=False)

    end_time = time.time()
    print("time consuming：" + str((end_time - start_time) / 60) + "min")
    print('End testing...')

if  __name__ == "__main__" :
    genome = Fasta("%s/union.peaks.pad1k.fa"%input_dir) 
    label_file='%s/union.peaks.labels.regress.txt'%input_dir
    tf_gexp_file = '%s/final_exp_matrix.txt'%input_dir
    tf_gexp_norm = pd.read_csv(tf_gexp_file ,sep='\t',header=0,index_col=[0])
    label = pd.read_csv(label_file,sep='\t',header=0,index_col=[0])
    fold_idx = int(sys.argv[1])
    ratio = 1
    if ratio<0 or ratio>1: 
        print('Input ratio between 0 and 1')
        sys.exit()    
    cellid_info = open('%s/Leave_one_out_cross_validation.txt'%input_dir).readlines()[fold_idx] 
    #train_cell_idx = [int(cellid) for cellid in cellid_info.split('\t')[1].split(' ')]
    test_cell_idx = [int(cellid) for cellid in cellid_info.strip().split('\t')[2].split(' ')]
    random.seed(1234) 
    torch.manual_seed(1234)
    region_idx = random.sample(list(range(label.shape[0])),int(label.shape[0]*ratio))
    region_idx.sort()
    model = torch.load('%s/regress_model.pth' % (input_dir))
    print('Start testing...')
    model_evaluation(model,batch_size=12000)