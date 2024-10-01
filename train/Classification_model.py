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
from torch.utils.data.sampler import WeightedRandomSampler
from torch.utils.tensorboard import SummaryWriter
from torch.optim.lr_scheduler import ExponentialLR,ReduceLROnPlateau
from sklearn.model_selection import train_test_split
from sklearn import metrics
from sklearn.metrics import roc_auc_score,precision_recall_curve,auc,average_precision_score
from sklearn import preprocessing

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
        self.fc4 = nn.Linear(64, 2)
        
        self.BN1 = nn.BatchNorm2d(160)
        self.BN4 = nn.BatchNorm1d(64)
        
        self.dropout1 = nn.Dropout(0.2)
        self.dropout2 = nn.Dropout(0.5)
        self.sigmoid = nn.Sigmoid()

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
        fc3_output = self.fc4(fc2_output)
        output = self.sigmoid(fc3_output)
        return output

class MyDataset(Dataset): 
    def __init__(self, training_queue, bootstrap=False):
        self.training_queue = training_queue
        self.bootstrap = bootstrap

    def __len__(self):
        return len(self.training_queue)

    def __getitem__(self, idx):
        region_idx, cellid = self.training_queue[idx]
        info = label.index[region_idx]
        sequence = str(genome[info])
        X_mat = seq_to_mat(sequence)
        X_gexp = tf_gexp_norm.loc[cellid].values.reshape((1, -1))
        Y = label.iloc[region_idx].loc[str(cellid)]
        return (X_mat, X_gexp), Y

################################## Model Training##################################
def model_training(model, batch_size, epochs):
    data_queue = []
    for idx in region_idx:
        for cellid in train_cell_idx:
            data_queue.append((idx, cellid))
    random.shuffle(data_queue)
    train_queue = data_queue[:int(0.95 * len(data_queue))]
    valid_queue = data_queue[int(0.95 * len(data_queue)):]
    print('{} training and {} validation examples'.format(len(train_queue), len(valid_queue)))

    train_dataset = MyDataset(train_queue, bootstrap=False)
    valid_dataset = MyDataset(valid_queue, bootstrap=False)

    train_samples = DataLoader(train_dataset, batch_size=batch_size, shuffle=False,num_workers=8,drop_last=False)
    
    num_classes = 2
    class_counts = torch.zeros(num_classes)
    for _, labels in train_samples:
        class_counts += torch.bincount(labels, minlength=num_classes)
    print(class_counts)

    total_samples = class_counts.sum()
    class_weights = total_samples / (num_classes * class_counts)
    print(class_weights)

    weights = []
    for data, label in train_samples:
        for cla in label:
            if cla ==0:
                weights.append(class_weights[0].item())
            else:
                weights.append(class_weights[1].item())
    print(len(weights))
    print(weights[:10])
    weights = pd.DataFrame(weights,columns=["weights"])
    weights = torch.tensor(weights["weights"].to_numpy(),dtype=torch.float)
    sampler = WeightedRandomSampler(weights,num_samples=len(train_dataset), replacement=True)
    
    train_dataloader = DataLoader(train_dataset, batch_size=batch_size, sampler = sampler,num_workers=8,drop_last=True)
    valid_dataloader = DataLoader(valid_dataset, batch_size=batch_size, shuffle=True,num_workers=8)

    X, y = next(iter(train_dataloader))
    print(y)
    criterion = nn.CrossEntropyLoss()
    criterion = criterion.cuda()
    optimizer = optim.Adam(model.parameters(), lr=0.001, weight_decay=0.005)

    if torch.cuda.is_available():
        model =model.cuda()

    patience = 12  # If there is no improvement in the validation set metrics for 12 consecutive epochs, stop training
    min_delta = 0.0005  
    early_stopping_counter = 0 
    best_validation_loss = np.Inf 
    
    total_train_step = 0
    total_test_step = 0

    start_time = time.time()
    start_time = time.localtime(start_time)
    start_time = time.strftime("%Y-%m-%d %H:%M:%S", start_time)
    print("start_time: ", start_time)
    for epoch in range(1,epochs+1):
        print("----------The {} training----------".format(epoch))
        train_loss = 0.0
        model.train()
        for batch_data, batch_labels in train_dataloader:
            seq_input, gexp_input = batch_data
            seq_input = seq_input.cuda()
            gexp_input = gexp_input.cuda()
            batch_labels = batch_labels.cuda()

            optimizer.zero_grad()
            outputs = model(seq_input, gexp_input)
            loss = criterion(outputs, batch_labels)
            loss.backward()
            optimizer.step()

            train_loss += loss.item()

        total_train_step += 1
        train_loss /= len(train_dataloader)

        valid_loss = 0.0
        model.eval()
        with torch.no_grad():
            for batch_data, batch_labels in valid_dataloader:
                seq_input, gexp_input = batch_data
                seq_input = seq_input.cuda()
                gexp_input = gexp_input.cuda()
                batch_labels = batch_labels.cuda()
                
                outputs = model(seq_input, gexp_input)
                loss = criterion(outputs, batch_labels)
                valid_loss += loss.item()
                
            total_test_step += 1
            valid_loss /= len(valid_dataloader)

        if valid_loss < best_validation_loss - min_delta:
            best_validation_loss = valid_loss
            early_stopping_counter = 0

            torch.save(model,'%s/class_model.pth'%(output_dir))

        else:
            early_stopping_counter += 1
            if early_stopping_counter >= patience:
                print("Early stopping triggered")
                break

        print('Epoch {}/{} - Train Loss: {:.4f} - Valid Loss: {:.4f}'.format(epoch, epochs, train_loss, valid_loss))
    
        end_time = time.time()
        end_time = time.localtime(end_time)
        end_time = time.strftime("%Y-%m-%d %H:%M:%S", end_time)
        print("current time: ", end_time,"\n")
    return model

if  __name__ == "__main__" :
    genome = Fasta("%s/union.peaks.pad1k.fa"%input_dir) 
    label_file='%s/union.peaks.labels.class.txt'%input_dir
    tf_gexp_file = '%s/final_exp_matrix.txt'%input_dir
    tf_gexp_norm = pd.read_csv(tf_gexp_file ,sep='\t',header=0,index_col=[0])
    label = pd.read_csv(label_file,sep='\t',header=0,index_col=[0])
    ratio = 1
    if ratio<0 or ratio>1: 
        print('Input ratio between 0 and 1')
        sys.exit()    
    cellid_info = open('%s/Leave_one_out_cross_validation.txt'%input_dir).readlines()[fold_idx] 
    train_cell_idx = [int(cellid) for cellid in cellid_info.split('\t')[1].split(' ')]
    #test_cell_idx = [int(cellid) for cellid in cellid_info.strip().split('\t')[2].split(' ')]
    random.seed(1234) 
    torch.manual_seed(1234)
    region_idx = random.sample(list(range(label.shape[0])),int(label.shape[0]*ratio))
    model = MyModel(NUM_GENE, SEQ_LEN)
    print('Start training...')
    model_training(model, batch_size=128*1, epochs=100)