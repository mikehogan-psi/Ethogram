# -*- coding: utf-8 -*-
"""
These functions load and save the data
Use with config.data['data_dir']
"""
from scipy.io import loadmat
import numpy as np
import torch
import os
from torch.utils.data import TensorDataset, Dataset, DataLoader, random_split

def load_data(datafiles):
    n = 0
    for name in datafiles:
        data = loadmat(name)
        if n == 0:
            X = data['X']
            y = data['y']
            cell_ids = data['cell_ids']
        else: 
            X = np.vstack((X,data['X']))
            y = np.vstack((y,data['y']))
            cell_ids = np.vstack((cell_ids,data['cell_ids']))
        n += 1    
        #print(f"{X.shape}")
    return X, cell_ids, y     

def prepare_data(X, cell_ids, y, batch_size):
    # Wrap into dataset
    X = torch.tensor(X, dtype=torch.float32)
    cell_ids = torch.tensor(cell_ids, dtype=torch.int16)
    y = torch.tensor(y, dtype=torch.float32)
    dataset = TensorDataset(X, cell_ids, y)
    # Load all data
    # all_data_loader = DataLoader(dataset, batch_size, shuffle=False)
    # Split train and test set
    train_size = int(0.5 * len(dataset))  
    test_size = len(dataset) - train_size  
    train_dataset, test_dataset = random_split(dataset,[train_size, test_size])
    train_loader = DataLoader(train_dataset, batch_size, shuffle=True)
    test_loader  = DataLoader(test_dataset, batch_size) 
    return train_loader, test_loader
    
        
        
        
    