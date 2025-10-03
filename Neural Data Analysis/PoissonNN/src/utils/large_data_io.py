# -*- coding: utf-8 -*-
"""
Use this to load big files
"""
import h5py #pip install h5py
import numpy as np
import os
import torch
from torch.utils.data import Dataset, DataLoader, random_split

def create_memmap_from_mat(datafiles, out_dir="memmap_data"):
    """
    Stacks multiple v7.3 .mat files into memory-mapped arrays.
    
    Parameters:
        datafiles: list of paths to .mat files
        out_dir: directory to store memmap files
    Returns:
        paths to memmap files (X_path, y_path, cell_ids_path)
    """
    os.makedirs(out_dir, exist_ok=True)

    # First, determine total number of rows and feature size
    total_rows = 0
    with h5py.File(datafiles[0], "r") as f:
        n_features = f['X'].shape[1]
    for fpath in datafiles:
        with h5py.File(fpath, "r") as f:
            total_rows += f['X'].shape[0]

    # Create memory-mapped arrays on disk
    X_path = os.path.join(out_dir, "X_stacked.dat")
    y_path = os.path.join(out_dir, "y_stacked.dat")
    cell_ids_path = os.path.join(out_dir, "cell_ids_stacked.dat")

    X_mm = np.memmap(X_path, dtype='float32', mode='w+', shape=(total_rows, n_features))
    y_mm = np.memmap(y_path, dtype='float32', mode='w+', shape=(total_rows, 1))
    cell_ids_mm = np.memmap(cell_ids_path, dtype='int64', mode='w+', shape=(total_rows, 1))

    # Fill memory-mapped arrays
    start = 0
    num_cells = 0
    for fpath in datafiles:
        with h5py.File(fpath, "r") as f:
            n_rows = f['X'].shape[0]
            X_mm[start:start+n_rows] = f['X'][:]
            y_mm[start:start+n_rows] = f['y'][:]
            cell_ids_mm[start:start+n_rows] = f['cell_ids'][:]
            start += n_rows
            num_cells = max(cell_ids_mm).item()

    # Flush to disk
    X_mm.flush()
    y_mm.flush()
    cell_ids_mm.flush()

    return X_path, y_path, cell_ids_path, n_features, total_rows, num_cells


# -------------------------
# PyTorch Dataset for memmap arrays
class MemmapDataset(Dataset):
    def __init__(self, X_path, y_path, cell_ids_path, n_features, total_rows):
        self.X = np.memmap(X_path, dtype='float32', mode='r', shape=(total_rows, n_features))
        self.y = np.memmap(y_path, dtype='float32', mode='r', shape=(total_rows, 1))
        self.cell_ids = np.memmap(cell_ids_path, dtype='int64', mode='r', shape=(total_rows, 1))

    def __len__(self):
        return self.X.shape[0]

    def __getitem__(self, idx):
        return (torch.tensor(self.X[idx], dtype=torch.float32),
                torch.tensor(self.cell_ids[idx], dtype=torch.int64),
                torch.tensor(self.y[idx], dtype=torch.float32))

# -------------------------
# Generate dataloaders for training, testing and prediction

def prepare_data(dataset, batch_size):
    # Load all data
    alldata_loader = DataLoader(dataset, batch_size, shuffle=False)
    # Split train and test set
    train_size = int(0.5 * len(dataset))  
    test_size = len(dataset) - train_size  
    train_dataset, test_dataset = random_split(dataset,[train_size, test_size])
    train_loader = DataLoader(train_dataset, batch_size, shuffle=False)
    test_loader  = DataLoader(test_dataset, batch_size, shuffle=False) 
    return alldata_loader, train_loader, test_loader
