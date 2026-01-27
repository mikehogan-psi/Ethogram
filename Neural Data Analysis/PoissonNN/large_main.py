# -*- coding: utf-8 -*-
"""
Main script to predict data. 
"""

import sys
import os
import torch
from scipy.io import savemat #install
#from pathlib import Path #install
import numpy as np
import importlib

# HERE CHANGE TO ADD YOUR DIRECTORY PATH
fpath = "C:\\Users\\G71044MH\\OneDrive - The University of Manchester\\Documents\\GitHub\\Ethogram\\Neural Data Analysis\\PoissonNN"
sys.path.append(fpath)
os.chdir(fpath)

from src.models import large_Poisson_NN_Model as PNN
from src.utils import config, large_data_io, show_results

# Reload the modules, in case of new edits
# ROOT = Path(os.getcwd())
# sys.path.append(str(ROOT))
importlib.reload(PNN)
importlib.reload(config)
importlib.reload(large_data_io) 

# Load configuration file
cfg = config.Config(fpath + "\\configs\\config.yaml")

print("Loading configuration file: ...")
print(repr(cfg))

# Load data and wrap into a Dataloader
print("Loading data: ...")
datafiles = cfg.data['data_dir']
batch_size = cfg.training['batch_size']

X_path, y_path, cell_ids_path, n_features, total_rows, num_cells = \
    large_data_io.create_memmap_from_mat(datafiles, out_dir=fpath)
    
num_cells += 1 

dataset = large_data_io.MemmapDataset(X_path, y_path, cell_ids_path, \
    n_features=n_features, total_rows=total_rows)
    
alldata_loader, train_loader, test_loader = \
    large_data_io.prepare_data(dataset, batch_size)

# Init Model
input_dim = n_features
hidden = cfg.model['hidden_dim_list']
drpout = cfg.training['dropout']
model = PNN.PoissonNN(input_dim, hidden, num_cells, drpout)

# Train model
mepchs = cfg.training['max_epochs']
lr = cfg.training['learning_rate']
pl = cfg.training['patience_length']
chkout = cfg.training['checkout_epochs']
nname = cfg.model['network_name']

#model initialised on cpu, here move the model to cuda
#device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
#model.to(device)

# train the model
model.train_model(train_loader, max_epochs=mepchs, learning_rate=lr,
                  patience_length=pl, checkpoints=chkout, net_name=nname)

# Test model
pR2 = np.zeros(num_cells)
for n in range(num_cells):
    pR2[n], _, _, _ = model.test_model(test_loader, n)

# Extract weights & biases
weights, biases = model.extract_weights_biases()

# Save model
# torch.save(model, str(ROOT) + "\\results\\full_model.pth")
torch.save(model, fpath + "\\results\\full_model.pth")

# Save cfg, predictions, pR2, model weights and biases
y_all = []
y_pred_all = []
for n in range(num_cells):
    _, y, y_pred, _ = model.test_model(alldata_loader, n)
    y_all.append(y)
    y_pred_all.append(y_pred)
results = {'config': cfg, 'pseudoR2': pR2, 'y_pred': y_pred_all,
           'y': y_all, 'weights': weights, 'biases': biases}
# savemat(str(ROOT) + "\\results\\results.mat", results) 
savemat(fpath + "\\results\\results.mat", results) 

# Example: Get single-cell predictions 
cell_id = 2
_, y, y_pred, _ = model.test_model(alldata_loader, cell_id)    
# Show single-cell results
show_results.plot_results(y, y_pred, cell_id, do_save=False)


