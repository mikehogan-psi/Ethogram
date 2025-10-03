# -*- coding: utf-8 -*-
"""
This is for a basic plot of results for individual cells
"""

import torch
import pandas as pd # Install
import matplotlib.pyplot as plt # Install

def plot_results(y, y_pred, cell_id, do_save=False):
    
    df = pd.DataFrame({'y_obs': y.numpy().flatten(), 
                       'y_pred': y_pred.detach().numpy().flatten()})
    
    plt.figure(dpi=300)
    plt.plot(df.index, df['y_obs'], 'k', label = 'observed')
    plt.plot(df.index, df['y_pred'], 'r', label = 'predicted')
    plt.xlabel("#Sample")
    plt.ylabel("#Counts")
    plt.title("Observed vs Predicted; cell#" + str(cell_id))
    plt.legend()  # adds the legend using the labels
    plt.show()
    if do_save == True:
        plt.savefig(["prediction cell#" + str(cell_id) + ".svg"], format="svg")