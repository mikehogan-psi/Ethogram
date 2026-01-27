# -*- coding: utf-8 -*-
"""
This is the core PoissonNN model
"""

import numpy as np
import torch
import torch.nn as nn

class PoissonNN(nn.Module):
    
    # Example: PoissonNN(5, [10,10], 10)
    def __init__(self, input_dim, hidden_dim_list, num_cells, dropout=0.):        
        super().__init__() 

        # Create a ModuleList of hidden layers
        self.hidden_layers = nn.ModuleList()
        self.dropouts = nn.ModuleList()
        
        # First hidden layer
        self.hidden_layers.append(nn.Linear(input_dim, hidden_dim_list[0]))
        self.dropouts.append(nn.Dropout(dropout))
        
        # Additional hidden layers
        for i in range(1, len(hidden_dim_list)):
            self.hidden_layers.append(nn.Linear(hidden_dim_list[i-1], hidden_dim_list[i]))
            self.dropouts.append(nn.Dropout(dropout))
        
        # Final hidden dimension = output size of last hidden layer
        final_hidden_dim = hidden_dim_list[-1]
        
        # Activation function
        self.activation = nn.ReLU()  
        
        # Cell-specific heads
        self.heads = nn.ModuleList([nn.Linear(final_hidden_dim, 1) for _ in range(num_cells)])
        
        # Add skip connections
        if input_dim != final_hidden_dim:
            self.skip_proj = nn.Linear(input_dim, final_hidden_dim)
        else:
            self.skip_proj = nn.Identity()
    
    def forward(self, x, cell_ids):
        """
        x: [batch_size, input_dim]
        cell_ids: [batch_size] tensor
        """
        # Shared network
        h = x
        for layer, dropout in zip(self.hidden_layers, self.dropouts):
            h = self.activation(layer(h))
            h = dropout(h)
        
        # Add skip connection from inputs
        h = h + self.skip_proj(x)
        
        # Individual cells outputs
        mu_list = []
        indices_list = []
        for cid in torch.unique(cell_ids):
            mask = (cell_ids == cid)
            indices = torch.where(mask)[0]
            h_c = h[indices]
            mu_c = torch.exp(self.heads[cid](h_c))
            mu_list.append(mu_c)
            indices_list.append(indices)
        
        # Preallocate on the same device & dtype as x
        # mu = x.new_empty(x.size(0), 1)
    
        # # Compute per-cell heads
        # for cid in torch.unique(cell_ids).tolist():        # -> Python ints
        #     mask = (cell_ids == cid)
        #     inds = torch.where(mask)[0]                    # same device as x
        #     h_c = h[inds]
        #     # Use exp to enforce positivity; softplus is more stable if needed
        #     mu_c = torch.exp(self.heads[cid](h_c))
        #     mu[inds] = mu_c
        
        # Reconstruct outputs in original order
        #mu = torch.empty(x.size(0), 1, dtype=h.dtype)
        mu = torch.empty(x.size(0), 1, dtype=h.dtype, device=x.device) #GPU
        for inds, mu_c in zip(indices_list, mu_list):
            mu[inds] = mu_c
        
        return mu
  
    def poisson_nll(self, y_pred, y_true):
        y_true = y_true.float().view_as(y_pred)
        #log_factorial = torch.lgamma(y_true + 1)
        #nll = -(y_true * torch.log(y_pred + 1e-8) - y_pred - log_factorial).mean()
        nll = -(y_true * torch.log(y_pred + 1e-8) - y_pred).mean()
        return nll

    def train_model(self, dataloader, max_epochs=1000, learning_rate=1e-3, 
                    patience_length=5, checkpoints=[], net_name='test'):
        # Select device
        # device = torch.device("cpu")
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu") #GPU
        print("Using device:", device)
        self.to(device)
        
        # Set training and optimising parameters
        optimizer = torch.optim.Adam(self.parameters(), lr=learning_rate)
        prev_loss = float('inf')
        tol = 1e-4
        patience = 0
        
        for epoch in range(max_epochs):
            epoch_loss = 0.0
            total_samples = 0
            for x_batch, cell_batch, y_batch in dataloader:
                # Move batches to device
                x_batch = x_batch.to(device)
                cell_batch = cell_batch.to(device)
                y_batch = y_batch.to(device)
                
                # Optimise
                optimizer.zero_grad()
                y_pred = self.forward(x_batch, cell_batch)
                loss = self.poisson_nll(y_pred, y_batch)
                loss.backward()
                optimizer.step()
                
                epoch_loss += loss.item() * x_batch.size(0)
                total_samples += x_batch.size(0)
             print(f"Epoch Loss: {epoch_loss}")
            
            epoch_loss /= total_samples
            
            if epoch in checkpoints:
                torch.save(self, f"{net_name}_checkpoint{epoch}.pth")
            
            if abs(prev_loss - epoch_loss) < tol:
                patience += 1
                if patience == patience_length:
                    print(f"Loss < {tol}, stopping at epoch {epoch}")
                    break
            else:
                prev_loss = epoch_loss
                patience = 0
            
            print(f"Epoch {epoch+1}/{max_epochs}, Loss: {epoch_loss:.4f}")

    def test_model(self, dataloader, cell_id):
        ################GPU START######################
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu") #GPU
        print("Using device:", device)
        self.to(device)
        ################GPU END######################
        self.eval()
        y = []
        y_pred = []
        X = []
        cell_id = torch.tensor(cell_id)
        
        # Retrieve the original values
        with torch.no_grad():
            for x_batch, cell_batch, y_batch in dataloader:
                mask = (cell_batch == cell_id)
                indices = torch.where(mask)[0]
                X.append(x_batch[indices,:])
                y.append(y_batch[indices])
                y_pred.append(self.forward(x_batch[indices,:], cell_batch[indices]))
        
        y = torch.cat(y, dim=0)
        y_pred = torch.cat(y_pred, dim=0)
        
        # Compute likelihood & pR2
        loss = self.poisson_nll(y_pred, y)
        mu_null = torch.full_like(y, y.mean())
        loss_null = self.poisson_nll(mu_null, y)
        loss_full = self.poisson_nll(y, y)
        pR2 = 1 - (loss - loss_full) / (loss_null - loss_full)
        print(f"cell: {cell_id}: pR2 = {pR2:.4f}")
        
        return pR2, y, y_pred, X

    def extract_weights_biases(self):
        weights = np.array([head.weight.detach().cpu().numpy().squeeze(0) for head in self.heads])
        biases = np.array([head.bias.detach().cpu().numpy().item() for head in self.heads])
        return weights, biases
