# -*- coding: utf-8 -*-
"""
Use LMMs to test difference in embedding dimensions
among groups (or brain areas)
This is also call a config.yaml to determine the hierarchical 
structure of the data 
"""

import pandas as pd
import statsmodels.formula.api as smf
import numpy as np

class embed_data():
   
    def __init__(self, subjects, groups, embeds):
        self.df = pd.DataFrame({
            "subject": subjects,
            "group": groups})
        Ndim = embeds.shape[1]
        for n in range(Ndim):
            name = "embed" + str(n+1)
            self.df[name] = embeds[:,n].flatten()
            
    def LMM(self):
        self.results = []
        for name, data in self.df.items():
            if "embed" in name:
                model = smf.mixedlm(f"{name} ~ group", self.df, groups=self.df["subject"])
                result = model.fit()
                # Extract coefficient + p-value for group effect
                coef = result.params.filter(like="group").values[0]
                pval = result.pvalues.filter(like="group").values[0]
    
                self.results.append({
                    "dimension": name,
                    "coef_group": coef,
                    "pval_group": pval
                })
            
#Generate some data
Nsubj = 10
Ngroup = 2
Ncell_subj = 10
Ncell = Nsubj*Ncell_subj
Ncell_group = Ncell/Ngroup
Ndim = 5
subj_bias = np.random.randn(Nsubj, Ndim)
#group_bias = np.random.randn(Ngroup, Ndim)
group_bias = np.repeat(np.transpose(np.matrix([0, 5])), Ndim, axis=1)
embeds = np.random.randn(Ncell, Ndim)
embeds += np.repeat(subj_bias, Ncell_subj, axis=0)
embeds += np.repeat(group_bias, Ncell_group, axis=0)
subjects = np.repeat(np.arange(Nsubj), Ncell_subj)
groups = np.repeat(np.arange(Ngroup), Ncell_group, axis=0)

#Generate data frame
data_embed = embed_data(subjects, groups, embeds)

#Test effect with LMM
data_embed.LMM()

