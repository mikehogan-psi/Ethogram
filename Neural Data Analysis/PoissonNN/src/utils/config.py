# -*- coding: utf-8 -*-
"""
Loads the configuration parameters for:
    - model structure 
    - training parameters
    - data locations
"""

import yaml

class Config:
    def __init__(self, config_path: str):
        with open(config_path, "r") as f:
            cfg = yaml.safe_load(f)
        
        # Split into sections
        self.model = cfg.get("model", {})
        self.training = cfg.get("training", {})
        self.data = cfg.get("data", {})

    def __repr__(self):
        repr_str = []
        repr_str.append(f"Config(model={self.model}")
        repr_str.append(f", training={self.training})")
        repr_str.append(f", data={self.data})")
        repr_str = "".join(repr_str)
        return repr_str

