#!/PHShome/ik936/anaconda3/envs/gpu/bin/python

## THIS LIBRARY MUST BE IMPORTED FIRST 
import scarches as sca

## And then again? 
import scarches as sca
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)
import scanpy as sc
import torch
from scarches.dataset.trvae.data_handling import remove_sparsity
import matplotlib.pyplot as plt
import numpy as np
import gdown
import time
import json
