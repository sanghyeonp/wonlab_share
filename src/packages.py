import pip

def import_or_install(package):
    try:
        __import__(package)
    except ImportError:
        pip.main(['install', package])       

packages_list = ['os',
                'pandas', 'csv', 'numpy',
                'argparse', 'code', 'tqdm',
                'subprocess', 'multiprocessing', 
                'statsmodels']

[import_or_install(p) for p in packages_list]    

import pandas as pd
pd.set_option('mode.chained_assignment',  None)
import os
from datetime import datetime
import csv
import argparse
from tqdm import tqdm
tqdm.pandas(leave=False)
import io
import numpy as np
import subprocess
from multiprocessing import Pool, Manager

import code
# code.interact(local=dict(globals(), **locals()))