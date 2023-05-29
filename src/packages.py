import pip

def import_or_install(package):
    try:
        __import__(package)
    except ImportError:
        pip.main(['install', package])       

packages_list = ['pandas', 'os', 'numpy',
                'csv', 'argparse', 'code', 'tqdm',
                'subprocess', 'multiprocessing']

[import_or_install(p) for p in packages_list]    

import pandas as pd
pd.set_option('mode.chained_assignment',  None)
import os
from datetime import datetime
import csv
import argparse
from tqdm import tqdm
import io
import numpy as np
import subprocess
from multiprocessing import Pool, Manager

import code
# code.interact(local=dict(globals(), **locals()))