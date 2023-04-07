import pip

def import_or_install(package):
    try:
        __import__(package)
    except ImportError:
        pip.main(['install', package])       

packages_list = ['pandas', 'os', 'subprocess', 'multiprocessing', 'csv', 'argparse', 'code', 'tqdm', 'statsmodels']

[import_or_install(p) for p in packages_list]    

import pandas as pd
pd.set_option('mode.chained_assignment',  None)
import os
import subprocess
from multiprocessing import Pool, Manager
from datetime import datetime
import csv
import argparse
from tqdm import tqdm
from statsmodels.sandbox.stats.multicomp import multipletests
import io

import code
# code.interact(local=dict(globals(), **locals()))