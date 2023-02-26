import pip

def import_or_install(package):
    try:
        __import__(package)
    except ImportError:
        pip.main(['install', package])       

packages_list = ['pandas', 'os', 'subprocess', 'csv', 'argparse', 'code']

[import_or_install(p) for p in packages_list]    

import pandas as pd
pd.set_option('mode.chained_assignment',  None)
import os
import subprocess
import csv
import argparse
import code