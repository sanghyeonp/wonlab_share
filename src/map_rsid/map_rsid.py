import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

from util import *
from packages import *


def parse_args():
    parser = argparse.ArgumentParser(description=":: Map rsID using chromosome and base position ::")
    parser.add_argument('--file', required=True,
                        help='Path to the input file.')
    