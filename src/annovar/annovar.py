import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

from util import *
from packages import *
from shared_data import ANNOVAR_software, ANNOVAR_humandb


def generate_cmd(annov_input):
    prefix = os.path.split(annov_input)[-1].replace(".annovin", "")
    cmd = "{} {} {} -buildver hg19 -protocol refGene,avsnp150 -operation g,f -remove -nastring . -out {}.1".format(
        ANNOVAR_software, annov_input, ANNOVAR_humandb, prefix)
    
    # print(cmd)

    return cmd


def annovar(annov_input1, annov_input2):
    print("## Running ANNOVAR with {}...".format(annov_input1))
    cmd1 = generate_cmd(annov_input1)
    run_bash(cmd1)

    print("## Running ANNOVAR with {}...".format(annov_input2))
    cmd2 = generate_cmd(annov_input2)
    run_bash(cmd2)
