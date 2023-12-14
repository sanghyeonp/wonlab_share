import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

from util import *
from packages import *
from shared_data import ANNOVAR_software, ANNOVAR_humandb


def generate_cmd(annov_input, outd):
    prefix = os.path.split(annov_input)[-1].replace(".annovin", "")
    cmd = "{} {} {} -buildver hg19 -protocol refGene,avsnp150 -operation g,f -remove -nastring . -out {}".format(
        ANNOVAR_software, os.path.join(outd, annov_input), ANNOVAR_humandb, os.path.join(outd, prefix))
    
    # print(cmd)

    return cmd


def annovar(annov_input1, annov_input2, outd, log_list=[]):
    log_list = logger(log_list, log="## Running ANNOVAR with {}...".format(annov_input1))
    cmd1 = generate_cmd(annov_input1, outd)
    run_bash(cmd1)

    log_list = logger(log_list, log="## Running ANNOVAR with {}...".format(annov_input2))
    cmd2 = generate_cmd(annov_input2, outd)
    run_bash(cmd2)
    
    return log_list
