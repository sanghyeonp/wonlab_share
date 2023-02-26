from src.packages import *
from src.util import run_bash

LIFTOVER_DIR = "/data1/software/liftOver" 

liftover_software = os.path.join(LIFTOVER_DIR, "liftOver")

chain_mapping = {(37, 38): os.path.join(LIFTOVER_DIR, "chain", "GRCh37_to_GRCh38.chain.gz"),
                (38, 37): os.path.join(LIFTOVER_DIR, "chain", "GRCh38_to_GRCh37.chain.gz"),
                (18, 19): os.path.join(LIFTOVER_DIR,"chain", "hg18ToHg19.over.chain.gz")
                }


def run_liftover(input_bed, build_from: int, build_to: int, outd):
    global liftover_software, chain_mapping
    _, filename = os.path.split(input_bed)
    if filename[-4:] == '.bed':
        filename = filename[:-4]
    bash_cmd = "{} {} {} {} {}".format(liftover_software, 
                                    input_bed, 
                                    chain_mapping[(build_from, build_to)],
                                    os.path.join(outd, filename+".lifted"),
                                    os.path.join(outd, filename+".unlifted")
                                    )
    # print(bash_cmd)
    stdout = run_bash(bash_cmd)
