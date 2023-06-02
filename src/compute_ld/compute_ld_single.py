import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

from util import *
from packages import *

from shared_data import PLINK_bfile, PLNIK1_9


def parse_args():
    parser = argparse.ArgumentParser(description=":: Compute LD between two SNPs ::")
    
    # Required arguments
    parser.add_argument('--snp1', required=True,
                        help='SNP in either rsID or chr:pos format.')
    parser.add_argument('--snp2', required=True,
                        help='SNP in either rsID or chr:pos format.')
    parser.add_argument('--variant-identifier', dest="variant_identifier", required=True,
                        help='Type of variant identifier. Choices = ["rsid", "chr:pos"]')
    parser.add_argument('--reference', required=True,
                        help='Type of reference. Choices = ["1kg", "ukb"]')
    parser.add_argument('--ancestry', required=True,
                        help='Ancestry of the input file. \
                            Choices from 1kg = ["European", "East-asian"] \
                            Choices from ukb = ["European"]')
    args = parser.parse_args()
    return args


def bfile_selection(reference, variant_identifier, ancestry):
    """
    reference: choices from ["1kg", "ukb"]
    variant_identifier: choices from ["rsid", "chr:pos"]
    ancestry: choices from 1kg ["European", "East-asian"]
            choices from ukb ["European"]
    """
    global PLINK_bfile
    
    assert reference in ["1kg", "ukb"], "Allowed choices for reference are ['1kg', 'ukb']"
    assert variant_identifier in ["rsid", "chr:pos"], "Allowed choices for variant_identifier are ['rsid', 'chr:pos']"
    if reference == '1kg':
        assert ancestry in ["European", "East-asian"], "Allowed choices for ancestry are ['European', 'East-asian']"
    elif reference == 'ukb':
        assert ancestry in ["European"], "Allowed choices for ancestry are ['European']"
    
    return PLINK_bfile[reference][variant_identifier][ancestry]


def extract_R2(stdout):
    for row in stdout:
        if 'R-sq' in row:
            _, r2 = row.split(sep="R-sq = ")
            value, _ = r2.split(sep="D' = ")
            return float(value.strip())
    return -9


def compute_ld_single(snp1, snp2, bfile, shared_list=[]):
    """
    snp1 and snp2 in format CHR:POS:Allele1:Allele2 where Allele 1 and 2 are in alphabetical order.
    """

    if snp1 == snp2:
        shared_list.append([snp1, snp2, 1])
        return 1

    cmd = "{} --bfile {} --ld {} {}".format(PLNIK1_9, bfile, snp1, snp2)

    stdout = run_bash(cmd)
    r2 = extract_R2(stdout)
    shared_list.append([snp1, snp2, r2])
    return r2


if __name__ == "__main__":
    args = parse_args()

    bfile = bfile_selection(args.reference, args.variant_identifier, args.ancestry)

    r2 = compute_ld_single(snp1=args.snp1, 
                    snp2=args.snp2, 
                    bfile=bfile
                    )
    
    print(r2)
    run_bash("rm plink.log plink.nosex")

    """
    UKB reference로 사용 권장 하지 않기.
    """