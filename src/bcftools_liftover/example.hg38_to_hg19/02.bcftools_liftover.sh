#!/bin/bash
#SBATCH -J 02.bcftools_liftover
#SBATCH -p cpu
#SBATCH -o ./%x_%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=14-00:00:00
#SBATCH --mail-user=sh.austin.park@gmail.com 
#SBATCH --mail-type=END,FAIL 

module purge
set -euo pipefail

CONDA_PATH=/data1/software/anaconda3 
source $CONDA_PATH/bin/activate /data1/software/anaconda3/envs/R_SP_4.4 

##################################################
# Required user defined variables
##################################################


# Input GWAS summary statistics file
INPUT="GCST90624699.ID.tsv"

# Output prefix
OUT_PREFIX="BMI_EUR.Abner2024"

# Trait name
TRAIT="BMI"

# Column names in the input file
COL_SNP="ID"
COL_CHROM="chromosome"
COL_POS="base_pair_location"
COL_A1="effect_allele"
COL_A2="other_allele"
COL_FRQ="effect_allele_frequency"
COL_BETA="beta"
COL_SE="standard_error"
COL_P="p_value"

# Keep only autosomes 1 to 22
KEEP_AUTOSOMES="true"

##################################################
# Checks
##################################################
source /data1/sanghyeon/wonlab_contribute/combined/src/bcftools_liftover/global_variable.sh

require_file() {
    local f="$1"
    [[ -f "$f" ]] || { echo "ERROR: file not found: $f" >&2; exit 1; }
}

for f in "$INPUT" "$HG38_REF" "$HG19_REF" "$CHAIN_hg38ToHg19"; do
    require_file "$f"
done

require_file "$bcftools_exe"

##################################################
# Output files
##################################################

HG38_BCF="${OUT_PREFIX}.hg38.bcf"
HG19_BCF="${OUT_PREFIX}.hg19.bcf"
HG19_AUTO_BCF="${OUT_PREFIX}.hg19.auto.bcf"
HG19_CHR_BCF="${OUT_PREFIX}.hg19.auto.chr.bcf"
HG19_NOJUMP_BCF="${OUT_PREFIX}.hg19.auto.chr.nojump.bcf"
FINAL_TSV="${OUT_PREFIX}.hg19.tsv"
COLHDR="${OUT_PREFIX}.colheaders.tsv"
REJECT_BCF="${OUT_PREFIX}.lift.reject.bcf"

SWAP_FLIP_TSV="${OUT_PREFIX}.allele_swap.strand_flip.tsv"
REJECT_TSV="${OUT_PREFIX}.lift.reject.tsv"
DROP_NONAUTO_TSV="${OUT_PREFIX}.drop.nonautosome.tsv"
DROP_CHRJUMP_TSV="${OUT_PREFIX}.drop.chromjump.tsv"
DROP_AMBIGUOUS_TSV="${OUT_PREFIX}.drop.ambiguous.tsv"

LOGFILE="${OUT_PREFIX}.liftover.$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -i "$LOGFILE") 2>&1

echo "START: $(date)"
echo "INPUT: $INPUT"
echo "OUT PREFIX: $OUT_PREFIX"
echo "TRAIT: $TRAIT"
echo "LOG FILE: $LOGFILE"


##################################################
# Output header labels for dropped SNP tables
##################################################

HDR_CHR="$COL_CHROM"
HDR_POS="$COL_POS"
HDR_SNP="$COL_SNP"
HDR_A1="$COL_A1"
HDR_A2="$COL_A2"

##################################################
# Build bcftools +munge column mapping
##################################################

{
    printf "%s\tCHR\n" "$COL_CHROM"
    printf "%s\tBP\n" "$COL_POS"
    printf "%s\tSNP\n" "$COL_SNP"
    printf "%s\tA1\n" "$COL_A1"
    printf "%s\tA2\n" "$COL_A2"
    printf "%s\tFRQ\n" "$COL_FRQ"
    printf "%s\tBETA\n" "$COL_BETA"
    printf "%s\tSE\n" "$COL_SE"
    printf "%s\tP\n" "$COL_P"
} > "$COLHDR"

echo
echo "Generated colheaders file:"
cat "$COLHDR"

##################################################
# Step 1. Convert input summary statistics to hg38 BCF
##################################################

"${bcftools_exe}" +munge --no-version -Ou \
    -C "$COLHDR" \
    -f "$HG38_REF" \
    -s "$TRAIT" \
    "$INPUT" \
    | "${bcftools_exe}" sort -Ob -o "$HG38_BCF" --write-index

echo
echo "Step 1 done:"
echo "  $HG38_BCF"
echo "  ${HG38_BCF}.csi"

[[ -s "$HG38_BCF" ]] || { echo "ERROR: Step 1 failed, $HG38_BCF not created properly"; exit 1; }

##################################################
# Step 2. Liftover hg38 to hg19 and normalize
##################################################

"${bcftools_exe}" +liftover --no-version -Ou "$HG38_BCF" -- \
    -s "$HG38_REF" \
    -f "$HG19_REF" \
    -c "$CHAIN_hg38ToHg19" \
    --write-src \
    --swap-tag SWAP \
    --flip-tag FLIP \
    --reject "$REJECT_BCF" \
    --reject-type b \
    | "${bcftools_exe}" norm -Ou -f "$HG19_REF" -m -any \
    | "${bcftools_exe}" sort -Ob -o "$HG19_BCF" --write-index

echo
echo "Step 2 done:"
echo "  $REJECT_BCF"
echo "  $HG19_BCF"
echo "  ${HG19_BCF}.csi"

[[ -s "$HG19_BCF" ]] || { echo "ERROR: Step 2 failed, $HG19_BCF not created properly"; exit 1; }

##################################################
# Step 2a. Save liftover-rejected variants
##################################################

{
    printf "%s\t%s\t%s\t%s\t%s\tINFO\n" \
        "$HDR_CHR" "$HDR_POS" "$HDR_SNP" "$HDR_A2" "$HDR_A1"
    "${bcftools_exe}" query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO\n' "$REJECT_BCF"
} > "$REJECT_TSV"

echo
echo "Saved liftover-rejected variants:"
echo "  $REJECT_TSV"

##################################################
# Step 3. Keep only autosomes if requested
##################################################

if [[ "$KEEP_AUTOSOMES" == "true" ]]; then
    {
        printf "%s\t%s\t%s\t%s\t%s\n" \
            "$HDR_CHR" "$HDR_POS" "$HDR_SNP" "$HDR_A2" "$HDR_A1"
        "${bcftools_exe}" query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' "$HG19_BCF" \
            | awk 'BEGIN{OFS="\t"}
                   $1 !~ /^(1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22)$/ {print}'
    } > "$DROP_NONAUTO_TSV"

    echo
    echo "Keeping only chromosomes 1 to 22"
    "${bcftools_exe}" view -r "$(seq -s, 1 22)" -Ob -o "$HG19_AUTO_BCF" "$HG19_BCF"
    "${bcftools_exe}" index -f "$HG19_AUTO_BCF"
else
    cp "$HG19_BCF" "$HG19_AUTO_BCF"
    cp "${HG19_BCF}.csi" "${HG19_AUTO_BCF}.csi"

    printf "%s\t%s\t%s\t%s\t%s\n" \
        "$HDR_CHR" "$HDR_POS" "$HDR_SNP" "$HDR_A2" "$HDR_A1" > "$DROP_NONAUTO_TSV"
fi

echo
echo "Step 3 done:"
echo "  $HG19_AUTO_BCF"
echo "  ${HG19_AUTO_BCF}.csi"
echo "Saved non-autosomal variants:"
echo "  $DROP_NONAUTO_TSV"

##################################################
# Step 4. Rename chromosomes if map exists
##################################################

if [[ -f "$CHRMAP" ]]; then
    echo
    echo "Renaming chromosomes using: $CHRMAP"
    "${bcftools_exe}" annotate --rename-chrs "$CHRMAP" -Ob -o "$HG19_CHR_BCF" "$HG19_AUTO_BCF"
    "${bcftools_exe}" index -f "$HG19_CHR_BCF"
else
    echo
    echo "Chromosome rename map not found. Skipping rename step."
    cp "$HG19_AUTO_BCF" "$HG19_CHR_BCF"
    cp "${HG19_AUTO_BCF}.csi" "${HG19_CHR_BCF}.csi"
fi

##################################################
# Step 4a. Save chromosome-jump variants
##################################################

{
    printf "%s\t%s\t%s\t%s\t%s\tSRC_%s\n" \
        "$HDR_CHR" "$HDR_POS" "$HDR_SNP" "$HDR_A2" "$HDR_A1" "$HDR_CHR"
    "${bcftools_exe}" query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/SRC_CHROM\n' "$HG19_CHR_BCF" \
        | awk 'BEGIN{OFS="\t"} $1 != $6 {print}'
} > "$DROP_CHRJUMP_TSV"

##################################################
# Step 4b. Remove chromosome-jump variants
##################################################

"${bcftools_exe}" view \
    -i 'INFO/SRC_CHROM==CHROM' \
    -Ob -o "$HG19_NOJUMP_BCF" "$HG19_CHR_BCF"

"${bcftools_exe}" index -f "$HG19_NOJUMP_BCF"

echo
echo "Step 4 done:"
echo "  $HG19_CHR_BCF"
echo "  ${HG19_CHR_BCF}.csi"
echo "  $HG19_NOJUMP_BCF"
echo "  ${HG19_NOJUMP_BCF}.csi"
echo "Saved chromosome-jump variants:"
echo "  $DROP_CHRJUMP_TSV"


##################################################
# Step 4c. Save allele-swapped variants
##################################################

{
    printf "%s\t%s\t%s\t%s\t%s\tSRC_%s\tSRC_%s\tSRC_REF_ALT%s\tSWAP\tFLIP\n" \
        "$HDR_CHR" "$HDR_POS" "$HDR_SNP" "$HDR_A2" "$HDR_A1" \
        "$HDR_CHR" "$HDR_POS"

    "${bcftools_exe}" query \
        -i 'INFO/SWAP=1 || INFO/FLIP=1' \
        -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/SRC_CHROM\t%INFO/SRC_POS\t%INFO/SRC_REF_ALT\t%INFO/SWAP\t%INFO/FLIP\n' \
        "$HG19_NOJUMP_BCF"
} > "$SWAP_FLIP_TSV"

echo
echo "Saved allele-swapped and strand-flipped variants:"
echo "  $SWAP_FLIP_TSV"

##################################################
# Step 5. Export final table
##################################################
echo
echo "Exporting final TSV"
FINAL_BODY_TMP="tmp.tsv"

{
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$COL_CHROM" "$COL_POS" "$COL_SNP" "$COL_A1" "$COL_A2" \
        "$COL_BETA" "$COL_SE" "$COL_P" "$COL_FRQ"

    "${bcftools_exe}" query -f '%CHROM\t%POS\t%ID\t%ALT\t%REF\t[%ES]\t[%SE]\t[%LP]\t[%AF]\n' "$HG19_NOJUMP_BCF" \
        | awk 'BEGIN{OFS="\t"}
               {
                   p = 10^(-$8)
                   print $1, $2, $3, $4, $5, $6, $7, p, $9
               }'
} | awk 'BEGIN{OFS="\t"} NR==1{print; next} {sub(/^chr/,"",$1); print}' > "$FINAL_BODY_TMP"

# Step 5b. Save variants with missing beta
{
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$COL_CHROM" "$COL_POS" "$COL_SNP" "$COL_A1" "$COL_A2" \
        "$COL_BETA" "$COL_SE" "$COL_P" "$COL_FRQ"

    awk -F'\t' 'BEGIN{OFS="\t"} $6=="." || $6==""' "$FINAL_BODY_TMP"
} > "$DROP_AMBIGUOUS_TSV"

# Step 5c. Save final table after removing variants with missing beta
{
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$COL_CHROM" "$COL_POS" "$COL_SNP" "$COL_A1" "$COL_A2" \
        "$COL_BETA" "$COL_SE" "$COL_P" "$COL_FRQ"

    awk -F'\t' 'BEGIN{OFS="\t"} !($6=="." || $6=="")' "$FINAL_BODY_TMP"
} > "$FINAL_TSV"

echo
echo "Step 5 done:"
echo "  $FINAL_TSV"
echo "Saved variants with missing beta:"
echo "  $DROP_AMBIGUOUS_TSV"

rm -f "$FINAL_BODY_TMP"


##################################################
# Cleanup intermediate BCF/CSI files
##################################################

echo
echo "Cleaning up intermediate BCF and CSI files..."

rm -f \
    "$HG38_BCF" "${HG38_BCF}.csi" \
    "$HG19_BCF" "${HG19_BCF}.csi" \
    "$HG19_AUTO_BCF" "${HG19_AUTO_BCF}.csi" \
    "$HG19_CHR_BCF" "${HG19_CHR_BCF}.csi" \
    "$HG19_NOJUMP_BCF" "${HG19_NOJUMP_BCF}.csi" \
    "$REJECT_BCF"

echo "Cleanup done."

##################################################
# Done
##################################################

echo
echo "Finished: $(date)"
echo "Outputs:"
echo "  $FINAL_TSV"
echo "  $REJECT_TSV"
echo "  $SWAP_FLIP_TSV"
echo "  $DROP_NONAUTO_TSV"
echo "  $DROP_CHRJUMP_TSV"
echo "  $DROP_AMBIGUOUS_TSV"
