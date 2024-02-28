#!/bin/bash -l
#SBATCH --job-name=soup
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=2-0
#SBATCH --mem=256G
#SBATCH --array=3
#SBATCH --mail-type=ALL
#SBATCH --output=logs/soup_%A_%a.err

set -euo pipefail

# load modules
module purge
module load container/singularity/3.10.4/

echo "Starting at `date`"

echo "Running on hosts: $SLURM_JOB_NODELIST"
echo "Running on $SLURM_JOB_NUM_NODES nodes."
echo "Running $SLURM_NTASKS tasks."
echo "Account: $SLURM_JOB_ACCOUNT"
echo "Job ID: $SLURM_JOB_ID"
echo "Job name: $SLURM_JOB_NAME"
echo "Node running script: $SLURMD_NODENAME"
echo "Submit host: $SLURM_SUBMIT_HOST"

# root directory
cd /exports/archive/hg-funcgenom-research/

# specify path to metadata file
PROJECT_ROOT="mdmanurung/nih-malaria/"
CRPATH="${PROJECT_ROOT}/data/metadata/path_files.tsv"

# get variables from table
ROWID=$SLURM_ARRAY_TASK_ID
read BATCH PATH_RAW PATH_BAM N <<<$(awk -F, "NR == $ROWID" ${CRPATH} | cut -f 1,2,3,4)
PATH_BAM="${PATH_BAM}/sample_alignments.bam"
BARCODES="${PROJECT_ROOT}results/barcodes_filtered/barcodes_${BATCH}.tsv.gz"
SOUPORCELL_OUTDIR="${PROJECT_ROOT}results/souporcell/soup_${BATCH}"
VCF="${PROJECT_ROOT}data/snp_vcf/GRCh38_1000G_MAF0.01_GeneFiltered_ChrEncoding.vcf"
FASTA="${PROJECT_ROOT}/data/refgenome/refdata-gex-GRCh38-2020-A/fasta/genome.fa"

if [ ! -d "$SOUPORCELL_OUTDIR" ]; then

    mkdir $SOUPORCELL_OUTDIR

    echo "Running souporcell for capture ${BATCH}"
    singularity exec \
        --bind /exports/archive/hg-funcgenom-research/ \
        mdmanurung/nih-malaria/containers/Demuxafy.sif souporcell_pipeline.py -i $PATH_BAM \
        -b $BARCODES \
        -f $FASTA \
        -t 12 \
        -o $SOUPORCELL_OUTDIR \
        -k $N \
        --common_variants $VCF
fi

echo "Program finished with exit code $? at: `date`"
