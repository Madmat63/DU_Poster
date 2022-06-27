#!/bin/bash

#SBATCH --mem=32G
#SBATCH --cpus-per-task=12
#SBATCH --array=0-2            # limit to 3 samples. Use --array=0-46 for all 47 samples.

# le script va s'arrêter
# - à la première erreur
# - si une variable n'est pas définie
# - si une erreur est recontrée dans un pipe
set -euo pipefail


# chargement des modules nécessaires
module load fastqc/0.11.9
module load minimap2/2.18
module load samtools/1.9
module load bedtools/2.30.0

# répertoire de base
base_dir="/shared/projects/form_2021_29/mchicard/Nanopore"
# répertoire contenant les fichiers du génome de référence
# (séquence et annotations)
genome_dir="${base_dir}/genome"
# répertoire contenant les fichiers .fastq
fastq_dir="${base_dir}/raw"
# liste de tous les fichiers .fastq
fastq_files=(${fastq_dir}/*fastq)
# extraction de l'identifiant de l'échantillon
# à partir du nom de fichier : /shared/ifbstor1/projects/form_2021_29/mchicard/Nanopore/raw/barcode01.fastq
# on extrait : barcode01
sample=$(basename -s .fastq "${fastq_files[$SLURM_ARRAY_TASK_ID]}}")



echo "=============================================================="
echo "Contrôle qualité - échantillon ${sample}"
echo "=============================================================="
mkdir -p reads_qc
#srun fastqc "${fastq_dir}/${sample}.fastq" --outdir reads_qc

echo "=============================================================="
echo "Alignement des reads sur le génome de référence - échantillon ${sample}"
echo "=============================================================="
mkdir -p map
srun minimap2  "${genome_dir}/pRMCE_CD4_CTG" "${fastq_dir}/${sample}.fastq" > "map/minimap2-${sample}.sam" -ax map-ont

echo "=============================================================="
echo "Conversion en binaire, tri et indexation des reads alignés - échantillon ${sample}"
echo "=============================================================="
srun samtools view -b "map/minimap2-${sample}.sam" -o "map/minimap2-${sample}.bam"
srun samtools sort "map/minimap2-${sample}.bam" -o "minimap2-${sample}.sorted.bam"
srun samtools index "map/minimap2-${sample}.sorted.bam"

echo "=============================================================="
echo "Couverture - échantillon ${sample}"
echo "=============================================================="
mkdir -p coverage
srun bedtools genomecov -ibam "map/bowtie-${sample}.sorted.bam" "${annotations}" > "count/count-${sample}.bedgraph"

echo "=============================================================="
echo "Nettoyage des fichiers inutiles - échantillon ${sample}"
echo "=============================================================="
rm -f "map/minimap2-${sample}.sam" "map/minimap2-${sample}.bam"
