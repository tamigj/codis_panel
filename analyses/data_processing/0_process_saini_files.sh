#!/bin/bash

#SBATCH --job-name=0_process_saini
#SBATCH --time=08:00:00
#SBATCH --mem=30G
#SBATCH --output=logs/0_process_saini.log
#SBATCH --error=logs/0_process_saini.err
#SBATCH -p hns,owners,normal
#SBATCH -c 1

source /scratch/groups/noahr/tami/codis_panel/config.sh

cd ${DIR_DATA_RAW}

# Identify chromosomes that have CODIS STRs
CHROMOSOMES=($(tail -n +2 marker_positions.txt | cut -d' ' -f2 | sort -n | uniq))

for CHR in "${CHROMOSOMES[@]}"; do

    # Download and unzip VCF file
    if [[ ! -f "1kg.snp.str.chr${CHR}.vcf" ]]; then
      wget "https://s3.amazonaws.com/snp-str-imputation/1000genomes/1kg.snp.str.chr${CHR}.vcf.gz"
      gunzip "1kg.snp.str.chr${CHR}.vcf.gz"
    fi

    # Process each STR on this chromosome
    while read -r str chr pos; do

      if [[ "$chr" == "$CHR" ]]; then

        start=$((pos - 500000))
        end=$((pos + 500000))

        # Extract 1MB region around STR using vcftools
        $DIR_VCFTOOLS/vcftools --vcf "1kg.snp.str.chr${CHR}.vcf" \
                               --chr ${CHR} \
                               --from-bp ${start} --to-bp ${end} \
                               --recode \
                               --out "${str}_halfwindow500000WithSTR"

        # Rename output file to drop recode
        mv "${str}_halfwindow500000WithSTR.recode.vcf" "${str}_halfwindow500000WithSTR.vcf"
        
      fi
  done < <(tail -n +2 marker_positions.txt)
done

rm *log
rm *chr*vcf
