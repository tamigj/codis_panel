#/!bin/bash

DIR_DATA_NRUN=$1
DIR_OUTPUT_NRUN=$2
str=$3

source /scratch/groups/noahr/tami/codis_panel/config.sh

mkdir -p $DIR_OUTPUT_NRUN

chr=$(grep -nF ${str} \
    $DIR_DATA_NRUN/reference/${str}.recode.vcf | cut -f1 | awk -F ":" '{print $2}')

rm -f $DIR_OUTPUT_NRUN/${str}*

ml java

java -Xmx1g -jar $DIR_SOFTWARE/beagle.22Jul22.46e.jar \
    ref=$DIR_DATA_NRUN/reference/${str}.recode.vcf \
    gt=$DIR_DATA_NRUN/SNP/${str}_SNP.recode.vcf \
    out=$DIR_OUTPUT_NRUN/${str} \
    map=$DIR_PLINK_MAPS/plink.chr${chr}.GRCh37.map \
    ap=true \
    gp=true \
    impute=true \
    nthreads=1

gunzip $DIR_OUTPUT_NRUN/${str}.vcf.gz

"$DIR_VCFTOOLS/vcftools" \
    --vcf $DIR_OUTPUT_NRUN/${str}.vcf \
    --snp ${str} \
    --extract-FORMAT-info GT \
    --out $DIR_OUTPUT_NRUN/${str} &&

"$DIR_VCFTOOLS/vcftools" \
    --vcf $DIR_OUTPUT_NRUN/${str}.vcf \
    --snp ${str} \
    --extract-FORMAT-info GP \
    --out $DIR_OUTPUT_NRUN/${str}
