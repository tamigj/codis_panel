# Minimal SNP sets for record-matching with CODIS STR profiles

This repository contains all the code necessary to replicate all analyses described in paper [LINK]. Below, we describe the key dependencies, raw data and relevant references.

## Dependencies

- BEAGLE 5.4, specifically [beagle.22Jul22.46e.jar](https://faculty.washington.edu/browning/beagle/)
- Java version 8 - required by BEAGLE. 
- [VCFtools](https://vcftools.github.io/)
- [PLINK 1.9](https://www.cog-genomics.org/plink/)
- RecordMatchin R [package](https://github.com/jk2236/RecordMatching)

## Dataset

- A phased reference SNP-STR haplotype panel from [Saini et al.](https://www.nature.com/articles/s41467-018-06694-0) from the 1000 Genomes Project phase 3 containing 2504 individuals can be downloaded from [here](https://gymreklab.com/2018/03/05/snpstr_imputation.html). Processed data containing 1-Mb SNP windows extending 500-Kb in each direction from each CODIS STR midpoint can be found [here](https://github.com/jk2236/RM_WGS/tree/main/data/1KGP).
- HapMap GrCh37 genetic maps in PLINK format can be downloaded from the [BEAGLE page](https://example.com/beagle-page).

## Code description

- The full description of the code, inluding the input, output, and description of each script can be found [here](https://github.com/tamigj/codis_panel/blob/main/README_CODE_DESCRIPTION)

## Reference

(paper)

Kim J, Rosenberg NA (2022). Record-matching of STR profiles with fragmentary genomic SNP data. bioRxiv, 2022.09.01.505545. 10.1101/2022.09.01.505545.

Kim J, Edge MD, Algee-Hewitt BFB, Li JZ, Rosenberg NA (2018). Statistical detection of relatives typed with disjoint forensic and biomedical loci. Cell, 175(3):848-858.e6. 10.1016/j.cell.2018.09.008.

Edge MD, Algee-Hewitt BFB, Pemberton TJ, Li JA, Rosenberg NA (2017). Linkage disequilibrium matches forensic genetic records to disjoint genomic marker sets. PNAS, 114(22):5671-5676. 10.1073/pnas.1619944114.
