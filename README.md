# Toward minimal SNP sets for record-matching with CODIS STR profiles


This repository contains all the code necessary to replicate all analyses described in paper [LINK]. Below, we describe the key dependencies, raw data and relevant references.

## Dependencies
All software required for this project are provided in this repository under [software](https://github.com/tamigj/codis_panel/tree/main/software). These include:
- BEAGLE 5.4, specifically [beagle.22Jul22.46e.jar](https://faculty.washington.edu/browning/beagle/)
- Java version 8 - required by BEAGLE. 
- [VCFtools](https://vcftools.github.io/)
- [PLINK 1.9](https://www.cog-genomics.org/plink/)
- RecordMatching R [package](https://github.com/jk2236/RecordMatching)

## Dataset

- A phased reference SNP-STR haplotype panel from [Saini et al.](https://www.nature.com/articles/s41467-018-06694-0) from the 1000 Genomes Project phase 3 containing 2504 individuals can be downloaded from [here](https://gymreklab.com/2018/03/05/snpstr_imputation.html). For our work, we download and process this data to keep 1-Mb SNP windows extending 500-Kb in each direction from the midpoint of each of 18 CODIS STRs. This resulting data was used as a starting point for all of our analyses. In this repository, we provide both the [code](https://github.com/tamigj/codis_panel/blob/main/analyses/data_processing/0_process_saini_files.sh) to download data from Saini et al. and process it, as well as the resulting processed [dataset](https://github.com/tamigj/codis_panel/tree/main/data/raw) that was used for all of our analyses. 
- HapMap GrCh37 genetic maps in PLINK format can be downloaded from the [BEAGLE page](https://example.com/beagle-page).

## Code description

- Below is a schematic representation describing the overall structure of our code base. Each grey box corresponds to a directory under [/analyses](https://github.com/tamigj/codis_panel/tree/main/analyses).
- More detailed information about scripts that make up each analysis is included in this [README](https://github.com/tamigj/codis_panel/blob/main/README_CODE_DESCRIPTION).

![CODIS Codebase Schematic](/analyses/CODIS_codebase_schematic_v2.png)

## References

(our paper, when published)

Kim, J., & Rosenberg, N. A. (2023). Record-matching of STR profiles with fragmentary genomic SNP data. European Journal of Human Genetics, 31(11), 1283-1290.

Kim J, Edge MD, Algee-Hewitt BFB, Li JZ, Rosenberg NA (2018). Statistical detection of relatives typed with disjoint forensic and biomedical loci. Cell, 175(3):848-858.e6. 10.1016/j.cell.2018.09.008.

Edge MD, Algee-Hewitt BFB, Pemberton TJ, Li JA, Rosenberg NA (2017). Linkage disequilibrium matches forensic genetic records to disjoint genomic marker sets. PNAS, 114(22):5671-5676. 10.1073/pnas.1619944114.

## Note: data in related manuscripts

The phased haplotypes in this repository and accompanying publication trace to Saini et al. (2018)[https://www.nature.com/articles/s41467-018-06694-0]; Kim & Rosenberg (2023)[https://www.nature.com/articles/s41431-023-01430-9] and Lappo & Rosenberg (2024)[https://www.sciencedirect.com/science/article/pii/S258900422400052X?via%3Dihub] have also used phased haplotypes based on Saini et al (2018). Kim & Rosenberg (2023) used both phased data from Saini et al. (2018) and unphased data that they produced from the phased data, and they posted the unphased data in a github repository. Lappo & Rosenberg (2024) used the phased data from Kim & Rosenberg (2023) in their analysis and included the unphased data from Kim & Rosenberg (2023) for download in their electronic supplementary material. For clarity here, the code in this repository obtains phased haplotypes directly from Saini et al. (2018).








