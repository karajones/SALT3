## Scripts used to process data for *Environmental DNA biodiversity assessment of mesophotic reefs in the northern Gulf of America*.

> [!CAUTION]
> This repository is currently under development.

These scripts are intended to provide additional details on the methods used for bioinformatic analyses outlined in the manuscript and to provide additional documentation to improve reproducibility of the results.

#### Data accessibility 
- Sequencing files (fastqs) are deposited in the NCBI Sequence Read Archive under BioProject [PRJNA1198004](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1198004)
- The final processed data are deposited in the USGS ScienceBase repository: [Detections of marine fauna using environmental DNA metabarcoding of seawater samples from the northern Gulf of America](https://doi.org/10.5066/P1Z4JWKF)

## Table of contents

#### Sequence processing scripts
- `01_fastq_processing_MDBC.sh`: Shell script used to prepare raw sequencing (fastq) files for import into the `dada2` denoising pipeline.
- `02_dada2_MDBC.R`: R script pipeline used to infer amplicon sequence variants (ASVs) from clean fastq files.
- `03_database_prep_MDBC.R`: R script used to create the 12S MiFish reference sequence database with data downloaded from NCBI.
- `04_parsing_BLAST_results_MDBC.R`: R script used to assign taxonomy to ASVs.

#### Plot scripts
- `12_plots_fish_depth.R`: R script used to create a plot of fish detections by depth.
- `13_plots_coral_dotplot.R`: R script used to create a plot comparing coral detections with eDNA and video.
- `14_heatmap.R`: R script used to create a plot of read abundance and number of ASVs found for each taxonomic class by depth class.

#### Other files
- `additional_data`: Directory of files used to help create some of the plots.
- `DISCLAIMER_provisional.md`: Required disclaimer for scripts.
- `LICENSE.md`: License for use.

## Disclaimer

This software is preliminary or provisional and is subject to revision. It is being provided to meet the need for timely best science. The software has not received final approval by the U.S. Geological Survey (USGS). No warranty, expressed or implied, is made by the USGS or the U.S. Government as to the functionality of the software and related material nor shall the fact of release constitute any such warranty. The software is provided on the condition that neither the USGS nor the U.S. Government shall be held liable for any damages resulting from the authorized or unauthorized use of the software.

Any use of trade, firm, or product names is for descriptive purposes only and does not imply endorsement by the U.S. Government. Although this information product, for the most part, is in the public domain, it also may contain copyrighted materials as noted in the text. Permission to reproduce copyrighted items must be secured from the copyright owner.

Please see the license and disclaimer files for additional details.
