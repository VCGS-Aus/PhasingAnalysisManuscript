# Phasing Analysis Code and Results
This repository contains the analysis scripts, data-processing code, and results associated with the paper:
**Genome-wide phasing in routine diagnostic samples using long-read sequencing**.


## Overview
This repository provides the code and supporting materials required to reproduce the analyses and results presented in the paper. It's organized to separate data processing, statistical analysis, figure/table generation, and resulting outputs.

```
.
├── README.md
├── LICENSE
├── data/
│   ├── mendeliome_incidentalome_no_overlap_autosomes.bed    # Bed file with autosomal gene coordinates from Mendeliome and Incidentalome gene panel from PanelApp Australia
│   └── clinvar_p_lp.bed                                     # Data generated during preprocessing
│
├── scripts/
│   ├── extract_phase_sets_to_calculate_phasing_accuracy.sh  # Script to extract phase sets from vcf used to calculate phasing accuracy
│   ├── calculate_phasing_accuracy.py                        # Script to calculate phasing accuracy using phase set data extracted using extract_phase_sets_to_calculate_phasing_accuracy.sh
│   ├── calculate_pairwise_phasing_accuracy.py               # Script to calculate pairwise phasing accuracy using phase set data from  extracted using extract_phase_sets_to_calculate_phasing_accuracy.sh
│   └── extract_phase_blocks.groovy                          # Script to extract phase blocks from bam
│   ├── pairwise_comparison_clinvar_variants.groovy          # Script to calculate pairwise phasing probability of ClinVar variants using phase blocks extracted using extract_phase_blocks.groovy  
│   └── calculate_heterozygosity.sh                          # Calculates heterozygosity from vcf files
│
├── manuscript_figures.html                                  # HTML with code used to generate figures, tables and statistical output for the manuscript
```
