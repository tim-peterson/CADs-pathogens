# sphingolipids-pathogens-aging
Relationship between sphingolipids, pathogens, and aging. Uses cationic amphipathic drugs (CADs) to unveil the relationship of sphingolipids, which are natural cationic amphipaths, to pathogens and aging.

[Biorxiv paper - v1](https://www.biorxiv.org/content/10.1101/2020.04.10.035683v1.full.pdf).

# Figure 1C

- CRISPRi_combine_pvals.py. Code to combining p-values for the six CADs screens.


# Figure 1D
- [pipeline to generate CADs and non-CAD screens all in one file](https://github.com/tim-peterson/morpheome/blob/master/app/Console/Commands/Morpheome/pipeline/top_hits_small_mol_screens_with_external_morpheome.php)

-  [CRISPRi_CADs_greater_than2.5.sql](https://github.com/tim-peterson/CADs-pathogens/blob/master/CRISPRi_CADs_greater_than2.5.sql): get all genes with greater than 2.5x standard deviations from the mean for all six CADs.

# Figure 2B

- [pathogens_screen_count_graphing.R](https://github.com/tim-peterson/CADs-pathogens/blob/master/pathogens_screen_count_graphing.R): code to visualize the common gene hits identified in
genome-wide screens with the 15 pathogens.

# Figure 2C

- [infectious_agent_intersect_nonCADs_CADs.sql](https://github.com/tim-peterson/morpheome/blob/master/sql_commands/ssri/infectious_agent_intersect_nonCADs_CADs.sql)
# Figure 3A

- RNAseq_CADs_nonCADs.sql: get all genes with greater than 3x fold up or downregulated with all CADs where dmso control has a counts per million of greater than 0.1.


# Figure 3B

- [pathogens_transcriptional_profiling_count_graphing.R](https://github.com/tim-peterson/CADs-pathogens/blob/master/pathogens_transcriptional_profiling_count_graphing.R): code to visualize the common gene hits identified in transcriptional profiling of mammalian tissues/cells in response to 13 diverse pathogens. Most code to generate the raw data is SQL code, namely in pathogens_transcriptional_profiling.sql. 

[//]: # (- pathogens_transcriptomics.py: gets gene occurrence counts across transcriptional profiling of 13 pathogens.)
- Equivalent code to the SQL queries is available via a Juptyer notebook: https://colab.research.google.com/drive/1CoYDm8LvrK9tRUsFbhweioUD1yF5Ozdy?usp=sharing

# Figure 4A

[aging_gwas_count_graphing.R](https://github.com/tim-peterson/CADs-pathogens/blob/master/aging_gwas_count_graphing.R): code to visualize the genes in the top 200 by p-value counted by aging-related disease. P-values for each gene for each disease were aggregated using this Juptyer notebook script: https://colab.research.google.com/drive/1WpuKph4o8hZgMC7mJU_SpBRipWJy7A8C?usp=sharing 

# Figure 4

Gene expression to identify drugs with broad anti-pathogen activity: https://colab.research.google.com/drive/1_wuYJcJgLK3MFl33f5IdhFF2U5tBmbum?usp=sharing 

# additional underlying data

- Additional underlying data was obtained from analysis from https://github.com/tim-peterson/morpheome. Search for CADs, SSRIs, infectious disease (IDs).


