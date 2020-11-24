# CADs-pathogens
Relationship between cationic amphipathic drugs (CADs) and pathogens.

[Biorxiv paper - v1](https://www.biorxiv.org/content/10.1101/2020.04.10.035683v1.full.pdf).

# Figure 1C

- CRISPRi_combine_pvals.py. Code to combining p-values for the six CADs screens.


# Figure 1D

-  CRISPRi_CADs_greater_than2.5.sql: get all genes with greater than 2.5x standard deviations from the mean for all six CADs.

# Figure 3B

- pathogens_screen_count.R: code to visualize the common gene hits identified in
genome-wide screens with the 15 pathogens.

# Figure 4A

- RNAseq_CADs_nonCADs.sql: get all genes with greater than 3x fold up or downregulation change with all CADs where dmso control has a counts per million of greater than 0.1.


# Figure 4C

- pathogens_transcriptional_profiling_count.R: code to visualize the common gene hits identified in
transcriptional profiling of mammalian tissues/cells in response to 13 diverse pathogens. Most code to generate the raw data is SQL code, namely in pathogens_transcriptional_profiling.sql. 

[//]: # (- pathogens_transcriptomics.py: gets gene occurrence counts across transcriptional profiling of 13 pathogens.)
- Equivalent code to the SQL queries is available via a Juptyer notebook: https://colab.research.google.com/drive/1CoYDm8LvrK9tRUsFbhweioUD1yF5Ozdy?usp=sharing

# underlying data

- Many of the underlying data was obtained from analysis that can be obtained via https://github.com/tim-peterson/morpheome. Search for CADs, SSRIs, infectious disease (IDs).


