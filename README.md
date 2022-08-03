## Trajectories of Life Expectancy by Commuting Zone in the US: 1990-2019.

Under Review at American Journal of Public Health.

**NCHS MORTALITY DATA folder contains:**

Originally, the NCHS vital registration data. However, since this data cannot be shared publicly, we just put a link to the data request. 

**POPULATION DATA folder contains:**

Bridged-race post-censal population estimates obtained from Census Bureau (https://www.cdc.gov/nchs/nvss/bridged_race/data_documentation.htm)

**LE Trajectories folder contains:**
- Analysis
  - CZ_LE_prep.R
    - uses NCHS vital registration data and population denomoinators to generate _1990_2019_cbsa_mortality.rdata_, master data for overall analysis that contains fips, year (1990-2019), age 5 yr group (0,1,5,10...), sex (M, F), race (NHW, NHB, H, NHAPI, NHAIAN), death counts, population denominators, CBSA delineation, and CZ delineation. 
    - see CZ_LE_prep Documentation.docx for further detail.
    - Please note that we cannot share the rdata files publicly since they are summaries of the NCHS data (see above) with cells with < 10 counts. 
  - CZ_LE_analysis.R 
    - main R script for analyses conducted in manuscript. Please enable R script document outline pane for an overview and navigation to sections.
  - LifeTableFUN.R
    - Helper code/function to create lifetables in CZ_LE_analysis.R 
  - Crosswalks
    - contains crosswalks from FIPS to CBSA.
    - contains crosswalks from FIPS to CZ.
- Tables & Figures
  - contains all ggplot figures generated from CZ_LE_analysis.R
