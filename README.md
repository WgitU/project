# Integrative single-cell clustering analysis

# Wu

## Data

### Abstract 

The data used for the analyses consist of 61 mouse serum embryonic stem cells. This dataset is publicly available with the accession code GSE74535 on the website of NCBI . The data sets consists of 61 ESCs that were measured by scRNA-seq expression of 7356 genes and DNA methylation levels at more than three million positions. 

### Availability 

The data are publicly available for download via the online data portal at (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74535). No registration is required.

### Description 

All output from this project published online is available according to the conditions of the Creative Commons License
(https://creativecommons.org/licenses/by-nc-sa/2.0/).

In the preprocessing procedure, the same DNA methylation positions on the promoters for all ESCs were selected since methylation in the vicinity of the promoter is associated with the absence of transcription. For RNA expression, 500 most variable genes were used for clustering analysis for a good cell cluster separation and the scRNA-seq data were normalized to correct for technical factors: the library size for each cell, the sum of read counts across all genes, and the median of all library sizes were calculated and the original counts were divided by its corresponding library size and multiplied the ratio by the median library size. We took the floor of normalized data. Therefore, both the 500 gene expression and 89 DNA methylation positions data were presented as counts and used in the clustering analysis.

However, this process is extremely time-consuming; therefore, we
also provide a preprocessed .RData file.

## Code

### Abstract

All of the data processing and analysis for this report were done in R. The corresponding code is provided to take exploratory data analysis on the raw data; conduct various preprocessing steps; fit a set of Bayesian models to the preprocessed data via Markov chain Monte Carlo (MCMC) methods; and generate descriptive plots used in the
report.

### Description

All of the R scripts used in the report are available in a public repository on GitHub [https://github.com/WgitU/project]. The MIT license applies to all code, and no permissions are required to access the code.

### Optional Information

R version 3.5.1 (2018-07-02, “Feather Spray”) was used for the analyses in this report. The
necessary R libraries for the code used for data processing and analysis are:

- readr, version 1.3.1 (https://CRAN.R-project.org/package=readr)
- dplyr, version 0.8.0.1 (https://CRAN.R-project.org/package=dplyr)
- plotly, version 4.8.0 (https://CRAN.R-project.org/package=plotly)
- ggplot2, version 3.1.0 (http://ggplot2.org)
- Rtsne, version 0.15 (https://CRAN.R-project.org/package=Rtsne)
- factoextra, version 1.0.5 (https://CRAN.R-project.org/package=factoextra)
- gridExtra, version 2.3 (https://CRAN.R-project.org/package=gridExtra)
- Seurat, version 2.3.4 (https://CRAN.R-project.org/package=Seurat)

## Instructions for Use

### Reproducibility

All data preparation and analyses are reproduced, as well as all Figures in the
report.

All workflow information is contained in the MASTER_reproducibility.R script. The general steps
are:

1. Take exploratory data analysis on the raw data.
2. Conduct data processing/preparation for the analyses.
3. Fit a set of Bayesian models to the preprocessed data via Markov chain Monte Carlo (MCMC) methods
4. Generate all plots in the report.