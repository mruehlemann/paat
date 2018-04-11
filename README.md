# PAAT (Phylogeny-aware Abundance Testing)

PAAT is a framework to test for differential abundances in 16S rRNA gene sequencing data.

By using phylogenetic information of (Z)OTUs/ASVs PAAT can identify signals exhibited by sub-trees of the initial phylogenetic tree.

The method does not depend on classification using databases, which is error-prone due to wrong annotations of references sequences or uncertain assignments.

## Prerequisites

R Packages:
* MASS (glm.nb function)
* phyloseq 

PAAT is created to work with phyloseq objects that contain:

1. The OTU/ASV abundance table
2. The Phylogenetic Tree of the OTU/ASV Sequences
3. Sample grouping information (and additional covariates to be corrected for)

## Example steps

The test dataset is a subset of the dataset from the article "The treatment-naive microbiome in new-onset Crohn's disease" by Gevers *et al.*. Description on how this dataset was created and filtered is found [here](../blob/master/Rscripts/testdataset_gevers.R).

Assuming you have your phyloseq object (e.g. from DADA2 output) named `ps`

1. load `phyloseq` and PAAT functions
```
library(phyloseq)
source("https://github.com/mruehlemann/paat/raw/master/Rscripts/paat_functions.R")
```

2. 
