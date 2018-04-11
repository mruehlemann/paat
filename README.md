# PAAT (Phylogeny-aware Abundance Testing)

PAAT is a framework to test for differential abundances in 16S rRNA gene sequencing data.

By using phylogenetic information of (Z)OTUs/ASVs PAAT can identify signals exhibited by sub-trees of the initial phylogenetic tree.

The method does not depend on classification using databases, which is error-prone due to wrong annotations of references sequences or uncertain assignments.

## Prerequisites

R Packages:
* MASS (glm.nb function)
* phyloseq 
* vegan (distance calculation)
* ape (tree methods)

PAAT is created to work with phyloseq objects that contain:

1. The OTU/ASV abundance table (taxa_are_rows==T)
2. The Phylogenetic Tree of the OTU/ASV Sequences
3. Sample grouping information (and additional covariates to be corrected for)
4. Taxonomic information about the OTUs to annotate clusters

## Example steps

The test dataset is a subset of the dataset from the article "The treatment-naive microbiome in new-onset Crohn's disease" by Gevers *et al.*. Description on how this dataset was created and filtered is found [here](../master/Rscripts/testdataset_gevers.R).

1. load `phyloseq` and PAAT functions
```
> library(phyloseq)
> library(vegan)
> library(ape)
> source("https://github.com/mruehlemann/paat/raw/master/Rscripts/paat_functions.R")
```

2. load testdataset `ps`
```
> load("gevers.Robj")
> ps
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 9965 taxa and 321 samples ]
sample_data() Sample Data:       [ 321 samples by 56 sample variables ]
tax_table()   Taxonomy Table:    [ 9965 taxa by 7 taxonomic ranks ]
phy_tree()    Phylogenetic Tree: [ 9965 tips and 9964 internal nodes ]
```

3. We set a threshold to subsample to 3000 sequences per sample; some samples will be excluded because the fall below this threshold
```
> subset_depth=3000 
> ps.even<-rarefy_even_depth(ps,sample.size=subset_depth,replace=F)
```

4. Now we define which variable includes the information on group membership which will be used for testing. Additionally we define the sets we use for contrasting (no IBD vs. CD) and which covariates are used for correction in the models (age and sex).
```
> testvar="diagnosis"
> set1=c("no")
> set2=c("CD")
> covar=c("sex","age")
```

5. Now we filter out all samples not belonging to our target groups and design data frame which includes all variables for testing.
```
> testdataset<-prune_samples(unlist(sample_data(ps.even)[,testvar]) %in% c(set1,set2), ps.even)
> model.specs<-data.frame(sample_data(testdataset))[,c(covar,testvar)]
> model.specs$Group<-factor(ifelse(model.specs[,testvar]==set1,"G1","G2"),levels=c("G1","G2"))
```

6. Convert Phylogenetic Tree to Branch-Matrix
```
> phylomat<-branchNodeAdjacency(phy_tree(testdataset))
```

7. Matrix multiplication of the abundance matrix with the Branch-Matrix gives branch abundances (in counts)
```
> phylocount<-as(otu_table(testdataset),"matrix") %*% phylomat
```

8. To reduce multiple testing burden and exclude low abundant branches, we define filtering thresholds. We want to remove all branches that are present in less than 25% of the samples in both groups. 
...Additionally, we filter out all branches with less than 0.1% mean abundance across all samples. However, since we also do not want branches that are too broad, we also define an upper threshold of 50% maximum mean abundance.

```
> presence_thresh=0.25
> abu_thresh_upper=0.5
> abu_thresh_lower=0.001
> phylocount.filter<-phylocount[,colMeans(phylocount/subset_depth)<abu_thresh_upper & colMeans(phylocount/subset_depth)>abu_thresh_lower & apply(aggregate(. ~ model.specs$Group,data=data.frame(ifelse(phylocount>0,1,0)),mean)[,-1],2,min)>presence_thresh] /subset_depth
```

