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
* ggtree (plotting)

PAAT is created to work with phyloseq objects that contain:

1. The OTU/ASV abundance table (taxa_are_rows==T)
2. The Phylogenetic Tree of the OTU/ASV Sequences
3. Sample grouping information (and additional covariates to be corrected for)
4. Taxonomic information about the OTUs to annotate clusters

## Example steps

The test dataset is a subset of the dataset from the article "The treatment-naive microbiome in new-onset Crohn's disease" by Gevers *et al.*. Description on how this dataset was created and filtered is found [here](../master/examples/gevers_generate_testdataset.R).

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

3. We set a threshold to subsample to 3,000 sequences per sample; some samples will be excluded because the fall below this threshold
```
> subset_depth=3000 
> ps.even<-rarefy_even_depth(ps,sample.size=subset_depth,replace=F,rngseed=666)
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
> phylocount<-t(as(otu_table(testdataset),"matrix")) %*% phylomat
```

8. To reduce multiple testing burden and exclude low abundant branches, we define filtering thresholds. We want to remove all branches that are present in less than 25% of the samples in both groups.  
   Additionally, we filter out all branches with less than 0.1% mean abundance across all samples. However, since we also do not want branches that are too broad, we also define an upper threshold of 50% maximum mean abundance.
```
> presence_thresh=0.25
> abu_thresh_upper=0.5
> abu_thresh_lower=0.001
> phylocount.filter<-filter_abundance(abutab=phylocount, groups=model.specs$Group, min.mean=abu_thresh_lower, 
   max.mean=abu_thresh_upper, group.min.presence=presence_thresh)
```

10. To futher reduce redundant branches, parent branches that are more than 95% similar (Bray-Curtis) to one of their sub-branches are excluded.  
   For this we first create a matrix that expresses parent-child relationship, that is 1 if the branch (rows) is a sub-branch of the parent branch (columns). 
```
> phylomat.filter<-phylomat[,colnames(phylocount.filter)]
> phylomat.filter.cor<-cor(phylomat.filter)
> diag(phylomat.filter.cor)<-0
> phylomat.filter.cor[upper.tri(phylomat.filter.cor)]<-0
> phylomat.filter.cor<-ifelse(phylomat.filter.cor>0,1,0)
```
... Then we calculate the branch similarity. By multplying these to matrices field by field, we only get similarity for parent branches to all their sub-branches. We filter out all parent branches that are > 95% similar to one of their child branches.
```
> phylocount.bc<-1-as.matrix(vegdist(t(phylocount.filter)))
> diag(phylocount.bc)<-0
> phylocount.filter<-phylocount.filter[,apply((phylocount.bc * phylomat.filter.cor),2,function(x) any(x>0.95))==F]
```

11. As we are trying to find meaningful clusters that are not too broad, we exclude all branches and their parents, that cover more than 25% of all leaves in the tree and arise from merging two larger subtrees (each of the two subtree consist of > 20% of the total leaves in this branch).
```
> phylomat.jc<-1-as.matrix(vegdist(t(phylomat.filter),"jaccard"))
> excl<-names(which(apply(phylomat.jc*phylomat.filter.cor,2,max)<0.8 & colSums(phylomat.filter)>(nrow(phylomat.filter)*0.25)))
> phylocount.filter<-phylocount.filter[,!(colnames(phylocount.filter) %in% unique(c(excl,names(which(apply(phylomat.filter.cor[c(excl),,drop=F]==1,2,any))))))]
```

12. Now the pre-filtering is completed. The next step is the differential abundance testing. Function `test_abundance` is a wrapper for the glm.nb function from the MASS package. Branches with more than 5% zeroes are only tested for
differential abundance in the non-zero samples. For branches with < 5% zeroes, all samples are included. The GLMs are adjusted for the previously defined covariates.  
```
> outmat<-test_abundance(phylocount.filter,model.specs,testingvar="Group",covars=covar,zerothresh=0.05)
```
   Alternatively, there is also a wrapper which uses relative abundances as input, arcsin-squareroot transforms them and uses a linear model for differential abundance testing. Generally, it is up to the user which test they prefer.
For the purpose of this example, we stick to the arcsin-squareroot method as this is the method used in Gevers *et al.*
```
> phylocount.filter.rel<-phylocount.filter/subset_depth
> outmat<-test_abundance_asinsqrt(phylocount.filter.rel,model.specs,testingvar="Group",covars=covar)
> outmat$p.adj<-p.adjust(outmat$p,"fdr")
> outmat.phylo<-outmat[outmat$p.adj<0.05 & is.na(outmat$p.adj)==F,]
> nrow(outmat.phylo)
[1] 129
```

13. After correction for multiple testing using the Benjamini-Hochberg method, a total of 135 branches are significantly associated with Crohn's disease. However, these results also include nested signals, as signals from parts
closer to the tips influence the results at broader parent-branches. Before we can filter these out, we need to identify clusters from the same part of the phylogenetic tree, that show the same effect directions.
```
> phylomat.filter.cor2<-phylomat.filter.cor
> diag(phylomat.filter.cor2)<-1
> phylomat.filter.cor2.sub<-phylomat.filter.cor2[rownames(outmat.phylo),rownames(outmat.phylo)]
> outmat.phylo$cluster<-factor(apply(phylomat.filter.cor2.sub,1,function(x) min(which(x==1)))*sign(outmat.phylo$beta))
> outmat.phylo$cluster<-as.numeric(factor(outmat.phylo$cluster,levels=unique(outmat.phylo$cluster)))
> unique(outmat.phylo$cluster)
```

14. The 129 signals could be assigned to 21 clusters. To refine the signals and filter out weak results, we take each cluster and order them by absolute effect size. The strongest signal of the cluster is retained and all parent and child branches of this signal are removed from the results list. Signals from the same cluster,
but from neighbouring branches are retained and again ordered by absolute effect size. This itereation is repeated until no independent signals are left in the cluster. Then the procedure is performed in the same manner for the
other clusters. Through this, only the strongest independent signals are retained. 
```
> outmat.phylo.refined<-outmat.phylo[0,]
> for(cl in unique(outmat.phylo$cluster)){
    branch_in_cluster<-outmat.phylo[outmat.phylo$cluster==cl,]
    branch_in_cluster<-branch_in_cluster[order(abs(branch_in_cluster$beta),decreasing=T),]
    while(nrow(branch_in_cluster)>0){
      outmat.phylo.refined[nrow(outmat.phylo.refined)+1,]<-branch_in_cluster[1,]
      branch_in_cluster<-branch_in_cluster[-1,]
      branch_in_cluster<-branch_in_cluster[apply(cor(phylomat.filter[,rownames(branch_in_cluster),drop=F],phylomat.filter[,rownames(outmat.phylo.refined),drop=F]),1,function(x) any(x>0))==F,]
    }
  }
> nrow(outmat.phylo.refined)
[1] 28
```

15. Out of the 129 Signal, 29 are retained as independent. We annotate these branches by weighing the tip labels at the different taxonomic levels by the mean abundances in the dataset. Only if the cumulative weight at a levels
surpasses 0.5 the label is reported at this level. Addtionally, a tag-tips are defined, as the highest abundant tips in the respective branches.
```
> phylocount.final<-phylocount.filter[,rownames(outmat.phylo.refined)]
> phylomat.final<-phylomat.filter[,rownames(outmat.phylo.refined)]

> outmat.annotated<-annotate_branches(phyloseq=testdataset,results=outmat.phylo.refined, phylomat=phylomat.final)
```

16. We will now create a plot of the tree with subtrees highlighted using the results. However, since the complete dataset has ~ 3,500 tips, we first need to filter these before plotting. We use the same theshold as before for the
branch abundances. Then we create the plot and save it.
```
> otu.mat<-as.data.frame(t(otu_table(testdataset)))
> otu.mat.sub<-filter_abundance(abutab=otu.mat, groups=model.specs$Group, min.mean=abu_thresh_lower, 
   group.min.presence=presence_thresh)
> tips_to_keep<-colnames(otu.mat.sub)
> treep<-plot_annotated_tree(phyloseq=testdataset, results=outmat.annotated, 
   phymat=phylomat.final, tips=tips_to_keep)
> ggsave(treep, file=paste0("examples/gevers_",paste0(c(set1),collapse=""),".vs.",paste0(c(set2),collapse=""),".bc.new.pdf"),
   height=16,width=20)
```
   The resulting image (after only slight movements of overlapping/truncated labels):

![gevers-paat](https://github.com/mruehlemann/paat/raw/master/examples/gevers_no.vs.CD.paat.clean.png)
