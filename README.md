# PhAAT (Phylogeny-aware Abundance Testing)

PhAAT is a framework to test for differential abundances in 16S rRNA gene sequencing data.

By using phylogenetic information of (Z)OTUs/ASVs PhAAT can identify signals exhibited by sub-trees of the initial phylogenetic tree.

The method does not depend on classification using databases, which is error-prone due to wrong annotations of references sequences or uncertain assignments. It can be used with denovo generated phylogenetic OTUs/trees or with pre computed reference trees from closed-reference OTU picking (e.g. greengenes).

PhAAT is based on a matrix representation of the bifurcated phylogenetic relationship of the representative sequences. All branches of the phylogenetic tree are in a first step converted to this Sequence-to-branch matrix which is 1 if a sequence is a tip in the respective branch and 0 if not. Through matrix multiplication of the sequence/OTU abundance table
with this phylogenetic matrix, we obtain a branch-abundance matrix. To reduce multiple testing burden, we first filter out branches that are too low abundant/prevalent, addtionally also too high abundant branches can be reduced, as they like represent too broad signals. A second filtering step removes branches that differt too little from their child branches. 
This similarity is defined as the Bray-Curtis similarity of the parent branch to its child branches. We use a default value of > 95% similarity for a branch to be removed. A last step to remove broad signals is implemented by calculating the Jaccard index of a parent branch to its child-branches. If a branche is large, defined as containing > 25% of all tips in 
the tree, and both of its direct child branches have a Jaccard index of > 0.2 compared to the parent branch (meaning they both contain at least 20% of the tips in the parent branch), this branch and all it's parent branches are removed.

Once the filtering is completed, any statistical framework to test for differential abundance can be applied. We implemented two default methods: linear models using arcsin-sqareroot transformed relative abundances (`test_abundance_asinsqrt`) and a negative-binomial GLM approach (`test_abundance`), both with the possibility to include covariates
to be corrected for in the model. After abundance testing, we are likely left with redundant and/or
overlapping results, as a strong differential abundance in a sub-branch will likely lead to a (though less-strong) signal in parent branches. To refine signals, we first cluster sub-branches with the same effect direction and belonging to the same part of the tree. For each of these clusters, we select the best signal (defined as the absolute beta value 
in the statistical model) and remove all parent- and  child-branches to this signal, of the remaining signals we again pick the best and remove all it's related signals. This iteration is carried on until only independent signals remain in the cluster. After performing this for all clusters we are left with the final set of differentially abundant phylogenetic groups.

Additionally we implemented an annotation algorithm for these final signals. The annotation is based on the taxa_table() of the phyloseq object used in the analysis. By weighing the taxonomic annotations of the sequences/OTUs based on their abundance, we try to define the major groups affected by the differential abundance. This is performed by the `annotate_branches` function.

Finally, we implemented a function that visualizes the final results on a plotted tree depicting the phylogenetic signals in diffrential abundance and the annotation.

## Pre-requisites

R Packages:
* MASS (glm.nb function)
* phyloseq 
* vegan (distance calculation)
* ape (tree methods)
* ggtree (plotting)

PhAAT is created to work with phyloseq objects that contain:

1. The OTU/ASV abundance table (taxa_are_rows==T)
2. The Phylogenetic Tree of the OTU/ASV Sequences
3. Sample grouping information (and additional covariates to be corrected for)
4. Taxonomic information about the OTUs to annotate clusters

## Example: Treatment-naive microbiota in inflammatory bowel disease

The test dataset is a subset of the dataset from the article "The treatment-naive microbiome in new-onset Crohn's disease" by Gevers *et al.*. Description on how this dataset was created and filtered is found [here](../master/examples/gevers_generate_testdataset.R).

Another dataset used for demonstration is the data from [Giloteaux et al. (2016)](https://doi.org/10.1186/s40168-016-0171-4). Here we generated the dataset using the DADA2 pipeline, the scripts can be found [here](../master/examples/giloteaux_generate_testdataset.R).
The analysis was performed in analogy to the following step, the code can be found [here](../master/examples/giloteaux_code.R).

1. load `phyloseq` and PhAAT functions
```
> library(phyloseq)
> library(vegan)
> library(ape)
> source("https://github.com/mruehlemann/phaat/raw/master/Rscripts/phaat_functions.R")
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

![gevers-phaat](https://github.com/mruehlemann/phaat/raw/master/examples/gevers_no.vs.CD.paat.clean.png)

17. We also apply the same linear model framework to the OTU abundance table and plot the tree.
```
> res.otu<-test_abundance_asinsqrt(abundance_data=otu.mat.sub/subset_depth, sample_specs=model.specs, testingvar="Group", covars=covar)
> res.otu$p.adj<-p.adjust(res.otu$p,"fdr")
> outmat.otu<-res.otu[res.otu$p.adj<0.05 & is.na(res.otu$p.adj)==F,]
> phylomat.otu<-ifelse(cor(otu.mat.sub)==1,1,0)
> outmat.otu.annotated<-annotate_branches(phyloseq=testdataset,results=outmat.otu, phylomat=phylomat.otu)
> treep.otu<-plot_annotated_tree(phyloseq=testdataset, results=outmat.otu.annotated, phymat=phylomat.otu, tips=c(tips_to_keep,outmat.annotated$tag))
> ggsave(treep.otu, file=paste0("examples/gevers_",paste0(c(set1),collapse=""),".vs.",paste0(c(set2),collapse=""),".otu.pdf"), height=16,width=20)
```
   It shows, that the PhAAT picks up all signals also seen in the OTU based analysis, plus additional signals not seen when analyzing OTUs (e.g. Prevotella copri)

![gevers-otu](https://github.com/mruehlemann/phaat/raw/master/examples/gevers_no.vs.CD.otu.clean.png)

18. We can also perform the calculations on groups using taxonomic assigments. For this we collapse OTUs with the same annotation down to genus level into clusters and use them in differential abundance calculation.
```
> tax.gen<-unique(tax_table(testdataset)[,1:6])
> tt.gen<-tax_table(testdataset)[,1:6]
> alltax<-apply(tax.gen,1,paste,collapse=";")
> cl<-data.frame(id=rownames(tt.gen), cluster=sapply(apply(tt.gen,1,paste,collapse=";"),function(x) match(x, alltax)))
> library(reshape2)
> ii<-data.frame(dcast(id ~ cluster,data=cl,fill=0),row.names=1)
> ii[ii>0]<-1
> ii<-ii[colnames(otu.mat),]
> colnames(ii)<-as.character(alltax)

> gen.mat<-as.matrix(otu.mat) %*% as.matrix(ii)
> gen.mat.sub<-filter_abundance(abutab=gen.mat, groups=model.specs$Group, min.mean=abu_thresh_lower, group.min.presence=presence_thresh)
> res.gen<-test_abundance_asinsqrt(abundance_data=gen.mat.sub/subset_depth, sample_specs=model.specs, testingvar="Group", covars=covar)
> res.gen$p.adj<-p.adjust(res.gen$p,"fdr")
> outmat.gen<-res.gen[res.gen$p.adj<0.05 & is.na(res.gen$p.adj)==F,]

> otu_to_siggen<-names(which(rowSums(ii[,rownames(outmat.gen)])>0))
> res.gen<-res.otu[rownames(res.otu) %in% otu_to_siggen,]

> outmat.gen.annotated<-res.gen
> outmat.gen.annotated$tag<-rownames(outmat.gen.annotated)
> outmat.gen.annotated$tax<-apply(tt.gen[rownames(outmat.gen.annotated),],1,function(x) x[max(which(is.na(x)==F))])
> outmat.gen.annotated$ntips<-colSums(ii)[match(apply(tt.gen[rownames(outmat.gen.annotated),],1,paste, collapse=";"),alltax)]

> treep.gen<-plot_annotated_tree(phyloseq=testdataset, results=outmat.gen.annotated, phymat=phylomat.otu, tips=c(tips_to_keep,outmat.annotated$tag))
> ggsave(treep.gen, file=paste0("examples/gevers_",paste0(c(set1),collapse=""),".vs.",paste0(c(set2),collapse=""),".gen.pdf"), height=16,width=20)
```

   We see, that the overall signals are in some regions much broader than for the other methods. Faecalibacterium signal is not picked up, likely because of two different signals in different direction within the Faecalibacterium genus.  

![gevers-otu](https://github.com/mruehlemann/phaat/raw/master/examples/gevers_no.vs.CD.gen.clean.png)

```
> sessionInfo()
R version 3.4.1 (2017-06-30)
Platform: x86_64-apple-darwin16.6.0 (64-bit)
Running under: macOS High Sierra 10.13

Matrix products: default
BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libLAPACK.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] reshape2_1.4.2  ggtree_1.8.2    treeio_1.0.2    ggplot2_2.2.1  
 [5] MASS_7.3-47     ape_4.1         vegan_2.4-4     lattice_0.20-35
 [9] permute_0.9-4   phyloseq_1.20.0

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.15        compiler_3.4.1      plyr_1.8.4         
 [4] XVector_0.16.0      iterators_1.0.8     tools_3.4.1        
 [7] zlibbioc_1.22.0     jsonlite_1.5        tibble_1.3.4       
[10] nlme_3.1-131        rhdf5_2.20.0        gtable_0.2.0       
[13] mgcv_1.8-17         pkgconfig_2.0.1     rlang_0.1.2        
[16] Matrix_1.2-10       foreach_1.4.4       igraph_1.1.2       
[19] rvcheck_0.0.9       parallel_3.4.1      stringr_1.2.0      
[22] cluster_2.0.6       Biostrings_2.44.2   S4Vectors_0.14.4   
[25] IRanges_2.10.3      multtest_2.32.0     stats4_3.4.1       
[28] ade4_1.7-8          grid_3.4.1          glue_1.1.1         
[31] Biobase_2.36.2      data.table_1.10.4   survival_2.41-3    
[34] purrr_0.2.3         tidyr_0.7.1         magrittr_1.5       
[37] splines_3.4.1       scales_0.5.0        codetools_0.2-15   
[40] BiocGenerics_0.22.0 biomformat_1.4.0    colorspace_1.3-2   
[43] labeling_0.3        stringi_1.1.5       lazyeval_0.2.0     
[46] munsell_0.4.3      
```
