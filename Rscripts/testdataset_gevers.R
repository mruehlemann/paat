### import phyloseq
library(phyloseq); packageVersion("phyloseq")

### path to the greengenes 97% OTU reference tree
path_to_greengenes_tree="<enter path here>"

### the Gevers et al biom file and sample metadata was downloaded from Qiita (https://qiita.ucsd.edu/study/description/1939)
path_to_biom="<enter path here>"
path_to_metadata="<enter path here>"

### load data
x = import_biom(path_to_biom,path_to_greengenes_tree,parseFunction=parse_taxonomy_greengenes)
specs<-read.table(path_to_metadata,head=T,sep="\t",row.names=1,quote="")

### we only keep samples with age > 10 and < 25, without drug intake (antibiotics, biologics, mesalamine or steroids), terminal ileum as biopsy location and only first timepoints
specs<-with(specs,specs[diagnosis %in% c("CD","no") & age>10 & age<25 & is.na(age)==F & antibiotics=="false" & biologics=="false" & biopsy_location=="Terminal ileum" & mesalamine=="false" & steroids=="false" & duplicated(anonymized_name)==F,] )

### create phyloseq object from this
ps<-merge_phyloseq(x,sample_data(specs))

