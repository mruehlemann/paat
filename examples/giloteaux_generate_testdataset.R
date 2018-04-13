#### all data were downloaded from the EBI-ENA Servers, Accession PRJEB13092; sample metadata was provided as supplement to the original article by Giloteaux et al https://doi.org/10.1186/s40168-016-0171-4

library(dada2); packageVersion("dada2")

path <- getwd()
list.files(path)

fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(150,150),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

pdf()
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
plotErrors(errF, nominalQ=TRUE)
dev.off()

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

dadaFs[[1]]

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

seqtab.all <- makeSequenceTable(mergers)
dim(seqtab.all)

seqtab <- removeBimeraDenovo(seqtab.all, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab)

#### taxonomic assignment

taxa.silva <- assignTaxonomy(seqtab, "/ifs/data/nfs_share/sukmb276/references/dada2_ref/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
taxa.silva.spe <- addSpecies(taxa.silva, "/ifs/data/nfs_share/sukmb276/references/dada2_ref/silva_species_assignment_v132.fa.gz")

taxa.rdp <- assignTaxonomy(seqtab, "/ifs/data/nfs_share/sukmb276/references/dada2_ref/rdp_train_set_16.fa.gz", multithread=TRUE)
taxa.rdp.spe <- addSpecies(taxa.rdp, "/ifs/data/nfs_share/sukmb276/references/dada2_ref/rdp_species_assignment_16.fa.gz")


#### phylogenetic tree construction

library(dada2); packageVersion("dada2")
library(phangorn); packageVersion("phangorn")
library(DECIPHER); packageVersion("DECIPHER")

seqs <- getSequences(seqtab)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)


phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)

## negative edges length changed to 0!

fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                      rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)

### phyloseq
library(phyloseq); packageVersion("phyloseq")

specs<-read.table("specs.csv",head=T,row.names=1,dec=",")
specs<-specs[rownames(seqtab),]

colnames(seqtab)<-paste("ASV",seq(1:ncol(seqtab)),sep="_")
rownames(taxa.silva)<-rownames(taxa.silva.spe)<-rownames(taxa.rdp)<-rownames(taxa.rdp.spe)<-paste("ASV",seq(1:ncol(seqtab)),sep="_")
fitGTR$tree$tip.label<-paste("ASV",seq(1:ncol(seqtab)),sep="_")

ps <- phyloseq(tax_table(taxa.rdp), sample_data(specs),otu_table(seqtab, taxa_are_rows = FALSE),phy_tree(fitGTR$tree))

