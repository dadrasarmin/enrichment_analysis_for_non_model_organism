# Enrichment analysis in R for non-model organism(s)
In this tutorial, I show how to perform enrichment analysis using two packages called `Tidyverse` and `clusterProfiler` for a non-model organism. The tutorial is designed to be a demo. Feel free to adjust the protocol accroding to your need.
## Set up the environment
You need at least two packages to do this analysis: `Tidyverse` and `clusterProfiler`. If you do not have them, please install them as follow:
```
install.packages(c("BiocManager", "devtools"))
library("devtools")
library("BiocManager")
BiocManager::install("tidyverse")
devtools::install_github("GuangchuangYu/GOSemSim")
devtools::install_github("GuangchuangYu/clusterProfiler")
library(tidyverse)
library(clusterProfiler)
```
If there is any error, please read the error message and try to install the missing package. For some OS, you need to install `Rcpp` or other packages before you start the installation for these packages.

## Similar but different: Over-Representation Analysis (ORA) vs. Gene Set Enrichment Analysis (GSEA)
If you are completely new to the topic, I highly suggest that read about these two topics elsewhere:
1. [A short video on the subject](https://www.youtube.com/watch?v=B7F7a9NcGS0)
2. [clusterProfiler webpage](https://yulab-smu.top/biomedical-knowledge-mining-book/enrichment-overview.html)

## What do you need to do ORA?
We basically need 4 tables to do OAR:
* All genes that are present in our study
* A list of genes that for some reason are interesting for us
* A table that matches Term to Gene ID
* A table that matches Term to names

## GO enrichment analysis
First, I cover GO enrichment analysis.
### Preparing Term to Gene table
Let's use an example. Let's say I am interested in an alga called *Meostaenium endlicherianum*. I can perform functional annotation for its proteome (or genes) using different bioinformatic tools such as [eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper), [InterProScan](https://interproscan-docs.readthedocs.io/en/latest/HowToRun.html), [Funannotate](https://funannotate.readthedocs.io/en/latest/), [BLAST](https://www.ncbi.nlm.nih.gov/books/NBK569856/), etc. Let's assume that I used eggNOG-mapper to predict functional annotation for *M. endlicherianum* proteins and I have a table as `eggNOG_mapper_m.endlicherianum_protein.emapper.annotations.tsv`. Let's assume that I want to perform analysis at gene level (instead of transcript level) and I assume that I can collapse attributes of all isoforms of a gene to that gene. Therefore, we can use the following piece of code to create term to gene table.
```
# prepare the term to gene table
eggNOG <- read_tsv("eggNOG_mapper_m.endlicherianum_protein.emapper.annotations.tsv") %>%
    dplyr::select(GOs, `#query`) %>%
    dplyr::filter(GOs != "-") %>%
    separate_rows(GOs, sep = ",") %>%
    mutate(gene = gsub("\\..*", "", `#query`)) %>%
    select(GOs, gene) %>%
    distinct() %>%
    drop_na()
colnames(eggNOG) <- c("term", "gene")
```
### Preparing Term to name table
There are different ways to map a GO term to a GO term. I used core ontology from [gene ontology](https://geneontology.org/docs/download-ontology/) database. You can also look at different GO-term tables on different databases such as [OLS](https://www.ebi.ac.uk/ols/ontologies/go).
```
# prepare the term to name table
ontology <- get_ontology(file = "go.obo",
                         propagate_relationships = "is_a",
                         extract_tags = "everything",
                         merge_equivalent_terms = TRUE)
eggNOG_term <- eggNOG %>%
    mutate(name = ontology$name[term]) %>%
    select(c(term, name)) %>%
    distinct() %>%
    drop_na() %>%
    filter(!grepl("obsolete", name))

eggNOG <- eggNOG %>%
    filter(term %in% eggNOG_term$term)
```
### Save both tables
I can save both tables as follows:
```
# save the results
write_tsv(x = eggNOG, file = "term2gene_GO.tsv")
write_tsv(x = eggNOG_term, file = "term2name_GO.tsv")
```
## The background gene list
We call all genes that are present in a study background. For example, if you are analyzing a RNA-Seq dataset, all genes that are expressed form your background (and not all annotated genes in the genome your are investigating). Let's assume I have a table that contains all genes that I believe they are expressed in my samples called `counts_filetered_cpm.tsv`. I can create a background list as follows:
```
background_genes <- read_tsv("counts_filetered_cpm.tsv") %>%
    dplyr::select("geneID") %>%
    unlist() %>%
    as.vector()
```
## The gene set of interest
Let's assume that I have a table that contains log2(fold change) values of genes and adjusted p-values for different genes called `comparison_6.tsv` (data is from Dadras et al. 2023)[https://www.nature.com/articles/s41477-023-01491-0]. We can calculate these tables based on RNA-Seq tables via edgeR, limma, DESeq2, etc. Please keep in mind that this is just an example. You can perform enrichment analysis with any kind of input data and it is not limited to RNA-Seq.
Let's assume that I am only instersted in genes that have abosulte fold change of 2 or more with adjusted p-value of at least 0.05. Then, I can create the list of interesting genes via:
```
# read the gene list of interest
interesting_set <- read_tsv("comparison_6.tsv") %>%
    dplyr::filter(abs(logFC) >= 1 & adj.P.Val <= 0.05) %>%
    dplyr::select("geneID") %>%
    unlist() %>%
    as.vector()
```

## ORA
Now, I have prepared everything that I need. Therefore, I can perform ORA via clusterProfiler package with specific p- and q-values. For details of function input, please use the documentation of clusterProfiler. I saved the results of enrichment analysis as `enrichmet_results.csv`. I also visualized the enrichment results via a dotplot. For more information about visualization please read the package [documantation](https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html).
```
# perform ORA
term2gene <- read_tsv("term2gene_GO.tsv")
term2name <- read_tsv("term2name_GO.tsv")

enrichment <- enricher(interesting_set,
                       TERM2GENE = term2gene,
                       TERM2NAME = term2name,
                       pvalueCutoff = 0.05,
                       universe = background_genes,
                       qvalueCutoff = 0.05)
#save the enrichment result
write.csv(file = paste0("enrichment_results.csv"),                 # EDIT THIS
          x = enrichment@result)

if (any(enrichment@result$p.adjust <= 0.05)){
    p <- dotplot(enrichment,
                 x= "geneRatio", # Options: GeneRatio, BgRatio, pvalue, p.adjust, qvalue
                 color="p.adjust",
                 orderBy = "x", # Options: GeneRatio, BgRatio, pvalue, p.adjust, qvalue
                 showCategory=100,
                 font.size=8) +
        ggtitle("dotplot for GO ORA")
    
    ggsave(filename = paste0("enrichment_dotplot.pdf"),                # EDIT THIS
           plot =  p,  dpi = 300, width = 21, height = 42, units = "cm")
}
```

## What about KEGG?
If you are interested to perform KEGG enrichment analysis, you have to create a table that annotate genes with KEGG ortholog IDs. There are different ways to do that including [KAAS](https://www.genome.jp/kegg/kaas/) and [eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper). Let's assume I want to use eggNOG-mapper for this tutorial. We have to prepare similar term to gene and term to name tables for KEGG. However, since clusterProfiler can work with KEGG ortholog under the hood, the task is a little bit easier. I saved the enrichment results as a table called `enrichment_KEGG_restuls.csv` and visualized it via dotplot as `enrichment_KEGG_dotplot.pdf`. The code is as follows:
```
# KEGG ORA
eggNOG_kegg <- read_tsv("eggNOG_mapper_m.endlicherianum_protein.emapper.annotations.tsv") %>%
    dplyr::select(KEGG_ko, `#query`) %>%
    dplyr::filter(KEGG_ko != "-") %>%
    separate_rows(KEGG_ko, sep = ",") %>%
    dplyr::mutate(gene = gsub("\\..*", "", `#query`)) %>%
    dplyr::mutate(term = gsub("ko:", "", KEGG_ko)) %>%
    dplyr::select(term, gene) %>%
    distinct() %>%
    drop_na()
# create a list of kegg ortholog that I am interested in
interesting_set_kegg <- eggNOG_kegg %>%
    dplyr::filter(gene %in% interesting_set) %>%
    unlist() %>%
    as.vector()
# create a list of kegg ortholog that includes all kegg orthologs which form my background
background_kegg <- eggNOG_kegg %>%
    dplyr::filter(gene %in% background_genes) %>%
    unlist() %>%
    as.vector()

enrichment_kegg <- enrichKEGG(interesting_set_kegg,
           organism = "ko",
           keyType = "kegg",
           pvalueCutoff = 0.05,
           pAdjustMethod = "BH",
           universe = background_kegg,
           minGSSize = 10,
           maxGSSize = 500,
           qvalueCutoff = 0.05,
           use_internal_data = FALSE)
#save the enrichment result
write.csv(file = paste0("enrichment_KEGG_results.csv"),                 # EDIT THIS
          x = enrichment_kegg@result)

if (any(enrichment_kegg@result$p.adjust <= 0.05)){
    p <- dotplot(enrichment_kegg,
                 x= "geneRatio", # Options: GeneRatio, BgRatio, pvalue, p.adjust, qvalue
                 color="p.adjust",
                 orderBy = "x", # Options: GeneRatio, BgRatio, pvalue, p.adjust, qvalue
                 showCategory=100,
                 font.size=8) +
        ggtitle("dotplot for KEGG ORA")
    
    ggsave(filename = paste0("enrichment_KEGG_dotplot.pdf"),                # EDIT THIS
           plot =  p,  dpi = 300, width = 21, height = 42, units = "cm")
}

```
## What about GSEA?
If you want to perform GSEA, instead of a list of genes that are filetered based on some criterias, we need a list of genes and their changes (let's say fold changes). We can do this analysis similar to ORA as follows. I saved the results in a table called `enrichment_GSEA_results.csv` and visualized it via a ridge plot called `enrichment_GSEA_ridgeplot.pdf`.
```
# GSEA
# This time I am not filtering the changes
gsea_list <- read_tsv("comparison_6.tsv") %>%
    dplyr::arrange(desc(logFC))
gsea_input <- gsea_list %>%
    dplyr::select(logFC) %>%
    unlist() %>%
    as.vector()
names(gsea_input) <- gsea_list$geneID
# do the analysis below
enrichment_gsea <- GSEA(geneList = gsea_input,
                        TERM2GENE = term2gene,
                        TERM2NAME = term2name,
                        minGSSize = 10,
                        maxGSSize = 500,
                        eps = 1e-10,
                        pvalueCutoff = 0.05,
                        pAdjustMethod = "BH")
#save the enrichment result
write.csv(file = paste0("enrichment_GSEA_results.csv"),                 # EDIT THIS
          x = enrichment_gsea@result)

if (any(enrichment_gsea@result$p.adjust <= 0.05)){
    p <- ridgeplot(enrichment_gsea,
                   core_enrichment= FALSE,
                   fill="p.adjust",
                   orderBy = "NES",
                   showCategory=100) +
        ggtitle("Ridge plot for GSEA")
    
    ggsave(filename = paste0("enrichment_GSEA_ridgeplot.pdf"),                # EDIT THIS
           plot =  p,  dpi = 300, width = 21, height = 42, units = "cm")
}
```

## Final advice
There are many more ways to perform analysis and visualize your data using [clusterProfiler](https://yulab-smu.top/biomedical-knowledge-mining-book/index.html) or other packages such as [g:Profiler](https://biit.cs.ut.ee/gprofiler/gost). Always try to read the documents and try new things for yourself. Enjoy.