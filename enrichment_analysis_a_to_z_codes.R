# set up the environment
install.packages(c("BiocManager", "devtools", "ontologyIndex"))
library("devtools")
library("BiocManager")
BiocManager::install("tidyverse")
devtools::install_github("GuangchuangYu/GOSemSim")
devtools::install_github("GuangchuangYu/clusterProfiler")
library(tidyverse)
library(clusterProfiler)
library(ontologyIndex)
# prepare the term to gene table
eggNOG <- read_tsv("eggNOG_mapper_m.endlicherianum_protein.emapper.annotations.tsv") %>%
    dplyr::select(GOs, `#query`) %>%
    dplyr::filter(GOs != "-") %>%
    separate_rows(GOs, sep = ",") %>%
    dplyr::mutate(gene = gsub("\\..*", "", `#query`)) %>%
    dplyr::select(GOs, gene) %>%
    distinct() %>%
    drop_na()
colnames(eggNOG) <- c("term", "gene")

# prepare the term to name table
ontology <- get_ontology(file = "go.obo",
                         propagate_relationships = "is_a",
                         extract_tags = "everything",
                         merge_equivalent_terms = TRUE)
eggNOG_term <- eggNOG %>%
    mutate(name = ontology$name[term]) %>%
    dplyr::select(c(term, name)) %>%
    distinct() %>%
    drop_na() %>%
    dplyr::filter(!grepl("obsolete", name))

eggNOG <- eggNOG %>%
    filter(term %in% eggNOG_term$term)
# save the results
write_tsv(x = eggNOG, file = "term2gene_GO.tsv")
write_tsv(x = eggNOG_term, file = "term2name_GO.tsv")
# read the background gene
background_genes <- read_tsv("counts_filetered_cpm.tsv") %>%
    dplyr::select("geneID") %>%
    unlist() %>%
    as.vector()
# read the gene list of interest
interesting_set <- read_tsv("comparison_6.tsv") %>%
    dplyr::filter(abs(logFC) >= 1 & adj.P.Val <= 0.05) %>%
    dplyr::select("geneID") %>%
    unlist() %>%
    as.vector()
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