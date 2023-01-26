#!/usr/bin/env Rscript

library(clusterProfiler)
library(org.Hs.eg.db)
library(data.table)
library(tidyverse)

#INPUT
args = commandArgs()

print(args)

#   tuple val(EP_list), path (annotated_ORs)
#    each path(full_GWAS_hg19)

(EP_list = args[8])
annotated_ORs = args[9]
full_GWAS_hg19 = args[10]
EP_annotations = args[11]



#OUTPUT
GO_figure = paste0(EP_list, "_", Sys.Date(), "_GO_figure.png")


(all_genes_with_significant_TS_EPs <- fread(EP_annotations) %>% 
    dplyr::select(ensembl_gene_id) %>% distinct() )


(significant_recessive_annotated_ORs <- fread(annotated_ORs) %>% 
    dplyr::filter(P<0.01)%>% 
    dplyr::select(ensembl_gene_id) )
genes_with_significant_EPs_RECsignificant_for_EP_list<- unique(
    as.character(
        unlist(
            lapply(significant_recessive_annotated_ORs$ensembl_gene_id, as.character) %>% 
            str_split(pattern = ",")))) %>% as_tibble()


(genes_with_significant_EP_eQTLs <- fread(EP_annotations) %>% 
    filter(tidytype=="ESpos_Contact_overlap_eQTLs") %>% 
    dplyr::select(ensembl_gene_id) %>% distinct() )

(genes_with_significant_facet_EPs_with_contact <- 
    fread(EP_annotations) %>%
    filter(tidytype=="neuron_facet_not_eQTLs") %>% 
    dplyr::select(ensembl_gene_id) %>% distinct())

# GO term analysis and plotting
gene_lists <- c("all_genes_with_significant_TS_EPs", "genes_with_significant_EPs_RECsignificant_for_EP_list",
                "genes_with_significant_EP_eQTLs","genes_with_significant_facet_EPs_with_contact")
gene_list_N <- length(gene_lists)
  
Enrich_go_function <- function (gene_list) {
  enrichGO(gene = gene_list,
                OrgDb = "org.Hs.eg.db",
                keyType = "ENSEMBL",
                ont = "ALL",
                qvalueCutoff = 0.05,
                readable = TRUE)
}


i=1
for (gene_list in gene_lists) { 
  print (gene_list ) 
  ego <- Enrich_go_function (unlist(get(gene_list)))
  (nam <- paste("fig", i, sep = ""))
  print(nam)
  assign(nam, 
         dotplot(ego, showCategory=25, title=paste("GO term analysis",gene_list,"genes, N=",nrow(data.frame(get(gene_list)))))
  )
  i<-i+1
  
  }


png(GO_figure, width = 25, height = 15, units = 'in', res = 300)
fig1 + fig2 +fig3+fig4
dev.off()