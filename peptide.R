library(ggplot2)
library(clusterProfiler)
library(org.Ce.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(celegans.db)
library(enrichplot)
library(DOSE)


df = read.csv('EV0-EVP_seg_0.05.csv',header = TRUE, sep = '\t')
df <- df[order(-df$log2FoldChange_EVPS_EV0S),]


fc <- df$log2FoldChange_EVPS_EV0S

names(fc) <- df$Gene_stable_ID

#Biological processes
gse_bp<- gseGO(fc, ont = 'BP',
               keyType = 'WORMBASE',
               OrgDb = 'org.Ce.eg.db',
               pvalueCutoff = 0.05)

results_bp <- as.data.frame(gse_bp)


pm_bp <- pairwise_termsim(gse_bp)



p2 <- plot(dotplot(gse_bp,showCategory = 15)+ggtitle('Biological Processes'))

png('bp_dotplot.png', res = 250, width = 1500, height = 2000)
print(p2)
dev.off()

p1 <- plot(emapplot(pm_bp, showCategory = 15, layout = 'nicely')+ggtitle('Biological Processes'))
png('bp_map.png', res = 250, width = 2000, height = 2000)
print(p1)
dev.off()


p3 <- plot(treeplot(pm_bp))
png('bp_tree.png', res = 200, width = 2800, height = 1500)
print(p3)
dev.off()

p2.1 <- goplot(gse_bp, showCategory = 5)
png('bp_go.png', res = 200, width = 2800, height = 1500)
print(p2.1)
dev.off()

p2.2 <- plot(dotplot(gse_bp,showCategory =20, split = '.sign')+facet_grid(.~.sign)+ggtitle('Biological Processes'))
png('bp_dotplot_ad.png', res = 200, width = 1500, height = 2500)
print(p2.2)
dev.off()

#Cellular component
gse_cc<- gseGO(fc, ont = 'CC',
               keyType = 'WORMBASE',
               OrgDb = 'org.Ce.eg.db',
               pvalueCutoff = 0.05)

results_cc <- as.data.frame(gse_cc)


pm_cc <- pairwise_termsim(gse_cc)



p4 <- plot(dotplot(gse_cc,showCategory = 15)+ggtitle('Cellular Component'))

png('cc_dotplot.png', res = 250, width = 1500, height = 2000)
print(p4)
dev.off()

p5 <- plot(emapplot(pm_cc, showCategory = 15, layout = 'nicely')+ggtitle('Cellular Component'))
png('cc_map.png', res = 250, width = 2000, height = 2000)
print(p5)
dev.off()


p6 <- plot(treeplot(pm_cc))
png('cc_tree.png', res = 200, width = 2800, height = 1500)
print(p6)
dev.off()

#Molecular function

gse_mf<- gseGO(fc, ont = 'MF',
               keyType = 'WORMBASE',
               OrgDb = 'org.Ce.eg.db',
               pvalueCutoff = 0.05)

results_mf <- as.data.frame(gse_mf)


pm_mf <- pairwise_termsim(gse_mf)



p7 <- plot(dotplot(gse_mf,showCategory = 15)+ggtitle('Molecular Function'))

png('mf_dotplot.png', res = 250, width = 1500, height = 2000)
print(p7)
dev.off()

p8 <- plot(emapplot(pm_mf, showCategory = 15, layout = 'nicely')+ggtitle('Molecular Function'))
png('mf_map.png', res = 250, width = 2000, height = 2000)
print(p8)
dev.off()


p9 <- plot(treeplot(pm_mf))
png('mf_tree.png', res = 200, width = 2800, height = 1500)
print(p9)
dev.off()


#KEGG

entrez_list <- df$log2FoldChange_EVPS_EV0S
names(entrez_list) <- df$NCBI_gene_.formerly_Entrezgene._ID

gse_kegg<- gseKEGG(entrez_list, organism = 'cel',
                   keyType = 'ncbi-geneid',
                   pvalueCutoff = 0.05)


p10 <-plot(dotplot(gse_kegg)+ggtitle('KEGG'))

png('kegg_dotplot.png', res = 150, width = 1000, height = 1000)
print(p10)
dev.off()



enr_kegg <- enrichKEGG(as.list(df$NCBI_gene_.formerly_Entrezgene._ID), organism = 'cel',
                       keyType = 'ncbi-geneid',
                       pvalueCutoff = 0.05)


p10.1 <- plot(dotplot(enr_kegg))
png('kegg_dotplot_1.1.png', res = 150, width = 1000, height = 1000)
print(p10.1)
dev.off()

