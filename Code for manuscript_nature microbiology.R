
#
require(ALDEx2)
require(ggplot2)
require(phyloseq)


### load Rel data
genus_r <- read.csv(file = 'metaG-200918/Abun_S/Abun_prop.tsv',
                    header = TRUE, sep = '\t', row.names = 1)
colnames(genus_r) <- gsub(pattern = '_prop', '', colnames(genus_r))
### load Count data
genus <- read.csv(file = 'metaG-200918/Abun_S/Abun_cnt.abun',
                  header = TRUE, sep = '\t', row.names = 1)
colnames(genus) <- gsub(pattern = '_cnt', '', colnames(genus))
genus <- genus[rowSums(genus) > 10, ]

genus_anno <- read.csv(file = 'metadata.tsv', header = TRUE, sep = '\t', stringsAsFactors = FALSE, strip.white = TRUE)
rownames(genus_anno) <- genus_anno$sample

### filtering
isRM_g1 <- rowSums(genus[, genus_anno$sample[genus_anno$group == 'NC']]) < 30
isRM_g2 <- rowSums(genus[, genus_anno$sample[genus_anno$group == 'CDHFD']]) < 30
isRM_g3 <- rowSums(genus[, genus_anno$sample[genus_anno$group == 'CDHFD-I']]) < 30
isRM_g4 <- rowSums(genus[, genus_anno$sample[genus_anno$group == 'CDHFD-C']]) < 30
isRM_g5 <- rowSums(genus[, genus_anno$sample[genus_anno$group == 'CDHFD-I-12C']]) < 30
isRM_g6 <- rowSums(genus[, genus_anno$sample[genus_anno$group == 'CDHFD-I-13C']]) < 30
isRM <- isRM_g1 & isRM_g2 & isRM_g3 & isRM_g4 & isRM_g5 & isRM_g6
genus <- genus[!isRM, ]
genus_mark <- genus[, genus_anno$sample[genus_anno$group %in% c('CDHFD-I-12C', 'CDHFD-I-13C')]]


### beta diversity
require(vegan)
# c('NC', 'CDHFD', 'CDHFD-I', 'CDHFD-C', 'CDHFD-I-12C', 'CDHFD-I-13C')
used <- c('NC', 'CDHFD-C') # c('NC', 'CDHFD', 'CDHFD-I') 
genus_anno.used <- genus_anno[genus_anno$group %in% used, ]
ID <- 'common'

genus_anno.used$group <- factor(genus_anno.used$group, used)
genus_anno.used$col <- col[genus_anno.used$group]
genus_anno.used$pch <- pch[genus_anno.used$group]
genus_used <- genus[, rownames(genus_anno.used)]


#### NMDS
phyl <- phyloseq(otu_table(genus_used, taxa_are_rows = TRUE),
                 sample_data(genus_anno.used))

phyl.ord <- ordinate(physeq = phyl, method = "PCoA", distance = "bray") # DCA, CCA, RDA, NMDS, unifrac
phyl.ord_m <- cbind(phyl.ord$vectors[, 1:2], 
                    genus_anno.used[rownames(phyl.ord$vectors), ])


####
betad <- betadiver(t(genus_used), "z") # Arrhenius
R1Q2 <- adonis2(betad ~ group, genus_anno.used, perm=10000);  R1Q2 
R1Q2 <- data.frame(R1Q2, check.rows = FALSE, check.names = FALSE);  R1Q2
write.csv(R1Q2, file = 'Result/beta/R1Q2.csv')

#### Draw
if (TRUE) {
  cent <- aggregate(cbind(Axis.1, Axis.2) ~ group, data = phyl.ord_m, FUN = mean)
  segs <- merge(phyl.ord_m, setNames(cent, c('group','oAxis.1','oAxis.2')),
                by = 'group', sort = FALSE)
  
  g <- ggplot(phyl.ord_m, aes(x = Axis.1, y = Axis.2, fill = group, col = group)) +
    geom_segment(data = segs,
                 mapping = aes(xend = oAxis.1, yend = oAxis.2),
                 size = 0.1) + # spiders
    geom_point(shape = 21, alpha = 1, col = '#333333', size = 1.5) +
    geom_point(data = cent, size = 3, col = 'black', pch = 3) +
    # geom_path(data = cent, size = 2, col = 'black') +
    labs(x = sprintf('PCoA 1 (%s)', round(phyl.ord$values$Relative_eig[1] * 100, 2)), 
         y = sprintf('PCoA 2 (%s)', round(phyl.ord$values$Relative_eig[2] * 100, 2)))
  g <- g + stat_ellipse(aes(col = group), level = 0.7, size = 1, lty = 'solid')
  # g
  g <- g + theme_bw() + 
    theme(panel.grid = element_line(color = 'white', linetype = 1), 
          legend.position = 'right', #, c(1.0, 1.0), 
          legend.spacing.x = unit(x = 6, units = 'pt'),
          legend.title = element_blank(),
          legend.text = element_text(size = 8),
          legend.key.size = unit(x = 8, units = 'pt'),
          legend.background = element_rect(fill="white", size=.5, linetype="solid", colour = 'black'), 
          axis.title = element_text(size = 18), 
          axis.text = element_text(size = 16, color = 'black'),
          axis.ticks.length = unit(.25, "cm")) 
  g <- g +
    scale_fill_manual(values = col) +
    scale_color_manual(values = col)
  g
  ggsave('Result/beta/beta.PCoA-nc_C.pdf', plot = g, width = 5, height = 3.6)
}


### differential analysis
# c('NC', 'CDHFD', 'CDHFD-I', 'CDHFD-C', 'CDHFD-I-12C', 'CDHFD-I-13C')
grp_list <- list(c('NC', 'CDHFD'), c('NC', 'CDHFD-I'), c('NC', 'CDHFD-C'),
                 c('CDHFD', 'CDHFD-I'), c('CDHFD', 'CDHFD-C'),
                 c('CDHFD-I', 'CDHFD-C'))
dat_list <- list()
for (grp in grp_list) {
  print(grp)
  genus_anno.used <- genus_anno[genus_anno$group %in% grp, ]
  
  genus_anno.used$group <- factor(genus_anno.used$group, grp)
  genus_anno.used$col <- col[genus_anno.used$group]
  genus_anno.used$pch <- pch[genus_anno.used$group]
  genus_used <- genus[, rownames(genus_anno.used)]
  
  
  dat <- aldex(reads = genus_used, conditions = factor(genus_anno[colnames(genus_used), 'group'], grp))
  if (FALSE) {
    
    mm <- model.matrix(~ Genotype + TimePoint, genus_anno.used)
    x <- aldex.clr(reads = genus_used, conds = mm, 
                   mc.samples = 128, denom="all")
    dat <- aldex.glm(x, mm)
  }
  
  dat_list[[sprintf('%s-%s', grp[2], grp[1])]] <- dat
}


genus_prop <- read.csv(file = 'metaG-200918/Abun_S/Abun_prop.tsv', 
                       header = TRUE, sep = '\t', row.names = 1)
colnames(genus_prop) <- gsub(pattern = '_prop', '', colnames(genus_prop))


dat_list_prop <- lapply(names(dat_list), FUN = function(name) {
  dat <- dat_list[[name]]
  dat <- cbind(genus_prop[rownames(dat), ], dat)
  write.csv(x = dat, file = sprintf('Result/diff_with expr.%s.csv', name))
  return(dat)
})
names(dat_list_prop) <- names(dat_list)


##### Heat map
### Microbes
require(factoextra)
require(dendextend)
require(ComplexHeatmap)
require(circlize)
# load(file = 'Result/diff.RData')
load(file = 'Result/Diff/diff_with expr.RData')

Cor <- list()

DEGene <- read.csv(file = 'Meta information for heatmap_2020.10.15.tsv', header = TRUE, sep = '\t', 
                   strip.white = TRUE, stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)
mapping <- read.csv(file = 'Meta information for id.csv', row.names = 1)
DEGene <- DEGene[rownames(mapping)[mapping$MARK != ""], ]

all(rownames(DEGene) %in% rownames(dat_list_prop$`CDHFD-I-CDHFD`))

pheno <- genus_anno[colnames(DEGene), ]
pheno$group <- factor(pheno$group, c('NC', 'CDHFD', 'CDHFD-I', 'CDHFD-C'))

ds <- c('euclidean', 'pearson', 'spearman')
Colv  <- lapply(X = ds, FUN = function(d) {
  colv <- as.matrix(DEGene) %>% 
    t %>% 
    factoextra::get_dist(method = d) %>% 
    hclust(method = 'ward.D2') %>% 
    as.dendrogram %>% rev %>%
    set("branches_lwd", 2) %>% 
    set("leaves_pch", 19)
  colv
})
names(Colv) <- ds

mapping <- data.frame(NAME = gsub('^(.*?)\\|(.*?)\\|(.*?)$', '\\1', rownames(DEGene)), 
                      row.names = rownames(DEGene))

DEGene <- DEGene[!(rownames(DEGene) %in% c('Bacteroides vulgatus|821|S', 'Anaerostipes caccae|105841|S')), ]
Cor$microbe <- DEGene

DEGene <- t(scale(t(DEGene)))

C12 <- dat_list_prop$`CDHFD-I-CDHFD`[rownames(DEGene), genus_anno$sample[genus_anno$group == 'CDHFD-I-12C' ] ]
C13 <- dat_list_prop$`CDHFD-I-CDHFD`[rownames(DEGene), genus_anno$sample[genus_anno$group == 'CDHFD-I-13C' ] ]
DEGene_mark <- data.frame(C12_mean = rowMeans(C12),
                          C13_mean = rowMeans(C13),
                          C12_median = apply(C12, 1, median),
                          C13_median = apply(C13, 1, median),
                          row.names = rownames(DEGene))

col_fun = colorRamp2(c(0, 0.025, 0.3), c("white", '#be29ec', "#660066"))

### draw heatmap
if (TRUE) {
  if (FALSE) {
    Fun_row_anno <- rowAnnotation(
      # C12 = DEGene_mark[rownames(DEGene), 'C12_mean'],
      C13 = DEGene_mark[rownames(DEGene), 'C13_mean'],
      col = list(C12 = col_fun, 
                 C13 = col_fun),
      Lab = anno_text(x = mapping[rownames(DEGene), 'NAME'], # gsub(pattern = '[gp]__', replacement = '', x = row_label),
                      location = 0,
                      rot = 0,
                      just = 'left',
                      gp = gpar(fontsize = 6)),
      # border = TRUE,
      # gap = unit(3, "points"),
      simple_anno_size = unit(0.4, "cm"),
      gp = gpar(col = "grey", width = 0.1),
      annotation_legend_param = list(legend_direction = "vertical",
                                     nrow = 9,
                                     border = 'black',
                                     grid_width = unit(5, "mm")))
    
  } else {
    Fun_row_anno <- NULL
  }
  
  # top_col <- rainbow(nlevels(pheno$group)); names(top_col) <- levels(pheno$group)
  top_col <- c('NC' = '#00ff00', 'CDHFD' = '#fffb96', 'CDHFD-I' = '#d62d20')
  Fun_top_anno <- HeatmapAnnotation(` ` = pheno$group,
                                    col = list(` ` = top_col),
                                    gp = gpar(col = "black"),
                                    row_title = NULL)
  
  RowV <- as.matrix(DEGene) %>% 
    factoextra::get_dist(method = 'spearman') %>% 
    hclust(method = 'ward.D2') %>% rev %>%
    as.dendrogram
  
  ordering_row <- seq(dim(DEGene)[1]); names(ordering_row) <- rownames(DEGene)
  ordering_row <- ordering_row[order.dendrogram(RowV)]
  
  idx_l <- which(grepl('Beta', names(ordering_row))) # symbiosum
  idx <- which(grepl('Gluco', names(ordering_row))) # symbiosum
  
  ordering_row <- c(ordering_row[(idx_l+1):idx], ordering_row[1:idx_l], ordering_row[(idx+1):dim(DEGene)[1]])
  
  
  for (d in ds) {
    pdf(file = 'Result/heatmap_plot_metabolic---.pdf', width = 6, height = 6)
    hp <- Heatmap(matrix = as.matrix(DEGene), 
                  # col = colorRamp2(breaks = c(-2, -0.0, 2), colors = c("#03396c", "#000000", "#e0301e"), space = 'sRGB'), 
                  # col = colorRamp2(breaks = c(-2, 0, 2), colors = c("#1b85b8", "white", "#ae5a41"), space = 'sRGB'),
                  col = colorRamp2(breaks = c(-2.5, -0.5, 2.5), colors = c("#028900", "#000000", "#ff0000"), space = 'sRGB'),
                  # col = colorRamp2(breaks = c(-2, 0, 2), colors = c("#028900", "white", "#a200ff"), space = 'sRGB'),
                  # rect_gp = gpar(col = NA, lwd = 0.1),
                  border = 'black',
                  
                  # cluster_rows = RowV,
                  row_order = ordering_row,
                  cluster_columns = Colv$spearman,
                  # column_order = ordering,
                  show_row_dend = TRUE,
                  
                  show_row_names = FALSE,
                  show_column_names = TRUE,
                  column_names_gp = gpar(cex = 0.6),
                  
                  right_annotation = Fun_row_anno,
                  top_annotation = Fun_top_anno,
                  
                  show_heatmap_legend = TRUE, 
                  heatmap_legend_param = list(legend_direction = "ver", title = ''),
                  column_title = paste('', '', sep = ''))
    draw(hp, show_heatmap_legend = FALSE, show_annotation_legend = FALSE,
         heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
    dev.off()
  }
}


### Metablic
pheno <- genus_anno[colnames(DEGene), ]
pheno$group <- factor(pheno$group, c('NC', 'CDHFD', 'CDHFD-I', 'CDHFD-C'))

MetaB <- read.csv(file = 'metabolic/Stool metabolites for heatmap_2.tsv', header = TRUE, sep = '\t', 
                  strip.white = TRUE, stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)
MetaB <- MetaB[!is.na(MetaB$M1), ]

DEGene <- MetaB[, colnames(MetaB) != "13C %"]
DEGene <- DEGene[, pheno$sample]
mapping <- data.frame(NAME = rownames(DEGene), 
                      row.names = rownames(DEGene))

ds <- c('euclidean', 'pearson', 'spearman')
Colv  <- lapply(X = ds, FUN = function(d) {
  colv <- as.matrix(DEGene) %>% 
    t %>% 
    factoextra::get_dist(method = d) %>% 
    hclust(method = 'ward.D2') %>% 
    as.dendrogram %>% rev %>%
    set("branches_lwd", 2) %>% 
    set("leaves_pch", 19)
  colv
})
names(Colv) <- ds

Cor$metabolic <- DEGene

DEGene <- t(scale(t(DEGene)))
DEGene_mark <- data.frame(C13_mean = MetaB$`13C %`,
                          row.names = rownames(MetaB))

col_fun = colorRamp2(c(0.2, 0.8, 1), c("white", '#be29ec', "#660066"))


### draw heatmap
if (TRUE) {
  if (FALSE) {
    Fun_row_anno <- rowAnnotation(
      # C12 = DEGene_mark[rownames(DEGene), 'C12_mean'],
      C13 = DEGene_mark[rownames(DEGene), 'C13_mean'],
      col = list(C12 = col_fun, 
                 C13 = col_fun),
      Lab = anno_text(x = mapping[rownames(DEGene), 'NAME'], # gsub(pattern = '[gp]__', replacement = '', x = row_label),
                      location = 0,
                      rot = 0,
                      just = 'left',
                      gp = gpar(fontsize = 6)),
      # border = TRUE,
      # gap = unit(3, "points"),
      simple_anno_size = unit(0.4, "cm"),
      gp = gpar(col = "grey", width = 0.1),
      annotation_legend_param = list(legend_direction = "vertical",
                                     nrow = 9,
                                     border = 'black',
                                     grid_width = unit(5, "mm")))
    
  } else {
    Fun_row_anno <- NULL
  }
  
  # top_col <- rainbow(nlevels(pheno$group)); names(top_col) <- levels(pheno$group)
  top_col <- c('NC' = '#00ff00', 'CDHFD' = '#fffb96', 'CDHFD-I' = '#d62d20')
  Fun_top_anno <- HeatmapAnnotation(` ` = pheno$group,
                                    col = list(` ` = top_col),
                                    gp = gpar(col = "black"),
                                    row_title = NULL)
  
  RowV <- as.matrix(DEGene) %>% 
    factoextra::get_dist(method = 'spearman') %>% 
    hclust(method = 'ward.D2') %>% rev %>%
    as.dendrogram
  
  ordering_row <- seq(dim(DEGene)[1]); names(ordering_row) <- rownames(DEGene)
  ordering_row <- ordering_row[order.dendrogram(RowV)]
  
  idx_l <- which(grepl('Beta', names(ordering_row))) # symbiosum
  idx <- which(grepl('Gluco', names(ordering_row))) # symbiosum
  
  ordering_row <- c(ordering_row[(idx_l+1):idx], ordering_row[1:idx_l], ordering_row[(idx+1):dim(DEGene)[1]])
  
  
  for (d in ds) {
    pdf(file = 'Result/heatmap_plot_metabolic---.pdf', width = 6, height = 6)
    hp <- Heatmap(matrix = as.matrix(DEGene), 
                  # col = colorRamp2(breaks = c(-2, -0.0, 2), colors = c("#03396c", "#000000", "#e0301e"), space = 'sRGB'), 
                  # col = colorRamp2(breaks = c(-2, 0, 2), colors = c("#1b85b8", "white", "#ae5a41"), space = 'sRGB'),
                  col = colorRamp2(breaks = c(-2.5, -0.5, 2.5), colors = c("#028900", "#000000", "#ff0000"), space = 'sRGB'),
                  # col = colorRamp2(breaks = c(-2, 0, 2), colors = c("#028900", "white", "#a200ff"), space = 'sRGB'),
                  # rect_gp = gpar(col = NA, lwd = 0.1),
                  border = 'black',
                  
                  # cluster_rows = RowV,
                  row_order = ordering_row,
                  cluster_columns = Colv$spearman,
                  # column_order = ordering,
                  show_row_dend = TRUE,
                  
                  show_row_names = FALSE,
                  show_column_names = TRUE,
                  column_names_gp = gpar(cex = 0.6),
                  
                  right_annotation = Fun_row_anno,
                  top_annotation = Fun_top_anno,
                  
                  show_heatmap_legend = TRUE, 
                  heatmap_legend_param = list(legend_direction = "ver", title = ''),
                  column_title = paste('', '', sep = ''))
    draw(hp, show_heatmap_legend = FALSE, show_annotation_legend = FALSE,
         heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
    dev.off()
  }
}



### correlation
load(file = 'Result/Diff/diff_with expr.RData')
microbe <-  read.table(file = 'Meta information_114 species for correlation analysis.tsv', header = FALSE, sep = '\n', strip.white = TRUE, stringsAsFactors = FALSE)$V1

all(microbe %in% rownames(dat_list_prop$`CDHFD-I-CDHFD`))
DEGene <- dat_list_prop$`CDHFD-I-CDHFD`[microbe, ]


## --- microbe --- ##
microbe <- genus[microbe, pheno$sample]

C12 <- dat_list_prop$`CDHFD-I-CDHFD`[rownames(microbe), genus_anno$sample[genus_anno$group == 'CDHFD-I-12C' ] ]
C13 <- dat_list_prop$`CDHFD-I-CDHFD`[rownames(microbe), genus_anno$sample[genus_anno$group == 'CDHFD-I-13C' ] ]
microbe_mark <- data.frame(C12_mean = rowMeans(C12),
                           C13_mean = rowMeans(C13),
                           C12_median = apply(C12, 1, median),
                           C13_median = apply(C13, 1, median),
                           row.names = rownames(microbe))

microbe <- rbind(other = colSums(genus[, pheno$sample]) - colSums(microbe), microbe)
# write.table(microbe, file = 'mmvec-2nd/microbe.tsv', quote = FALSE, sep = '\t')
microbe <- as.data.frame(t(compositions::clr(t(microbe))))

microbe_mark <- rbind(other = 0, microbe_mark)


## --- metabolic --- ##
# --- for stool
metabolic <- read.csv(file = 'metabolic/stool neg & pos_normalized to QC_unfilterd.tsv', header = TRUE, sep = '\t',
                      row.names = 1, stringsAsFactors = FALSE, strip.white = TRUE)
metabolic <- metabolic[, pheno$sample]
metabolic <- log(metabolic)

if (FALSE) {
  # --- for liver
  metabolic <- read.csv(file = 'metabolic/liver neg & pos_normalized to QC.csv', header = TRUE,
                        row.names = 1, stringsAsFactors = FALSE, strip.white = TRUE)
  metabolic <- metabolic[, pheno$sample]
  metabolic <- log(metabolic)
  
  # --- for serum
  metabolic <- read.csv(file = 'metabolic/serum neg & Pos_normalized to QC.csv', header = TRUE,
                        row.names = 1, stringsAsFactors = FALSE, strip.white = TRUE)
  metabolic <- metabolic[, pheno$sample]
  metabolic <- log(metabolic)
}

#### find differential metabolic
if (TRUE) {
  g1 <- pheno$sample[pheno$group == 'CDHFD-I']
  g2 <- pheno$sample[pheno$group == 'CDHFD']
  g_other <- pheno$sample[pheno$group == 'NC']
  
  mat <- NULL
  apply(X = metabolic[, c(g1, g2)], MARGIN = 1, FUN = function(l) {
    tab <- data.frame(wilcox = wilcox.test(x = l[g1], y = l[g2])$p.value, 
                      ttest = t.test(x = l[g1], y = l[g2])$p.value, 
                      FC = mean(l[g1]) / mean(l[g2]), 
                      FC.original = mean(exp(l[g1])) / mean(exp(l[g2])), 
                      check.names = FALSE)
    mat <<- rbind(mat, tab)
    return(NA)
  })
  rownames(mat) <- rownames(metabolic)
  
  mat$wilcox.BH <- p.adjust(mat$wilcox, "BH")
  mat$ttest.BH <- p.adjust(mat$ttest, "BH")
  
  mat <- cbind(metabolic[rownames(mat), ], 
               `mean.NC` = rowMeans(metabolic[, g_other]), 
               `mean.CDHFD` = rowMeans(metabolic[, g2]), 
               `mean.CDHFD-I` = rowMeans(metabolic[, g1]), 
               `median.NC` = apply(metabolic[, g_other], 1, median), 
               `median.CDHFD` = apply(metabolic[, g2], 1, median), 
               `median.CDHFD-I` = apply(metabolic[, g1], 1, median), 
               mat)
  
  # write.csv(mat, file = 'metabolic/diff-metabolic-serum.csv')
  
  mat <- mat[mat$ttest < 0.05, ]
} 

metabolic <- metabolic[rownames(mat), ]

metabolic_mark <- read.csv(file = 'metabolic/Stool metabolites for heatmap_2.tsv', header = TRUE, sep = '\t', 
                           strip.white = TRUE, stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)
metabolic_mark <- metabolic_mark[rownames(metabolic_mark) %in% rownames(metabolic), ]
metabolic_mark <- data.frame(C13_mean = metabolic_mark$`13C %`,
                             row.names = rownames(metabolic_mark))



Comb <- Comb.p <- as.data.frame(matrix(data = NA, nrow = dim(microbe)[1], ncol = dim(metabolic)[1], 
                                       dimnames = list(rownames(microbe), rownames(metabolic))), 
                                check.rows = FALSE)
idx <- pheno$sample[pheno$group == 'CDHFD-I']


tmp <- sapply(X = rownames(microbe), FUN = function(b) {
  sapply(X = rownames(metabolic), FUN = function(m) {
    # Comb[b, m] <<- cor(as.numeric(microbe[b, idx]), 
    #                    as.numeric(metabolic[m, idx]), 
    #                    method = 'pearson')
    corTest <- cor.test(as.numeric(microbe[b, idx]), 
                        as.numeric(metabolic[m, idx]), 
                        method = 'pearson')
    Comb[b, m] <<- corTest$estimate
    Comb.p[b, m] <<- corTest$p.value
    
  })
  return(NA)
})
idx <- which(apply(X = Comb, MARGIN = 1, FUN = function(b) {
  all(is.na(b))
}))
Comb <- Comb[!(rownames(Comb) %in% c('other', names(idx))), ]


#### Filtered
slected_metabolic <- c("1-Palmitoyl-2-hydroxy-sn-glycero-3-phosphoethanolamine", 
                       "5-Hydroxymethylcytidine", 
                       "Adenosine", 
                       "Deoxyadenosine", 
                       "Deoxythymidine 5'-phosphate (dTMP)", 
                       "Dimethylaminopurine", 
                       "Dodecanoic acid", 
                       "Hydroxyisocaproic acid", 
                       "L-Glutamate", 
                       "Myristic acid", 
                       "N-Acetyl-L-glutamate", 
                       "Pantothenate", 
                       "Pentadecanoic Acid", 
                       "3-(3-Hydroxyphenyl)propanoic acid", 
                       "S-Methyl-5'-thioadenosine", 
                       "sn-Glycerol 3-phosphoethanolamine", 
                       "Undecanoic Acid", 
                       "Xanthine", 
                       "Pyridoxal (Vitamin B6)")
Comb <- Comb[c(#'Anaerostipes sp. BG01|2025494|S', 
  'Bacteroides uniformis|820|S', 
  'Bacteroides acidifaciens|85831|S', 
  'Parabacteroides distasonis|823|S'), 
  slected_metabolic]


#### Draw
col_fun_microbe = colorRamp2(c(0, 0.025, 0.3), c("white", '#be29ec', "#660066"))
related <- read.table(file = 'Meta information selected.tsv', header = FALSE, sep = '\n', 
                      stringsAsFactors = FALSE, strip.white = TRUE)$V1
microbe_mark <- microbe_mark[rownames(Comb), ]
tab <- data.frame(id = related, 
                  idx = match(related, rownames(Comb)), 
                  stringsAsFactors = FALSE)
tab <- tab[!is.na(tab$idx), ]
tab <- cbind(tab, Name = gsub('^(.*?)\\|(.*?)\\|(.*?)$', '\\1', tab$id))
Fun_row_anno <- rowAnnotation(
  # C12 = microbe_mark[rownames(Comb), 'C12_mean'],
  C13 = microbe_mark[rownames(Comb), 'C13_mean'],
  col = list(C12 = col_fun_microbe, 
             C13 = col_fun_microbe),
  # Lab = anno_mark(at = tab$idx, labels = tab$Name, labels_gp = gpar(cex = 0.8)),
  Text = anno_text(gsub('^(.*?)\\|(.*?)\\|(.*?)$', '\\1', rownames(Comb)), gp = gpar(cex = 0.5) ),
  border = TRUE,
  show_annotation_name = FALSE,
  # gap = unit(3, "points"),
  simple_anno_size = unit(0.3, "cm"),
  # gp = gpar(col = "grey", width = 0.1),
  annotation_legend_param = list(legend_direction = "vertical",
                                 nrow = 9,
                                 border = 'black',
                                 grid_width = unit(5, "mm"))
)


col_fun_metabolic = colorRamp2(c(0.2, 0.8, 1), c("white", '#be29ec', "#660066"))
tmp <- data.frame(C13_mean = rep(0, dim(metabolic)[1]), 
                  row.names = rownames(metabolic))
tmp[rownames(metabolic_mark), 'C13_mean'] <- metabolic_mark$C13_mean
metabolic_mark <- tmp
related <- rownames(metabolic_mark)[metabolic_mark$C13_mean > 0]
related <- c(related, "PS(16:0/16:0)")
tab_col <- data.frame(id = related, 
                      idx = match(related, colnames(Comb)), 
                      stringsAsFactors = FALSE)
tab_col <- tab_col[!is.na(tab_col$idx), ]
Fun_col_anno <- HeatmapAnnotation(
  which = 'column',
  C13 = metabolic_mark[colnames(Comb), 'C13_mean'],
  col = list(C13 = col_fun_metabolic),
  # Lab = anno_mark(at = tab_col$idx, labels = tab_col$id, labels_gp = gpar(cex = 0.8), 
  #                 which = 'column', side = 'bottom', labels_rot = 60),
  Text = anno_text(colnames(Comb), rot = 60, gp = gpar(cex = 0.5)),
  border = TRUE,
  show_annotation_name = FALSE,
  # gap = unit(3, "points"),
  simple_anno_size = unit(0.3, "cm"),
  # gp = gpar(col = "grey", width = 0.1),
  annotation_legend_param = list(legend_direction = "vertical",
                                 nrow = 9,
                                 border = 'black',
                                 grid_width = unit(5, "mm"))
)

colv <- as.matrix(Comb) %>% 
  t %>% 
  factoextra::get_dist(method = "euclidean") %>% 
  hclust(method = 'complete') %>% 
  as.dendrogram #%>% rev 
plot(colv)

ordering_col <- seq(dim(Comb)[2]); names(ordering_col) <- colnames(Comb)
ordering_col <- ordering_col[order.dendrogram(colv)]
if (FALSE) {
  write.csv(data.frame(taxa = names(ordering_col), 
                       idx = ordering_col),
            file = 'Cor/cor-selected.col.csv', 
            row.names = FALSE, quote = FALSE)
  
  # after process
  tmp <- read.csv(file = 'Cor/cor-selected.col.csv', header = TRUE, stringsAsFactors = FALSE)
  ordering_col <- tmp$idx; names(ordering_col) <- tmp$taxa
  # ordering_col <- rev(ordering_col)
  save(ordering_col, file = 'Cor/cor-selected.col..Rds')
  
  # load me
  load(file = 'Cor/cor-selected.col..Rds')
}

backup <- Comb
# Comb <- backup

Comb <- (scale((Comb)))

hp <- Heatmap(matrix = as.matrix(Comb), 
              # col = colorRamp2(breaks = c(-1, -0.0, 1), colors = c("#03396c", "#000000", "#e0301e"), space = 'sRGB'),
              # col = colorRamp2(breaks = c(-1, 0, 1), colors = c("#1b85b8", "white", "#ae5a41"), space = 'sRGB'),
              col = colorRamp2(breaks = c(-1, 0, 1), colors = c("#0000ff", "#000000", "#ff0000"), space = 'sRGB'),
              # col = colorRamp2(breaks = c(-1, 0, 1), colors = c("#028900", "white", "#a200ff"), space = 'sRGB'),
              # rect_gp = gpar(col = NA, lwd = 0.1),
              border = 'black', rect_gp = gpar(color = '#555555'),
              
              cluster_rows = TRUE,
              # row_order = ordering_row,
              # cluster_columns = colv,
              column_order = ordering_col,
              show_row_dend = TRUE,
              show_column_dend = TRUE,
              
              show_row_names = FALSE,
              show_column_names = FALSE,
              column_names_gp = gpar(cex = 0.6),
              
              right_annotation = Fun_row_anno,
              bottom_annotation = Fun_col_anno,
              
              show_heatmap_legend = TRUE, 
              heatmap_legend_param = list(legend_direction = "ver", title = ''),
              column_title = paste('', '', sep = ''))
pdf('Cor/cor-selected-.pdf', width = 5.49, height = 3.06)
draw(hp, show_heatmap_legend = FALSE, show_annotation_legend = FALSE,
     heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()





















#DNA-SIP
library("phyloseq")
library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library("dplyr")        # filter and reformat data frames
library("tibble")
library(HTSSIP)
library(RColorBrewer)

setwd("E:/9.Q-SIP/5.Re-analysis_2021.03.23_2")

tax_mat<- read_excel("1.Taxonomy.xlsx", sheet = "Sheet1")
samples_df<- read_excel("1.Sample.xlsx", sheet = "Sheet1")
otu_mat<- read_excel("1.OTU.xlsx", sheet = "Sheet1")

otu_mat <- otu_mat %>%
  tibble::column_to_rownames("OTU")
tax_mat <- tax_mat %>% 
  tibble::column_to_rownames("OTU")
samples_df <- samples_df %>% 
  tibble::column_to_rownames("Sample")

otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)

physeq_rep5 <- phyloseq(OTU, TAX, samples)

xls_data = read_excel('Summary_qPCR_4.xlsx',sheet = 'Sheet1')

qPCR_df=data.frame(
  IS_CONTROL = xls_data$IS_CONTROL,
  Sample = xls_data$Sample,
  Buoyant_density = xls_data$Buoyant_density,
  qPCR_tech_rep_mean = xls_data$qPCR_tech_rep_mean,
  qPCR_tech_rep_sd= xls_data$qPCR_tech_rep_sd,
  Gradient= xls_data$Gradient,
  Fraction= xls_data$Fraction,
  Treatment= xls_data$Treatment,
  Replicate=xls_data$Replicate)

qPCR <- list(list(),qPCR_df)
names(qPCR) <- c('raw','summary')
physeq_rep3_t = OTU_qPCR_trans(physeq_rep5, qPCR)

data.frame(sample_data(physeq_rep5)) %>%
  dplyr::select(Treatment, Replicate) %>%
  distinct
atomX = qSIP_atom_excess(physeq_rep3_t,
                         control_expr='Treatment=="12C_Con"',
                         treatment_rep='Replicate')
atomX %>% names

df_atomX_boot = qSIP_bootstrap(atomX, n_boot=100)
df_atomX_boot %>% head

CI_threshold = 0
df_atomX_boot = df_atomX_boot %>%
  mutate(Incorporator = A_CI_low > CI_threshold,
         OTU = reorder(OTU, -A))

n_incorp = df_atomX_boot %>%
  filter(Incorporator == TRUE) %>%
  nrow
cat('Number of incorporators:', n_incorp, '\n')

ggplot(df_atomX_boot, aes(OTU, A, ymin=A_CI_low, ymax=A_CI_high, color=Incorporator)) +
  geom_pointrange(size=0.25) +
  geom_linerange() +
  geom_hline(yintercept=0, linetype='dashed', alpha=0.5) +
  labs(x='OTU', y='Atom fraction excess') +
  theme_bw() +
  theme(
    axis.text.x = element_blank()
  )

write.table(df_atomX_boot, "E:/9.Q-SIP/5.Re-analysis_2021.03.23_2.csv", sep="\t")

write.table(df_atomX_boot, "E:/9.Q-SIP/5.Re-analysis_2021.03.23_2.txt", sep="\t")


xls_data<- read_excel("A fraction for plot.xlsx", sheet = "Sheet1")

picked_colors=c("red", "#FFD700", "#ADFF2F", "#05c77d", "#FF00FF", "#0000FF", "#00e5ff", "#bbbbff", "#FF1493", "#FFFF00", "#FFC0CB", "#ff6200", "#40E0D0", "#09ff00", "#1E90FF", "#bf00ff", "#ff9100", "#A0522D", "#8A2BE2")

#picked_colors=sample(picked_colors,size=19)
> #getPalette = colorRampPalette(brewer.pal(9, "Set1"))
  > #spec_colors=getPalette(19)
  
  ggplot(xls_data, aes( A, Rank,xmin=A_CI_low, xmax=A_CI_high, color=Phylum)) +
  geom_linerange() +
  geom_vline(xintercept=0, linetype='dashed', size=0.6, colour="#4f4a4a", alpha=1)  +
  labs(x='Excess atom fraction 13C', y='') +
  theme_bw() +
  geom_point(size = 5,shape=21,colour="#292828",stroke=1,aes(fill=Phylum)) +
  scale_color_manual(values =picked_colors ) +
  scale_fill_manual(values =picked_colors) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(face="bold", colour="#333232", size=20), axis.text.x = element_text(angle=0, vjust=0.5, size=16,colour="#000000"),
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank(),
  )+
  scale_x_continuous(breaks=seq(-0.8, 1, 0.2))
# device = "target format" , one of "eps", "ps", "tex" (pictex), "pdf", "jpeg", "tiff", "png", "bmp", "svg" or "wmf" (windows only).
# first parameter is filename, it could be "anything you want", but it won't decide the file format, the "device" parameter will. 
ggsave("plot.pdf", width = 6.18, height = 10, device="pdf") 






