library(VGAM)

# for gene: get counts
# fit VGAM:glm with and without celltype 
# LRT test

plot_colours <- function(colours){
  # mod from https://stackoverflow.com/questions/51867716/plot-colors-with-hex-values-in-r
  # take colors, name the column, then cbind
  
  color_df = data.frame(
    ColorName = names(colours),
    ColorHex = colours
  )
  
  # plot it, using scale_color_identity to interpret the hex as a color
  ggplot(color_df, aes(x = ColorName, y = 1, fill = ColorHex)) +
    geom_tile() +
    geom_text(aes(label = ColorName)) +
    scale_fill_identity() + 
    coord_flip() + theme_void()
}

#colour_vector <- c("SC" = "#7E2678",
#                   "RS" = "#6167AF",
#                   "ES" = "grey50")

colour_vector <- c("SG" = "#046735",
                   "SC" = "#D31281",
                   "RS" = "#282B69",
                   "ES" = "grey60",
                   "Sertoli" = "#643F18",
                   "Leydig" = "#985D25",
                   "Immune" = "#FFDFB2")

plot_colours(colour_vector)

annotate_chromosome <- function(sce){
  library(biomaRt)
  ensembl <- useMart("ensembl")
  ensembl <- useDataset(dataset = "mmusculus_gene_ensembl", 
                        mart = ensembl)
  genenames <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name"),
                     mart = ensembl)
  rownames(genenames) <- make.unique(genenames$external_gene_name)
  genenames <- genenames[!duplicated(genenames),]
  rowData(sce) <- genenames[rownames(sce),]
  return(sce)
}

annotate_chromosome_list <- function(genes){
  library(biomaRt)
  ensembl <- useMart("ensembl")
  ensembl <- useDataset(dataset = "mmusculus_gene_ensembl", 
                        mart = ensembl)
  genenames <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name"),
                     mart = ensembl)
  rownames(genenames) <- make.unique(genenames$external_gene_name)
  genenames <- genenames[!duplicated(genenames),]
  return(genenames[genes,])
}

counts_reference <- function(sce){
  assays(sce)[["counts_reference"]]
}

counts_alternative <- function(sce){
  assays(sce)[["counts_alternative"]]
}

theme_paper <- function(){
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_line(color = "grey"),
        plot.background=element_blank(), 
        theme(text = element_text(size=15)))
}

theme_tsne <- function(){
  
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_line(color = "grey"),
        plot.background=element_blank())
}

setwd("~/Desktop/Projects/ASE_Spermatogenesis_Paper/")
sce <- readRDS("./Data/processed/sce_merged_new_sparse.rds")
data_f1 <- sce[,sce$Library %in% c("Sample5", "Sample6")]
sce <- data_f1

model.selection.results.sc <- readRDS("~/Desktop/Projects/ASE_Spermatogenesis_Paper/Data/processed/model_selection_sc.rds")
model.selection.results.rs <- readRDS("~/Desktop/Projects/ASE_Spermatogenesis_Paper/Data/processed/model_selection_rs.rds")
model.selection.results.es <- readRDS("~/Desktop/Projects/ASE_Spermatogenesis_Paper/Data/processed/model_selection_es.rds")

celltype1 = "SC"
celltype2 = "RS"

genes_check <- intersect(model.selection.results.sc$Gene, model.selection.results.rs$Gene)

ds_ref_ct1 <- counts_reference(sce[genes_check, data_f1$CellType == celltype1])
ds_alt_ct1 <- counts_alternative(sce[genes_check, data_f1$CellType == celltype1])
ds_ref_ct2 <- counts_reference(sce[genes_check, data_f1$CellType == celltype2])
ds_alt_ct2 <- counts_alternative(sce[genes_check, data_f1$CellType == celltype2])
library_ct1 <- data_f1[,data_f1$CellType == celltype1]$Library
library_ct2 <- data_f1[,data_f1$CellType == celltype2]$Library

# pseudobulk - split our replicates into technical replicates across cells

data_test_ref <- data.frame(
  CT1_REP1 = rowSums(counts_reference(sce[genes_check, data_f1$CellType == "SC" & data_f1$Library == "Sample5"][,1:150])), 
  CT1_REP2 = rowSums(counts_reference(sce[genes_check, data_f1$CellType == "SC" & data_f1$Library == "Sample5"][,150:300])), 
  CT1_REP3 = rowSums(counts_reference(sce[genes_check, data_f1$CellType == "SC" & data_f1$Library == "Sample5"][,300:450])), 
  CT1_REP4 = rowSums(counts_reference(sce[genes_check, data_f1$CellType == "SC" & data_f1$Library == "Sample6"][,1:150])), 
  CT1_REP5 = rowSums(counts_reference(sce[genes_check, data_f1$CellType == "SC" & data_f1$Library == "Sample6"][,150:300])), 
  CT1_REP6 = rowSums(counts_reference(sce[genes_check, data_f1$CellType == "SC" & data_f1$Library == "Sample6"][,300:450])), 
  CT2_REP1 = rowSums(counts_reference(sce[genes_check, data_f1$CellType == "RS" & data_f1$Library == "Sample5"][,1:300])), 
  CT2_REP2 = rowSums(counts_reference(sce[genes_check, data_f1$CellType == "RS" & data_f1$Library == "Sample5"][,300:600])), 
  CT2_REP3 = rowSums(counts_reference(sce[genes_check, data_f1$CellType == "RS" & data_f1$Library == "Sample5"][,600:900])), 
  CT2_REP4 = rowSums(counts_reference(sce[genes_check, data_f1$CellType == "RS" & data_f1$Library == "Sample6"][,1:300])), 
  CT2_REP5 = rowSums(counts_reference(sce[genes_check, data_f1$CellType == "RS" & data_f1$Library == "Sample6"][,300:600])), 
  CT2_REP6 = rowSums(counts_reference(sce[genes_check, data_f1$CellType == "RS" & data_f1$Library == "Sample6"][,600:889]))
)

data_test_alt <- data.frame(
  CT1_REP1 = rowSums(counts_alternative(sce[genes_check, data_f1$CellType == "SC" & data_f1$Library == "Sample5"][,1:150])), 
  CT1_REP2 = rowSums(counts_alternative(sce[genes_check, data_f1$CellType == "SC" & data_f1$Library == "Sample5"][,150:300])), 
  CT1_REP3 = rowSums(counts_alternative(sce[genes_check, data_f1$CellType == "SC" & data_f1$Library == "Sample5"][,300:450])), 
  CT1_REP4 = rowSums(counts_alternative(sce[genes_check, data_f1$CellType == "SC" & data_f1$Library == "Sample6"][,1:150])), 
  CT1_REP5 = rowSums(counts_alternative(sce[genes_check, data_f1$CellType == "SC" & data_f1$Library == "Sample6"][,150:300])), 
  CT1_REP6 = rowSums(counts_alternative(sce[genes_check, data_f1$CellType == "SC" & data_f1$Library == "Sample6"][,300:450])), 
  CT2_REP1 = rowSums(counts_alternative(sce[genes_check, data_f1$CellType == "RS" & data_f1$Library == "Sample5"][,1:300])), 
  CT2_REP2 = rowSums(counts_alternative(sce[genes_check, data_f1$CellType == "RS" & data_f1$Library == "Sample5"][,300:600])), 
  CT2_REP3 = rowSums(counts_alternative(sce[genes_check, data_f1$CellType == "RS" & data_f1$Library == "Sample5"][,600:900])), 
  CT2_REP4 = rowSums(counts_alternative(sce[genes_check, data_f1$CellType == "RS" & data_f1$Library == "Sample6"][,1:300])), 
  CT2_REP5 = rowSums(counts_alternative(sce[genes_check, data_f1$CellType == "RS" & data_f1$Library == "Sample6"][,300:600])), 
  CT2_REP6 = rowSums(counts_alternative(sce[genes_check, data_f1$CellType == "RS" & data_f1$Library == "Sample6"][,600:889]))
)

# add metadata

metadata_bulks <- data.frame(
  CT = c(rep("SC", 6), rep("RS", 6)), 
  LIB = rep(c(rep("LIB1", 3), rep("LIB2", 3)), 2), 
  SAMP = rep(paste0("S", 1:6), 2)
)

# testing
# estimate mean and overdispersions
means_emp <-  rowMeans(as.matrix(data_test_ref / (data_test_ref + data_test_alt)), na.rm = T)
variances_emp <- rowVars(as.matrix(data_test_ref / (data_test_ref + data_test_alt)), na.rm = T)

ggplot(data.frame(x = means_emp, y = variances_emp), aes(abs(x - 0.5), y)) + geom_point(size = 0.1)


#
library(VGAM)

# Example 1

all_coefs <- data.frame(do.call("rbind", lapply(1:nrow(data_test_ref), function(i){
  tryCatch({
    y = as.numeric(data_test_ref[i, ])
    N = as.numeric(data_test_ref[i, ] + data_test_alt[i, ])
    fit <- vglm(cbind(y, N-y) ~ 1, betabinomialff, trace = TRUE)
    coef(fit)
  },error=function(cond) {
    return(c(NA, NA))
  })
})))

colnames(all_coefs) <- c("alpha", "beta")
all_coefs_save <- all_coefs
all_coefs$ks = rowSums(data_test_ref)
all_coefs$Ns = rowSums(data_test_alt) + rowSums(data_test_alt)

rownames(all_coefs) <- rownames(data_test_alt)

# how many fits failed?
table(is.na(all_coefs$alpha))

all_coefs$mean = all_coefs$alpha / ( all_coefs$alpha + all_coefs$beta )
all_coefs$rho = 1 / ( all_coefs$alpha + all_coefs$beta + all_coefs$beta )

all_coefs <- all_coefs[all_coefs$mean < 1 & all_coefs$mean > 0, ]

ggplot(all_coefs, aes(Ns, ks)) + geom_point() + scale_x_log10() + scale_y_log10()
ggplot(all_coefs, aes(Ns, mean)) + geom_point() + scale_x_log10()
ggplot(all_coefs, aes(ks, mean)) + geom_point() + scale_x_log10()
ggplot(all_coefs, aes(Ns, rho)) + geom_point() + scale_x_log10() + scale_y_log10()
ggplot(all_coefs, aes(abs(mean - 0.5), rho)) + geom_point() + scale_y_log10() + scale_x_log10() + geom_smooth()

# 

       