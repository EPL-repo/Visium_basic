library(Seurat)
library(dplyr)
library(tidyverse)
library(purrr)
library(janitor)
library(magrittr)
library(patchwork)
library(shiny)
library(stringr)
library(R.utils)
library(Giotto)
library(ggplot2)
library(viridis)

#-------- Setup

source("helperFunctions.R")

# Data path to h5 output object
data_path <- "path/to/SpaceRangerOutput/SlideName/"

results_folder <- "figures"

# Load scRNAseq reference
scRNA <- readRDS(file = "NSCLC_reference/refquery_final.rds")

#-------- Integration in SEURAT

# Load Visium sample using preProcessVisiumSeurat() from helperFunctions
# This loader will create a Seurat object from a h5, calculate mitochondrial gene % per 
# spot, log-normalize the data, find variable features and perform PCA on the re-scaled data
obj.Seurat <- preProcessSeuratVisium(data_path, normalization = "LogNormalize")

# Integrate the pre-processed Visium data with the scRNA-seq referencedata
# Anchor-gene approach from Seurat compacted into  anchorMapping() from helperFunctions
features <- SelectIntegrationFeatures(list(scRNA,obj.Seurat))
obj.Seurat <- anchorMapping(scRNA, 
						    obj.Seurat, 
						    feats = features, 
						    query.dims=30, 
						    anchor.labels = levels(as.factor(scRNA$Cell_Cluster_level2)))

# Plot the predicted signatures on the Visium sample on top of the H&E image of the tissue
scplots <- purrr::map(levels(as.factor(scRNA$Cell_Cluster_level2)), function(x) 
					  SpatialFeaturePlot(obj.Seurat, x) +
		   theme(legend.key.size = unit(10, "mm"),
		   legend.text = element_text(size = 15),
		   legend.title = element_text(size = 20)))
patchwork::wrap_plots(scplots, ncol=4) %T>% 
	ggsave(filename = paste0(results_folder,"Ref_annot_NSCLC_lvl2_Seurat.pdf"), 
		   width = 25, height = 25, 
		   units = "in", dpi = 300)
			   
# Shiny for visualisation
# Code repurposed from Seurat to produce an interactive Shiny app
# Modified from Ting Lab, Angela Ting's laboratory at Cleveland Clinic
# https://github.com/tingalab
runApp('shiny-app')

#-------- Integration in GIOTTO

# Note 1: Giotto seems to be more effective than other methods at detecting sparse cell 
# types, such as immune cells
# Note 2: Giotto objects not supported for Seurat conversion yet!!!

# Configure workspace with Giotto

instrs <- createGiottoInstructions(save_dir = results_folder,
								   save_plot = TRUE,
								   show_plot = FALSE)

# Create Giotto object
# Note: adjust x-y coordinates to align image (iterative approach - unless image scaling
# factor is known)				   
obj <- createGiottoVisiumObject(visium_dir = data_path, 
							    expr_data = 'filter',
							    h5_visium_path = paste0(data_path,'filtered_feature_bc_matrix.h5'),
							    h5_tissue_positions_path = paste0(data_path,'spatial/tissue_positions.csv'),
							    h5_image_png_path = paste0(data_path,'spatial/tissue_lowres_image.png'),
							    gene_column_index = 2, 
							    instructions = instrs, 
							    xmax_adj = 1090, 
							    ymin_adj = 760, 
							    ymax_adj = 680, 
							    xmin_adj = 290)

# Pre-process in the Giotto object format 
# - Log-normalize (same scaling factor as with Seurat)
# - Calculate highly variable genes
# - Filter the main count matrix for genes that are all highly variable, expressed in 
#   more than 3% of cells, and have a mean expression > 0.4 in cells expressing the 
#   gene in question
# - Adjust expression matrix for technical covariates (%mito/nr genes)
# - PCA is then run on those genes
# - Run UMAP on PCA space 
# - Leiden clustering
# - Giotto creates a spatial network
# - Plot the PCA results and the spatial network formed by the Visium spots
obj.Giotto <- preProcessGiotto(obj, "ZZ269007")

# Create signature matrix using the scRNA-seq data (rank method)
sc_sign_matrix <- makeSignMatrixRank(sc_matrix = as.matrix(scRNA@assays$RNA@data),
									 sc_cluster_ids = scRNA$Cell_Cluster_level2,
									 ties_method = c("random"),
									 gobject = NULL)

obj.Giotto <- runSpatialEnrich(obj.Giotto,
							   enrich_method = c("rank"),
							   sign_matrix = sc_sign_matrix,
							   expression_values = c("normalized"))

# Plot spatial enrichment signatures
scplots <- purrr::map(levels(as.factor(scRNA$Cell_Cluster_level2)), 
					  function(i) spatPlot(gobject = obj.Giotto, 
					  cell_color=unlist(c(obj.Giotto@spatial_enrichment$rank[,..i])), 
					  point_size = 2) +
	theme(title = element_text(size=18),
	legend.text = element_text(size = 15),
	legend.title = element_text(size = 15),
	axis.title = element_text(size=15),
	axis.text = element_text(size=15)) +
	ggtitle(i) + scale_fill_distiller(palette = “Spectral”))
patchwork::wrap_plots(scplots, ncol=4) %T>% ggsave(filename = paste0(results_folder,
					  "Ref_annot_NSCLC_lvl2_Giotto.pdf", width = 25, height = 20, 
					  units = "in", dpi = 300)






