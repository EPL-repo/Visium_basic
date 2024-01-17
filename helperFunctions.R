library(Seurat)
library(dplyr)
library(tidyverse)
library(janitor)
library(magrittr)
library(purrr)
library(patchwork)

#------ Seurat

preProcessSeuratVisium<- function(file, normalization = "SCT"){
    if(normalization == "SCT"){
      Load10X_Spatial(file) %>% 
        PercentageFeatureSet(pattern="^MT-", col.name="percent.mt") %>%
        SCTransform(assay="Spatial", vars.to.regress = "percent.mt", 
        			return.only.var.genes = FALSE) %>%
        RunPCA() %>%
        return()
    }else if(normalization == "LogNormalize"){
      Load10X_Spatial(file) %>% 
        PercentageFeatureSet(pattern="^MT-", col.name="percent.mt") %>%
        NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
        FindVariableFeatures(nfeatures = 3000) %>%
        ScaleData(vars.to.regress = "percent.mt") %>%
        RunPCA() %>%
        return()
    } else {
      warning("Please pick either SCT or LogNormalize for normalization.")
    }
}

UMAP <- function(sample, n_dims=15, res=0.5){
  sample <- FindNeighbors(sample,reduction="pca", dims=1:n_dims)
  sample <- FindClusters(sample, resolution=res)
  sample <- RunUMAP(sample,reduction="pca", dims=1:n_dims)
  return(sample)
}

anchorMapping <- function(reference, query, feats, query.dims=15, 
						  normalization.method = "LogNormalize", anchor.labels,
						  save.loc=FALSE){
  DefaultAssay(query) <- "Spatial"
  anchors = FindTransferAnchors(reference, 
  								query = query,  
  								normalization.method = normalization.method, 
  								features = feats)
  predictions.assay <- TransferData(anchorset = anchors, refdata = reference$Cell_Cluster_level2, 
  									prediction.assay=TRUE, weight.reduction=query[["pca"]],
  									dims=1:query.dims)
  non.mapping <- c()
  for(i in 1:dim(predictions.assay)[1]){ if(sum(predictions.assay@data[i,])==0) non.mapping <- c(non.mapping, rownames(predictions.assay[i]))}
  #print("Non-integrated features: ", non.mapping)
  predictions.assay@misc$non.mapping <- non.mapping
  predictions.assay@misc$mapping <- setdiff(anchor.labels, non.mapping)
  query[["predictions"]] <- predictions.assay
  DefaultAssay(query) <- 'predictions'
  return(query)
}

#------ Giotto

preProcessGiotto <- function(gobject, name){
  metadata <- pDataDT(gobject)
  in_tissue_barcodes <- metadata[in_tissue == 1]$cell_ID
  spatPlot(gobject = gobject, cell_color = 'in_tissue', point_size = 2,
    	   cell_color_code = c('0' = 'lightgrey', '1' = 'blue'),
           save_param = list(save_name = paste0(name,'_1_spatplot_image_inTissue')))
  ## subset on spots that were covered by tissue
  gobject <- subsetGiotto(gobject, cell_ids = in_tissue_barcodes)
  ## define filters
  filterCombinations(gobject = gobject,
                     expression_thresholds = c(1, 2),
                     gene_det_in_min_cells = c(10, 20, 50),
                     min_det_genes_per_cell = c(200, 500, 1000),
                     save_param = list(save_name = paste0(name,'_2_filterCombinations')))
  ## normalize
  gobject <- normalizeGiotto(gobject, scalefactor = 6000, verbose = T)
  ## add gene & cell statistics
  gobject <- addStatistics(gobject)
  ## visualize
  spatPlot2D(gobject = gobject, show_image = T, point_alpha = 0.7,
             cell_color = 'nr_feats', color_as_factor = F,
             save_param = list(save_name = paste0(name,'_3_nr_genes')))
  ## batch effects and technical covariates
  cell_metadata    <- pDataDT(gobject)
  feature_metadata <- fDataDT(gobject)
  mitochondrial_genes <- grep('MT-', x = feature_metadata$feat_ID, value = TRUE)
  gobject <- addFeatsPerc(gobject,
                          expression_values = 'normalized',
                          feats = mitochondrial_genes,
                          vector_name = "mito")
  ## visualize number of genes and mitochondrial content per spot
  spatPlot2D(gobject = gobject,
             show_image = F,
             point_alpha = 0.7,
             cell_color = 'nr_feats', color_as_factor = F,
             coord_fix_ratio = 1,
             save_param = list(save_name = paste0(name,'_4_spatPlot2D_nr_feats'), 
             				   base_width = 6,base_height = 6, save_format = 'pdf'))
  spatPlot2D(gobject = gobject,
             show_image = F,
             point_alpha = 0.7,
             cell_color = 'mito', color_as_factor = F,
             coord_fix_ratio = 1,
             save_param = list(save_name = paste0(name,'_5_spatPlot2D_mito'),
             				   base_width = 6, base_height = 6, save_format = 'pdf'))
  ## adjust expression matrix
  gobject <- adjustGiottoMatrix(gobject = gobject,
                                covariate_columns = c('nr_feats', 'mito'),
                                update_slot = 'custom')
  ## dimension reduction
  ## highly variable features (HVF)
  gobject <- calculateHVF(gobject, save_param = list(save_name = paste0(name, '_6_HVFplot')))
  gene_metadata <- fDataDT(gobject)
  featgenes <- gene_metadata[hvf == 'yes' & perc_cells > 3 & mean_expr_det > 0.4]$feat_ID
  gobject <- runPCA(gobject, 
            		genes_to_use = featgenes, 
        		    scale_unit = F, center = T, 
          		    method="factominer")
  plotPCA(gobject = gobject,
          save_param = list(save_name = paste0(name,'_7_PCA_reduction')))
  ## run UMAP on PCA space (default)
  gobject <- runUMAP(gobject, dimensions_to_use = 1:10)
  plotUMAP(gobject = gobject,
           save_param = list(save_name = paste0(name,'_8_UMAP_reduction')))
  gobject <- createNearestNetwork(gobject = gobject, dimensions_to_use = 1:10, k = 15)
  ## Leiden clustering
  gobject <- doLeidenCluster(gobject = gobject, resolution = 0.4, n_iterations = 1000)
  plotUMAP(gobject = gobject,
           cell_color = 'leiden_clus', show_NN_network = T, point_size = 2.5,
           save_param = list(save_name = paste0(name,'_9_UMAP_leiden')))
  spatDimPlot(gobject = gobject, cell_color = 'leiden_clus',
              dim_point_size = 2, spat_point_size = 2.5,
              save_param = list(save_name = paste0(name,'_10_coVis_leiden')))
  spatDimPlot(gobject = gobject, cell_color = 'nr_feats', color_as_factor = F,
              dim_point_size = 2, spat_point_size = 2.5,
              save_param = list(save_name = paste0(name,'_11_coVis_nr_genes')))
  gobject <- createSpatialGrid(gobject,
                         sdimx_stepsize = 400,
                         sdimy_stepsize = 400,
                         minimum_padding = 0)
  
  gobject <- createSpatialNetwork(gobject = gobject, 
                            method = 'kNN', k = 5, 
                            maximum_distance_knn = 400,
                            minimum_k = 1,
                            name = paste0(name, '_12_spatial_network'))
  return(gobject)
}
