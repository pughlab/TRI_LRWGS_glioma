library(Seurat)
library(limma)
library(EnhancedVolcano)
library(pheatmap)

DATA_DIR <- "/Users/ryanyutian/Desktop/NucSeq_data"       

dat <- readRDS(file.path(DATA_DIR,'GBM_Cohort_NoDoublets_Annotated_July152021_Seurat.rds'))

normalized_data <- Seurat::GetAssayData(dat, assay = "RNA", slot = "data")

sample_ids <- c('A_R_GBM607_180329', 'A_RR_GBM809_180413', 
                'B_P_GBM593.2_181205','B_R_GBM898_181205', 
                'F_P_GBM620_200206', 'F_R_GBM691_181026', 
                'G_P_GBM454_181102', 'G_R_GBM833_181102', 
                'H_P_GBM460_181205', 'H_R_GBM492_181026', 
                'J_P_GBM401_181129', 'J_R_GBM498_190703', 'J_RR_GBM551_181026', 
                'K_P_GBM529_190712', 'K_R_GBM832_200218', 
                'N_P_BT2013110_190703', 'N_R_GBM745_200218', 
                'L_P_GBM618_190703', 'C_P_GBM577.2_181205', 
                'I_P_GBM440_190712', 'M_P_GBM672_190620', 
                'X20_R_SMTB302_190923', 'X23_R_SMTB814_190923')


plot_colnames <- c('A_R_GBM607', 'A_RR_GBM809',
                   'B_P_GBM593','B_R_GBM898',
                   'D_P_GBM620', 'D_R_GBM691', 
                   'E_P_GBM454', 'E_R_GBM833', 
                   'F_P_GBM460', 'F_R_GBM492', 
                   'G_P_GBM401', 'G_R_GBM498', 'G_RR_GBM551', 
                   'H_P_GBM529', 'H_R_GBM832', 
                   'I_P_BT2013110', 'I_R_GBM745', 
                   'X_P_GBM618', 'X_P_GBM577', 
                   'X_P_GBM440', 'X_P_GBM672', 
                   'X_R_SMTB302', 'X_R_SMTB814')

plot_col_order <- c('D_P_GBM620', 'G_P_GBM401', 'G_R_GBM498', 'G_RR_GBM551', 'X_P_GBM577', 'X_P_GBM672', 'D_R_GBM691', 
                    'H_P_GBM529', 'H_R_GBM832', 'X_R_SMTB814', 'E_R_GBM833', 'I_P_BT2013110', 'I_R_GBM745', 'A_R_GBM607', 
                    'A_RR_GBM809', 'B_P_GBM593', 'B_R_GBM898', 'E_P_GBM454', 'F_P_GBM460', 'F_R_GBM492', 'X_P_GBM618', 
                    'X_P_GBM440', 'X_R_SMTB302')

# Function to calculate the average expression for malignant cells for genes of interest
calculate_avg_malignant_expression <- function(sample_id, genes_of_interest) {
  
  all_cells_in_sample <- colnames(dat)[which(dat@meta.data[, 1] == sample_id)]
  
  sample_seurat <- subset(dat, cells = all_cells_in_sample)
  
  malignant_cells <- colnames(sample_seurat)[sample_seurat@meta.data$Final_Annotation == "Malignant"]
  
  avg_expression <- sapply(genes_of_interest, function(gene) {
    mean(sample_seurat@assays$RNA@data[gene, malignant_cells], na.rm = TRUE)
  })
  
  return(avg_expression)
  
}

# Function to calculate the average expression for non-malignant cells for genes of interest
calculate_avg_non_malignant_expression <- function(sample_id, genes_of_interest) {
  
  all_cells_in_sample <- colnames(dat)[which(dat@meta.data[, 1] == sample_id)]
  
  sample_seurat <- subset(dat, cells = all_cells_in_sample)
  
  non_malignant_cells <- colnames(sample_seurat)[sample_seurat@meta.data$Final_Annotation != "Malignant"]
  
  avg_expression <- sapply(genes_of_interest, function(gene) {
    mean(sample_seurat@assays$RNA@data[gene, non_malignant_cells], na.rm = TRUE)
  })
  
  return(avg_expression)
  
}
################################################################################

# Heatmap of Malignant Cell Expressions

genes_of_interest = c('EGFR', 'CDKN2A', 'CDKN2B', 'TERT', 'PTEN', 'NF1', 'ARID1B')

# Loop through all the sample ids and store the result in a matrix
avg_expression_matrix_malig <- matrix(0, nrow=length(genes_of_interest), ncol=length(sample_ids))
rownames(avg_expression_matrix_malig) <- genes_of_interest
colnames(avg_expression_matrix_malig) <- plot_colnames

for (i in 1:length(sample_ids)) {
  avg_expression_matrix_malig[, i] <- calculate_avg_malignant_expression(sample_ids[i], genes_of_interest)
}

# Plot the heatmap
colors <- colorRampPalette(c("blue", "white", "red"))(400)
breaks_scale <- seq(0, 4, length.out = 401)

pheatmap(avg_expression_matrix_malig[],
         scale = "none",
         color = colors,
         breaks = breaks_scale,
         main = "Average Malignant Gene Expression",
         show_colnames = TRUE, 
         show_rownames = TRUE,
         cluster_rows = FALSE,  # Disable row clustering
         cluster_cols = FALSE,  # Disable column clustering
)

################################################################################

# Heatmap of Non-Malignant Cell Expressions

genes_of_interest = c('EGFR', 'CDKN2A', 'CDKN2B', 'TERT', 'PTEN', 'NF1', 'ARID1B')

# Loop through all the sample ids and store the result in a matrix
avg_expression_matrix_nonmalig <- matrix(0, nrow=length(genes_of_interest), ncol=length(sample_ids))
rownames(avg_expression_matrix_nonmalig) <- genes_of_interest
colnames(avg_expression_matrix_nonmalig) <- plot_colnames

for (i in 1:length(sample_ids)) {
  avg_expression_matrix_nonmalig[, i] <- calculate_avg_non_malignant_expression(sample_ids[i], genes_of_interest)
}

# Plot the heatmap
colors <- colorRampPalette(c("blue", "white", "red"))(400)
breaks_scale <- seq(0, 4, length.out = 401)

pheatmap(avg_expression_matrix_nonmalig[],
         scale = "none",
         color = colors,
         breaks = breaks_scale,
         main = "Average Non-Malignant Gene Expression",
         show_colnames = TRUE, 
         show_rownames = TRUE,
         cluster_rows = FALSE,  # Disable row clustering
         cluster_cols = FALSE,  # Disable column clustering
)

