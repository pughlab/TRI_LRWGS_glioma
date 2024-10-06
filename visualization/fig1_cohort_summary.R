setwd("~/Desktop/Manuscript/cohort/")
library(reticulate)
library(Seurat)
library(ComplexHeatmap)
library(readxl)

file <- "~/Desktop/Manuscript/cohort/202404_clinical_info.csv"
meta <- read.csv(file)
head(meta)

meta$NewSampleID <- factor(meta$NewSampleID, levels = unique(meta$NewSampleID))

### Tumour Stage

plot <- as.data.frame.matrix(xtabs(~NewSampleID + TumourStage, data=meta))
plot$Recurrence <- gsub(1, 2, plot$Recurrence)
plot$SecondRecurrence <- gsub(1, 3, plot$SecondRecurrence)

stage <- Heatmap(t(plot), 
                 name = "Stage",
                 col = c("white", "darkgreen", "darkred", "darkblue"),
                 border_gp = gpar(col = "black", lty = 1, lwd = 1),
                 show_row_names = TRUE,
                 rect_gp = gpar(col = "white", lwd = 2),
                 cluster_rows = FALSE,
                 cluster_columns = FALSE
)
stage

### Sex

var <- "Sex"
plot <- data.frame(meta[ ,var])
colnames(plot) <- var
rownames(plot) <- meta$NewSampleID

sex <- Heatmap(t(plot), 
               name = "Sex",
               col = c("orchid3", "lightskyblue1"),
               border_gp = gpar(col = "black", lty = 1, lwd = 1),
               show_row_names = FALSE,
               rect_gp = gpar(col = "white", lwd = 2)
)
sex

### Age Cohort

var <- "AgeCohort"
plot <- data.frame(meta[ ,var])
colnames(plot) <- var
rownames(plot) <- meta$NewSampleID

age <- Heatmap(t(plot), 
               name = var,
               col = c("darkblue", "gold"),
               border_gp = gpar(col = "black", lty = 1),
               show_row_names = FALSE,
               cluster_columns = FALSE,
               rect_gp = gpar(col = "white", lwd = 2)
)
age

### Tumour Grade

plot <- as.data.frame.matrix(xtabs(~NewSampleID + TumourGrade, data=meta))
head(plot)

var <- "TumourGrade"

grade <- Heatmap(t(plot), 
                 name = var,
                 col = c("white", "black"),
                 border_gp = gpar(col = "black", lty = 1, lwd = 1),
                 show_row_names = TRUE,
                 rect_gp = gpar(col = "white", lwd = 2),
                 cluster_rows = FALSE,
                 cluster_columns = FALSE
)
grade

### Pathology

plot <- as.data.frame.matrix(xtabs(~NewSampleID + Pathology, data=meta))
head(plot)

var <- "Pathology"


path <- Heatmap(t(plot), 
                name = var,
                col = c("white", "black"),
                border_gp = gpar(col = "black", lty = 1, lwd = 1),
                show_row_names = TRUE,
                rect_gp = gpar(col = "white", lwd = 2),
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                row_order = c(1, 3, 4, 2)
                
)
path

### Treatment

plot <- as.data.frame.matrix(xtabs(~NewSampleID + Treatment, data=meta))
head(plot)

var <- "Treatment"
plot$N.D. <- gsub(1, 2, plot$N.D. )

treat <- Heatmap(t(plot), 
                 name = var,
                 col = c("white", "black", "grey"),
                 border_gp = gpar(col = "black", lty = 1, lwd = 1),
                 show_row_names = TRUE,
                 rect_gp = gpar(col = "white", lwd = 2),
                 cluster_rows = FALSE,
                 cluster_columns = FALSE,
                 row_order = c(1, 4, 3, 2)
)
treat

### Matched Normal

plot <- as.data.frame.matrix(xtabs(~ NewSampleID + MatchedNormal, data=meta))
plot$Yes[plot$N.D. > 0] <- 2
head(plot)

var <- "MatchedNormal"

normal <- Heatmap(t(plot$Yes), 
                  name = var,
                  col = c("white", "black"),
                  border_gp = gpar(col = "black", lty = 1, lwd = 1),
                  show_row_names = TRUE,
                  rect_gp = gpar(col = "white", lwd = 2),
                  cluster_rows = FALSE,
                  cluster_columns = FALSE
)
normal


### NucleiSeq

plot <- as.data.frame.matrix(xtabs(~ NewSampleID + NucleiSeq, data=meta))
plot$Yes[plot$N.D. > 0] <- 2
head(plot)

var <- "NucleiSeq"

nuclei <- Heatmap(t(plot$Yes), 
               name = var,
               col = c("white", "black"),
               border_gp = gpar(col = "black", lty = 1, lwd = 1),
               show_row_names = TRUE,
               rect_gp = gpar(col = "white", lwd = 2),
               cluster_rows = FALSE,
               cluster_columns = FALSE
)
nuclei

### IDH Status

plot <- as.data.frame.matrix(xtabs(~ NewSampleID + IDH1StatusFromLongRanger, data=meta))
plot$R132H <- gsub(1, 2, plot$R132H)
plot$R132S <- gsub(1, 3, plot$R132S)
plot$R132C <- gsub(1, 4, plot$R132C)
plot$WT <- gsub(1, 5, plot$WT)

head(plot)

var <- "IDH1StatusFromLongRanger"

idh <- Heatmap(t(plot), 
                  name = var,
                  col = c("white", "palevioletred1", "royalblue1", "purple3", "darkorange", "black"),
                  border_gp = gpar(col = "black", lty = 1, lwd = 1),
                  show_row_names = TRUE,
                  rect_gp = gpar(col = "white", lwd = 2),
                  cluster_rows = FALSE,
                  cluster_columns = FALSE
)
idh

##

ht_list = age %v% sex %v% path %v% grade %v% idh %v% treat %v% normal %v% nuclei %v% stage
draw(ht_list)