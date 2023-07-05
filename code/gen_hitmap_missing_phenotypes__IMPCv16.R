library(data.table)
library(tidyverse)

KOMPv16.file = "/Volumes/GoogleDrive/My Drive/Miami_IMPC/data/v16/statistical-results-ALL.csv.gz"

KOMPv16 = fread(KOMPv16.file, header=TRUE, sep=",")
KOMPv16$parameter_name <- trimws(KOMPv16$parameter_name) #remove white spaces
KOMPv16$proc_param_name <- paste0(KOMPv16$procedure_name,"_",KOMPv16$parameter_name)

table(KOMPv16$procedure_name, KOMPv16$data_type)
table(KOMPv16$procedure_name, KOMPv16$statistical_method) # ABR method use reference range plus method instead of linear mixed model
#dat <- KOMPv10.1 %>% select(procedure_name=="Gross Pathology and Tissue Collection")

length(unique(KOMPv16$parameter_name)) # 1541
length(unique(KOMPv16$marker_symbol)) # 8348
table(KOMPv16$data_type)#

KOMPv16.uni <- KOMPv16 %>% filter(data_type=="unidimensional")
dim(KOMPv16.uni)
param.count <- table(KOMPv16.uni$parameter_name)
param.list <- names(param.count)[param.count>500] # use phenotypes with # of tested genes > 500
KOMPv16.uni <- KOMPv16.uni%>%filter(parameter_name %in% param.list)
dim(KOMPv16.uni)
length(unique(KOMPv16.uni$parameter_name)) # 303
length(unique(KOMPv16.uni$marker_symbol)) # 8348

#proc.count <- table(KOMPv16.uni$procedure_name)
#proc.list <- names(proc.count)[proc.count>5000]
#KOMPv16.uni <- KOMPv16.uni%>%filter(procedure_name %in% proc.list)
#dim(KOMPv16.uni)
#length(unique(KOMPv16.uni$parameter_name)) # 465
#length(unique(KOMPv16.uni$marker_symbol)) # 8348


mtest <- table(KOMPv16.uni$parameter_name, KOMPv16.uni$marker_symbol)
pheno.list <- rownames(mtest)
head(pheno.list)

BC.list <- c("BMC/Body weight",
             "Body length",
             "Bone Area",
             "Bone Mineral Content (excluding skull)",
             "Bone Mineral Density (excluding skull)",
             "Fat mass",
             "Fat/Body weight",
             "Lean mass",
             "Lean/Body weight")
BC.list%in%pheno.list
CC.list <- c("Alanine aminotransferase", "Albumin", "Alkaline phosphatase", "Alpha-amylase", "Aspartate aminotransferase",
                           "Calcium", "Chloride", "Creatinine", "Glucose", "HDL-cholesterol", "Iron", "Phosphorus", "Potassium",
                           "Sodium", "Total bilirubin", "Total cholesterol", "Total protein", "Triglycerides",
                           "Urea (Blood Urea Nitrogen - BUN)")

CC.list%in%pheno.list
OF.list <- c("Center average speed", "Center distance travelled", "Center permanence time", "Center resting time",
                           "Distance travelled - total", "Latency to center entry", "Number of center entries", "Percentage center time",
                           "Periphery average speed", "Periphery distance travelled", "Periphery permanence time",
                           "Periphery resting time", "Whole arena average speed", "Whole arena resting time")
OF.list%in%pheno.list

new.pheno.list <- pheno.list
new.pheno.list[pheno.list %in% BC.list] <- "BC"
new.pheno.list[pheno.list %in% CC.list] <- "CC"
new.pheno.list[pheno.list %in% OF.list] <- "OF"
new.pheno.list[!pheno.list %in% c(BC.list, CC.list, OF.list)] <- "Others"
table(new.pheno.list)

mtest <- as.data.frame.matrix(mtest)
#rownames(mtest) <- new.pheno.list
max(mtest)
mtest <- mtest>0
class(mtest) <- "numeric"
100*sum(mtest)/(nrow(mtest)*ncol(mtest)) # only 24.4% tested
dim(mtest)

library(viridis)
library(ComplexHeatmap)

## create row annotation
ann_data <- data.frame(Domain = new.pheno.list)
ann <- HeatmapAnnotation(df = ann_data, which="row", show_annotation_name = FALSE,
                         col = list(Domain = c("BC" = "blue", "CC" = "red", "OF" = "green", "Others"="grey")))

## heatmap of gene - phenotype (red: tested, white: untested)
png("docs/figure/figures.Rmd/missing_heatmap_viridis.png", width=8, height=5, units="in", res=300)
ht = Heatmap(mtest,
             cluster_rows = T, clustering_distance_rows ="binary",
             cluster_columns = T, clustering_distance_columns = "binary",
             show_row_dend = F, show_column_dend = F,  # do not show dendrogram
             show_row_names = F, show_column_names = F, col = c("#440154","#fde725"),
             row_names_gp = gpar(fontsize = 10),
             show_heatmap_legend = FALSE,
             right_annotation = ann,
             row_title = "Phenotypes",
             column_title = "Genes")
draw(ht)
dev.off()


