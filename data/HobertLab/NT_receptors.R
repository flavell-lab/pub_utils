##### Analysis of Neurotransmitter receptor genes.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
N2_neu = qs::qread("../data/cell_datasets/240214_N2_neurons.qs") 
AF16_neu = qs::qread("../data/cell_datasets/240214_AF16_neurons.qs")
NIC203_neu = qs::qread("../data/cell_datasets/240214_NIC203_neurons.qs") 
### Gene and orthologies tables
N2_families_PFAM = qs::qread("../data/N2_familiesPFAM.qs")
AF16_familiesPFAM = qs::qread("../data/AF16_familiesPFAM.qs")
NIC203_familiesPFAM = qs::qread("../data/NIC203_familiesPFAM.qs")

# library(monocle3)
library(pheatmap)
library(RColorBrewer)
library(readxl)
library(writexl)
library(ggplot2)
library(viridis)
library(dplyr)
library(ggbeeswarm)
library(dunn.test)
library(scico)
require(tidyr)
require(scales)

# Correcting N2 gene IDs to match name format in the other species.
test <- N2_neu
neuCDS_info <- data.frame(
  external_gene_id = rowData(test)$external_gene_id,
  wormbase_gene = rowData(test)$wormbase_gene
)
neuCDS_info$updated_gene_name <- N2_families_PFAM[match(neuCDS_info$wormbase_gene, N2_families_PFAM$gene_id),]$gene_name
rowData(test)$updated_gene_id <- NA
rowData(test)$updated_gene_id[match(rowData(test)$wormbase_gene, rowData(test)$wormbase_gene)] <- neuCDS_info$updated_gene_name
test <- subset(test, !is.na(rowData(test)$updated_gene_id))
rownames(test) <- rowData(test)$updated_gene_id
N2_neu <- test
test2 <- N2_families_PFAM[complete.cases(N2_families_PFAM$gene_name), ]
N2_families_PFAM <- test2
rm(test,test2,neuCDS_info)
#####

N2_binary = qs::qread("../data/20240215_N2_RF_binary.qs")
AF16_binary = qs::qread("../data/20240215_AF16_RF_binary.qs")
NIC203_binary = qs::qread("../data/20240215_NIC203_RF_binary.qs")
identical(rownames(N2_binary),rownames(N2_neu))
identical(rownames(AF16_binary),rownames(AF16_neu))
identical(rownames(NIC203_binary),rownames(NIC203_neu))

#### Correcting gene annotations based on manual examination of gene models. 
cbrbas1_index <- which(AF16_familiesPFAM$gene_short_name == "QX1410_007145")
AF16_familiesPFAM$cb_short_name[cbrbas1_index] <- "cb_bas-1"
AF16_familiesPFAM$gene_short_name[cbrbas1_index] <- "bas-1"
ctrbas1_index <- which(NIC203_familiesPFAM$ctrop_short_name == "ORF008196")
NIC203_familiesPFAM$ctrop_short_name[ctrbas1_index] <- "ct_bas-1"
NIC203_familiesPFAM$gene_short_name[ctrbas1_index] <- "bas-1"
rownames(AF16_binary)[cbrbas1_index] <- "bas-1"
rownames(NIC203_binary)[ctrbas1_index] <- "bas-1"
rownames(AF16_neu)[cbrbas1_index] <- "bas-1"
rownames(NIC203_neu)[ctrbas1_index] <- "bas-1"

cbrcha1_index <- which(AF16_familiesPFAM$gene_short_name == "QX1410_011872")
cbrcha1.2_index <- which(AF16_familiesPFAM$gene_short_name == "cha-1")
AF16_familiesPFAM$cb_short_name[cbrcha1_index] <- "cb_cha-1"
AF16_familiesPFAM$gene_short_name[cbrcha1_index] <- "cha-1"
AF16_familiesPFAM$cb_short_name[cbrcha1.2_index] <- "cb_cha-1.2"
AF16_familiesPFAM$gene_short_name[cbrcha1.2_index] <- "cha-1.2"
rownames(AF16_binary)[cbrcha1_index] <- "cha-1"
rownames(AF16_neu)[cbrcha1_index] <- "cha-1"
rownames(AF16_binary)[cbrcha1.2_index] <- "cha-1.2"
rownames(AF16_neu)[cbrcha1.2_index] <- "cha-1.2"

cbrceh37_index <- which(AF16_familiesPFAM$gene_short_name == "QX1410_017391")
AF16_familiesPFAM$cb_short_name[cbrceh37_index] <- "cb_ceh-37"
AF16_familiesPFAM$gene_short_name[cbrceh37_index] <- "ceh-37"
rownames(AF16_binary)[cbrceh37_index] <- "ceh-37"
rownames(AF16_neu)[cbrceh37_index] <- "ceh-37"

ctrmab5_index <- which(NIC203_familiesPFAM$ID == "ORF010213")
NIC203_familiesPFAM$ctrop_short_name[ctrmab5_index] <- "ct_mab-5"
NIC203_familiesPFAM$gene_short_name[ctrmab5_index] <- "mab-5"
rownames(NIC203_binary)[ctrmab5_index] <- "mab-5"
rownames(NIC203_neu)[ctrmab5_index] <- "mab-5"

identical(rownames(N2_binary),rownames(N2_neu))
identical(rownames(AF16_binary),rownames(AF16_neu))
identical(rownames(NIC203_binary),rownames(NIC203_neu))

### Keep only neuron classes that were detected in all three species. ASE and I3 are missing in the C. tropicalis datasets. AWA have very few cells sequenced in the C. briggsae dataset.
N2_classes_toremove = c("unknown_5","unknown_5b","unknown_6","ASEs","I3", "AWA")
AF16_classes_toremove = c("ASEs","I3","AWA")
NIC203_classes_toremove = c("unknown_2","unknown_3","AWA")

N2_neu_common = N2_neu[,!N2_neu@colData$unified_Aug2023 %in% N2_classes_toremove]
AF16_neu_common = AF16_neu[,!AF16_neu@colData$unified_Aug2023 %in% AF16_classes_toremove]
NIC203_neu_common = NIC203_neu[,!NIC203_neu@colData$unified_Aug2023 %in% NIC203_classes_toremove]


### Keep only 1-to-1-to-1 orthologs (11,280 genes). 3 are going to be removed later on because their gene models are bad. Total of 11,277 1-1-1 orthologs in our datasets.
mat_list <- list(N2_neu_common, AF16_neu_common, NIC203_neu_common)
common_rows <- Reduce(intersect, lapply(mat_list, rownames)) 
Only_common_orthologs = common_rows ### Only_common_orthologs is the list of 1-to-1-to-1 orthologs. Will be used to filter gene families to keep only the common genes.
###

amphidsensillaN = c("ADF","ADL","AFD","ASE","ASER","ASEL","ASEs","ASG","ASH_PHA_PHB","ASH","PHA","PHB","PHA_PHB","ASI","ASJ","ASK","AWA","AWB","AWC","AWCs","AWC_ON","AWC_OFF")
other_sensoryN = c("ADE_PDE_CEP","PDE","CEP","ALM_PLM_AVM_PVM","ALM","PLM","AVM","PVM","ALN_PLN","ALN","PLN","AQR_PQR_URX","AQR","PQR","URX","BAG","FLP_PVD","FLP","PVD","IL1","IL2_DV","IL2_LR","OLL","OLQ","PHC","SDQ","URB","URY")
sensoryN = c("ADE_PDE_CEP","PDE","CEP","ADF","ADL","AFD","ALM_PLM_AVM_PVM","ALM","PLM","AVM","PVM","ALN_PLN","ALN","PLN","AQR_PQR_URX","AQR","PQR","URX","ASE","ASER","ASEL","ASEs","ASG","ASH_PHA_PHB","ASH","PHA","PHB","PHA_PHB","ASI","ASJ","ASK","AWA","AWB","AWC","AWCs","AWC_ON","AWC_OFF","BAG","FLP_PVD","FLP","PVD","IL1","IL2_DV","IL2_LR","OLL","OLQ","PHC","SDQ","URB","URY")
interN = c("ADA","AIA","AIB","AIM","AIN","AIY","AIZ","ALA","AUA","AVA","AVB","AVD","AVE","AVF","AVG","AVH","AVJ","AVK","BDU","DVA","DVC","LUA","PVC","PVN","PVP","PVQ","PVR","PVT","PVW","PVW?","RIA","RIB","RIC","RID","RIF","RIG","RIH","RIM","RIP","RIR","RIS","RIV","RMF","RMG","SAA")
motorN = c("AS","AS_SAB","SAB","AVL","DVB","PDA","PDB","RMD_LR","RMD_DV","RME_LR","RME_DV","RMH", "SIA","SIB","SMB","SMD","VC","VA","VB","DA","DB","URA","VD","DD","VA_VB_DA_DB","VD_DD")
pharynxN = c("I1","I2","I3","I4","I5","I6","M1","M2","M3","M4","M5","MC","MI","NSM")
CAN = "CAN"
allNeuronTypes <- c(amphidsensillaN,other_sensoryN, interN, motorN, pharynxN, CAN)
categorize_neuron <- function(cell_group) {
  if (cell_group %in% amphidsensillaN) {
    return("sensilla")
  } else if (cell_group %in% other_sensoryN) {
    return("other_sensory") 
  } else if (cell_group %in% interN) {
    return("Inter")
  } else if (cell_group %in% motorN) {
    return("Motor")
  } else if (cell_group %in% pharynxN) {
    return("Pharynx")
  } else if (cell_group %in% CAN) {
    return("CAN")
  } else {
    return("other")
  }
}
color_vector <- c(sensilla = "magenta", other_sensory = "#ab164c", Inter = "#1E88E5", Motor = "#FFC107", Pharynx = "#004D40", CAN = "black") # from David Nichol's blog.
neuron_cat <- sapply(colnames(N2_binary), categorize_neuron)


# Organize NT receptor df. 
Divscore_df = qs::qread("../data/20241113_Divscore_df.qs")
NT_recDF <- read_excel(path = "../data/gene families/20240110 NT receptors.xlsx", sheet = 1)
NT_recDF = NT_recDF[, 3:7]
colnames(NT_recDF)[1] <- "Gene"

NT_recDF_Div <- merge(NT_recDF, Divscore_df, by = "Gene", all.x = TRUE) #79 genes
NT_recDF_Div <- NT_recDF_Div[!(NT_recDF_Div$Gene %in% c("dop-4","ser-5","glc-1")), ] # remove genes with no 1-to-1-to-1 orthologs. 3 genes dop-4 ser-5 and glc-1. Total remmaining 76 genes.
NT_recDF_Div <- NT_recDF_Div[!is.nan(NT_recDF_Div$DivScore_basic), , drop = FALSE] # Remove genes not expressed in any of the nervous systems (NA in Divscore). 1 gene, 75 genes remaining.
NT_recLabels <- NT_recDF_Div[, c(2,4)] # Gene labels for pheatmaps
rownames(NT_recLabels) <- NT_recDF_Div[[1]] # Gene labels for pheatmaps
neuron_catdf <- as.data.frame(neuron_cat, stringsAsFactors = FALSE) #Cell labels for pheatmaps
colnames(neuron_catdf) <- "neuron_category" #Cell labels for pheatmaps
pheatmap::pheatmap(N2_binary[NT_recDF_Div$Gene,], cluster_cols = F, cluster_rows = F,border_color = "grey40", color = c("grey95","#56B4E9"), legend_breaks = c(0,1), legend_labels = c("off","expressed"), annotation_row = NT_recLabels)
pheatmap::pheatmap(AF16_binary[NT_recDF_Div$Gene,], cluster_cols = F, cluster_rows = F,border_color = "grey40", color = c("grey95","#E69F00"), legend_breaks = c(0,1), legend_labels = c("off","expressed"), annotation_row = NT_recLabels)
pheatmap::pheatmap(NIC203_binary[NT_recDF_Div$Gene,], cluster_cols = F, cluster_rows = F,border_color = "grey40", color = c("grey95","#009E73"), legend_breaks = c(0,1), legend_labels = c("off","expressed"), annotation_row = NT_recLabels)


NT_recFAST <- NT_recDF_Div[grepl("Ach|GABA|Glu", NT_recDF_Div$Confirmed_Ligand), ] ### Only rows that bind Ach Glu or GABA. This include the bimodal receptor lgc-39. 
NT_recFAST[NT_recFAST == "Ach_octopamine_tyramine"] <- "Ach" # In the context of the Fast acting NTs, I consider lgc-39 as an Ach receptor.
ordered_NTrecFAST<- NT_recFAST$Gene[order(NT_recFAST$DivScore_basic, decreasing = TRUE)]
orderedcells = intersect(allNeuronTypes, colnames(N2_binary))

### Analyze Jaccard distances of Neurotransmitter genes grouped by their ligand. 
plot_GeneFamilies_Divscore =function(DivDF, by_Category = Gene_Family, whichScore = DivScore, group_labels = "", xLabel = deparse(substitute(whichScore)), bee = FALSE, yLabel = "Gene Category")
{
  set.seed(123)
  ggplot(DivDF, aes(x = {{ by_Category }}, y = {{ whichScore }})) +
    geom_boxplot(outlier.shape = NA, width= 0.4) +  
    scale_x_discrete(limits = group_labels) + 
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +  
    labs(x = yLabel, y = xLabel) +
    theme_classic() + 
    coord_flip() + 
    theme(
      axis.title.x=element_text(size=20, vjust=-1),
      axis.title.y=element_text(size=20, vjust=2),
      axis.line.x=element_line(lineend="square",size=1),
      axis.line.y = element_blank(),  
      axis.text.x = element_text(size=16),
      axis.text.y = element_text(size=16),
      axis.ticks.y=element_blank(),
      axis.ticks.length=unit(.25, "cm"),
      plot.background = element_rect(fill = "transparent"),
      panel.background = element_rect(fill = "transparent")
    ) +
    if(bee == T){ geom_beeswarm(cex = 1.2,color = "#E69F00", size = 0.7, alpha = 1)
    } else {geom_point(position = position_jitter(width = 0.1), color = "#E69F00", size =0.6, alpha = 0.7)
    }
}

DivScore_NT = qs::qread("../data/20241113_Divscore_df.qs")
colnames(DivScore_NT)[colnames(DivScore_NT) == "DivScore_basic"] <- "Jacc_dis"
DivScore_NT_LGICsGrouped = DivScore_NT
LGICs_subgroups = c("Ionotropic_Glu_rec","nAch_channels_LGIC","Other_LGICs")
DivScore_NT_LGICsGrouped[which(DivScore_NT_LGICsGrouped$Gene_family %in% LGICs_subgroups),]$Gene_family = "NT_LGICs"
head(DivScore_NT_LGICsGrouped)
DivScore_NT_LGICsGrouped <- merge(DivScore_NT_LGICsGrouped, NT_recFAST[, c("Gene", "Confirmed_Ligand")], by = "Gene", all.x = TRUE)
DivScore_NT_LGICsGrouped$Confirmed_Ligand[is.na(DivScore_NT_LGICsGrouped$Confirmed_Ligand)] <- "leftout"
ligand_colors = c("Ach" = "#ef476f", "Glu" = "#ffa600", "GABA" = "#26547c","monoamines" = "#06d6a0", "other" = "grey")

# Reorganize with NTs only. 
NT_only_geneFamilies = c("NT_Synthesis_and_transport","NT_GPCRs","NT_LGICs","NT_Reuptake")
DivScore_NT_only = subset(DivScore_NT_LGICsGrouped, Gene_family %in% NT_only_geneFamilies)
DivScore_NT_only <- DivScore_NT_only[, !colnames(DivScore_NT_only) %in% "Confirmed_Ligand"]
DivScore_NT_only <- merge(DivScore_NT_only, NT_recDF_Div[, c("Gene", "Confirmed_Ligand")], by = "Gene", all.x = TRUE)
DivScore_NT_only <- DivScore_NT_only[!DivScore_NT_only$Gene %in% c("eat-2","F59E12.8"), ] # Removed because not expressed in the NS (Jaccard score is NA). 
DivScore_NT_only <- DivScore_NT_only[,c(1,2,12,3:11)]
cat_other = c("Betaine","adenosine")
cat_monoamine = c("DA","DA_Tyr","DA_Tyr_serotonin","octopamine","serotonin","Tyr","unknown","Ach_octopamine_tyramine") # Included lgc-39 in monoamines because it appears as such in other shown analyses of NT receptors.
DivScore_NT_only[DivScore_NT_only$Confirmed_Ligand %in% cat_other,]$Confirmed_Ligand <- "other"
DivScore_NT_only[DivScore_NT_only$Confirmed_Ligand %in% cat_monoamine,]$Confirmed_Ligand <- "monoamines"
DivScore_NT_only[DivScore_NT_only$Gene %in% c("snf-3"),]$Confirmed_Ligand <- "other"
DivScore_NT_only[DivScore_NT_only$Gene %in% c("bas-1","cat-1","cat-2","cat-4","dat-1", "mod-5","tbh-1","tdc-1","tph-1"),]$Confirmed_Ligand <- "monoamines" 
DivScore_NT_only[DivScore_NT_only$Gene %in% c("cha-1","cho-1","unc-17","snf-6"),]$Confirmed_Ligand <- "Ach"
DivScore_NT_only[DivScore_NT_only$Gene %in% c("eat-4","glt-1","glt-3","glt-4","glt-5","glt-6","glt-7","C08B6.5","T25E4.2","W02A2.5","ZK867.2"),]$Confirmed_Ligand <- "Glu"
DivScore_NT_only[DivScore_NT_only$Gene %in% c("unc-25","unc-46","unc-47","snf-11"),]$Confirmed_Ligand <- "GABA"

DivScore_NT_only_noReup <- subset(DivScore_NT_only, DivScore_NT_only$Gene_family != "NT_Reuptake")

### Plot Jaccard distances of Neurotransmitter-related genes, colored by their ligand.
### 

set.seed(123)
pPerNT = ggplot(DivScore_NT_only_noReup, aes(x=Gene_family, y=Jacc_dis, color = Confirmed_Ligand))+
  geom_boxplot(outlier.shape = NA, width = 0.4, color = "black", fill="transparent", size = 0.8) +
  geom_beeswarm(cex = 2, size = 1.4, alpha = 1) +
  scale_x_discrete(limits = rev(c("NT_Synthesis_and_transport","NT_LGICs","NT_GPCRs")), labels = rev(c("NT_Syn","NT_LGICs","NT_GPCRs"))) + 
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) + 
  scale_color_manual(values = ligand_colors) + 
  labs(y = "Jaccard distance") +
  theme_classic() + 
  coord_flip() + 
  theme(
    axis.title.y=element_blank(),
    axis.title.x=element_text(size=20, vjust=2),
    axis.line.y= element_line(lineend="square",size=1),
    axis.line.x = element_line(lineend="square",size=1),
    axis.text.y = element_text(size=16),
    axis.text.x = element_text(size=16),
    axis.ticks = element_line(size = 1),
    axis.ticks.length=unit(.25, "cm"),
    plot.background = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent"),
    panel.border = element_rect(fill = "transparent", color = "black",size = 1)
  )
pPerNT
ggsave(plot = pPerNT,width = 8,height = 4, path="outputs/", filename = "Jacc_NTrec_per_ligand.pdf", device = "pdf", dpi = 300,units = "in")

dunn_NT_only <- dunn.test(DivScore_NT_only_noReup$Jacc_dis, DivScore_NT_only_noReup$Confirmed_Ligand, method = "bh") 
aggregate(Jacc_dis ~ Confirmed_Ligand, data = DivScore_NT_only_noReup, FUN = function(x) c(Median = median(x), Mean = mean(x), SD = sd(x), N = length(x)))
dunn.test(DivScore_NT_only_noReup[which(DivScore_NT_only_noReup$Gene_family == "NT_LGICs"),]$Jacc_dis, DivScore_NT_only_noReup[which(DivScore_NT_only_noReup$Gene_family == "NT_LGICs"),]$Confirmed_Ligand, method = "bh") 
aggregate(Jacc_dis ~ Confirmed_Ligand, data = DivScore_NT_only_noReup[which(DivScore_NT_only_noReup$Gene_family == "NT_LGICs"),], FUN = function(x) c(Median = median(x), Mean = mean(x), SD = sd(x), N = length(x)))
dunn.test(DivScore_NT_only_noReup[which(DivScore_NT_only_noReup$Gene_family == "NT_GPCRs"),]$Jacc_dis, DivScore_NT_only_noReup[which(DivScore_NT_only_noReup$Gene_family == "NT_GPCRs"),]$Confirmed_Ligand, method = "bh") 
aggregate(Jacc_dis ~ Confirmed_Ligand, data = DivScore_NT_only_noReup[which(DivScore_NT_only_noReup$Gene_family == "NT_GPCRs"),], FUN = function(x) c(Median = median(x), Mean = mean(x), SD = sd(x), N = length(x)))



# For each neuron, calculate how much receptiveness it has to every NT in every species. 
# Meaning, for every neuron in each species, how many Glu receptors. Then subdivide => Glu_exc, Glu_inh, Glu_mod.
# Then aminergic receptors too.

# Classify the NT receptor genes into their ligand and function group.
GluRs_total <- NT_recFAST$Gene[NT_recFAST$Confirmed_Ligand == "Glu"]
AchRs_total <- NT_recFAST$Gene[NT_recFAST$Confirmed_Ligand == "Ach"]
GABARs_total <- NT_recFAST$Gene[NT_recFAST$Confirmed_Ligand == "GABA"]
GluRs_exc <- NT_recFAST$Gene[NT_recFAST$Confirmed_Ligand == "Glu" & NT_recFAST$Effect == "excitatory"]
GluRs_inh <- NT_recFAST$Gene[NT_recFAST$Confirmed_Ligand == "Glu" & NT_recFAST$Effect == "inhibitory"]
GluRs_mod <- NT_recFAST$Gene[NT_recFAST$Confirmed_Ligand == "Glu" & NT_recFAST$Effect == "modulatory"]
AchRs_exc <- NT_recFAST$Gene[NT_recFAST$Confirmed_Ligand == "Ach" & NT_recFAST$Effect == "excitatory"]
AchRs_inh <- NT_recFAST$Gene[NT_recFAST$Confirmed_Ligand == "Ach" & NT_recFAST$Effect == "inhibitory"]
AchRs_mod <- NT_recFAST$Gene[NT_recFAST$Confirmed_Ligand == "Ach" & NT_recFAST$Effect == "modulatory"]
GABARs_exc <- NT_recFAST$Gene[NT_recFAST$Confirmed_Ligand == "GABA" & NT_recFAST$Effect == "excitatory"]
GABARs_inh <- NT_recFAST$Gene[NT_recFAST$Confirmed_Ligand == "GABA" & NT_recFAST$Effect == "inhibitory"]
GABARs_mod <- NT_recFAST$Gene[NT_recFAST$Confirmed_Ligand == "GABA" & NT_recFAST$Effect == "modulatory"]
Receptors_cat = c(GluRs_total,AchRs_total,GABARs_total,GluRs_exc,GluRs_inh,GluRs_mod,AchRs_exc,AchRs_inh,AchRs_mod,GABARs_exc,GABARs_inh,GABARs_mod)

# For every neuron, counts how many receptors in the category are expressed (column sum).
quantify_NT_receptivity = function(species_name = "", binary_mat){
  newDF <- data.frame(
    Species = rep(species_name, length(orderedcells)),
    Neuron_category = neuron_catdf$neuron_category[match(orderedcells, rownames(neuron_catdf))],
    cell_group = orderedcells
  )
  newDF$GluRs_total <- colSums(binary_mat[GluRs_total,])[match(newDF$cell_group, colnames(binary_mat))]
  newDF$AchRs_total <- colSums(binary_mat[AchRs_total,])[match(newDF$cell_group, colnames(binary_mat))]
  newDF$GABARs_total <- colSums(binary_mat[GABARs_total,])[match(newDF$cell_group, colnames(binary_mat))]
  newDF$GluRs_mod <- colSums(binary_mat[GluRs_mod,])[match(newDF$cell_group, colnames(binary_mat))]
  newDF$GluRs_exc <- colSums(binary_mat[GluRs_exc,])[match(newDF$cell_group, colnames(binary_mat))]
  newDF$GluRs_inh <- colSums(binary_mat[GluRs_inh,])[match(newDF$cell_group, colnames(binary_mat))]
  newDF$AchRs_mod <- colSums(binary_mat[AchRs_mod,])[match(newDF$cell_group, colnames(binary_mat))]
  newDF$AchRs_exc <- colSums(binary_mat[AchRs_exc,])[match(newDF$cell_group, colnames(binary_mat))]
  newDF$AchRs_inh <- colSums(binary_mat[AchRs_inh,])[match(newDF$cell_group, colnames(binary_mat))]
  newDF$GABARs_mod <- colSums(binary_mat[GABARs_mod,])[match(newDF$cell_group, colnames(binary_mat))]
  newDF$GABARs_exc <- colSums(binary_mat[GABARs_exc,])[match(newDF$cell_group, colnames(binary_mat))]
  newDF$GABARs_inh <- colSums(binary_mat[GABARs_inh,])[match(newDF$cell_group, colnames(binary_mat))]
  return(newDF)
}

N2_Rec = quantify_NT_receptivity("elegans",N2_binary)
AF16_Rec = quantify_NT_receptivity("briggsae",AF16_binary)
NIC203_Rec = quantify_NT_receptivity("tropicalis",NIC203_binary)
CrossRec = rbind(N2_Rec, rbind(AF16_Rec,NIC203_Rec))
CrossRec$Species <- factor(CrossRec$Species, levels = c("elegans", "briggsae", "tropicalis"))

### Monoamines
# Generate NTamine_rec df. 
NT_recAmine <- NT_recDF_Div[grepl("octopamine|Tyr|DA|serotonin", NT_recDF_Div$Confirmed_Ligand), ] ### This includes the bimodal receptor lgc-39. 
ordered_NT_recAmine<- NT_recAmine$Gene[order(NT_recAmine$DivScore_basic, decreasing = TRUE)]
NT_recAmineLabels <- NT_recAmine[, c("Gene","Confirmed_Ligand", "Effect")] # Only include the necessary annotation
rownames(NT_recAmineLabels) <- NT_recAmineLabels$Gene # Set row names
NT_recAmineLabels <- NT_recAmineLabels[ , c("Confirmed_Ligand","Effect"), drop = FALSE]

AmineRs_mod <- NT_recAmine$Gene[NT_recAmine$Effect == "modulatory"]
AmineRs_exc <- NT_recAmine$Gene[NT_recAmine$Effect == "excitatory"]
AmineRs_inh <- NT_recAmine$Gene[NT_recAmine$Effect == "inhibitory"]
N2_Rec$AmineRs_total <- colSums(N2_binary[NT_recAmine$Gene,])[match(N2_Rec$cell_group, colnames(N2_binary))]
N2_Rec$AmineRs_mod <- colSums(N2_binary[AmineRs_mod,])[match(N2_Rec$cell_group, colnames(N2_binary))]
N2_Rec$AmineRs_exc <- colSums(t(as.matrix(N2_binary[AmineRs_exc, ])))[match(N2_Rec$cell_group, colnames(N2_binary))]
N2_Rec$AmineRs_inh <- colSums(N2_binary[AmineRs_inh,])[match(N2_Rec$cell_group, colnames(N2_binary))]
AF16_Rec$AmineRs_total <- colSums(AF16_binary[NT_recAmine$Gene,])[match(AF16_Rec$cell_group, colnames(AF16_binary))]
AF16_Rec$AmineRs_mod <- colSums(AF16_binary[AmineRs_mod,])[match(AF16_Rec$cell_group, colnames(AF16_binary))]
AF16_Rec$AmineRs_exc <- colSums(t(as.matrix(AF16_binary[AmineRs_exc, ])))[match(AF16_Rec$cell_group, colnames(AF16_binary))]
AF16_Rec$AmineRs_inh <- colSums(AF16_binary[AmineRs_inh,])[match(AF16_Rec$cell_group, colnames(AF16_binary))]
NIC203_Rec$AmineRs_total <- colSums(NIC203_binary[NT_recAmine$Gene,])[match(NIC203_Rec$cell_group, colnames(NIC203_binary))]
NIC203_Rec$AmineRs_mod <- colSums(NIC203_binary[AmineRs_mod,])[match(NIC203_Rec$cell_group, colnames(NIC203_binary))]
NIC203_Rec$AmineRs_exc <- colSums(t(as.matrix(NIC203_binary[AmineRs_exc, ])))[match(NIC203_Rec$cell_group, colnames(NIC203_binary))]
NIC203_Rec$AmineRs_inh <- colSums(NIC203_binary[AmineRs_inh,])[match(NIC203_Rec$cell_group, colnames(NIC203_binary))]

testtt = CrossRec
testtt = rbind(N2_Rec, rbind(AF16_Rec,NIC203_Rec))
testtt$Species <- factor(testtt$Species, levels = c("elegans", "briggsae", "tropicalis"))
unique_columns <- setdiff(colnames(testtt), colnames(CrossRec))
for (col in unique_columns) {CrossRec[[col]] <- NA}
matching_rows <- intersect(rownames(CrossRec), rownames(testtt))
CrossRec[matching_rows, unique_columns] <- testtt[matching_rows, unique_columns]


# hm_GluR_df = CrossRec[,c(1,3,8,9,7,4)] 
hm_GluR_df = CrossRec[,c("Species","cell_group","GluRs_exc","GluRs_inh","GluRs_mod","GluRs_total")]
# hm_AchR_df = CrossRec[,c(1,3,11,12,10,5)]
hm_AchR_df = CrossRec[,c("Species","cell_group","AchRs_exc","AchRs_inh","AchRs_mod","AchRs_total")]
# hm_GABAR_df = CrossRec[,c(1,3,14,15,13,6)]
hm_GABAR_df = CrossRec[,c("Species","cell_group","GABARs_exc","GABARs_inh","GABARs_mod","GABARs_total")]
# hm_Amine_df = CrossRec[,c(1,3,24,25,23,22)]
hm_Amine_df = CrossRec[,c("Species","cell_group","AmineRs_exc","AmineRs_inh","AmineRs_mod","AmineRs_total")]

formatDFforHM = function(df, receptor = "GluRs"){
  long = df %>% pivot_longer(cols = starts_with(receptor), names_to = "category", values_to = "count") %>%
    unite("Species_Receptor", Species, category, sep = "_")
  wide <- long %>% pivot_wider(names_from = cell_group, values_from = count)
  # widedf <- wide[c(1,6,11,2,7,12,5,10,15,3,8,13,4,9,14),]
  widedf <- wide[c(1,5,9,2,6,10,3,7,11,4,8,12),]
  matrix_data <- as.matrix(widedf[,-1])  # Exclude the first column which is the row names
  rownames(matrix_data) <- widedf[[1]]   # Set the row name
  return(matrix_data)
}
hm_GluR = formatDFforHM(hm_GluR_df, receptor = "GluRs")
hm_AchR = formatDFforHM(hm_AchR_df, receptor = "AchRs")
hm_GABAR = formatDFforHM(hm_GABAR_df, receptor = "GABARs")
hm_Amine = formatDFforHM(hm_Amine_df, receptor = "AmineRs")

pal_zero = rev(scico(14, palette = 'lapaz'))
pal_zero[pal_zero == pal_zero[1]] <- "white"
pal_zero[pal_zero == pal_zero[2]] <- "#fce7e2"

#### Heatmap for the # of expressed NT receptor per cell type. 
#### 
p_hm_GluR = pheatmap(hm_GluR, cluster_rows = F, cluster_cols = F, color = pal_zero[1:(max(hm_GluR)+1)] ,gaps_row = seq(3,12,3), gaps_col = c(10,25,67,84,97),legend_breaks = 0:max(hm_GluR), cellwidth = 7, cellheight = 12, fontsize = 8)
ggsave(plot = p_hm_GluR, filename = "GluR_pheat_stringent.pdf", device = "pdf",width = 12,height = 4, dpi = 300, units = "in", path="outputs/")
p_hm_GABAR = pheatmap(hm_GABAR, cluster_rows = F, cluster_cols = F, color = pal_zero[1:(max(hm_GABAR)+1)], gaps_row = seq(3,12,3), gaps_col = c(10,25,67,84,97),legend_breaks = 0:max(hm_GABAR), cellwidth = 7, cellheight = 12, fontsize = 8)
ggsave(plot = p_hm_GABAR, filename = "GABAR_pheat_stringent.pdf", device = "pdf",width = 12,height = 4, dpi = 300, units = "in", path="outputs/")
p_hm_Amine = pheatmap(hm_Amine, cluster_rows = F, cluster_cols = F, color = pal_zero[1:(max(hm_Amine)+1)], gaps_row = seq(3,12,3), gaps_col = c(10,25,67,84,97),legend_breaks = 0:max(hm_Amine), cellwidth = 7, cellheight = 12, fontsize = 8)
ggsave(plot = p_hm_Amine, filename = "AmineR_pheat_stringent.pdf", device = "pdf",width = 12,height = 4, dpi = 300, units = "in", path="outputs/")



###### Thresholding of cases of change.
###### A receptivity is considered altered in a specific cell type if this cell type express 0 receptors from the category in at least one species and >2 receptors from the category in at least one species.
binary_alt_events = function(mat){
  nm <- matrix(0, nrow = 1, ncol = ncol(mat))
  colnames(nm) <- colnames(mat)
  for (col in 1:ncol(mat)) {
    column_values <- mat[, col]
    if (any(column_values == 0) && any(column_values >= 2)) {
      nm[1, col] <- 1
    } else {
      nm[1, col] <- 0
    }
  }
  return(nm)
}
hm_GluR_total = binary_alt_events(hm_GluR[c(10:12),])
hm_GluR_mod = binary_alt_events(hm_GluR[c(7:9),])
hm_GluR_exc = binary_alt_events(hm_GluR[c(1:3),])
hm_GluR_inh = binary_alt_events(hm_GluR[c(4:6),])
# p_hm_GluR_total = pheatmap(hm_GluR_total, cluster_rows = F, cluster_cols = F, color = c("grey95","black"), gaps_col = c(10,27,65,84,97),legend_breaks = c(0,1), legend_labels = c("consistent receptivity","altered receptivity") , cellwidth = 10, cellheight = 10, fontsize = 8)
# p_hm_GluR_mod = pheatmap(hm_GluR_mod, cluster_rows = F, cluster_cols = F, color = c("grey95","black"), gaps_col = c(10,27,65,84,97),legend_breaks = c(0,1), legend_labels = c("consistent mod","altered mod") , cellwidth = 10, cellheight = 10, fontsize = 8)
# p_hm_GluR_exc = pheatmap(hm_GluR_exc, cluster_rows = F, cluster_cols = F, color = c("grey95","black"), gaps_col = c(10,27,65,84,97),legend_breaks = c(0,1), legend_labels = c("consistent exc","altered exc") , cellwidth = 10, cellheight = 10, fontsize = 8)
# p_hm_GluR_inh = pheatmap(hm_GluR_inh, cluster_rows = F, cluster_cols = F, color = c("grey95","black"), gaps_col = c(10,27,65,84,97),legend_breaks = c(0,1), legend_labels = c("consistent inh","altered inh"), cellwidth = 10, cellheight = 10, fontsize = 8)

hm_GABAR_total = binary_alt_events(hm_GABAR[c(10:12),])
hm_GABAR_mod = binary_alt_events(hm_GABAR[c(7:9),])
hm_GABAR_exc = binary_alt_events(hm_GABAR[c(1:3),])
hm_GABAR_inh = binary_alt_events(hm_GABAR[c(4:6),])
# p_hm_GABAR_total = pheatmap(hm_GABAR_total, cluster_rows = F, cluster_cols = F, color = c("grey95","black"), gaps_col = c(10,27,65,84,97),legend_breaks = c(0,1), legend_labels = c("consistent receptivity","altered receptivity") , cellwidth = 10, cellheight = 10, fontsize = 8)
# p_hm_GABAR_mod = pheatmap(hm_GABAR_mod, cluster_rows = F, cluster_cols = F, color = c("grey95","black"), gaps_col = c(10,27,65,84,97),legend_breaks = c(0,1), legend_labels = c("consistent mod","altered mod") , cellwidth = 10, cellheight = 10, fontsize = 8)
# p_hm_GABAR_exc = pheatmap(hm_GABAR_exc, cluster_rows = F, cluster_cols = F, color = c("grey95","black"), gaps_col = c(10,27,65,84,97),legend_breaks = c(0,1), legend_labels = c("consistent exc","altered exc") , cellwidth = 10, cellheight = 10, fontsize = 8)
# p_hm_GABAR_inh = pheatmap(hm_GABAR_inh, cluster_rows = F, cluster_cols = F, color = c("grey95","black"), gaps_col = c(10,27,65,84,97),legend_breaks = c(0,1), legend_labels = c("consistent inh","altered inh"), cellwidth = 10, cellheight = 10, fontsize = 8)

hm_AmineR_total = binary_alt_events(hm_Amine[c(10:12),])
hm_AmineR_mod = binary_alt_events(hm_Amine[c(7:9),])
hm_AmineR_exc = binary_alt_events(hm_Amine[c(1:3),])
hm_AmineR_inh = binary_alt_events(hm_Amine[c(4:6),])
# p_hm_AmineR_total = pheatmap(hm_AmineR_total, cluster_rows = F, cluster_cols = F, color = c("grey95","black"), gaps_col = c(10,27,65,84,97),legend_breaks = c(0,1), legend_labels = c("consistent receptivity","altered receptivity") , cellwidth = 10, cellheight = 10, fontsize = 8)
# p_hm_AmineR_mod = pheatmap(hm_AmineR_mod, cluster_rows = F, cluster_cols = F, color = c("grey95","black"), gaps_col = c(10,27,65,84,97),legend_breaks = c(0,1), legend_labels = c("consistent mod","altered mod") , cellwidth = 10, cellheight = 10, fontsize = 8)
# p_hm_AmineR_exc = pheatmap(hm_AmineR_exc, cluster_rows = F, cluster_cols = F, color = c("grey95","black"), gaps_col = c(10,27,65,84,97),breaks=c(0,1),legend_breaks = c(0,1), legend_labels = c("consistent exc","altered exc") , cellwidth = 10, cellheight = 10, fontsize = 8)
# p_hm_AmineR_inh = pheatmap(hm_AmineR_inh, cluster_rows = F, cluster_cols = F, color = c("grey95","black"), gaps_col = c(10,27,65,84,97),legend_breaks = c(0,1), legend_labels = c("consistent inh","altered inh"), cellwidth = 10, cellheight = 10, fontsize = 8)


### Now a combined matrix of total receptors from the three NTs to calculate how many cell types display a receptivity loss event in any of the three NTs.
ThreeNTsChanges = rbind(hm_GluR_total,hm_GABAR_total,hm_AmineR_total)
rownames(ThreeNTsChanges) <- c("GluR_alterations", "GabaR_alterations", "AmineR_alterations")
howmanyalterations = colSums(ThreeNTsChanges)
ThreeNTsChanges <- rbind(ThreeNTsChanges, howmanyalterations) ## Create a fourth row with the sum result of the three existing rows. 
sum(ThreeNTsChanges["howmanyalterations", ] >= 1) # 32 cell groups out of 98 had an alteration (complete loss of receptivity, stringent) in at least one NT type. That's 33%. 

#Complete loss of receptivity (all kinds of receptors) to a neurotransmitter system in a cell type in at least one species. 
p_ThreeNTsChanges = pheatmap(ThreeNTsChanges, cluster_rows = F, cluster_cols = F, color = c("grey95","#fdc70c","#b01111"), gaps_col = c(10,25,67,84,97),legend_breaks = c(0,1,2), cellwidth = 7, cellheight = 12, fontsize = 8, labels_row = c("Glu alterations","GABA alterations","Monoamine alterations","Combined"), gaps_row = 3, )
ggsave(plot = p_ThreeNTsChanges, path="outputs/", filename = "NTcombined_ReceptivityLoss.eps", device = "eps",width = 18,height = 2, dpi = 600, units = "in", bg = "transparent")
ggsave(plot = p_ThreeNTsChanges, path="outputs/", filename = "NTcombined_ReceptivityLoss.pdf", device = "pdf",width = 18,height = 2, dpi = 600, units = "in", bg = "transparent")


### Now same but for alteration of receptivity.
### Not complete loss. Instead, loss of all excitatory OR all inhibitory OR all modulatory receptors for a given neurotransmitter ligand.

### A combined matrix for any type of change in any of the NTs.
Glu_anyChange = rbind(hm_GluR_total,hm_GluR_mod,hm_GluR_exc,hm_GluR_inh)
rownames(Glu_anyChange) <- c("GluR_total", "GluR_mod","GluR_exc","GluR_inh")
Glu_howmanyalterations = colSums(Glu_anyChange)
Glu_anyChange <- rbind(Glu_anyChange, Glu_howmanyalterations) 
binarized_Glu_anyChange <- ifelse(Glu_anyChange["Glu_howmanyalterations", ] > 0, 1, Glu_anyChange["Glu_howmanyalterations", ])
Glu_anyChange <- rbind(Glu_anyChange, binarized_Glu_anyChange) 

GABA_anyChange = rbind(hm_GABAR_total,hm_GABAR_mod,hm_GABAR_exc,hm_GABAR_inh)
rownames(GABA_anyChange) <- c("GABAR_total", "GABAR_mod","GABAR_exc","GABAR_inh")
GABA_howmanyalterations = colSums(GABA_anyChange)
GABA_anyChange <- rbind(GABA_anyChange, GABA_howmanyalterations) 
binarized_GABA_anyChange <- ifelse(GABA_anyChange["GABA_howmanyalterations", ] > 0, 1, GABA_anyChange["GABA_howmanyalterations", ])
GABA_anyChange <- rbind(GABA_anyChange, binarized_GABA_anyChange) 

Amine_anyChange = rbind(hm_AmineR_total,hm_AmineR_mod,hm_AmineR_exc,hm_AmineR_inh)
rownames(Amine_anyChange) <- c("AmineR_total", "AmineR_mod","AmineR_exc","AmineR_inh")
Amine_howmanyalterations = colSums(Amine_anyChange)
Amine_anyChange <- rbind(Amine_anyChange, Amine_howmanyalterations) 
binarized_Amine_anyChange <- ifelse(Amine_anyChange["Amine_howmanyalterations", ] > 0, 1, Amine_anyChange["Amine_howmanyalterations", ])
Amine_anyChange <- rbind(Amine_anyChange, binarized_Amine_anyChange) 

ThreeNTs_anyChangeBinarized = rbind(binarized_Glu_anyChange, binarized_GABA_anyChange, binarized_Amine_anyChange)
Sum_anyChange = colSums(ThreeNTs_anyChangeBinarized)
ThreeNTs_anyChangeBinarized <- rbind(ThreeNTs_anyChangeBinarized, Sum_anyChange) ## Create a fourth row with the sum result of the three existing rows. 
sum(ThreeNTs_anyChangeBinarized["Sum_anyChange", ] >= 1) # 53 cell groups out of 98 had an alteration (stringent) in at least one NT type. That's 54%. 
p_ThreeNTsanyChange = pheatmap(ThreeNTs_anyChangeBinarized, cluster_rows = F, cluster_cols = F, color = c("grey95","#fdc70c","#b01111"), gaps_col = c(10,25,67,84,97),legend_breaks = c(0,1,2), cellwidth = 7, cellheight = 12, fontsize = 8, labels_row = c("Glu alterations","GABA alterations","Monoamine alterations","Combined"), gaps_row = 3)
ggsave(plot = p_ThreeNTsanyChange, path="outputs/", filename = "NTcombined_anyAlteration.pdf", device = "pdf",width = 18,height = 2, dpi = 300, units = "in", bg = "transparent")
ggsave(plot = p_ThreeNTsanyChange, path="outputs/", filename = "NTcombined_anyAlteration.eps", device = "eps",width = 18,height = 2, dpi = 600, units = "in", bg = "transparent")


#### Part of whole plots (Right side of Figures 3h and 3i)
# Receptivity alteration summary
rowSums(ThreeNTs_anyChangeBinarized == 0)
rowSums(ThreeNTs_anyChangeBinarized == 1)
rowSums(ThreeNTs_anyChangeBinarized == 2)

rowSums(ThreeNTs_anyChangeBinarized > 0)
#53 out of 98 groups of neuron classes display an NT alteration (54%)

# Receptivity_alteration in any subtype category of NT receptor
NT <- c(rep("Glu",3), rep("GABA",3),rep("Amine",3), rep("Combined",3))
Change_value<- c(rep(c("consistent","1 altered NT","2 altered NT"),4))
counts <- c(66,32,0,88,10,0,75,23,0,45,41,12)
anyChangePoW <-  data.frame(NT,Change_value,counts)
anyChangePoW$NT <- factor(anyChangePoW$NT, levels = rev(c("Glu","GABA","Amine","Combined")))
anyChangePoW$Change_value <- factor(anyChangePoW$Change_value, levels = c("consistent","1 altered NT","2 altered NT"))
BP_anyChange = ggplot(anyChangePoW, aes(fill=Change_value, x=counts, y=NT)) + 
  geom_bar(position="fill", stat="identity", width = 0.8, color = "black", size=1.5)+
  # scale_y_continuous(limits = c(0,1),breaks = c(0,100), expand = c(0,0))+
  scale_x_continuous(position = "top", labels = percent_format(), breaks = c(0,0.5,1), expand = c(0.002,0.002))+
  scale_fill_manual(values = c("grey95","#fdc70c","#b01111"))+
  theme_classic()+
  theme(axis.line.x = element_line(size=1),
        axis.line.y = element_blank(),
        axis.ticks.x = element_line(size=1),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(hjust = 1),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
  )
BP_anyChange

rowSums(ThreeNTsChanges == 0)
rowSums(ThreeNTsChanges == 1)
rowSums(ThreeNTsChanges == 2)
rowSums(ThreeNTsChanges > 0)
#32 out of 98 groups of neuron classes display complete loss to at least one NT system (33%)

# Receptivity_loss of an entire NT.
NT <- c(rep("Glu",3), rep("GABA",3),rep("Amine",3), rep("Combined",3))
Change_value<- c(rep(c("consistent","1 altered NT","2 altered NT"),4))
counts <- c(84,14,0,89,9,0,86,12,0,66,29,3)
RecLossPoW <-  data.frame(NT,Change_value,counts)
RecLossPoW$NT <- factor(RecLossPoW$NT, levels = rev(c("Glu","GABA","Amine","Combined")))
RecLossPoW$Change_value <- factor(RecLossPoW$Change_value, levels = c("consistent","1 altered NT","2 altered NT"))
BP_RecLoss = ggplot(RecLossPoW, aes(fill=Change_value, x=counts, y=NT)) + 
  geom_bar(position="fill", stat="identity", width = 0.8, color = "black", size=1.5)+
  # scale_y_continuous(limits = c(0,1),breaks = c(0,100), expand = c(0,0))+
  scale_x_continuous(position = "top", labels = percent_format(), breaks = c(0,0.5,1), expand = c(0.002,0.002))+
  scale_fill_manual(values = c("grey95","#fdc70c","#b01111"))+
  theme_classic()+
  theme(axis.line.x = element_line(size=1),
        axis.line.y = element_blank(),
        axis.ticks.x = element_line(size=1),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(hjust = 1),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
  )
BP_RecLoss
BP_anyChange
ggsave(plot = BP_anyChange, path="outputs/", filename = "OOWBP_anyAlteration.pdf", device = "pdf",width = 8,height = 4.9, dpi = 600, units = "in", bg = "transparent")
ggsave(plot = BP_RecLoss, path="outputs/", filename = "OOWBP_ReceptivityLoss.pdf", device = "pdf",width = 8,height = 4.9, dpi = 600, units = "in", bg = "transparent")


### Further analysis of cases of "neurotransmitter deaf" in one species but not the others. 
hm_total_NTs = rbind(hm_GluR[c(10:12),],hm_GABAR[c(10:12),],hm_Amine[c(10:12),],ThreeNTsChanges)
hm_total_NTs_deaf = hm_total_NTs[,hm_total_NTs["howmanyalterations",] >0]
hm_total_NTs_deaf = as.data.frame(hm_total_NTs_deaf, row.names = rownames(hm_total_NTs_deaf))
hm_total_NTs_deaf_export <- cbind(Rowname = rownames(hm_total_NTs_deaf), hm_total_NTs_deaf)
#export to excel for manual inspection and analysis
# write.xlsx(hm_total_NTs_deaf_export, file = "outputs/NT_deaf.xlsx")






