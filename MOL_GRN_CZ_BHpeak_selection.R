
#title: "MOL_GRN_Chao_BHpeak_selection"
#Eneritz  
#date: "2023-11-20"


library(Pando)
library(Seurat)
library(tidyverse)
library(Signac)
library(BSgenome.Mmusculus.UCSC.mm10)
library(GenomicScores)
library(GenomeInfoDb)
library(rtracklayer)
library(GenomicRanges)
library(chromVAR)
library(TFBSTools)
library(chromVARmotifs)
library(JASPAR2022)
library(biomaRt)
library(janitor)

############################
HS_multi.broad_OLG.f <- readRDS("MORNAATACCOL.rds")


#down sample of MOL56

Idents(HS_multi.broad_OLG.f) <- "cellType_OL_merge"
unique(HS_multi.broad_OLG.f$cellType_OL_merge)
MOL56_CC <- subset( HS_multi.broad_OLG.f , ident=("MOL56"))
MOL2_CC <- subset( HS_multi.broad_OLG.f , ident=("MOL2"))
rm(HS_multi.broad_OLG.f)


Idents(MOL56_CC) <- "orig.ident_merge"
table(MOL56_CC$orig.ident_merge)

MOL56.down <- subset( MOL56_CC , downsample = 2000)
###############################################################33
table(MOL56.down$orig.ident_merge)
rm(MOL56_CC)


Idents(MOL2_CC) <- "orig.ident_merge"
table(MOL2_CC$orig.ident_merge)

MOL2.down <- subset( MOL2_CC , downsample = 2000)
###############################################################33
table(MOL2.down$orig.ident_merge)
rm(MOL2_CC)

#MOL2 (same for MOL56)

BH_peaks <- readRDS("DPA_OL.rds")
BH_peaks

MOL2_peaks <- unique(unlist(BH_peaks$MOL2))

length(MOL2_peaks )
rm(BH_peaks)

###########

DefaultAssay(MOL2.down) <- 'peaks'

Idents(MOL2.down) <- "orig.ident_merge"


length(MOL2_peaks)

 DefaultAssay(MOL2.down) <- 'peaks'
 
fragment.path <- 'atac_fragments.tsv.gz'
fragments <- CreateFragmentObject(fragment.path , cells = colnames(MOL2.down) )
# create a gene by cell matrix


macs2_counts <- FeatureMatrix(
 fragments = fragments ,
  features = MOL2_peaks,
  cells = colnames(MOL2.down)
)


DefaultAssay(MOL2.down) <- 'peaks'
annotations <- Annotation(MOL2.down)


# create a new assay using the MACS2 peak set and add it to the Seurat object
MOL2.down[["MOL2_peaks_BH"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = Fragments(MOL2.down),
  annotation = Annotation(MOL2.down) #, 
  #min.features = 1 
)

#################################################
##motifs mm

pfm <- getMatrixSet(
x = JASPAR2022,
opts = list(collection = "CORE",  tax_group = 'vertebrates', all_versions = FALSE))

x <- character()
for(i in 1:length(pfm@listData)){
x[i] <- pfm@listData[[i]]@name
}

motif_2_tf <- data.frame(motif = names(pfm@listData), tf = x, origin = "JASPAR", gene_id = NA, name = NA, family = NA, symbol = NA, motif_tf = NA)
data("motif2tf")
motif2tf <- subset(motif2tf, motif %in% motif_2_tf$motif)

firstup <- function(x) {
substr(x, 1, 1) <- toupper(substr(x, 1, 1))
substr(x, 2, nchar(x)) <- tolower(substr(x, 2, nchar(x)))
x
}

 
mouse_human_genes <- read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")

# separate human and mouse 
mouse <- split.data.frame(mouse_human_genes,mouse_human_genes$Common.Organism.Name)[[2]]
human <- split.data.frame(mouse_human_genes,mouse_human_genes$Common.Organism.Name)[[1]]

# remove some columns
mouse <- mouse[,c(1,4)]
human <- human[,c(1,4)]

# merge the 2 dataset  (note that the human list is longer than the mouse one)
mh_data <- merge.data.frame(mouse,human,by = "DB.Class.Key",all.y = TRUE) 
colnames(mh_data) <- c("DB.Class.Key" , "Symbol_mm"  ,   "tf" )


df1 <- left_join( motif2tf, mh_data, by = "tf") #%>% bind_cols( c(motif , Symbol_mm , origin , gene_id , family) )
df2 <- as_tibble(as.data.frame(cbind( df1$motif , df1$Symbol_mm , df1$origin , df1$gene_id , df1$family)))
colnames(df2) <- c( "motif" , "tf" , "origin" , "gene_id" , "family")

unique(df2$tf)
motif2tf <- df2

library(tidyr) 
motif2tf <- motif2tf %>% drop_na(tf)
###

genes = motif2tf$tf


ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")
gmi <- getBM(attributes=c('ensembl_gene_id', "mgi_symbol","chromosome_name", "start_position"), 
      filters = 'mgi_symbol', 
      values = (genes), 
      mart = ensembl)
dim(motif2tf)
dim(gmi)
motif2tf$tf 
colnames(gmi) <- c( "ensembl_gene_id" , "tf" ,       "chromosome_name", "start_position" )


df1 <- left_join( motif2tf, gmi, by = "tf") 
df2 <- as_tibble(as.data.frame(cbind( df1$motif , df1$tf , df1$origin , df1$ensembl_gene_id , df1$family)))
colnames(df2) <- c( "motif" , "tf" , "origin" , "gene_id" , "family")

motif2tf_mm <- df2

#################################

#GRN

#preprocess to prepare the grn for the inference

MOL2.down <- initiate_grn(MOL2.down, peak_assay = "MOL2_peaks_BH", rna_assay = "RNA", exclude_exons = FALSE)
MOL2.down <- FindVariableFeatures(MOL2.down, assay = "RNA", nfeatures = 2000 ) 
params <- Params(MOL2.down)
genes <- VariableFeatures(MOL2.down, assay = params$rna_a)
genes <- genes[!grepl("Rik", genes)]
length(genes)
VariableFeatures(MOL2.down, assay = params$rna_a) <- genes

MOL2.down <-   FindTopFeatures(MOL2.down , assay = "MOL2_peaks_BH", min.cutoff = 'q0') %>% 
  RunTFIDF(assay = "MOL2_peaks_BH") %>%
find_motifs(pfm = motifs_use, motif_tfs = motif2tf_mm, genome = BSgenome.Mmusculus.UCSC.mm10)

########################################

#infer grn

MOL2.down <- infer_grn(MOL2.down, peak_to_gene_method = "Signac", genes =NULL , parallel = F  ) %>%
find_modules(p_thresh = 1, nvar_thresh = 1, min_genes_per_module = 1, rsq_thresh = 0) %>%  get_network_graph(MOL2.down , graph_name = "full_graph")

########################################

MOL2.filter_mod <-  find_modules( MOL2.down , p_thresh = 0.1,
    nvar_thresh = 2, 
    min_genes_per_module = 1, 
    rsq_thresh = 0.05 ) 

MOL2.down <- MOL2.filter_mod 
#########################################
#tf coefficients

DefaultAssay(MOL2.down) <- 'RNA'
test1 <- coef(MOL2.down)
tf_peak_expr <- aggregate_matrix(t(GetAssayData(MOL2.down)[unique(test1$tf), ]), groups=MOL2.down$orig.ident_merge )
tf_expr_df <- tf_peak_expr %>% as_tibble(rownames='orig.ident_merge') %>% 
    pivot_longer(!orig.ident_merge, names_to='tf', values_to='expr')

#my selection of expressed tfs
transpose_df <- function(df) {
  t_df <- data.table::transpose(df)
  colnames(t_df) <- rownames(df)
  rownames(t_df) <- colnames(df)
  t_df <- t_df %>%
    tibble::rownames_to_column(.data = .) %>%
    tibble::as_tibble(.)
  return(t_df)
}


summary(tf_expr_df$expr )


tf_expr_ea <- tf_peak_expr %>% as_tibble(rownames='orig.ident_merge') %>% as.data.frame() %>% transpose_df() %>% column_to_rownames( var = "rowname") %>% row_to_names(1) %>%  mutate(
    expr_evidence = case_when(
      # Ctrl > 0.5  | Early > 0.5 | Peak > 0.5 | Late > 0.5 ~ "Expr" ,
       Ctrl < 0.01285  & Early < 0.01285 & Peak < 0.01285 & Late < 0.01285 ~ "No" ,
      Ctrl > 0.01285  & Early < Ctrl & Peak < Ctrl & Late < Ctrl ~ "Ctrl" ,
      Ctrl < Early  & Early > 0.01285 & Peak < Early & Late < Early ~ "Early" ,
      Ctrl < Peak  & Early < Peak & Peak > 0.01285 & Late < Peak ~ "Peak" , 
       Ctrl < Late  & Early < Late & Peak < Late & Late > 0.01285 ~ "Late" ,
      .default = "other"
    )
  )

table(tf_expr_ea$expr_evidence)

tf_exp <- tf_expr_ea %>% 
    filter(expr_evidence != 'No')

ct_grn <- test1 %>% dplyr::select(tf, target, corr,  estimate, padj, region) %>% 
    dplyr::filter(tf%in%rownames(tf_exp))
#############################################################

 tf_expr_ct_grn <- tf_expr_ea %>% 
   rownames_to_column( "tf") %>%
                 full_join( ct_grn , 
                           by = "tf") %>%
   dplyr::filter(tf%in%rownames(tf_exp)) %>%
  mutate(tf_act_Ctrl=as.numeric(Ctrl)*estimate) %>%
  mutate(tf_act_Early=as.numeric(Early)*estimate) %>%
  mutate(tf_act_Peak=as.numeric(Peak)*estimate) %>%
  mutate(tf_act_Late=as.numeric(Late)*estimate)  %>% 
  
   pivot_longer(Ctrl:Late, names_to='timepoint', values_to='is_active') %>% 
  group_by(tf, timepoint) %>% 
 mutate(tf_act=as.numeric(is_active)*estimate) %>%
    summarize(tf_act=mean(tf_act)) %>% 
    group_by(tf) %>% 
 
 mutate(
        max_est=max(abs(tf_act)),
        span=max(tf_act) - min(tf_act),
        crossing=sign(min(tf_act))!=sign(max(tf_act)),
        abs_span=abs(max(tf_act)) + abs(min(tf_act))
    ) %>% 
    ungroup() %>% 
    group_by(tf) %>% 
    mutate(high_Timepoint=timepoint[which.max(abs(tf_act))])
 
 
  tf_expr_ct_grn
##########################################################################

MOL2.down@grn
coef(MOL2.down)
modules <- NetworkModules(MOL2.down) 

####################################################

library(ggpubr)
Peak_high_tfs <-  tf_expr_ct_grn %>% filter(high_Timepoint=='Peak', timepoint=='Peak') %>% mutate(max_act=-max(tf_act)) %>% arrange(max_act) %>% pull(tf)
Late_high_tfs <-  tf_expr_ct_grn %>% filter(high_Timepoint=='Late', timepoint=='Late') %>% mutate(max_act=-max(tf_act)) %>% arrange(max_act) %>% pull(tf)
Early_high_tfs <-  tf_expr_ct_grn %>% filter(high_Timepoint=='Early', timepoint=='Early') %>% mutate(max_act=-max(tf_act)) %>% arrange(max_act) %>% pull(tf)
Ctrl_high_tfs <-  tf_expr_ct_grn %>% filter(high_Timepoint=='Ctrl', timepoint=='Ctrl') %>% mutate(max_act=-max(tf_act)) %>% arrange(max_act) %>% pull(tf)

Peak_high_tfs <-  tf_expr_ct_grn %>% filter(high_Timepoint=='Peak') %>% mutate(max_act=-max(tf_act)) %>% arrange(max_act) %>% pull(tf)

tf_order <- unique(c(Early_high_tfs, Peak_high_tfs , Late_high_tfs , Ctrl_high_tfs))
tf_order <- unique(c(Late_high_tfs))    
tf_order <- unique(c(Ctrl_high_tfs ,Early_high_tfs, Peak_high_tfs , Late_high_tfs  ))

# Default order
levels(factor(tf_expr_ct_grn$timepoint))

tf_expr_ct_grn$timepoint <- factor(tf_expr_ct_grn$timepoint, levels = c("Ctrl", "Early", "Peak", "Late"))

#order:
p1 <- ggplot(tf_expr_ct_grn, aes(factor(tf, levels=tf_order), tf_act)) +
    geom_line() +
    geom_point(mapping=aes(color=timepoint)) +
    geom_text(mapping=aes(label=tf), size=0, angle=90) +
    geom_hline(yintercept = 0) +
    scale_color_manual(values = c("#6E876E", "#F5AB00", "#E92A0C", "#97615E")) +
    labs(y= "TF activity (expr. x mean coef.)", x = " ") +
    theme(axis.text.x = element_text(size=5.5  , angle=90)) #+

p1

setEPS()
postscript("MOL2_TFact_GRN_bhpeaks.eps",
 width = 9, 
    height = 6) 
p1
dev.off()


pdf(file = "MOL2_TFact_GRN_bhpeaks.pdf",  
    width = 9,
    height = 6) 
p1
dev.off()

################################################################

##selected tfs


tf_order <- unique((Ctrl_high_tfs)) 
ctrl_grn <- tf_order %>%  as.tibble() %>% mutate(  timepoint= "Ctrl") %>% column_to_rownames("value") %>% rownames_to_column(var="tf") 
tf_order <- unique((Early_high_tfs)) 
early_grn <- tf_order  %>%  as.tibble() %>% mutate(  timepoint= "Early") %>% column_to_rownames("value") %>% rownames_to_column(var="tf") 
tf_order <- unique((Peak_high_tfs)) 
peak_grn <- tf_order  %>%  as.tibble() %>% mutate(  timepoint= "Peak") %>% column_to_rownames("value") %>% rownames_to_column(var="tf") 
tf_order <- unique((Late_high_tfs)) 
late_grn <- tf_order  %>%  as.tibble() %>% mutate(  timepoint= "Late") %>% column_to_rownames("value") %>% rownames_to_column(var="tf") 

MOL2_tf <- bind_rows(ctrl_grn , early_grn , peak_grn , late_grn)
#######################################################################
#figures:

#ctrl
tf_order <- unique((Ctrl_high_tfs)) 
ctrl_grn <- tf_expr_ct_grn %>%  dplyr::filter(high_Timepoint=='Ctrl') 
p1 <- ggplot(ctrl_grn, aes(factor(tf, levels=tf_order), tf_act)) +
    geom_line() +
    geom_point(mapping=aes(color=timepoint)) +
    geom_text(mapping=aes(label=tf), size=0, angle=90) +
    geom_hline(yintercept = 0) +
  
  scale_color_manual(values = c("#6E876E", "#F5AB00", "#E92A0C", "#97615E")) +
  
   labs(y= "TF activity (expr. x mean coef.)", x = " ") +
    theme(axis.text.x = element_text(size=8  , angle=90))

p1 

setEPS()
postscript("MOL2_Ctrl_TFact_GRN_bhpeaks.eps",
 width = 9, 
    height = 6) 
p1
dev.off()


pdf(file = "MOL2_Ctrl_TFact_GRN_bhpeaks.pdf", 
    width = 9, 
    height = 6) 
p1
dev.off()
#Early
tf_order <- unique((Early_high_tfs)) 
early_grn <- tf_expr_ct_grn %>%  dplyr::filter(high_Timepoint=='Early') 
p1 <- ggplot(early_grn, aes(factor(tf, levels=tf_order), tf_act)) +
    geom_line() +
    geom_point(mapping=aes(color=timepoint)) +
    geom_text(mapping=aes(label=tf), size=0, angle=90) +
    geom_hline(yintercept = 0) +
   scale_color_manual(values = c("#6E876E", "#F5AB00", "#E92A0C", "#97615E")) +  
  labs(y= "TF activity (expr. x mean coef.)", x = " ") +
  theme(axis.text.x = element_text(size=5  , angle=90))
p1

setEPS()
postscript("MOL2_Early_TFact_GRN_bhpeaks.eps",
 width = 9,
    height = 6) 
p1
dev.off()


pdf(file = "MOL2_Early_TFact_GRN_bhpeaks.pdf",  
    width = 9, 
    height = 6) 
p1
dev.off()

#Peak
tf_order <- unique((Peak_high_tfs)) 
peak_grn <- tf_expr_ct_grn %>%  dplyr::filter(high_Timepoint=='Peak') 
p1<- ggplot(peak_grn, aes(factor(tf, levels=tf_order), tf_act)) +
    geom_line() +
    geom_point(mapping=aes(color=timepoint)) +
    geom_text(mapping=aes(label=tf), size=0, angle=90) +
    geom_hline(yintercept = 0) +
   scale_color_manual(values = c("#6E876E", "#F5AB00", "#E92A0C", "#97615E")) +  
  labs(y= "TF activity (expr. x mean coef.)", x = " ") +
    theme(axis.text.x = element_text(size=5  , angle=90))

p1
setEPS()
postscript("MOL2_Peak_TFact_GRN_bhpeaks.eps",
 width = 9, 
    height = 6) 
dev.off()


pdf(file = "MOL2_Peak_TFact_GRN_bhpeaks.pdf",  
    width = 9,
    height = 6) 
p1
dev.off()

#Early
tf_order <- unique((Late_high_tfs)) 
late_grn <- tf_expr_ct_grn %>%  dplyr::filter(high_Timepoint=='Late') 
p1 <- ggplot(late_grn, aes(factor(tf, levels=tf_order), tf_act)) +
    geom_line() +
    geom_point(mapping=aes(color=timepoint)) +
    geom_text(mapping=aes(label=tf), size=0, angle=90) +
    geom_hline(yintercept = 0) +
      scale_color_manual(values = c("#6E876E", "#F5AB00", "#E92A0C", "#97615E")) +  
  labs(y= "TF activity (expr. x mean coef.)", x = " ") +
    theme(axis.text.x = element_text(size=5  , angle=90 )) 
 
p1

setEPS()
postscript("MOL2_Late_TFact_GRN_bhpeaks.eps",
 width = 9, 
    height = 6) 
p1
dev.off()


pdf(file = "MOL2_Late_TFact_GRN_bhpeaks.pdf", 
    width = 9,
    height = 6) 
p1
dev.off()
##########################################################
