##### Libraries  ---------------------------------

#First remove all loaded packages and clean the environment

rm(list=ls())
lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)

#Put a seed to have the same results (pGLS results can vary)
set.seed(2712)

#Load necessary packages
library("caper")
library("ape")
library(dplyr)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(data.table)
library(phytools)
library(ggtree)
library(RColorBrewer)
library("gt")
library(phangorn)


#### Load some useful functions  ---------------------------------

args = commandArgs(trailingOnly=TRUE)


PGLS_pvalue <- function(pgls_rslt) {
  sum_cor <- summary(pgls_rslt)
  pvalue =  formatC(sum_cor$coefficients[8], digits = 3)
  if (pvalue == "   0"){ pvalue = 2.2e-16}
  return(pvalue)
}

PGLS_R2 <- function(pgls_rslt) {
  sum_cor <- summary(pgls_rslt)
  r2 = formatC(sum_cor$r.squared, digits = 2)
  return(r2)
}

PGLS_lambda <- function(pgls_rslt) {
  sum_cor <- summary(pgls_rslt)
  PGLS_lambda = sum_cor$param[2]
  return(PGLS_lambda)
}


GLS_function <- function(gls_rslt) {
  GLS_cc <- coef(gls_rslt)
  GLS_f <- function(x) GLS_cc[1] + GLS_cc[2]*x
  return(GLS_f)
}


split_tibble <- function(tibble, column = 'col') {
  tibble %>% split(., .[,column]) %>% lapply(., function(x) x[,setdiff(names(x),column)])
}




#### Load Species table  ---------------------------------

#Import a species table with genomes and BUSCO informations
species_table <- 
  read.table("Species_table.tsv",
             sep="\t",
             header=TRUE)
species_table <-species_table %>% dplyr::select(-BUSCO_percentage)
species_table <- species_table %>% mutate(BUSCO_percentage = BUSCO_complete/BUSCO_searched)


#Rename Amblyopsis rosae to Troglichthys rosae
species_table[(species_table$Species == "Amblyopsis_rosae"),"Species"] <- "Troglichthys_rosae"


species_table %>%
  group_by(Habitat) %>%
  summarise(n())


#### Load BUSCO analysis and figure  ---------------------------------

species_table_fig <- 
  read.table("Species_table.tsv",
             sep="\t",
             header=TRUE)

species_table_fig[(species_table_fig$Species == "Amblyopsis_rosae"),"Species"] <- "Troglichthys_rosae"
species_table_fig$Species <- gsub("_", " ", species_table_fig$Species)


species_table_fig <- species_table_fig %>% 
  mutate(BUSCO_Complete_single_ps = BUSCO_Complete_single + BUSCO_Pseudogene)

species_table_fig <- 
  as.data.frame(
    species_table_fig %>%
      rowwise() %>%
      mutate(full_name = paste(Species, Assembly, sep=" ")))


species_table_for_fig_long <- 
  as.data.frame(species_table_fig %>%
                  dplyr::select(full_name, BUSCO_Complete_single_ps, BUSCO_Complete_duplicated, 
                                BUSCO_Fragmented, BUSCO_Missing) %>%
                  pivot_longer(!full_name, names_to = "category", values_to = "count"))



species_table_for_fig_long$category <-
  factor(species_table_for_fig_long$category,
         levels=rev(c("BUSCO_Complete_single_ps", "BUSCO_Complete_duplicated", 
                      "BUSCO_Fragmented", "BUSCO_Missing")))



pdf(file = "BUSCO_graphic.pdf",width = 8.34,  height = 4.61)

species_table_for_fig_long %>%
  ggplot(., aes(x=full_name, y=count, fill=category)) +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(axis.text=element_text(size=6),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") +
  scale_fill_manual(values =
                      c("BUSCO_Complete_single_ps"="#648FFF",
                        "BUSCO_Complete_duplicated"="#785EF0", 
                        "BUSCO_Fragmented"="#FFB000",
                        "BUSCO_Missing"="#DC267F")) +
  xlab("") +
  ylab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dev.off()



#Now lets look at the pseudogene percentage among BUSCO genes

species_table_for_fig_long_pseudo <- 
  as.data.frame(species_table_fig %>%
                  dplyr::select(full_name, BUSCO_complete, BUSCO_Pseudogene) %>%
                  pivot_longer(!full_name, names_to = "category", values_to = "count"))



species_table_for_fig_long_pseudo$category <-
  factor(species_table_for_fig_long_pseudo$category,
         levels=rev(c("BUSCO_complete", "BUSCO_Pseudogene")))


pdf(file = "BUSCO_graphic_pseudo.pdf",width = 8.34,  height = 4.61)

species_table_for_fig_long_pseudo %>%
  ggplot(., aes(x=full_name, y=count, fill=category)) +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(axis.text=element_text(size=6),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") +
  scale_fill_manual(values =
                      c("BUSCO_Pseudogene"="#000000",
                        "BUSCO_complete"="#009E73")) +
  xlab("") +
  ylab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dev.off()



#### Load Species tree  ---------------------------------

individual_tree_AA <- read.tree("AMAS_concatenated_alignment.prot.aln.treefile")
individual_tree_NT <- read.tree("AMAS_concatenated_alignment.cds.aln.treefile")


individual_tree_AA_rooted <- 
  root(individual_tree_AA, 
       outgroup=c("Percopsis_omiscomaycus.IRGN_VKOYSXN1ZJ",
                  "Percopsis_omiscomaycus.IRGN_B4MP9LQ3H0",
                  "Percopsis_omiscomaycus.TJN_237",
                  "Percopsis_transmontana.KU_KUI_29775"),
       resolve.root = TRUE)

individual_tree_NT_rooted <- 
  root(individual_tree_NT, 
       outgroup=c("Percopsis_omiscomaycus.IRGN_VKOYSXN1ZJ",
                  "Percopsis_omiscomaycus.IRGN_B4MP9LQ3H0",
                  "Percopsis_omiscomaycus.TJN_237",
                  "Percopsis_transmontana.KU_KUI_29775"),
       resolve.root = TRUE)



#change name of A.rosae
individual_tree_AA_rooted$tip.label[individual_tree_AA_rooted$tip.label=="Amblyopsis_rosae.C5"] <- "Troglichthys_rosae.C5"
individual_tree_AA_rooted$tip.label[individual_tree_AA_rooted$tip.label=="Amblyopsis_rosae.C3"] <- "Troglichthys_rosae.C3"
individual_tree_AA_rooted$tip.label[individual_tree_AA_rooted$tip.label=="Amblyopsis_rosae.L1"] <- "Troglichthys_rosae.L1"

individual_tree_NT_rooted$tip.label[individual_tree_NT_rooted$tip.label=="Amblyopsis_rosae.C5"] <- "Troglichthys_rosae.C5"
individual_tree_NT_rooted$tip.label[individual_tree_NT_rooted$tip.label=="Amblyopsis_rosae.C3"] <- "Troglichthys_rosae.C3"
individual_tree_NT_rooted$tip.label[individual_tree_NT_rooted$tip.label=="Amblyopsis_rosae.L1"] <- "Troglichthys_rosae.L1"

ggtree(individual_tree_AA_rooted, branch.length="none") + geom_tiplab() + xlim(0, 20)
ggtree(individual_tree_NT_rooted, branch.length="none") + geom_tiplab() + xlim(0, 20)


#Merge individuals to form a species tree

species_tree <- 
  keep.tip(individual_tree_NT_rooted,
           c("Amblyopsis_hoosieri.MLN_0242",
             "Troglichthys_rosae.C3",
             "Aphredoderus_sayanus.ASAY_01",
             "Percopsis_omiscomaycus.TJN_237",
             "Percopsis_transmontana.KU_KUI_29775",
             "Chologaster_cornuta.CCOR_08",
             "Typhlichthys_eigenmanni.EN1",
             "Speoplatyrhinus_poulsoni.IRGN_2EJSDHOURK",
             "Typhlichthys_subterraneus.MLN_0051",
             "Forbesichthys_agassizii.UTTC_242",
             "Forbesichthys_papilliferus.FPAP_01",
             "Amblyopsis_spelaea.YFTC_23868"))
species_tree$tip.label[species_tree$tip.label=="Amblyopsis_hoosieri.MLN_0242"] <- "Amblyopsis_hoosieri"
species_tree$tip.label[species_tree$tip.label=="Troglichthys_rosae.C3"] <- "Troglichthys_rosae"
species_tree$tip.label[species_tree$tip.label=="Aphredoderus_sayanus.ASAY_01"] <- "Aphredoderus_sayanus"
species_tree$tip.label[species_tree$tip.label=="Percopsis_omiscomaycus.TJN_237"] <- "Percopsis_omiscomaycus"
species_tree$tip.label[species_tree$tip.label=="Percopsis_transmontana.KU_KUI_29775"] <- "Percopsis_transmontana"
species_tree$tip.label[species_tree$tip.label=="Chologaster_cornuta.CCOR_08"] <- "Chologaster_cornuta"
species_tree$tip.label[species_tree$tip.label=="Typhlichthys_eigenmanni.EN1"] <- "Typhlichthys_eigenmanni"
species_tree$tip.label[species_tree$tip.label=="Speoplatyrhinus_poulsoni.IRGN_2EJSDHOURK"] <- "Speoplatyrhinus_poulsoni"
species_tree$tip.label[species_tree$tip.label=="Typhlichthys_subterraneus.MLN_0051"] <- "Typhlichthys_subterraneus"
species_tree$tip.label[species_tree$tip.label=="Forbesichthys_agassizii.UTTC_242"] <- "Forbesichthys_agassizii"
species_tree$tip.label[species_tree$tip.label=="Forbesichthys_papilliferus.FPAP_01"] <- "Forbesichthys_papilliferus"
species_tree$tip.label[species_tree$tip.label=="Amblyopsis_spelaea.YFTC_23868"] <- "Amblyopsis_spelaea"


ggtree(species_tree, branch.length="none") + geom_tiplab() + xlim(0, 20)


species_tree <- makeNodeLabel(species_tree, method="number", prefix="Node")

species_tree_chase <- 
  read.nexus("Chase_CAVE_NODE_123SUM.tre")
species_tree_chase <- drop.tip(species_tree_chase, "Typhlichthys_subterraneus_SAMN14308618")
species_tree_chase$tip.label[species_tree_chase$tip.label=="Typhlichthys_eigenmanni_SAMN14308592"] <- "Typhlichthys_eigenmanni"
species_tree_chase$tip.label[species_tree_chase$tip.label=="Typhlichthys_subterraneus_NCBI_96345"] <- "Typhlichthys_subterraneus"
species_tree_chase$tip.label[species_tree_chase$tip.label=="Percopsis_transmontana_NCBI_7149"] <- "Percopsis_transmontana"
species_tree_chase$tip.label[species_tree_chase$tip.label=="Forbesichthys_papilliferus_SAMN14308584"] <- "Forbesichthys_papilliferus"
species_tree_chase$tip.label[species_tree_chase$tip.label=="Forbesichthys_agassizii_SAMN14308580"] <- "Forbesichthys_agassizii"
species_tree_chase$tip.label[species_tree_chase$tip.label=="Chologaster_cornuta_SAMN14308572"] <- "Chologaster_cornuta"
species_tree_chase$tip.label[species_tree_chase$tip.label=="Amblyopsis_spelaea_SAMN14308568"] <- "Amblyopsis_spelaea"
species_tree_chase$tip.label[species_tree_chase$tip.label=="Amblyopsis_hoosieri_SAMN20209475"] <- "Amblyopsis_hoosieri"
species_tree_chase$tip.label[species_tree_chase$tip.label=="Amblyopsis_rosae_SAMN14308599"] <- "Troglichthys_rosae"
species_tree_chase$tip.label[species_tree_chase$tip.label=="Speoplatyrhinus_poulsoni_SAMN14308586"] <- "Speoplatyrhinus_poulsoni"
species_tree_chase$tip.label[species_tree_chase$tip.label=="Percopsis_omiscomycus"] <- "Percopsis_omiscomaycus"

#### Load Vision sequence table  ---------------------------------

vision_seq_df <- 
  read.table("Vision_sequences_table.tsv",
             sep="\t",
             header=FALSE)
colnames(vision_seq_df) <- 
  c("Gene_clade","Species","gene_name","Genomic_position","CDS_length","Gene_Type","Note",
    "Exon_Count","Sequence")

vision_seq_df <- 
  vision_seq_df %>% filter(Species != "") %>%
  dplyr::select(Gene_clade, Species, gene_name, Genomic_position,
                CDS_length, Gene_Type, Note, Exon_Count, Sequence)


vision_seq_df <- 
  vision_seq_df %>%
  mutate(Gene_Type_Simp = case_when(
    Gene_Type == "Complete" ~ "Complete_Incomplete",
    Gene_Type == "Incomplete" ~ "Complete_Incomplete",
    Gene_Type == "Pseudogene" ~ "Pseudogene"
  ))



vision_seq_df[(vision_seq_df$Species == "Amblyopsis_rosae"),"Species"] <- "Troglichthys_rosae"

#### Load LoF Tables  ---------------------------------

Lof_stop <- 
  read.table("Stop_codon_table.tsv",
             sep="\t",
             header=TRUE)
Lof_stop <- Lof_stop %>% filter(Species != "")
Lof_stop <- Lof_stop %>%
  mutate(Position_CDS = Position_AA*3)

Lof_FS <- 
  read.table("Frameshift_table.tsv",
             sep="\t",
             header=TRUE)

Lof_other <- 
  read.table("OtherLoF_table.tsv",
             sep="\t",
             header=TRUE)
colnames(Lof_other) <- c("Species", "Pseudogene", "LoF_Type", "PositionAA")

Lof_other <- Lof_other %>%
  mutate(Position_CDS = PositionAA*3)




Table_all_LoF <- 
  rbind(
    Lof_stop %>% dplyr::select(Species, Pseudogene, LoF_Type, Position_CDS),
    Lof_FS %>% dplyr::select(Species, Pseudogene, LoF_Type, Position_CDS),
    Lof_other %>% dplyr::select(Species, Pseudogene, LoF_Type, Position_CDS)
  )


Table_all_LoF <- 
  Table_all_LoF %>%
  mutate(LoF_name_temp = paste(Pseudogene,LoF_Type, sep="_")) %>%
  mutate(LoF_name = paste(LoF_name_temp,Position_CDS, sep="_")) %>%
  dplyr::select(-LoF_name_temp)





Table_all_LoF[(Table_all_LoF$Species == "Amblyopsis_rosae"),"Species"] <- "Troglichthys_rosae"
Lof_other[(Lof_other$Species == "Amblyopsis_rosae"),"Species"] <- "Troglichthys_rosae"
Lof_FS[(Lof_FS$Species == "Amblyopsis_rosae"),"Species"] <- "Troglichthys_rosae"
Lof_stop[(Lof_stop$Species == "Amblyopsis_rosae"),"Species"] <- "Troglichthys_rosae"



#Verify that the LoF table and Sequence table are similar in term of pseudogene nb


as.data.frame(
  vision_seq_df %>%
    filter(Gene_Type_Simp == "Pseudogene") %>%
    group_by(Species) %>%
    summarise(countP = n())
)


as.data.frame(
  Table_all_LoF %>%
    dplyr::select(Species, Pseudogene) %>%
    distinct() %>%
    group_by(Species) %>%
    summarise(countP = n())
)



#Count the number of LoF mutation per species

Table_all_LoF %>%
  group_by(Species) %>%
  summarise(count = n()) %>%
  arrange(desc(count))


#### Load IGV Table ---------------------------------

IGV_df <- 
  read.table("IGV_verif_table.tsv",
             sep="\t",
             header=TRUE)

IGV_df[(IGV_df$Species == "Amblyopsis_rosae"),"Species"] <- "Troglichthys_rosae"
IGV_df[(IGV_df$Individual == "UAIC_14148"),"Individual"] <- "UAIC_14148_01"
IGV_df[(IGV_df$Individual == "UTC_243"),"Individual"] <- "UTTC_243"
IGV_df[(IGV_df$Individual == "UTC_242"),"Individual"] <- "UTTC_242"
IGV_df[(IGV_df$Individual == "UTC_649"),"Individual"] <- "UTTC_649"




IGV_df <- 
  IGV_df %>%
  mutate(LoF_name = paste(LoF_Type, Pos_CDS, sep="_")) %>%
  mutate(LoF_gene_name = paste(Gene, LoF_name, sep="_"))


IGV_df_long <- 
  as.data.frame(
    IGV_df %>%
      pivot_longer(!c(Species, Gene, LoF_Type, Pos_CDS,LoF_Pos, Note, LoF_gene_name, GeneSet,
                      Individual, LoF_name, LoF_gene_name, Detection.method, Description.of.LoF.detected.on.IGV), 
                   names_to = "Support", values_to = "count")
  )



IGV_df_long$Support <-
  factor(IGV_df_long$Support ,
         levels=c("Support_Reads", "Non_support_Reads"))



# Assign a genotype to each individual, following a binomial probability (p = 0.5) 

IGV_df <- 
  IGV_df %>% 
  rowwise() %>%
  mutate(Total_Reads = Non_support_Reads + Support_Reads) %>%
  mutate(Non_support_Reads_prop = Non_support_Reads/Total_Reads) %>%
  mutate(Support_Reads_prop = Support_Reads/Total_Reads) %>% 
  mutate(Minimum_read_number = min(Support_Reads, Non_support_Reads)) %>%
  mutate(binom_proba = pbinom(Minimum_read_number, Total_Reads , p = 0.5)) %>%
  mutate(Genotype = case_when(
    binom_proba >= 0.05 & binom_proba <= 0.95 & Non_support_Reads != 0 & Support_Reads != 0 ~ "Heterozygous",
    binom_proba >= 0.05 & binom_proba <= 0.95 & Non_support_Reads == 0  & Support_Reads != 0 ~ "LoF_Homozygous",
    binom_proba >= 0.05 & binom_proba <= 0.95 & Non_support_Reads != 0  & Support_Reads == 0 ~ "no_LoF_Homozygous",
    binom_proba < 0.05 & (Support_Reads > Non_support_Reads) ~ "LoF_Homozygous",
    binom_proba > 0.95 & (Support_Reads > Non_support_Reads) ~ "LoF_Homozygous",
    binom_proba < 0.05 & (Support_Reads < Non_support_Reads) ~ "no_LoF_Homozygous",
    binom_proba > 0.95 & (Support_Reads < Non_support_Reads) ~ "no_LoF_Homozygous",
  )) %>%
  mutate(LoF_genotype = case_when(
    Genotype == "LoF_Homozygous" ~ 2,
    Genotype == "Heterozygous"  ~ 1,
    Genotype == "no_LoF_Homozygous" ~ 0
  )) %>%
  mutate(non_LoF_genotype = case_when(
    Genotype == "LoF_Homozygous" ~ 0,
    Genotype == "Heterozygous"  ~ 1,
    Genotype == "no_LoF_Homozygous" ~ 2
  ))  %>% 
  mutate(total_allele_nb = LoF_genotype + non_LoF_genotype)  

IGV_df <- as.data.frame(IGV_df)

#rename splice mutations to add the intron number

IGV_df[(IGV_df$LoF_gene_name == "Cryba4_Splice site_NA"),"LoF_gene_name"] <- "Cryba4_Splice_site_2"
IGV_df[(IGV_df$LoF_gene_name == "Cryba1b_Splice site_NA"),"LoF_gene_name"] <- "Cryba1b_Splice_site_4"
IGV_df[(IGV_df$LoF_gene_name == "Opn6a-2_Splice site_NA"),"LoF_gene_name"] <- "Opn6a-2_Splice_site_5"
IGV_df[(IGV_df$LoF_gene_name == "pde6hb_Splice site_NA"),"LoF_gene_name"] <- "pde6hb_Splice_site_1"
IGV_df[(IGV_df$LoF_gene_name == "gnat2_Splice site_NA"),"LoF_gene_name"] <- "gnat2_Splice_site_6"

IGV_df[(IGV_df$LoF_gene_name == "Rh2.1_Start loss_NA"),"LoF_gene_name"] <- "Rh2.1_Start_loss"
IGV_df[(IGV_df$LoF_gene_name == "Rh2.1_Start loss_3"),"LoF_gene_name"] <- "Rh2.1_Start_loss"
IGV_df[(IGV_df$LoF_gene_name == "Crygn2_Start loss_3"),"LoF_gene_name"] <- "Crygn2_Start_loss"
IGV_df[(IGV_df$LoF_gene_name == "gcap3_Start loss_3"),"LoF_gene_name"] <- "gcap3_Start_loss"


IGV_df <- IGV_df %>% distinct()

#put the most pasimonious state for genotype equal to NA

IGV_df %>% filter(is.na(Genotype))

IGV_df <- 
  IGV_df %>%
  mutate(Genotype = ifelse(Individual == "MLN_0296" & LoF_gene_name == "Lws1_STOP_150", "no_LoF_Homozygous", Genotype))

IGV_df <- 
  IGV_df %>%
  mutate(Genotype = ifelse(Individual == "MLN_0296" & LoF_gene_name == "gcap3_Start_loss", "LoF_Homozygous", Genotype))

IGV_df <- 
  IGV_df %>%
  mutate(Genotype = ifelse(Individual == "MLN_0296" & LoF_gene_name == "Lws1_Del_316", "LoF_Homozygous", Genotype))


#write table to file
write.table(IGV_df,
            "IGV_df.tsv",
            sep="\t",
            quote=FALSE,
            col.names=TRUE,
            row.names=FALSE)



# Compute LoF allele frequency


IGV_df_freq <- 
  as.data.frame(IGV_df %>%
                  group_by(Species, GeneSet, LoF_gene_name) %>%
                  summarise(sum_tot_allele = sum(total_allele_nb, na.rm = TRUE),
                            sum_LoF_allele = sum(LoF_genotype, na.rm = TRUE))) %>%
  mutate(LoF_freq = sum_LoF_allele/sum_tot_allele)


write.table(IGV_df_freq,
            "IGV_df_freq.tsv",
            sep="\t",
            quote=FALSE,
            col.names=TRUE,
            row.names=FALSE)




#Verify that all LoF present in the LoF tables are present on the IGV table

igv_count <- 
  IGV_df %>%
  filter(Detection.method == "Detected on Assembly") %>%
  dplyr::select(Species, LoF_gene_name) %>%
  distinct() %>%
  group_by(Species) %>%
  summarise(count_igv_table = n())


table_lof_count <- 
  Table_all_LoF %>%
  dplyr::select(Species, LoF_name) %>%
  distinct() %>%
  group_by(Species) %>%
  summarise(count_lof_table = n())


left_join(igv_count, table_lof_count, by='Species')



#Add LoF identified on IGV to the LoF Table

IGV_df %>%
  filter(Detection.method != "Detected on Assembly")  %>%
  dplyr::select(Species, Gene, Pos_CDS, Description.of.LoF.detected.on.IGV) %>%
  distinct()


temp_df <- 
  rbind(
    cbind("Typhlichthys_subterraneus", "Opn6a-1", "Del", "664", "Opn6a-1_Del_664"),
    cbind("Typhlichthys_subterraneus", "Rpe65b", "Ins", "723", "Rpe65b_Ins_723"),
    cbind("Typhlichthys_subterraneus", "arr3a", "STOP", "879", "arr3a_STOP_879"),
    cbind("Typhlichthys_subterraneus", "gnat2", "STOP", "213", "gnat2_STOP_213"),
    cbind("Typhlichthys_subterraneus", "gnat2", "Ins", "375", "gnat2_Ins_375"),
    cbind("Typhlichthys_subterraneus", "grk7a", "Del", "183", "grk7a_Del_183"),
    cbind("Typhlichthys_subterraneus", "grk1b", "STOP", "1512", "grk1b_STOP_1512")
  )

temp_df <- as.data.frame(temp_df)
colnames(temp_df) <- colnames(Table_all_LoF)
temp_df$Position_CDS <- as.numeric(temp_df$Position_CDS)

Table_all_LoF_igv <- rbind(Table_all_LoF, temp_df)

temp_df <- 
  rbind(
    cbind("Typhlichthys_subterraneus", "Opn6a-1", "Del", "2","664"),
    cbind("Typhlichthys_subterraneus", "Rpe65b", "Ins", "2","723"),
    cbind("Typhlichthys_subterraneus", "gnat2", "Ins", "5","375"),
    cbind("Typhlichthys_subterraneus", "grk7a", "Del", "5","183")
  )
temp_df <- as.data.frame(temp_df)
colnames(temp_df) <- colnames(Lof_FS)
temp_df$Position_CDS <- as.numeric(temp_df$Position_CDS)
temp_df$Nb_nuc <- as.numeric(temp_df$Nb_nuc)


Lof_FS_igv <- rbind(Lof_FS, temp_df)



#Remove LoF not retrieved in any specimen from the LoF Table

IGV_df %>%
  group_by(LoF_gene_name) %>%
  summarize(never_found = all(! Genotype %in% c("LoF_Homozygous", "Heterozygous"))) %>%
  filter(never_found) %>%
  pull(LoF_gene_name)


Table_all_LoF_igv <- 
  Table_all_LoF_igv %>%
  filter(! LoF_name %in% c("Exorh_STOP_708", "Opn6a-2_STOP_372", 
                           "gngt2b_STOP_192", "Cryba4_Splice site_6", "grk7a_Del_388"))


#Also remove the deletion in grk7a never observed from the table listing frameshifts
Lof_FS_igv <- 
  Lof_FS_igv %>%
  filter(! Position_CDS == "388")



#### Compute stats on observed LoF mutations (corrected with IGV analysis) ---------------------------------


### What are the number and type of observed LoF ? 

Table_all_LoF_igv_uniq <- 
  Table_all_LoF_igv %>%
  dplyr::select(Pseudogene, LoF_Type, Position_CDS) %>%
  distinct()


LoF_type_summary <- 
  as.data.frame(Table_all_LoF_igv_uniq %>%
                  group_by(LoF_Type) %>%
                  summarise(count = n()))

pdf(file = "LoF_Type_Pie.pdf",width = 8.34,  height = 4.61)

LoF_type_summary %>%
  ggplot(., aes(x="", y=count, fill=LoF_Type)) +
  geom_bar(width = 1, stat = "identity") + 
  coord_polar("y", start=0) +
  scale_fill_manual(values = c(
    "Ins" = "#000000",
    "Del" = "gray",
    "STOP" = "#D55E00",
    "Splice site" = "#CC79A7",
    "Start loss" = "#009E73",
    "Stop loss" = "#0072B2"
  )) + 
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    axis.text.x=element_blank(),
    plot.title=element_text(size=14, face="bold")
  ) +
  theme(legend.position="none") 

dev.off()


### What is the distribution of indel sizes ? 


Lof_FS_uniq <- 
  Lof_FS_igv %>%
  dplyr::select(Pseudogene, LoF_Type, Nb_nuc, Position_CDS) %>%
  distinct()


pdf(file = "Frameshift_sizes.pdf",width = 11.34,  height = 4.61)

Lof_FS_uniq %>%
  group_by(LoF_Type, Nb_nuc) %>%
  summarise(count = n()) %>%
  ggplot(., aes(x=Nb_nuc, y=count, fill=LoF_Type)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c(
    "Ins" = "#000000",
    "Del" = "gray"
  )) + 
  theme_classic() +
  theme(legend.position="none") +
  xlab("Number of nucleotides") +
  ylab("count") +
  scale_x_continuous(breaks=c(0,1,2,4,5,7,8,13,19))


dev.off()



#### Make a summary table of vision genes in genome assemblies ---------------------------------

summary_vision_long <- 
  vision_seq_df %>%
  group_by(Species, Gene_Type) %>%
  summarise(count = n())

summary_vision_wide <- 
  as.data.frame(
    summary_vision_long %>%
      pivot_wider(names_from = Gene_Type, values_from = count, values_fill = 0)
  )



summary_vision_wide <- 
  summary_vision_wide %>%
  mutate(Total = Complete + Incomplete + Pseudogene ) %>%
  mutate(Prop_pseudo = Pseudogene / Total)

gt(summary_vision_wide %>% arrange(Prop_pseudo))



#add species info

summary_vision_wide_info <- 
  left_join(summary_vision_wide, species_table, by="Species")

## look at the table

summary_vision_wide_info %>%
  dplyr::select(Species, Complete, Incomplete, Pseudogene, Prop_pseudo, Habitat) %>%
  arrange(Prop_pseudo)



## Add the number of nucleotides and the number of LoF  

LoF_per_sp <- 
  as.data.frame(Table_all_LoF %>%
                  group_by(Species) %>%
                  summarise(LoF_nb = n()) %>%
                  arrange(desc(LoF_nb)))

nt_per_sp <- 
  as.data.frame(vision_seq_df %>%
                  group_by(Species) %>%
                  summarise(nt_nb = sum(CDS_length)))

summary_vision_wide_info_NT <- 
  left_join(summary_vision_wide_info, nt_per_sp, by="Species") 
summary_vision_wide_info_NT_LOF <- 
  left_join(summary_vision_wide_info_NT, LoF_per_sp, by="Species") 


summary_vision_wide_info_NT_LOF <- 
  summary_vision_wide_info_NT_LOF %>%
  mutate(LoF_per_nt = LoF_nb/nt_nb)

summary_vision_wide_info_NT_LOF[is.na(summary_vision_wide_info_NT_LOF)] <- 0


write.table(summary_vision_wide_info_NT_LOF,
            "summary_vision_wide_info_NT_LOF.tsv",
            sep="\t",
            quote=FALSE,
            col.names=TRUE,
            row.names=FALSE)


gt(summary_vision_wide_info_NT_LOF %>%
     dplyr::select(Species, Complete, Incomplete, Pseudogene, Total, Prop_pseudo, nt_nb, LoF_nb, LoF_per_nt) %>%
     arrange(Prop_pseudo))




#### Load dn/ds values -- At the species level and alignments without the outgroup G. morhua  ---------------------------------

vision_dNdS <- 
  read.table("Vision_genes_dNdS.csv",
             sep=",",
             header=FALSE)
colnames(vision_dNdS) <- 
  c("Gene_clade", "species_gene_name", "LB", "MLE", "UB", "dN", "dS")

#Remove bad estimates of omega. dS or dN > 1 => saturation; ds < 0.01 = unreliable estimation of omega
vision_dNdS <- 
  vision_dNdS %>%
  filter(dN < 1) %>%
  filter(dS < 1) %>%
  filter(dS > 0.01)


list_species <-
  vision_seq_df %>%
  pull(Species)  %>%
  unique()


#Add info about complete or incomplete 


all_gene_name <- vision_dNdS %>% pull(species_gene_name) %>% unique()

names_Pseudo <- all_gene_name[grepl("_Pseudogene", all_gene_name)]
names_Complete <- all_gene_name[grepl("_Complete", all_gene_name)]
names_Incomplete <- all_gene_name[grepl("_Incomplete", all_gene_name)]

vision_dNdS <- 
  vision_dNdS %>%
  filter(species_gene_name %in% c(names_Pseudo, names_Complete, names_Incomplete)) %>%
  mutate(Gene_Type = case_when(
    species_gene_name %in% names_Pseudo ~ "Pseudogene",
    species_gene_name %in% names_Complete ~ "Complete",
    species_gene_name %in% names_Incomplete ~ "Incomplete"
  ))


#Assign a species to every gene name

vision_dNdS <- 
  vision_dNdS %>% 
  mutate(across('species_gene_name', str_replace, 'Amblyopsis_rosae', 'Troglichthys_rosae'))


all_gene_name <- vision_dNdS %>% pull(species_gene_name) %>% unique()

vision_dNdS_df <- as.data.frame(NULL)
for(curr_species in list_species){
  
  curr_sp_gene_names <- 
    all_gene_name[grepl(curr_species, all_gene_name)]
  
  curr_sp_table <- 
    vision_dNdS %>%
    filter(species_gene_name %in% curr_sp_gene_names) %>%
    mutate(Species = curr_species)
  
  vision_dNdS_df <- 
    rbind(vision_dNdS_df,
          curr_sp_table)
  
  
}


#### Load dn/ds values -- At the species level and alignments with the outgroup G. morhua ---------------------------------

vision_dNdS_gadus <- 
  read.table("Vision_genes_dNdS.Gadus.csv",
             sep=",",
             header=FALSE)
colnames(vision_dNdS_gadus) <- 
  c("Gene_clade", "species_gene_name", "LB", "MLE", "UB", "dN", "dS")

#Remove bad estimates of omega
vision_dNdS_gadus <- 
  vision_dNdS_gadus %>%
  filter(dN < 1) %>%
  filter(dS < 1) %>%
  filter(dS > 0.01)


list_species <-
  vision_seq_df %>%
  pull(Species)  %>%
  unique()


#Add info about complete or incomplete 


all_gene_name <- vision_dNdS_gadus %>% pull(species_gene_name) %>% unique()

names_Pseudo <- all_gene_name[grepl("_Pseudogene", all_gene_name)]
names_Complete <- all_gene_name[grepl("_Complete", all_gene_name)]
names_Incomplete <- all_gene_name[grepl("_Incomplete", all_gene_name)]

vision_dNdS_gadus <- 
  vision_dNdS_gadus %>%
  filter(species_gene_name %in% c(names_Pseudo, names_Complete, names_Incomplete)) %>%
  mutate(Gene_Type = case_when(
    species_gene_name %in% names_Pseudo ~ "Pseudogene",
    species_gene_name %in% names_Complete ~ "Complete",
    species_gene_name %in% names_Incomplete ~ "Incomplete"
  ))


#Assign a species to every gene name

vision_dNdS_gadus <- 
  vision_dNdS_gadus %>% 
  mutate(across('species_gene_name', str_replace, 'Amblyopsis_rosae', 'Troglichthys_rosae'))


all_gene_name <- vision_dNdS_gadus %>% pull(species_gene_name) %>% unique()

list_species <- c(list_species, "Gadus_morhua")

vision_dNdS_gadus_df <- as.data.frame(NULL)
for(curr_species in list_species){
  
  curr_sp_gene_names <- 
    all_gene_name[grepl(curr_species, all_gene_name)]
  
  curr_sp_table <- 
    vision_dNdS_gadus %>%
    filter(species_gene_name %in% curr_sp_gene_names) %>%
    mutate(Species = curr_species)
  
  vision_dNdS_gadus_df <- 
    rbind(vision_dNdS_gadus_df,
          curr_sp_table)
  
  
}



#### Analysis of dN/dS values - Sequences from assembly (species level) ---------------------------------


vision_dNdS_df <- 
  left_join(vision_dNdS_df, 
            species_table %>% dplyr::select(Species, Habitat) , 
            by="Species")


pdf(file = "dnds_Vision_boxplot.pdf",width = 8.34,  height = 4.61)

vision_dNdS_df %>%
  ggplot(., aes(x=reorder(Species, MLE, mean), y=MLE, fill=Habitat)) +
  geom_boxplot() +
  scale_fill_manual(values = c("troglomorphic_subterranean" = "#FD914E", 
                               "Surface" = "#7DA0FB")) +
  ylab("dN/dS") +
  xlab("Species") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  theme(legend.position="none")

dev.off()



pdf(file = "dnds_Vision_violinplot.pdf",width = 8.34,  height = 4.61)

vision_dNdS_df %>%
  ggplot(., aes(x=reorder(Species, MLE, mean), y=MLE, fill=Habitat)) +
  geom_violin() +
  geom_jitter(shape=16, position=position_jitter(0.2), color="gray") +
  geom_boxplot(width=0.1) + 
  scale_fill_manual(values = c("troglomorphic_subterranean" = "#FD914E", 
                               "Surface" = "#7DA0FB")) +
  ylab("dN/dS") +
  xlab("Species") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  theme(legend.position="none")

dev.off()





vision_dNdS_df_mean <- 
  as.data.frame(
    vision_dNdS_df %>%
      group_by(Species) %>%
      summarise(mean_MLE = mean(MLE),
                nb_genes = n())
  ) %>%
  arrange(mean_MLE)


#pGLS between mean dN/dS and habitat


vision_dNdS_df_mean_info <- 
  left_join(vision_dNdS_df_mean,
            summary_vision_wide_info)

caper_dNdS_pseudo <- 
  comparative.data(phy = species_tree, 
                   data = vision_dNdS_df_mean_info,
                   names.col = Species, vcv = TRUE,
                   na.omit = FALSE, warn.dropped = TRUE)
caper_dNdS_pseudo_chase <- 
  comparative.data(phy = species_tree_chase, 
                   data = vision_dNdS_df_mean_info,
                   names.col = Species, vcv = TRUE,
                   na.omit = FALSE, warn.dropped = TRUE)



#Launch pGLS
dNdS_vs_Habitat <-
  pgls(mean_MLE ~ Habitat, 
       data = caper_dNdS_pseudo, 
       lambda = "ML")
summary(dNdS_vs_Habitat)

#Launch pGLS with an alternative tree
dNdS_vs_Habitat_chase <-
  pgls(mean_MLE ~ Habitat, 
       data = caper_dNdS_pseudo_chase, 
       lambda = "ML")
summary(dNdS_vs_Habitat_chase)

#### Analysis of dN/dS values - Same analysis but with alignments containing Gadus morhua ---------------------------------


vision_dNdS_gadus_df <- 
  left_join(vision_dNdS_gadus_df, 
            species_table %>% dplyr::select(Species, Habitat) , 
            by="Species")

vision_dNdS_gadus_df[(vision_dNdS_gadus_df$Species == "Gadus_morhua"),"Habitat"] <- "Surface"


#Check differences with the alignment which dont contain Gadus morhua

new_values <- vision_dNdS_gadus_df %>% dplyr::select(species_gene_name, MLE) 
colnames(new_values) <- c("gene_name", "MLE_new")
old_values <- vision_dNdS_df %>% dplyr::select(species_gene_name, MLE) 
colnames(old_values) <- c("gene_name", "MLE_old")
Old_vs_new_dn_ds <- left_join(old_values, new_values, by="gene_name")

Old_vs_new_dn_ds %>%
  ggplot(., aes(x=MLE_new, y=MLE_old)) + 
  geom_point() +
  xlab("new dN/dS") +
  ylab("old dN/dS")




pdf(file = "dnds_Vision_boxplot.gadus.pdf",width = 8.34,  height = 4.61)

vision_dNdS_gadus_df %>%
  ggplot(., aes(x=reorder(Species, MLE, mean), y=MLE, fill=Habitat)) +
  geom_boxplot() +
  scale_fill_manual(values = c("troglomorphic_subterranean" = "#FD914E", 
                               "Surface" = "#7DA0FB")) +
  ylab("dN/dS") +
  xlab("Species") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  theme(legend.position="none")

dev.off()



pdf(file = "dnds_Vision_violinplot.gadus.pdf",width = 8.34,  height = 4.61)

vision_dNdS_gadus_df %>%
  ggplot(., aes(x=reorder(Species, MLE, mean), y=MLE, fill=Habitat)) +
  geom_violin() +
  #geom_jitter(shape=16, position=position_jitter(0.2), color="gray") +
  geom_jitter(shape=16, position=position_jitter(0.2), aes(color= Gene_Type)) +
  scale_color_manual(values = c("Complete" = "gray", 
                                "Incomplete" = "gray",
                                "Pseudogene" = "black")) +
  geom_boxplot(width=0.1, outlier.shape = NA) + 
  scale_fill_manual(values = c("troglomorphic_subterranean" = "#FD914E", 
                               "Surface" = "#7DA0FB")) +
  ylab("dN/dS") +
  xlab("Species") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  theme(legend.position="none")

dev.off()





vision_dNdS_df_mean <- 
  as.data.frame(
    vision_dNdS_df %>%
      group_by(Species) %>%
      summarise(mean_MLE = mean(MLE),
                nb_genes = n())
  ) %>%
  arrange(mean_MLE)


#pGLS between mean dN/dS and habitat


vision_dNdS_df_mean_info <- 
  left_join(vision_dNdS_df_mean,
            summary_vision_wide_info)

caper_dNdS_pseudo <- 
  comparative.data(phy = species_tree, 
                   data = vision_dNdS_df_mean_info,
                   names.col = Species, vcv = TRUE,
                   na.omit = FALSE, warn.dropped = TRUE)
caper_dNdS_pseudo_chase <- 
  comparative.data(phy = species_tree_chase, 
                   data = vision_dNdS_df_mean_info,
                   names.col = Species, vcv = TRUE,
                   na.omit = FALSE, warn.dropped = TRUE)





dNdS_vs_Habitat <-
  pgls(mean_MLE ~ Habitat, 
       data = caper_dNdS_pseudo, 
       lambda = "ML")
summary(dNdS_vs_Habitat)

dNdS_vs_Habitat_chase <-
  pgls(mean_MLE ~ Habitat, 
       data = caper_dNdS_pseudo_chase, 
       lambda = "ML")
summary(dNdS_vs_Habitat_chase)


#### pGLS - Habitat vs pseudo proportion ---------------------------------

caper_vision <- 
  comparative.data(phy = species_tree, 
                   data = summary_vision_wide_info_NT_LOF,
                   names.col = Species, vcv = TRUE,
                   na.omit = FALSE, warn.dropped = TRUE)

caper_vision_chase <- 
  comparative.data(phy = species_tree_chase, 
                   data = summary_vision_wide_info_NT_LOF,
                   names.col = Species, vcv = TRUE,
                   na.omit = FALSE, warn.dropped = TRUE)


Pseudoprop_vs_Habitat <-
  pgls(Prop_pseudo ~ Habitat, 
       data = caper_vision, 
       lambda = "ML")
summary(Pseudoprop_vs_Habitat)


Pseudoprop_vs_Habitat_chase <-
  pgls(Prop_pseudo ~ Habitat, 
       data = caper_vision_chase, 
       lambda = "ML")
summary(Pseudoprop_vs_Habitat)


pdf(file = "PseudoProp_Violinplot_diploid.pdf",width = 8.34,  height = 4.61)

summary_vision_wide_info_NT_LOF %>%
  ggplot(., aes(x=Habitat, y=Prop_pseudo, fill=Habitat)) +
  geom_violin() +
  geom_jitter(size=2, shape=16, color="black", alpha=0.8, position=position_jitter(0.2)) +
  theme_classic() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") +
  ylab("Proportion of pseudogenes") +
  scale_fill_manual(values = c("troglomorphic_subterranean" = "#FD914E", 
                               "Surface" = "#7DA0FB")) 


dev.off()


#### Species tree + barplots for BUSCO and VISION ---------------------------------

#Prepare BUSCO table for figure

summary_vision_wide_info_NT_LOF_busco_fig <- 
  summary_vision_wide_info_NT_LOF %>%
  filter(Species != "Danio_rerio") %>%
  dplyr::select(Species, BUSCO_complete, BUSCO_Fragmented, BUSCO_Pseudogene)

summary_vision_wide_info_NT_LOF_busco_fig <- 
  as.data.frame(summary_vision_wide_info_NT_LOF_busco_fig %>%
                  pivot_longer(!Species, names_to = "category", values_to = "count"))
summary_vision_wide_info_NT_LOF_busco_fig$category <- 
  factor(summary_vision_wide_info_NT_LOF_busco_fig$category, 
         levels=c("BUSCO_complete", "BUSCO_Fragmented", "BUSCO_Pseudogene"))


#Prepare Vision table for figure

summary_vision_wide_info_NT_LOF_vision_fig <- 
  summary_vision_wide_info_NT_LOF %>%
  filter(Species != "Danio_rerio") %>%
  dplyr::select(Species, Complete, Incomplete, Pseudogene)

summary_vision_wide_info_NT_LOF_vision_fig <- 
  as.data.frame(summary_vision_wide_info_NT_LOF_vision_fig %>%
                  pivot_longer(!Species, names_to = "category", values_to = "count"))
summary_vision_wide_info_NT_LOF_vision_fig$category <- 
  factor(summary_vision_wide_info_NT_LOF_vision_fig$category, 
         levels=c("Complete", "Incomplete", "Pseudogene"))


#Draw the backbone first
p1 <- 
  ggtree(species_tree, branch.length="none", size=1)  %<+% species_table +
  aes(color=Habitat) + 
  scale_color_manual(values = c("troglomorphic_subterranean" = "#FD914E", 
                                "Surface" = "#7DA0FB")) +
  geom_tiplab(size=6, offset=0.1, fontface=3) + 
  theme(legend.position='none') +
  ggnewscale::new_scale_color() +
  geom_facet(panel = 'BUSCO', 
             data = summary_vision_wide_info_NT_LOF_busco_fig, geom = geom_bar, 
             mapping = aes(x = as.numeric(count), fill = as.factor(category)), 
             orientation = 'y', width = 0.8, stat='identity', 
             position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = c("BUSCO_complete" = "#08A278", 
                               "BUSCO_Fragmented" = "gray",
                               "BUSCO_Pseudogene" = "black")) +
  ggnewscale::new_scale_fill() +
  geom_facet(panel = 'Vision', 
             data = summary_vision_wide_info_NT_LOF_vision_fig, geom = geom_bar, 
             mapping = aes(x = as.numeric(count), fill = as.factor(category)), 
             orientation = 'y', width = 0.8, stat='identity', 
             position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = c("Complete" = "#08A278", 
                               "Incomplete" = "gray",
                               "Pseudogene" = "black")) +
  theme_tree2(legend.position = 'none') +
  xlim_tree(20)





pdf(file = "SpeciesTree_BUSCO_Vision.pdf",width = 13.34,  height = 4.61)

facet_widths(p1, c(Tree = 1))


dev.off()





#### IGV analysis and draw LoF tile plot ---------------------------------

#Generate a habitat table at the individual level
Indiv_df <- 
  as.data.frame(
    individual_tree_AA_rooted$tip.label)
colnames(Indiv_df) <- c("Individual")

Indiv_df <- 
  Indiv_df %>%
  mutate(Species = gsub("\\..*", "", Individual))

species_table_habitat <- species_table %>% dplyr::select(Species, Habitat)
individual_table_habitat <- left_join(Indiv_df, species_table_habitat, by="Species")


#Now lets make a IGV simplified table for graph

IGV_df <- 
  IGV_df %>%
  mutate(Species_Ind = 
           paste(Species, Individual, sep=".")) 

IGV_df_simp <- 
  IGV_df %>%
  dplyr::select(Species_Ind, Gene, LoF_gene_name, Genotype) %>%
  arrange(Gene) %>%
  dplyr::select(-Gene)

Indiv_df <- 
  as.data.frame(
    individual_tree_AA_rooted$tip.label)
colnames(Indiv_df) <- c("Individual")

Individuals <- Indiv_df %>% pull(Individual) %>% unique()
LoF_list <- IGV_df$LoF_gene_name %>% unique()

all_combinations <- expand.grid(Species_Ind = Individuals, LoF_gene_name = LoF_list)


IGV_df_expanded <- all_combinations %>%
  left_join(IGV_df_simp, by = c("Species_Ind", "LoF_gene_name")) %>%
  mutate(Genotype = ifelse(is.na(Genotype), "no_LoF_Homozygous", Genotype))


IGV_df_expanded <- 
  IGV_df_expanded %>% dplyr::select(Species_Ind, LoF_gene_name, Genotype)


IGV_df_expanded$LoF_gene_name <- factor(IGV_df_expanded$LoF_gene_name,
                                        levels = IGV_df_expanded$LoF_gene_name %>% unique()) 
IGV_df_expanded$Genotype <- factor(IGV_df_expanded$Genotype,
                                   levels = c("no_LoF_Homozygous", "LoF_Homozygous", "Heterozygous")) 



p1 <-
  ggtree(individual_tree_AA_rooted, size=1, branch.length="none" ) %<+% individual_table_habitat +
  aes(color = Habitat) +
  geom_tiplab(offset = 0.05, size=6, fontface = "italic") +
  scale_color_manual(values = c("troglomorphic_subterranean" = "#FD914E", 
                                "Surface" = "#7DA0FB")) +
  geom_facet(panel = "Gene", 
             data = IGV_df_expanded, 
             geom = geom_tile, 
             aes(x = as.integer(LoF_gene_name), 
                 fill = Genotype), 
             width = 1, 
             color="black") +
  scale_fill_manual(values = c("no_LoF_Homozygous" = "#009E73", 
                               "LoF_Homozygous" = "#E69F00",
                               "Heterozygous" = "#CC79A7")) +
  theme(legend.position = "none",
        strip.text = element_blank()) +
  xlim_tree(20) 



pdf(file = "SpeciesTree_TilePlot_LoF.pdf",width = 40.34,  height = 15.61)


facet_widths(p1, widths = c(2, 5))


dev.off()



#Now define which LoF are not covered (better define missing data)


uniq_LoF_clade_cov <- 
  read.table("LoF_assembly_covered.csv",
             sep=",",
             header=FALSE)

colnames(uniq_LoF_clade_cov) <- 
  c("Gene", "LoF_Type", "Pos_CDS", "Clade", "Species", "Coverage_Assembly", "Reason")


#Replace noLoF by Unknown for not covered positions


uniq_LoF_clade_cov <- 
  uniq_LoF_clade_cov %>%
  mutate(LoF_name = paste(LoF_Type, Pos_CDS, sep="_")) %>%
  mutate(LoF_gene_name = paste(Gene, LoF_name, sep="_"))
uniq_LoF_clade_cov[(uniq_LoF_clade_cov$Species == "Amblyopsis_rosae"),"Species"] <- "Troglichthys_rosae"

uniq_LoF_clade_cov[(uniq_LoF_clade_cov$LoF_gene_name == "Rh2.1_Start loss_NA"),"LoF_gene_name"] <- "Rh2.1_Start_loss"
uniq_LoF_clade_cov[(uniq_LoF_clade_cov$LoF_gene_name == "gnat2_Splice site_NA"),"LoF_gene_name"] <- "gnat2_Splice_site_6"
uniq_LoF_clade_cov[(uniq_LoF_clade_cov$LoF_gene_name == "Opn6a-2_Splice site_NA"),"LoF_gene_name"] <- "Opn6a-2_Splice_site_5"
uniq_LoF_clade_cov[(uniq_LoF_clade_cov$LoF_gene_name == "pde6hb_Splice site_NA"),"LoF_gene_name"] <- "pde6hb_Splice_site_1"
uniq_LoF_clade_cov[(uniq_LoF_clade_cov$LoF_gene_name == "Cryba1b_Splice site_NA"),"LoF_gene_name"] <- "Cryba1b_Splice_site_4"
uniq_LoF_clade_cov[(uniq_LoF_clade_cov$LoF_gene_name == "Crygn2_Start loss_3"),"LoF_gene_name"] <- "Crygn2_Start_loss"
uniq_LoF_clade_cov[(uniq_LoF_clade_cov$LoF_gene_name == "gcap3_Start loss_3"),"LoF_gene_name"] <- "gcap3_Start_loss"


IGV_df_expanded <- IGV_df_expanded %>% mutate(Species = gsub("\\..*", "", Species_Ind))

IGV_df_expanded_cov <- 
  left_join(IGV_df_expanded, uniq_LoF_clade_cov, by=c("Species", "LoF_gene_name"))

IGV_df_expanded_cov <- 
  IGV_df_expanded_cov %>%
  mutate(Genotye_Refined = if_else(
    Coverage_Assembly == "Covered",
    Genotype,
    "Unknown"
  ))




non_prez_lof <- 
  IGV_df_expanded_cov %>%
  group_by(LoF_gene_name, Genotype) %>%
  summarise(count = n()) %>%
  filter(count == 37) %>%
  filter(Genotype == "no_LoF_Homozygous") %>%
  pull(LoF_gene_name)


IGV_df_expanded_cov_prez <- 
  IGV_df_expanded_cov %>% filter(! LoF_gene_name %in% non_prez_lof)


IGV_df_expanded_cov_prez <- 
  IGV_df_expanded_cov_prez %>% dplyr::select(Species_Ind, LoF_gene_name, Genotye_Refined, Gene)

IGV_df_expanded_cov_prez$Genotye_Refined <- factor(IGV_df_expanded_cov_prez$Genotye_Refined,
                                                   levels = c("no_LoF_Homozygous", "LoF_Homozygous", "Heterozygous", "Unknown")) 


order_genes_plot <- 
  c("Rh1.1","Rh2.1","Lws1","Exorh","Opn3","Opn4m1","Opn4m3","Opn4x2",
    "Opn5","Opn6a-1","Opn6a-2","Opn7a","Opn7c","opn7e","Opn8b","Opn8c",
    "Parapinopsin-1","Parietopsin","Rgr1","Rgr2","Rrh","Va2","Va2-2","Tmt1b",
    "Tmt1b_2","Tmt3a","Tmt3b","Cryaa","Cryba1b","Cryba1l1","Cryba2b","Cryba4",
    "Cryba4_2","Crybb1","Crybb1-2","Crybb1l1","Crybb1l1-2","Crybb1l3","Crybgx",
    "Crygm5","Crygn2","Rpe65a","Rpe65b","arr3a","arr3b","saga","sagb","gcap1",
    "gcap2","gcap3","gcap7","grk1a","grk1b","grk7a","grk7b","pde6a","pde6b",
    "pde6c","pde6ga","pde6gb","pde6gb_2","pde6ha","pde6ha-2","pde6ha-3",
    "pde6hb","rcv1b","rcv2a","rcv2b","gc2","gc3","gucy2f","gucy2f-2",
    "gnat1","gnat2","gnb1a","gnb1a-2","gnb1b","gnb3a","gnb3b","gngt1","gngt2b")

IGV_df_expanded_cov_prez_order <- 
  left_join(data.frame(Gene=order_genes_plot), IGV_df_expanded_cov_prez,by="Gene") %>%
  dplyr::select(-Gene)

IGV_df_expanded_cov_prez_order <- IGV_df_expanded_cov_prez_order %>% filter(! is.na(LoF_gene_name))


IGV_df_expanded_cov_prez_order$LoF_gene_name <- factor(IGV_df_expanded_cov_prez_order$LoF_gene_name,
                                                       levels = IGV_df_expanded_cov_prez_order$LoF_gene_name %>% unique()) 

#opn7e = opn7d


p2 <-
  ggtree(individual_tree_AA_rooted, size=1, branch.length="none" ) %<+% individual_table_habitat +
  aes(color = Habitat) +
  geom_tiplab(offset = 0.05, size=6, fontface = "italic") +
  scale_color_manual(values = c("troglomorphic_subterranean" = "#FD914E", 
                                "Surface" = "#7DA0FB")) +
  geom_facet(panel = "Gene", 
             data = IGV_df_expanded_cov_prez_order, 
             geom = geom_tile, 
             aes(x = as.integer(LoF_gene_name), 
                 fill = Genotye_Refined), 
             width = 1, 
             color="black") +
  scale_fill_manual(values = c("no_LoF_Homozygous" = "#009E73", 
                               "LoF_Homozygous" = "#E69F00",
                               "Heterozygous" = "#CC79A7",
                               "Unknown"="#96AFA8")) +
  theme(legend.position = "none",
        strip.text = element_blank()) +
  xlim_tree(20) 


pdf(file = "SpeciesTree_TilePlot_LoF_Refined.pdf",width = 40.34,  height = 15.61)


facet_widths(p2, widths = c(2, 5))


dev.off()



IGV_df_expanded_cov_prez_order_wide <- 
  as.data.frame(IGV_df_expanded_cov_prez_order %>%
                  pivot_wider(names_from = Species_Ind, values_from = Genotye_Refined))


#### Count the number of LoFs and Pseudogenes in each branch ---------------------------------

individual_tree_AA_rooted <- makeNodeLabel(individual_tree_AA_rooted, method="number", prefix="Node")
write.tree(individual_tree_AA_rooted, "SpeciesTree_nodelabel.nwk")


species_test <- "Astyanax_mexicanus_Surface"
misc_table <- as.data.frame(species_test) %>% mutate(value = 1)
colnames(misc_table) <- c("node","value")
misc_table$node <- as.numeric(misc_table$node)
node_label_corresp <- 
  left_join(individual_tree_AA_rooted, misc_table, by = 'node') %>%
  dplyr::select(node, label)


#Count the number of LoF mutation which occured in each branch of the tree

uniq_LoF_vector <- 
  IGV_df_expanded %>% 
  filter(Genotype %in% c("LoF_Homozygous", "Heterozygous")) %>% 
  pull(LoF_gene_name) %>% unique()



LoF_event_branch_df <- as.data.frame(NULL)
for(curr_LoF in uniq_LoF_vector){
  
  curr_ind_mut <- 
    IGV_df_expanded %>%
    filter(LoF_gene_name == curr_LoF) %>%
    filter(Genotype %in% c("LoF_Homozygous", "Heterozygous")) %>%
    pull(Species_Ind)
  nb_ind_mut <- length(curr_ind_mut)
  
  if(nb_ind_mut > 1){
    curr_MRCA <- getMRCA(individual_tree_AA_rooted, curr_ind_mut)
    curr_desc <- Descendants(individual_tree_AA_rooted, curr_MRCA)[[1]]
    curr_desc_label <- node_label_corresp %>% filter(node %in% curr_desc) %>% pull(label)
    nb_inf_mut <- length(curr_desc_label)
    
    MRCA_label <- node_label_corresp %>% filter(node == curr_MRCA) %>% pull(label)
    merged_desc_label <- paste(curr_desc_label, collapse = ",")
    
    if(nb_inf_mut > nb_ind_mut) {
      curr_state="paraphyletic"
    } else {
      curr_state="monophyletic"
    }
    
    curr_df <- 
      as.data.frame(cbind(curr_LoF, curr_state, MRCA_label, merged_desc_label))
    
  } else {
    
    curr_state="terminal branch"
    
    curr_df <- 
      as.data.frame(cbind(curr_LoF, curr_state, curr_ind_mut, curr_ind_mut))
  }
  
  
  
  colnames(curr_df) <- c("LoF_gene_name", "state", "MRCA", "Species")
  
  LoF_event_branch_df <- rbind(LoF_event_branch_df, curr_df)
  
}

LoF_event_branch_df <- as.data.frame(LoF_event_branch_df)

LoF_event_branch_df %>%
  group_by(state) %>%
  summarise(n())


as.data.frame(
  LoF_event_branch_df %>%
    group_by(MRCA, Species) %>%
    summarise(count = n())
) %>%
  dplyr::select(MRCA, count, Species)




#Now lets do the same but with the number of Pseudogenes

IGV_df_expanded <- IGV_df_expanded %>%
  mutate(gene_name = gsub("_.*", "", LoF_gene_name))

IGV_df_expanded %>% arrange(gene_name) %>% pull(gene_name) %>% unique() 


LoF_event_branch_df <- LoF_event_branch_df %>% mutate(gene_name = gsub("_.*", "", LoF_gene_name))

uniq_pseudogenes <- LoF_event_branch_df %>% pull(gene_name) %>% unique()



as.data.frame(
  LoF_event_branch_df %>%
    group_by(gene_name, MRCA) %>%
    summarise(count = n())
) %>%
  arrange(gene_name)







#### Compute the gene repertoire at the individual level ---------------------------------

IGV_df_expanded <- IGV_df_expanded %>%
  mutate(gene_name = gsub("_.*", "", LoF_gene_name)) %>%
  mutate(Species = gsub("\\..*", "", Species_Ind))

list_species <- vision_seq_df %>% pull(Species) %>% unique()

Individual_repertoire_df <- as.data.frame(NULL)
for(curr_sp in list_species){
  
  curr_genes <- 
    vision_seq_df %>%
    filter(Species == curr_sp) %>%
    pull(gene_name)
  
  list_individuals <- 
    IGV_df_expanded %>%
    filter(Species == curr_sp) %>%
    pull(Species_Ind) %>%
    unique()
  
  for(curr_ind in list_individuals){
    
    
    curr_pseudos <- 
      IGV_df_expanded %>%
      filter(Species_Ind == curr_ind) %>%
      filter(Genotype != "no_LoF_Homozygous") %>%
      pull(gene_name) %>%
      unique()
    
    nb_pseudos_igv <- length(curr_pseudos)
    
    
    curr_df <- 
      vision_seq_df %>%
      filter(Species == curr_sp) %>%
      mutate(Species_Ind = curr_ind) %>%
      mutate(Gene_Type_ind = if_else(
        gene_name %in% curr_pseudos,
        "Pseudogene", 
        "Intact"
      )) %>%
      dplyr::select(Gene_clade, Species, gene_name, CDS_length, Exon_Count, Species_Ind, Gene_Type_ind)
    
    nb_pseudo_seq <- 
      length(curr_df %>% filter(Gene_Type_ind == "Pseudogene") %>% pull(gene_name) %>% unique())
    
    
    #print(paste(curr_ind,nb_pseudos_igv,nb_pseudo_seq, sep=","))
    
    Individual_repertoire_df <- rbind(Individual_repertoire_df, curr_df)
  }
  
  
}


#Now lets make a summary table :) 

Individual_repertoire_df_summary <- 
  Individual_repertoire_df %>%
  group_by(Species_Ind, Gene_Type_ind) %>%
  summarise(count = n())

Individual_repertoire_df_summary_wide <- 
  as.data.frame(
    Individual_repertoire_df_summary %>%
      pivot_wider(names_from = Gene_Type_ind, values_from = count, values_fill = 0)
  )



Individual_repertoire_df_summary_wide <- 
  Individual_repertoire_df_summary_wide %>%
  mutate(Total = Intact + Pseudogene ) %>%
  mutate(Prop_pseudo = Pseudogene / Total)

## Add the number of nucleotides and the number of LoF ... 

LoF_per_ind <- 
  as.data.frame(IGV_df_expanded %>%
                  filter(Genotype != "no_LoF_Homozygous") %>%
                  group_by(Species_Ind) %>%
                  summarise(LoF_nb = n()) %>%
                  arrange(desc(LoF_nb)))

nt_per_ind <- 
  as.data.frame(Individual_repertoire_df %>%
                  group_by(Species_Ind) %>%
                  summarise(nt_nb = sum(CDS_length)))

Individual_repertoire_df_summary_wide_NT <- 
  left_join(Individual_repertoire_df_summary_wide, nt_per_ind, by="Species_Ind") 
Individual_repertoire_df_summary_wide_NT_LOF <- 
  left_join(Individual_repertoire_df_summary_wide_NT, LoF_per_ind, by="Species_Ind") 


Individual_repertoire_df_summary_wide_NT_LOF <- 
  Individual_repertoire_df_summary_wide_NT_LOF %>%
  mutate(LoF_per_nt = LoF_nb/nt_nb)

Individual_repertoire_df_summary_wide_NT_LOF[is.na(Individual_repertoire_df_summary_wide_NT_LOF)] <- 0


write.table(Individual_repertoire_df_summary_wide_NT_LOF,
            "Individual_repertoire_df_summary_wide_NT_LOF.tsv",
            sep="\t",
            quote=FALSE,
            col.names=TRUE,
            row.names=FALSE)


gt(Individual_repertoire_df_summary_wide_NT_LOF %>%
     dplyr::select(Species_Ind, Intact, Pseudogene, Total, Prop_pseudo, nt_nb, LoF_nb, LoF_per_nt) %>%
     arrange(Prop_pseudo))




#### Draw the gene repertoire tile plot ---------------------------------

uniq_genes <- 
  Individual_repertoire_df %>%
  arrange(gene_name) %>%
  pull(gene_name) %>%
  unique()

uniq_individuals <- 
  Individual_repertoire_df %>%
  arrange(Species_Ind) %>%
  pull(Species_Ind) %>%
  unique()

all_combinations <- expand.grid(Species_Ind = uniq_individuals, gene_name = uniq_genes)


Individual_repertoire_df_expanded <- 
  all_combinations %>%
  left_join(Individual_repertoire_df, by = c("Species_Ind", "gene_name")) %>%
  mutate(Gene_Type_ind = ifelse(is.na(Gene_Type_ind), "Not_found", Gene_Type_ind))



order_genes_plot <- 
  c("Rh1.1","Rh2.1","Lws1","Exorh","Opn3","Opn4m1","Opn4m3","Opn4x2",
    "Opn5","Opn6a-1","Opn6a-2","Opn7a","Opn7c","opn7e","Opn8b","Opn8c",
    "Parapinopsin-1","Parietopsin","Rgr1","Rgr2","Rrh","Va2","Va2-2","Tmt1b",
    "Tmt1b_2","Tmt3a","Tmt3b","Cryaa","Cryba1b","Cryba1l1","Cryba2b","Cryba4",
    "Cryba4_2","Crybb1","Crybb1-2","Crybb1l1","Crybb1l1-2","Crybb1l3","Crybgx",
    "Crygm5","Crygn2","Rpe65a","Rpe65b","arr3a","arr3b","saga","sagb","gcap1",
    "gcap2","gcap3","gcap7","grk1a","grk1b","grk7a","grk7b","pde6a","pde6b",
    "pde6c","pde6ga","pde6gb","pde6gb_2","pde6ha","pde6ha-2","pde6ha-3",
    "pde6hb","rcv1b","rcv2a","rcv2b","gc2","gc3","gucy2f","gucy2f-2",
    "gnat1","gnat2","gnb1a","gnb1a-2","gnb1b","gnb3a","gnb3b","gngt1","gngt2b")


Individual_repertoire_df_expanded$gene_name <- factor(Individual_repertoire_df_expanded$gene_name,
                                                      levels = order_genes_plot)

Individual_repertoire_df_expanded$Gene_Type_ind <- factor(Individual_repertoire_df_expanded$Gene_Type_ind,
                                                          levels = c("Intact", "Pseudogene", "Not_found")) 



p1 <-
  ggtree(individual_tree_AA_rooted, size=1, branch.length="none" ) %<+% individual_table_habitat +
  aes(color = Habitat) +
  geom_tiplab(offset = 0.05, size=6, fontface = "italic") +
  scale_color_manual(values = c("troglomorphic_subterranean" = "#FD914E", 
                                "Surface" = "#7DA0FB")) +
  geom_facet(panel = "Gene", 
             data = Individual_repertoire_df_expanded, 
             geom = geom_tile, 
             aes(x = as.integer(gene_name), 
                 fill = Gene_Type_ind), 
             width = 1, 
             color="black") +
  scale_fill_manual(values = c("Intact" = "#08A278", 
                               "Pseudogene" = "#D55E00",
                               "Not_found" = "gray")) +
  theme(legend.position = "none",
        strip.text = element_blank()) +
  xlim_tree(20) 



pdf(file = "SpeciesTree_TilePlot_Pseudogene.pdf",width = 40.34,  height = 15.61)


facet_widths(p1, widths = c(2, 5))


dev.off()



#### Load dn/ds values -- Individual level ---------------------------------

vision_dNdS_ind <- 
  read.table("Vision_genes_dNdS.IND.csv",
             sep=",",
             header=FALSE)
colnames(vision_dNdS_ind) <- 
  c("Gene_clade", "species_gene_name", "LB", "MLE", "UB", "dN", "dS")

#Remove bad estimates of omega
vision_dNdS_ind <- 
  vision_dNdS_ind %>%
  filter(dN < 1) %>%
  filter(dS < 1) %>%
  filter(dS > 0.01)


list_species <-
  vision_seq_df %>%
  pull(Species)  %>%
  unique()

#Assign a species to every gene name

vision_dNdS_ind <- 
  vision_dNdS_ind %>% 
  mutate(across('species_gene_name', str_replace, 'Amblyopsis_rosae', 'Troglichthys_rosae'))


all_gene_name <- vision_dNdS_ind %>% pull(species_gene_name) %>% unique()

list_species <- c(list_species, "Gadus_morhua")

vision_dNdS_ind_df <- as.data.frame(NULL)
for(curr_species in list_species){
  
  curr_sp_gene_names <- 
    all_gene_name[grepl(curr_species, all_gene_name)]
  
  curr_sp_table <- 
    vision_dNdS_ind %>%
    filter(species_gene_name %in% curr_sp_gene_names) %>%
    mutate(Species = curr_species)
  
  vision_dNdS_ind_df <- 
    rbind(vision_dNdS_ind_df,
          curr_sp_table)
  
  
}


#Add info about complete or incomplete for each gene

vision_dNdS_ind_df_writtenState <- 
  vision_dNdS_ind_df %>%
  filter(Species %in% c("Danio_rerio", "Gadus_morhua", "Percopsis_transmontana", "Speoplatyrhinus_poulsoni"))


vision_dNdS_ind_df_unwrittenState <- 
  vision_dNdS_ind_df %>%
  filter(! Species %in% c("Danio_rerio", "Gadus_morhua", "Percopsis_transmontana", "Speoplatyrhinus_poulsoni"))

vision_dNdS_ind_df_unwrittenState <- 
  vision_dNdS_ind_df_unwrittenState %>%
  mutate(gene_name = 
           gsub(".*_", "", species_gene_name)) 


Pseudo_per_indiv <- 
  Individual_repertoire_df %>% dplyr::select(Species_Ind, gene_name, Gene_Type_ind) 


list_species <- Individual_repertoire_df %>% pull(Species_Ind) %>% unique()
all_gene_name <- vision_dNdS_ind_df_unwrittenState %>% pull(species_gene_name) %>% unique()
vision_dNdS_ind_df_unwrittenState_ind <- as.data.frame(NULL)
for(curr_species in list_species){
  
  curr_sp_gene_names <- 
    all_gene_name[grepl(curr_species, all_gene_name)]
  
  curr_sp_table <- 
    vision_dNdS_ind_df_unwrittenState %>%
    filter(species_gene_name %in% curr_sp_gene_names) %>%
    mutate(Species_Ind = curr_species)
  
  vision_dNdS_ind_df_unwrittenState_ind <- 
    rbind(vision_dNdS_ind_df_unwrittenState_ind,
          curr_sp_table)
  
  
}


#Rename some genes

vision_dNdS_ind_df_unwrittenState_ind$gene_name[vision_dNdS_ind_df_unwrittenState_ind$Species == "Troglichthys_rosae" & vision_dNdS_ind_df_unwrittenState_ind$gene_name == "Rh2"] <- "Rh2.1"
vision_dNdS_ind_df_unwrittenState_ind$gene_name[vision_dNdS_ind_df_unwrittenState_ind$Species == "Aphredoderus_sayanus" & vision_dNdS_ind_df_unwrittenState_ind$gene_name == "Rh2"] <- "Rh2.1"
vision_dNdS_ind_df_unwrittenState_ind$gene_name[vision_dNdS_ind_df_unwrittenState_ind$Species == "Chologaster_cornuta" & vision_dNdS_ind_df_unwrittenState_ind$gene_name == "Rh2"] <- "Rh2.1"
vision_dNdS_ind_df_unwrittenState_ind$gene_name[vision_dNdS_ind_df_unwrittenState_ind$Species == "Typhlichthys_subterraneus" & vision_dNdS_ind_df_unwrittenState_ind$gene_name == "Rh2"] <- "Rh2.1"
vision_dNdS_ind_df_unwrittenState_ind$gene_name[vision_dNdS_ind_df_unwrittenState_ind$Species == "Forbesichthys_papilliferus" & vision_dNdS_ind_df_unwrittenState_ind$gene_name == "Rh2"] <- "Rh2.1"

vision_dNdS_ind_df_unwrittenState_ind$gene_name[vision_dNdS_ind_df_unwrittenState_ind$Species == "Troglichthys_rosae" & vision_dNdS_ind_df_unwrittenState_ind$gene_name == "Rh1"] <- "Rh1.1"
vision_dNdS_ind_df_unwrittenState_ind$gene_name[vision_dNdS_ind_df_unwrittenState_ind$Species == "Aphredoderus_sayanus" & vision_dNdS_ind_df_unwrittenState_ind$gene_name == "Rh1"] <- "Rh1.1"
vision_dNdS_ind_df_unwrittenState_ind$gene_name[vision_dNdS_ind_df_unwrittenState_ind$Species == "Chologaster_cornuta" & vision_dNdS_ind_df_unwrittenState_ind$gene_name == "Rh1"] <- "Rh1.1"
vision_dNdS_ind_df_unwrittenState_ind$gene_name[vision_dNdS_ind_df_unwrittenState_ind$Species == "Typhlichthys_subterraneus" & vision_dNdS_ind_df_unwrittenState_ind$gene_name == "Rh1"] <- "Rh1.1"
vision_dNdS_ind_df_unwrittenState_ind$gene_name[vision_dNdS_ind_df_unwrittenState_ind$Species == "Forbesichthys_papilliferus" & vision_dNdS_ind_df_unwrittenState_ind$gene_name == "Rh1"] <- "Rh1.1"
vision_dNdS_ind_df_unwrittenState_ind$gene_name[vision_dNdS_ind_df_unwrittenState_ind$Species == "Percopsis_omiscomaycus" & vision_dNdS_ind_df_unwrittenState_ind$gene_name == "Rh1"] <- "Rh1.1"


vision_dNdS_ind_df_unwrittenState_ind$gene_name[vision_dNdS_ind_df_unwrittenState_ind$Species == "Troglichthys_rosae" & vision_dNdS_ind_df_unwrittenState_ind$gene_name == "Parapinopsin"] <- "Parapinopsin-1"
vision_dNdS_ind_df_unwrittenState_ind$gene_name[vision_dNdS_ind_df_unwrittenState_ind$Species == "Aphredoderus_sayanus" & vision_dNdS_ind_df_unwrittenState_ind$gene_name == "Parapinopsin"] <- "Parapinopsin-1"
vision_dNdS_ind_df_unwrittenState_ind$gene_name[vision_dNdS_ind_df_unwrittenState_ind$Species == "Chologaster_cornuta" & vision_dNdS_ind_df_unwrittenState_ind$gene_name == "Parapinopsin"] <- "Parapinopsin-1"
vision_dNdS_ind_df_unwrittenState_ind$gene_name[vision_dNdS_ind_df_unwrittenState_ind$Species == "Typhlichthys_subterraneus" & vision_dNdS_ind_df_unwrittenState_ind$gene_name == "Parapinopsin"] <- "Parapinopsin-1"
vision_dNdS_ind_df_unwrittenState_ind$gene_name[vision_dNdS_ind_df_unwrittenState_ind$Species == "Forbesichthys_papilliferus" & vision_dNdS_ind_df_unwrittenState_ind$gene_name == "Parapinopsin"] <- "Parapinopsin-1"
vision_dNdS_ind_df_unwrittenState_ind$gene_name[vision_dNdS_ind_df_unwrittenState_ind$Species == "Typhlichthys_eigenmanni" & vision_dNdS_ind_df_unwrittenState_ind$gene_name == "Parapinopsin"] <- "Parapinopsin-1"
vision_dNdS_ind_df_unwrittenState_ind$gene_name[vision_dNdS_ind_df_unwrittenState_ind$Species == "Forbesichthys_agassizii" & vision_dNdS_ind_df_unwrittenState_ind$gene_name == "Parapinopsin"] <- "Parapinopsin-1"
vision_dNdS_ind_df_unwrittenState_ind$gene_name[vision_dNdS_ind_df_unwrittenState_ind$Species == "Percopsis_omiscomaycus" & vision_dNdS_ind_df_unwrittenState_ind$gene_name == "Parapinopsin"] <- "Parapinopsin-1"

vision_dNdS_ind_df_unwrittenState_ind$gene_name[vision_dNdS_ind_df_unwrittenState_ind$Species == "Percopsis_omiscomaycus" & vision_dNdS_ind_df_unwrittenState_ind$gene_name == "2"] <- "Tmt1b_2"


vision_dNdS_ind_df_unwrittenState_ind_state <- 
  left_join(vision_dNdS_ind_df_unwrittenState_ind, Pseudo_per_indiv, by=c("Species_Ind", "gene_name"))






vision_dNdS_ind_df_writtenState <- 
  vision_dNdS_ind_df_writtenState %>%
  mutate(Species_Ind = case_when(
    Species == "Danio_rerio" ~ "Danio_rerio",
    Species == "Percopsis_transmontana" ~ "Percopsis_transmontana.KU_KUI_29775",
    Species == "Speoplatyrhinus_poulsoni" ~ "Speoplatyrhinus_poulsoni.IRGN_2EJSDHOURK",
    Species == "Gadus_morhua" ~ "Gadus_morhua"
  ))


vision_dNdS_ind_df_writtenState <- vision_dNdS_ind_df_writtenState %>%
  mutate(Gene_Type_ind = if_else(
    Gene_Type == "Pseudogene",
    "Pseudogene",
    "Intact"
  ))



#Combine dataframes

vision_dNdS_ind_df_writtenState_nogenetype <- vision_dNdS_ind_df_writtenState 
mycols <- colnames(vision_dNdS_ind_df_writtenState_nogenetype)
vision_dNdS_ind_df_unwrittenState_ind_state_cols <- vision_dNdS_ind_df_unwrittenState_ind_state %>% dplyr::select(mycols)


vision_dNdS_ind_df <- rbind(vision_dNdS_ind_df_unwrittenState_ind_state_cols, vision_dNdS_ind_df_writtenState_nogenetype)



#### Analysis of dN/dS values - Individual level ---------------------------------


vision_dNdS_ind_df <- 
  left_join(vision_dNdS_ind_df, 
            species_table %>% dplyr::select(Species, Habitat) , 
            by="Species")

vision_dNdS_ind_df[(vision_dNdS_ind_df$Species == "Gadus_morhua"),"Habitat"] <- "Surface"


#Lets plot the graphs


pdf(file = "dnds_Vision_boxplot.IND.pdf",width = 8.34,  height = 4.61)

vision_dNdS_ind_df %>%
  ggplot(., aes(x=reorder(Species_Ind, MLE, mean), y=MLE, fill=Habitat)) +
  geom_boxplot() +
  scale_fill_manual(values = c("troglomorphic_subterranean" = "#FD914E", 
                               "Surface" = "#7DA0FB")) +
  ylab("dN/dS") +
  xlab("Species") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  theme(legend.position="none")

dev.off()



pdf(file = "dnds_Vision_violinplot.IND.pdf",width = 8.34,  height = 4.61)

vision_dNdS_ind_df %>%
  ggplot(., aes(x=reorder(Species_Ind, MLE, mean), y=MLE, fill=Habitat)) +
  geom_violin() +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  geom_boxplot(width=0.1, outlier.shape = NA) + 
  scale_fill_manual(values = c("troglomorphic_subterranean" = "#FD914E", 
                               "Surface" = "#7DA0FB")) +
  ylab("dN/dS") +
  xlab("Species") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  theme(legend.position="none")

dev.off()





vision_dNdS_df_mean <- 
  as.data.frame(
    vision_dNdS_ind_df %>%
      group_by(Species) %>%
      summarise(mean_MLE = mean(MLE),
                nb_genes = n())
  ) %>%
  arrange(mean_MLE)


#pGLS between mean dN/dS and habitat


vision_dNdS_df_mean_info <- 
  left_join(vision_dNdS_df_mean,
            summary_vision_wide_info)

caper_dNdS_pseudo <- 
  comparative.data(phy = species_tree, 
                   data = vision_dNdS_df_mean_info,
                   names.col = Species, vcv = TRUE,
                   na.omit = FALSE, warn.dropped = TRUE)
caper_dNdS_pseudo_chase <- 
  comparative.data(phy = species_tree_chase, 
                   data = vision_dNdS_df_mean_info,
                   names.col = Species, vcv = TRUE,
                   na.omit = FALSE, warn.dropped = TRUE)





dNdS_vs_Habitat <-
  pgls(mean_MLE ~ Habitat, 
       data = caper_dNdS_pseudo, 
       lambda = "ML")
summary(dNdS_vs_Habitat)

dNdS_vs_Habitat_chase <-
  pgls(mean_MLE ~ Habitat, 
       data = caper_dNdS_pseudo_chase, 
       lambda = "ML")
summary(dNdS_vs_Habitat_chase)



#Now replot the prop of pseudo per individual + assign 0 to gadus and danio as no pseudogenes are observed 
#in these species

Intact_df <-  as.data.frame(Individual_repertoire_df %>%
                              filter(Gene_Type_ind == "Intact") %>%
                              group_by(Species_Ind, Species) %>%
                              summarise(Intact_nb = n()))

Pseudo_df <- as.data.frame(Individual_repertoire_df %>%
                             filter(Gene_Type_ind != "Intact") %>%
                             group_by(Species_Ind, Species) %>%
                             summarise(Pseudo_nb = n()))


Intact_Pseudo_df <- left_join(Intact_df, Pseudo_df, by=c("Species_Ind", "Species"))
Intact_Pseudo_df$Pseudo_nb[is.na(Intact_Pseudo_df$Pseudo_nb)] <- 0
Intact_Pseudo_df <- Intact_Pseudo_df %>% mutate(Prop_pseudo = Pseudo_nb / (Intact_nb + Pseudo_nb))


Intact_Pseudo_df_habitat <- 
  left_join(Intact_Pseudo_df, species_table %>% dplyr::select(Species, Habitat), by="Species")


dr_table <- cbind("Danio_rerio", "Danio_rerio", "90", "0", "0", "Surface")
colnames(dr_table) <- colnames(Intact_Pseudo_df_habitat)
gr_table <- cbind("Gadus_morhua", "Gadus_morhua", "90", "0", "0", "Surface")
colnames(gr_table) <- colnames(Intact_Pseudo_df_habitat)

Intact_Pseudo_df_habitat <- rbind(Intact_Pseudo_df_habitat, dr_table, gr_table)

Intact_Pseudo_df_habitat$Prop_pseudo <- as.numeric(Intact_Pseudo_df_habitat$Prop_pseudo)



pdf(file = "PseudoProp_Violinplot_diploid.IND.pdf",width = 8.34,  height = 4.61)

Intact_Pseudo_df_habitat %>%
  ggplot(., aes(x=Habitat, y=Prop_pseudo, fill=Habitat)) +
  geom_violin() +
  geom_jitter(size=2, shape=16, color="black", alpha=0.8, position=position_jitter(0.2)) +
  theme_classic() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") +
  ylab("Proportion of pseudogenes") +
  scale_fill_manual(values = c("troglomorphic_subterranean" = "#FD914E", 
                               "Surface" = "#7DA0FB")) 


dev.off()


#### Compute a LoF rate based on observed LoF mutations ---------------------------------


#Transition/transversion ratio computed with PAML on a concatenated alignment of vision genes

R_ratio = 1.75872 



#Stop codon probability

#extract sequences of P. transmontana (=reference)
Concat_vision_seq <- 
  vision_seq_df %>%
  filter(Species == "Percopsis_transmontana") %>%
  pull(Sequence) %>%
  paste(., collapse = '') %>%
  as.character(.)


n=0
x=0

for(codon_stop in seq(3, nchar(Concat_vision_seq), 3)){
  codon_start <- codon_stop-2
  codon_seq <- substr(Concat_vision_seq, codon_start, codon_stop)
  
  if(codon_seq == "CAA" | codon_seq == "CGA" | codon_seq == "CAG"){
    x <- x + R_ratio/((3*R_ratio)+6)
    n <- n + 3
  } else if(codon_seq == "GAA" | codon_seq == "AAA" | codon_seq == "AGA" | codon_seq == "TGT" | codon_seq == "TGC" | codon_seq == "AAG" | codon_seq == "GAG" | codon_seq == "TTG" | codon_seq == "TCG" | codon_seq == "GGA"){
    x <- x + 1.0/((3*R_ratio)+6)
    n <- n + 3
    
  } else if(codon_seq == "TTA" | codon_seq == "TCA" | codon_seq == "TAC" | codon_seq == "TAT" | codon_seq == "TTA" | codon_seq == "TAC"){
    x <- x + 2.0/((3*R_ratio)+6)
    n <- n + 3
    
  } else if(codon_seq == "TGG"){
    x <= x + (2.0*R_ratio)/((3*R_ratio)+6)
    n <- n + 3
    
  } else {
    
    x <- x + 0.0
    n <- n + 3
    
  }
  
  
  
}

proba_stop_gain <- x/(n/3)


#Probability of Frameshift
nb_stop <- 
  nrow(
    Table_all_LoF_igv_uniq %>%
      filter(LoF_Type == "STOP")
  )


nb_fs <- 
  nrow(
    Table_all_LoF_igv_uniq %>%
      filter(LoF_Type %in% c("Del", "Ins"))
  )


proba_frameshift <- proba_stop_gain * (nb_fs/nb_stop)


#Proba of splice site mutation (take the mean number of exons over species)

intron_number <-
  vision_seq_df %>%
  mutate(Intron_count = as.numeric(Exon_Count) - 1) %>%
  group_by(Species) %>%
  summarise(total_intron = sum(Intron_count)) %>%
  pull(total_intron) %>%
  sum()


base_number <- 
  vision_seq_df %>%
  group_by(Species) %>%
  summarise(sum_CDS_length = sum(CDS_length)) %>%
  pull(sum_CDS_length) %>%
  sum()

proba_splice <- 4 * (intron_number / base_number)

#Proba of Start and Stop loss

proba_start_loss <- 
  3 * (
    vision_seq_df %>% group_by(Species) %>% summarise(count = n()) %>% pull(count) %>% sum()
    / base_number
  )
proba_stop_loss <- 0.852 * proba_start_loss

# Sum of LoF probas


LoF_proba <- proba_stop_loss + proba_start_loss + proba_splice + proba_frameshift + proba_stop_gain


#Make a nice table

LoF_rates_df <- 
  as.data.frame(
    cbind(proba_stop_gain, proba_frameshift, 
          proba_start_loss, proba_stop_loss, proba_splice, LoF_proba))

colnames(LoF_rates_df) <- 
  c("Stop_gain", "Frameshift", "Start_loss", "Stop_loss", "Splice", "Total")



write.table(LoF_rates_df,
            "~/Amblyopsidae_Analysis/RawPlots/LoF_rates_df.tsv",
            sep="\t",
            quote=FALSE,
            col.names=TRUE,
            row.names=FALSE)



#### Get the observed distribution of LoF per gene per species  ---------------------------------

LoF_obs_df <- 
  as.data.frame(
    IGV_df_expanded %>%
      filter(Genotype != "no_LoF_Homozygous") %>%
      group_by(Species_Ind, gene_name) %>%
      summarise(LoF_nb = n())
  )

list_ind <- LoF_obs_df$Species_Ind %>% unique()

noLoF_obs_df <- as.data.frame(NULL)
for(curr_ind in list_ind){
  
  pseudo_prez <- 
    LoF_obs_df %>%
    filter(Species_Ind == curr_ind) %>%
    pull(gene_name)
  
  curr_df <- 
    Individual_repertoire_df_expanded %>%
    filter(Species_Ind == curr_ind) %>%
    filter(! gene_name %in% pseudo_prez) %>%
    filter(Gene_Type_ind != "Not_found") %>%
    mutate(LoF_nb = 0) %>%
    dplyr::select(Species_Ind, gene_name,LoF_nb)
  
  noLoF_obs_df <- rbind(noLoF_obs_df, curr_df)
}



LoF_distrib_obs_df <- rbind(LoF_obs_df, noLoF_obs_df)


LoF_distrib_obs_df_summary <-
  as.data.frame(LoF_distrib_obs_df %>%
                  group_by(Species_Ind, LoF_nb) %>%
                  summarise(count = n()))

sp_ind_df <- IGV_df %>% dplyr::select(Species, Species_Ind) %>% distinct()

LoF_distrib_obs_df_summary_sp <- 
  left_join(LoF_distrib_obs_df_summary, sp_ind_df, by="Species_Ind")


#Add the proportion

LoF_distrib_obs_df_summary_sp <- 
  as.data.frame(
    LoF_distrib_obs_df_summary_sp %>%
      group_by(Species_Ind) %>%
      mutate(total = sum(count))
  )

LoF_distrib_obs_df_summary_sp <- 
  LoF_distrib_obs_df_summary_sp %>%
  mutate(prop_count = count / total)

#Finally, extend for each category to be present in each species (from 0 to 6)

list_sp_ind <- LoF_distrib_obs_df_summary_sp$Species_Ind %>% unique
LoF_categories <- LoF_distrib_obs_df_summary_sp$LoF_nb %>% unique()

all_combinations <- expand.grid(Species_Ind = list_sp_ind, LoF_nb = LoF_categories)


LoF_distrib_obs_df_summary_sp_exp <- 
  all_combinations %>%
  left_join(LoF_distrib_obs_df_summary_sp, by = c("Species_Ind", "LoF_nb")) %>%
  mutate(across(where(is.numeric), ~ replace_na(., 0)),
         across(where(is.character), ~ replace_na(., "0")))


LoF_distrib_obs_df_summary_sp_exp <- 
  LoF_distrib_obs_df_summary_sp_exp %>%
  mutate(Species = gsub("\\..*", "", Species_Ind))

#lets plot

list_sp <- LoF_distrib_obs_df_summary_sp$Species %>% unique

pdf("MutationsNumber_profile_absolute.pdf",width = 6.34,  height = 4.61)

for(curr_sp in list_sp){
  
  
  p1 <- 
    LoF_distrib_obs_df_summary_sp_exp %>%
    filter(Species == curr_sp) %>%
    ggplot(., aes(x=as.numeric(LoF_nb), y=count, fill=Species_Ind)) +
    scale_fill_manual(values = c("#000000", "#E69F00", 
                                 "#56B4E9", "#009E73",
                                 "#F0E442", "#0072B2", 
                                 "#D55E00", "#CC79A7")) + 
    geom_bar(stat="identity", position = "dodge") +
    theme_classic() +
    scale_x_continuous("# of LoF", 
                       c(0,1, 2, 3, 4, 5, 6)) +
    #scale_x_continuous("# of LoF", 
    #                   x_axis_ticks) +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14)) +
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=18),
          plot.subtitle=element_text(size=16),
          axis.title.x=element_blank()) +
    theme(legend.position="none") 
  
  print(p1)
  
}

dev.off()




pdf("MutationsNumber_profile_proportion.pdf",width = 6.34,  height = 4.61)

for(curr_sp in list_sp){
  
  
  p1 <- 
    LoF_distrib_obs_df_summary_sp_exp %>%
    filter(Species == curr_sp) %>%
    ggplot(., aes(x=as.numeric(LoF_nb), y=prop_count, fill=Species_Ind)) +
    scale_fill_manual(values = c("#000000", "#E69F00", 
                                 "#56B4E9", "#009E73",
                                 "#F0E442", "#0072B2", 
                                 "#D55E00", "#CC79A7")) + 
    geom_bar(stat="identity", position = "dodge") +
    theme_classic() +
    scale_x_continuous("# of LoF", 
                       c(0,1, 2, 3, 4, 5, 6)) +
    ylim(0, 1) +
    #scale_x_continuous("# of LoF", 
    #                   x_axis_ticks) +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14)) +
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=18),
          plot.subtitle=element_text(size=16),
          axis.title.x=element_blank()) +
    theme(legend.position="none") 
  
  print(p1)
  
}

dev.off()





#### Prepare values for datations  ---------------------------------

species_to_date <- 
  species_table %>%
  filter(Habitat == "troglomorphic_subterranean") %>%
  pull(Species)

individual_to_date <- 
  Individual_repertoire_df %>%
  filter(Species %in% species_to_date) %>%
  pull(Species_Ind) %>%
  unique()


#Extract the mean length of genes per species


gene_length_df <- 
  as.data.frame(
    Individual_repertoire_df %>%
      filter(Species_Ind %in% individual_to_date) %>%
      group_by(Species_Ind) %>%
      summarise(mean_CDS_length = mean(CDS_length))
  )


Individual_repertoire_df_summary_wide_NT_LOF_date <-
  left_join(Individual_repertoire_df_summary_wide_NT_LOF, gene_length_df, by="Species_Ind") %>%
  filter(! is.na(mean_CDS_length))


All_stats_datation_df <- 
  Individual_repertoire_df_summary_wide_NT_LOF_date %>%
  dplyr::select(Species_Ind, Total, Pseudogene, mean_CDS_length)

#List of genes pseudogenized in at-least 1 cave specimen .. 


all_cave_pseudogene_list <- 
  Individual_repertoire_df %>%
  filter(Species_Ind %in% individual_to_date) %>%
  filter(Gene_Type_ind == "Pseudogene") %>%
  pull(gene_name) %>%
  unique()


#For each species X , make a list of genes pseudogenized removing pseudo unique to the species X
# and then compute mean CDS and gene number ... 

Subset_stats_datation_df <- as.data.frame(NULL)
for(curr_ind in individual_to_date){
  
  all_pseudogene_list_woCurrsp <- 
    Individual_repertoire_df %>%
    filter(Species_Ind %in% individual_to_date) %>%
    filter(Species_Ind != curr_ind) %>%
    filter(Gene_Type_ind == "Pseudogene") %>%
    pull(gene_name) %>%
    unique()
  
  
  curr_species_nb_subset_genes <- 
    nrow(Individual_repertoire_df %>%
           filter(Species_Ind == curr_ind) %>%
           filter(gene_name %in% all_pseudogene_list_woCurrsp))
  
  
  curr_species_nb_subset_pseudogenes <- 
    nrow(Individual_repertoire_df %>%
           filter(Species_Ind == curr_ind) %>%
           filter(gene_name %in% all_pseudogene_list_woCurrsp) %>%
           filter(Gene_Type_ind == "Pseudogene"))
  
  
  curr_species_subset_length <- 
    Individual_repertoire_df %>%
    filter(Species_Ind == curr_ind) %>%
    filter(gene_name %in% all_pseudogene_list_woCurrsp) %>%
    pull(CDS_length) %>%
    mean()
  
  
  curr_stats <- 
    as.data.frame(
      cbind(curr_ind, curr_species_nb_subset_genes, 
            curr_species_nb_subset_pseudogenes, curr_species_subset_length)
    )
  
  colnames(curr_stats) <- 
    c("Species_Ind", "Subset_GeneNb", "Subset_PseudogeneNb", "Subset_MeanCDSLength")
  
  
  Subset_stats_datation_df <- rbind(Subset_stats_datation_df, curr_stats)
}



#### Perform datation of cavefishes - Whole repertoire - mutation rate = 10E-8  ---------------------------------



Datation_results_df <- as.data.frame(NULL)
for(curr_ind in individual_to_date){
  
  curr_gene_number <- 
    All_stats_datation_df %>%
    filter(Species_Ind == curr_ind) %>%
    pull(Total)
  
  curr_meanCDSlength <- 
    All_stats_datation_df %>%
    filter(Species_Ind == curr_ind) %>%
    pull(mean_CDS_length)  
  
  curr_Dn <- 
    All_stats_datation_df %>%
    filter(Species_Ind == curr_ind) %>%
    pull(Pseudogene)  
  
  curr_LoF_rate <- 
    LoF_rates_df %>%
    pull(Total)
  
  
  curr_Datation_results_df <- as.data.frame(seq(1, 3000000, 10))
  colnames(curr_Datation_results_df) <- c("curr_time")
  
  
  curr_Datation_results_df <- 
    curr_Datation_results_df %>%
    mutate(gene_number = curr_gene_number) %>%
    mutate(meanCDSlength = curr_meanCDSlength) %>%
    mutate(LoF_rate = as.numeric(curr_LoF_rate)) %>%
    mutate(Dn = curr_Dn)
  
  
  curr_Datation_results_df <- 
    as.data.frame(
      curr_Datation_results_df %>%
        rowwise() %>%
        mutate(
          curr_proba=(((factorial(gene_number))/((factorial(Dn)*(factorial(gene_number-Dn))))))*((1-exp((-10^-8)*LoF_rate*curr_time*meanCDSlength))^Dn)*(exp((-10^-8)*LoF_rate*curr_time*meanCDSlength))^(gene_number-Dn))
    )
  
  
  curr_Datation_results_df <- 
    curr_Datation_results_df %>% mutate(Species_Ind = curr_ind)
  
  colnames(curr_Datation_results_df) <-
    c("GenerationNb","GeneNb","mean_cds_length", "LoF_rate", "PseudoNb","probability", "Species_Ind")
  
  Datation_results_df <- 
    rbind(Datation_results_df, curr_Datation_results_df)
  
}


individual_colors <- 
  c("Typhlichthys_subterraneus.MLN_0272" = "#E69F00",
    "Typhlichthys_subterraneus.MLN_0051" = "#E69F00",
    "Typhlichthys_subterraneus.MLN_18" = "#E69F00",
    "Typhlichthys_subterraneus.UAIC_14148_01"  = "#E69F00",
    "Typhlichthys_subterraneus.MLN_0296" = "#E69F00",
    "Typhlichthys_eigenmanni.EN1" = "#D55E00",
    "Typhlichthys_eigenmanni.MLN_0404" = "#D55E00",
    "Speoplatyrhinus_poulsoni.IRGN_2EJSDHOURK" = "#CC79A7",
    "Amblyopsis_hoosieri.MLN_0242" = "#0072B2",
    "Amblyopsis_hoosieri.MLN_0243" = "#0072B2",
    "Amblyopsis_hoosieri.MLN_0246" = "#0072B2",
    "Troglichthys_rosae.C3" = "#000000",
    "Troglichthys_rosae.C5" = "#000000",
    "Troglichthys_rosae.L1"   = "#000000",
    "Amblyopsis_spelaea.YFTC_23868" = "#56B4E9",
    "Amblyopsis_spelaea.YFTC_23869" = "#56B4E9",
    "Amblyopsis_spelaea.YFTC_23875" = "#56B4E9"
  )


pdf("Datation_all_repertoire.pdf",width = 6.34,  height = 4.61)

Datation_results_df %>%
  ggplot(., aes(x=GenerationNb, y=probability, color=Species_Ind)) +
  geom_line(linetype = "solid") +
  xlim(0, 1500000) +
  scale_color_manual(values = individual_colors) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  theme(legend.position="none")


dev.off()

#Get the max proba for each species

MaxProba_gen <- 
  as.data.frame(
    Datation_results_df %>%
      group_by(Species_Ind) %>%
      slice_max(probability, n=1) %>%
      dplyr::select(GenerationNb, Species_Ind, probability)
  )
colnames(MaxProba_gen) <- c("Best_GenNb", "Species_Ind", "MaxProba")


#Get the 0.05 confidence intervals for each species

LowProba05_gen <- 
  as.data.frame(
    Datation_results_df %>%
      filter(probability > 0.05) %>%
      group_by(Species_Ind) %>%
      slice_min(GenerationNb, n = 1) %>%
      dplyr::select(GenerationNb, Species_Ind)
  )
colnames(LowProba05_gen) <- c("MinGen_p_05", "Species_Ind")

HighProba05_gen <- 
  as.data.frame(
    Datation_results_df %>%
      filter(probability > 0.05) %>%
      group_by(Species_Ind) %>%
      slice_max(GenerationNb, n = 1) %>%
      dplyr::select(GenerationNb, Species_Ind)
  )
colnames(HighProba05_gen) <- c("MaxGen_p_05", "Species_Ind")

Summary_datation <- left_join(MaxProba_gen, LowProba05_gen, by="Species_Ind")
Summary_datation <- left_join(Summary_datation, HighProba05_gen, by="Species_Ind")


gt(Summary_datation %>% arrange(Best_GenNb))


### Perform datation of cavefishes - Subset repertoire - mutation rate = 10E-8  ---------------------------------

#Subset_stats_datation_df



Datation_results_subset_df <- as.data.frame(NULL)
for(curr_ind in individual_to_date){
  
  curr_gene_number <- 
    Subset_stats_datation_df %>%
    filter(Species_Ind == curr_ind) %>%
    pull(Subset_GeneNb)
  
  curr_meanCDSlength <- 
    Subset_stats_datation_df %>%
    filter(Species_Ind == curr_ind) %>%
    pull(Subset_MeanCDSLength)  
  
  curr_Dn <- 
    Subset_stats_datation_df %>%
    filter(Species_Ind == curr_ind) %>%
    pull(Subset_PseudogeneNb)  
  
  curr_LoF_rate <- 
    LoF_rates_df %>%
    pull(Total)
  
  
  curr_Datation_results_df <- as.data.frame(seq(1, 3000000, 10))
  colnames(curr_Datation_results_df) <- c("curr_time")
  
  
  curr_Datation_results_df <- 
    curr_Datation_results_df %>%
    mutate(gene_number = as.numeric(curr_gene_number)) %>%
    mutate(meanCDSlength = as.numeric(curr_meanCDSlength)) %>%
    mutate(LoF_rate = as.numeric(curr_LoF_rate)) %>%
    mutate(Dn = as.numeric(curr_Dn))
  
  
  curr_Datation_results_df <- 
    as.data.frame(
      curr_Datation_results_df %>%
        rowwise() %>%
        mutate(
          curr_proba=(((factorial(gene_number))/((factorial(Dn)*(factorial(gene_number-Dn))))))*((1-exp((-10^-8)*LoF_rate*curr_time*meanCDSlength))^Dn)*(exp((-10^-8)*LoF_rate*curr_time*meanCDSlength))^(gene_number-Dn))
    )
  
  
  curr_Datation_results_df <- 
    curr_Datation_results_df %>% mutate(Species_Ind = curr_ind)
  
  colnames(curr_Datation_results_df) <-
    c("GenerationNb","GeneNb","mean_cds_length", "LoF_rate", "PseudoNb","probability", "Species_Ind")
  
  Datation_results_subset_df <- 
    rbind(Datation_results_subset_df, curr_Datation_results_df)
  
}


individual_colors <- 
  c("Typhlichthys_subterraneus.MLN_0272" = "#E69F00",
    "Typhlichthys_subterraneus.MLN_0051" = "#E69F00",
    "Typhlichthys_subterraneus.MLN_18" = "#E69F00",
    "Typhlichthys_subterraneus.UAIC_14148_01"  = "#E69F00",
    "Typhlichthys_subterraneus.MLN_0296" = "#E69F00",
    "Typhlichthys_eigenmanni.EN1" = "#D55E00",
    "Typhlichthys_eigenmanni.MLN_0404" = "#D55E00",
    "Speoplatyrhinus_poulsoni.IRGN_2EJSDHOURK" = "#CC79A7",
    "Amblyopsis_hoosieri.MLN_0242" = "#0072B2",
    "Amblyopsis_hoosieri.MLN_0243" = "#0072B2",
    "Amblyopsis_hoosieri.MLN_0246" = "#0072B2",
    "Troglichthys_rosae.C3" = "#000000",
    "Troglichthys_rosae.C5" = "#000000",
    "Troglichthys_rosae.L1"   = "#000000",
    "Amblyopsis_spelaea.YFTC_23868" = "#56B4E9",
    "Amblyopsis_spelaea.YFTC_23869" = "#56B4E9",
    "Amblyopsis_spelaea.YFTC_23875" = "#56B4E9"
  )


pdf("Datation_subset_repertoire.pdf",width = 6.34,  height = 4.61)

Datation_results_subset_df %>%
  ggplot(., aes(x=GenerationNb, y=probability, color=Species_Ind)) +
  geom_line(linetype = "solid") +
  xlim(0, 2200000) +
  scale_color_manual(values = individual_colors) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  theme(legend.position="none")


dev.off()

#Get the max proba for each species

MaxProba_gen <- 
  as.data.frame(
    Datation_results_subset_df %>%
      group_by(Species_Ind) %>%
      slice_max(probability, n=1) %>%
      dplyr::select(GenerationNb, Species_Ind, probability)
  )
colnames(MaxProba_gen) <- c("Best_GenNb", "Species_Ind", "MaxProba")


#Get the 0.05 confidence intervals for each species

LowProba05_gen <- 
  as.data.frame(
    Datation_results_subset_df %>%
      filter(probability > 0.05) %>%
      group_by(Species_Ind) %>%
      slice_min(GenerationNb, n = 1) %>%
      dplyr::select(GenerationNb, Species_Ind)
  )
colnames(LowProba05_gen) <- c("MinGen_p_05", "Species_Ind")

HighProba05_gen <- 
  as.data.frame(
    Datation_results_subset_df %>%
      filter(probability > 0.05) %>%
      group_by(Species_Ind) %>%
      slice_max(GenerationNb, n = 1) %>%
      dplyr::select(GenerationNb, Species_Ind)
  )
colnames(HighProba05_gen) <- c("MaxGen_p_05", "Species_Ind")

Summary_datation_subset <- left_join(MaxProba_gen, LowProba05_gen, by="Species_Ind")
Summary_datation_subset <- left_join(Summary_datation_subset, HighProba05_gen, by="Species_Ind")


gt(Summary_datation_subset %>% arrange(Best_GenNb))

