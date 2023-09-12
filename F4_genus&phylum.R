setwd("D:/submission/others/energy/z25/Figure5/genus")
library(microeco)
#?microtable 
library(magrittr)
library(ggplot2)
library(randomcoloR)
theme_set(theme_bw())
#prepare data
otu_table_16S_1<-read.csv("genus_z25.csv",row.names = 1)
#data(sample_info_16S)
sample_info_16S_1<-read.csv("SampleIDZ25.csv",row.names = 1)
taxonomy_table_16S_1<-read.csv("tax.csv",row.names = 1)
dataset <- microtable$new(sample_table = sample_info_16S_1, 
                          otu_table = otu_table_16S_1,
                          tax_table = taxonomy_table_16S_1)
dataset$sample_table$Group <-  factor(dataset$sample_table$Group, levels = c("LF","HF","Z25"))

dataset$tidy_dataset()
dataset$cal_abund()
class(dataset$taxa_abund)
dir.create("taxa_abund")
dataset$save_abund(dirpath = "taxa_abund")


color_compaired <- c("#ed1299", "#09f9f5","#373bbf", "#246b93", "#cc8e12", "#d561dd", "#c93f00", "#ddd53e",
                     "#4aef7b", "#e86502", "#9ed84e","#03827f", "#99db27", "#e07233","#8249aa", "#ff523f", 
                     "#ce2523", "#f7aa5d", "#cebb10",  "#931635",  "#a1ce4c", "#ef3bb6", "#d66551",
                     "#1a918f", "#ff66fc", "#2927c4", "#7149af" ,"#57e559" ,"#8e3af4" ,"#f9a270" ,"#22547f", "#db5e92",
                     "#edd05e", "#6f25e8", "#0dbc21", "#280f7a", "#6373ed", "#5b910f" ,"#7b34c1" ,"#0cf29a" ,"#d80fc1",
                     "#dd27ce", "#07a301", "#167275", "#391c82", "#2baeb5","#925bea", "#63ff4f")

palette <- distinctColorPalette(30)

t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 10, groupmean = "Group")
t1$plot_bar(use_alluvium = FALSE, clustering = TRUE,  xtext_size = 20,color_values = c(color_compaired), ytitle_size = 30)
ggsave("phylumplotaverage.pdf", height = 12, width = 12)
t2 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 50, groupmean = "Group")
t2$plot_bar(use_alluvium = TRUE, clustering = TRUE,  xtext_size = 20,color_values = c(color_compaired), ytitle_size = 30)
ggsave("genusalluviumplotaverage.pdf", height = 12, width = 12)