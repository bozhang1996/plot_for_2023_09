setwd("D:/doctor1HF/16s/results/exported-table/randomForest")
######16s########
library(microeco)
# use pipe operator in magrittr package
library(magrittr)
# set.seed is used to fix the random number generation to make the results repeatable
set.seed(2023)
# make the plotting background same with the tutorial
library(ggplot2)
library(ggforce)
theme_set(theme_bw())
otu_table_16S_1<-read.csv("16s.csv",row.names = 1)
#data(sample_info_16S)
sample_info_16S_1<-read.csv("SampleID.csv",row.names = 1)
taxonomy_table_16S_1<-read.csv("tax.csv",row.names = 1)
dataset <- microtable$new(sample_table = sample_info_16S_1, 
                          otu_table = otu_table_16S_1,
                          tax_table = taxonomy_table_16S_1)
dataset$tidy_dataset()
dataset$cal_abund()
class(dataset$taxa_abund)
dir.create("taxa_abund")
dataset$save_abund(dirpath = "taxa_abund")

t1 <- trans_diff$new(dataset = dataset, method = "rf", group = "Group", taxa_level = "Genus")

# plot the MeanDecreaseGini bar
# group_order is designed to sort the groups
g1 <- t1$plot_diff_bar(use_number = 1:20)
# plot the abundance using same taxa in g1
g2 <- t1$plot_diff_abund(select_taxa = t1$plot_diff_bar_taxa)
# now the y axis in g1 and g2 is same, so we can merge them
# remove g1 legend; remove g2 y axis text and ticks
g1 <- g1 + theme(legend.position = "none")
g2 <- g2 + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
p1 <- gridExtra::grid.arrange(g1, g2, ncol = 2, nrow = 1, widths = c(2, 1.7))

ggsave("16s_g2.pdf",p1,path = "./",units = "in",width = 9.5,height = 7)
