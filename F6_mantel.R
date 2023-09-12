setwd("D:/submission/others/energy/z25/Figure6/MANTEL5")
library(ggcor)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(ggsci)
library(vegan)
library(linkET)
env<-read.csv("genus.csv",row.names = 1,header = T) #actual otu table
otu<-read.csv("met2.csv",row.names = 1,header = T) #actual metabolites table
otu <- as.data.frame(otu)
env <- as.data.frame(env)


head(otu)
head(env)

specmethod <-  dist_func(.FUN = "vegdist", method = "euclidean")
envmethod <- dist_func(.FUN = "vegdist", method = "bray")

# all_dist_method <- c(
#   "manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski",
#   "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup",
#   "binomial", "chao", "cao", "mahalanobis", "chisq", "chord", "hellinger",
#   "aitchison", "robust.aitchison", "maximum", "binary", "minkowski"
# )
mantel <- mantel_test(otu, env,
            spec_select = list(TRP  = 1:1, KYN = 2:5, INDOLE = 6:12, HT = 13:15),
            spec_dist = specmethod, #metabolites, method was set as euclidean
            env_dist = envmethod  #otu, method was set as bray
            )

mantel2 <- mantel %>%
  mutate(r = cut(r, breaks = c(-Inf, 0.1, 0.2, 0.4, Inf),
                    labels = c("< 0.1", "0.1 - 0.2", "0.2 - 0.4", ">= 0.4")),
         p = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                    labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
head(mantel2)

####graph
cor2 <- correlate(env)
corr2 <- cor2 %>% as_md_tbl()

colormap1 <- colorRampPalette(colors = c("purple","white","yellow"))(30)
colormap1 <- colorRampPalette(colors =c("#0dbc21","white", "#ff523f"))(30) 
p4 <- qcorrplot(cor2,
                grid_col = "grey50",
                grid_size = 0.2,
                type = "upper",
                diag = FALSE) +
  geom_star() +
  scale_fill_gradientn(colours = colormap1,
                       limits = c(-1, 1))
p4
p5 <- p4 +
  geom_mark(size = 6,
            only_mark = T,
            sig_level = c(0.05, 0.01, 0.001),
            sig_thres = 0.05,
            colour = '#391c82')
p5
p6 <- p5 +
  geom_couple(data = mantel2,
              aes(colour = p, size = r),
              curvature = nice_curvature())
p6
cols <-c("#ed1299", "#09f9f5", "#246b93", "#cc8e12", "#d561dd", "#c93f00", "#ddd53e",
"#4aef7b", "#e86502", "#9ed84e", "#8249aa", "#99db27", "#e07233", "#ff523f",
"#ce2523", "#f7aa5d", "#cebb10", "#03827f", "#931635", "#373bbf", "#a1ce4c", "#ef3bb6", "#d66551",
"#1a918f", "#ff66fc", "#2927c4", "#7149af" ,"#57e559" ,"#8e3af4" ,"#f9a270" ,"#22547f", "#db5e92",
"#edd05e", "#6f25e8", "#0dbc21", "#280f7a", "#6373ed", "#5b910f" ,"#7b34c1" ,"#0cf29a" ,"#d80fc1",
"#dd27ce", "#07a301", "#167275", "#391c82", "#2baeb5","#925bea", "#63ff4f")

p7 <- p6 +
  scale_size_manual(values = c(0.5, 1, 2, 3)) + 
  scale_colour_manual(values = c("#c93f00", "#a1ce4c","#a1ce4c","#d80fc1")) +
  guides(size = guide_legend(title = "Mantel r",
                             override.aes = list(colour = "grey35"),
                             order = 2),
         colour = guide_legend(title = "Mantel p",
                               override.aes = list(size = 3),
                               order = 1),
         fill = guide_colorbar(title = "Pearson r", order = 3))
  
p7

ggsave("p7.pdf", width = 15, height = 15)
