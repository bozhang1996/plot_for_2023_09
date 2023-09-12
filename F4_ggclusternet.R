#ggclusternet
setwd("D:/submission/others/energy/z25/Figure5/genus/inap/plothf")
library(phyloseq)
library(igraph)
library(network)
library(sna)
library(tidyverse)
library(ggClusterNet)
library(ggrepel)
metadata = read.csv("SampleIDlf.csv",row.names = 1, header = T)
otutab = read.csv("ggclusterlf.csv",row.names = 1, header = T)
taxonomy = read.csv("tax.csv",row.names = 1, header = T)
ps = phyloseq(sample_data(metadata),
              otu_table(as.matrix(otutab), taxa_are_rows=TRUE),
              tax_table(as.matrix(taxonomy))#,
              # phy_tree(tree),
              # refseq(rep)
)

library(randomcoloR)
palette <- distinctColorPalette(70)

tab.r = network.pip(
  ps = ps,
  N = 50,
  # ra = 0.05,
  big = FALSE,
  select_layout = FALSE,
  layout_net = "model_maptree2",#model_igraph2,model_igraph,model_maptree2,model_Gephi.2,fruchtermanreingold
  r.threshold = 0.6,
  p.threshold = 0.05,
  maxnode = 2,
  method = "pearson",
  label = T,
  lab = "elements",
  group = "group",
  fill = "Phylum",
  zipi = TRUE,
  ram.net = TRUE,
  clu_method = "cluster_fast_greedy",
  step = 100,
  R=10,
  ncpus = 1
)

plot = tab.r[[1]]
plot
p0 = plot[[1]]  + scale_fill_manual(values = c(palette))
p0

#-----others graph#-----
p0.1 = plot[[2]]
p0.1 
p0.2 = plot[[3]]
p0.2
dat = tab.r[[2]]
cortab = dat$net.cor.matrix$cortab
dat = module.compare.net.pip(
  ps = NULL,
  corg = cortab,
  degree = TRUE,
  zipi = FALSE,
  r.threshold= 0.8,
  p.threshold=0.05,
  method = "pearson",
  padj = F,
  n = 3)
res = dat[[1]]
head(res)
dat = tab.r[[2]]
node = dat$net.cor.matrix$node
edge = dat$net.cor.matrix$edge
head(edge)
head(node)

p <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2
),
data = edge, size = 0.03,alpha = 0.5,
color = "yellow") +
  geom_point(aes(X1, X2,
                 fill = Phylum,
                 size = igraph.degree),
             pch = 21, data = node,color = "gray40") +
  facet_wrap(.~ label,scales="free_y",nrow = 3) +
  # geom_text_repel(aes(X1, X2,label = elements),pch = 21, data = nodeG) +
  # geom_text(aes(X1, X2,label = elements),pch = 21, data = nodeG) +
  scale_colour_manual(values = c("#6D98B5","#D48852")) +
  scale_fill_hue()+
  scale_size(range = c(0.8, 5)) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5)
  ) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()
  ) +
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())

p

dat.z = dat$zipi.data
head(dat.z)
x1<- c(0, 0.62,0,0.62)
x2<- c( 0.62,1,0.62,1)
y1<- c(-Inf,2.5,2.5,-Inf)
y2 <- c(2.5,Inf,Inf,2.5)
lab <- c("peripheral",'Network hubs','Module hubs','Connectors')
roles.colors <- c("#E6E6FA","#DCDCDC","#F5FFFA", "#FAEBD7")
tab = data.frame(x1 = x1,y1 = y1,x2 = x2,y2 = y2,lab = lab)
tem = dat.z$group %>% unique() %>% length()
for ( i in 1:tem) {
  if (i == 1) {
    tab2 = tab
  } else{
    tab2 = rbind(tab2,tab)
  }
}

p <- ggplot() +
  geom_rect(data=tab2,
            mapping=aes(xmin=x1,
                        xmax=x2,
                        ymin=y1,
                        ymax=y2,
                        fill = lab))+
  guides(fill=guide_legend(title="Topological roles")) +
  scale_fill_manual(values = roles.colors)+
  geom_point(data=dat.z,aes(x=p, y=z,color=module)) + theme_bw()+
  guides(color= F) +
  ggrepel::geom_text_repel(data = dat.z,
                           aes(x = p, y = z,
                               color = module,label=label),size=4)+
  # facet_wrap(.~group) +
  facet_grid(.~ group, scale='free') +
  theme(strip.background = element_rect(fill = "white"))+
  xlab("Participation Coefficient")+ylab(" Within-module connectivity z-score")
p


result = corMicro (ps = ps,
                   N = 150,
                   method.scale = "TMM",
                   r.threshold=0.8,
                   p.threshold=0.05,
                   method = "pearson"
                   
                   
)

cor = result[[1]]
head(cor)

ps_net = result[[3]]

otu_table = ps_net %>% 
  vegan_otu() %>%
  t() %>%
  as.data.frame()

#based on ps of microbiota
result = corMicro (ps = ps,
                   N = 200,
                   method.scale = "TMM",
                   r.threshold=0.0,
                   p.threshold=0.05,
                   method = "pearson"
)

cor = result[[1]]
head(cor)

ps_net = result[[3]]

otu_table = ps_net %>% 
  vegan_otu() %>%
  t() %>%
  as.data.frame()
tax = ps_net %>% vegan_tax() %>%
  as.data.frame()
tax$filed = tax$Phylum
group2 <- data.frame(ID = row.names(tax),group = tax$Phylum)
group2$group  =as.factor(group2$group)
result2 = PolygonClusterG (cor = cor,nodeGroup =group2 )
node = result2[[1]]

nodes = nodeadd(plotcord =node,otu_table = otu_table,tax_table = taxonomy)

edge = edgeBuild(cor = cor,node = node)

pnet <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)),
                                data = edge, linewidth = 0.5) +
  geom_point(aes(X1, X2,fill = Phylum,size = mean),pch = 21, data = nodes) +
  scale_colour_brewer(palette = "Set1") +
  scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
  # labs( title = paste(layout,"network",sep = "_"))+
  # geom_text_repel(aes(X1, X2,label=Phylum),size=4, data = plotcord)+
  # discard default grid + titles in ggplot2
  theme(panel.background = element_blank()) +
  # theme(legend.position = "none") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
pnet

result2 = PolygonRrClusterG (cor = cor,nodeGroup =group2 )
node = result2[[1]]

nodes = nodeadd(plotcord =node,otu_table = otu_table,tax_table = taxonomy)

edge = edgeBuild(cor = cor,node = node)

pnet <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)),
                                data = edge, size = 0.5) +
  geom_point(aes(X1, X2,fill = Phylum,size = mean),pch = 21, data = nodes) +
  scale_colour_brewer(palette = "Set1") +
  scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
  # labs( title = paste(layout,"network",sep = "_"))+
  # geom_text_repel(aes(X1, X2,label=Phylum),size=4, data = plotcord)+
  # discard default grid + titles in ggplot2
  theme(panel.background = element_blank()) +
  # theme(legend.position = "none") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
pnet

#------model.gephi2#-----
cols <- c("#ed1299", "#09f9f5", "#6373ed", "#4aef7b", "#1a918f", "#c93f00", "#ddd53e",
          "#4aef7b", "#e86502", "#9ed84e", "#8249aa", "#99db27", "#e07233", "#ff523f",
          "#ce2523", "#f7aa5d", "#cebb10", "#03827f", "#931635", "#373bbf", "#a1ce4c", "#ef3bb6", "#d66551",
          "#1a918f", "#ff66fc", "#2927c4", "#7149af" ,"#57e559" ,"#8e3af4" ,"#f9a270" ,"#22547f", "#db5e92",
          "#edd05e", "#6f25e8", "#0dbc21", "#280f7a", "#6373ed", "#5b910f" ,"#7b34c1" ,"#0cf29a" ,"#d80fc1",
          "#dd27ce", "#07a301", "#167275", "#391c82", "#2baeb5","#925bea", "#63ff4f")
show_col(cols) # 以图形形式展示颜色

result = network(ps = ps,
                 N = 200,
                 layout_net = "randomClusterG",
                 r.threshold=0.5,
                 p.threshold=0.05,
                 label = F,  
                 group = "group",
                 fill = "Phylum",
                 zipi = TRUE)
p = result[[1]] + scale_fill_manual(values = c(cols))
p
data = result[[2]]
data
