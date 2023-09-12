setwd("D:/submission/others/energy/z25/Figure6/netcir")
library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)
library(reshape2)
library(ggplot2)
library(psych)
library(ggsci)
library(scales)

data<- read.csv("data.csv",header = TRUE,row.names = 1,check.names=FALSE,comment.char = "")

type = read.csv("type.csv",header = TRUE,check.names = FALSE)
head(type)

cor <- psych::corr.test(data, use = "pairwise",
                        method="spearman",
                        adjust="fdr", 
                        alpha=0.10,
                        ci=T) 
cor.r <- data.frame(cor$r) # r data
cor.p <- data.frame(cor$p) # p data
colnames(cor.r) = rownames(cor.r)
colnames(cor.p) = rownames(cor.p) 

# write.csv(cor.r,"cor.r.csv",quote = FALSE,col.names = NA,row.names = TRUE) # save data
# write.csv(cor.p,"cor.p.csv",quote = FALSE,col.names = NA,row.names = TRUE) # save data

head(cor.r)
head(cor.p)


cor.r$from = rownames(cor.r) 
cor.p$from = rownames(cor.p)
p = cor.p %>% 
  gather(key = "to", value = "p", -from) %>%
  data.frame()

cor.data = cor.r %>% 
  gather(key = "to", value = "r", -from) %>%
  data.frame() %>%
  left_join(p, by=c("from","to")) %>%
  filter(p <= 0.05, from != to) %>%
  mutate(
    linecolor = ifelse(r > 0,"positive","negative"), 
    linesize = abs(r) 
  ) 
head(cor.data)



# degree calculated

c(as.character(cor.data$from),as.character(cor.data$to)) %>%
  as_tibble() %>%
  group_by(value) %>%
  summarize(n=n()) -> vertices
colnames(vertices) <- c("name", "n")

## add plot cluster
vertices <- vertices %>%
  select(-n) %>% 
  left_join(type,by="name")

## sort plot cluster
vertices$type <- factor(vertices$type,
                        levels = c("Firmicutes","Bacteroidetes","Actinobacteria","Deferribacteres","Verrucomicrobia"
                        ,"Proteobacteria","TRP","5HT","Indole" ,"Kynurenine" ))
vertices <- vertices %>%
  arrange(type,group) 

dim(vertices)
head(vertices)


## prepare graph
graph <- graph_from_data_frame(cor.data, vertices = vertices, directed = FALSE )
graph 

## 
is.simple(graph) 
E(graph)$weight <- 1 
graph <- igraph::simplify(graph,edge.attr.comb = "first")
is.simple(graph)
E(graph)$weight <- 1
is.weighted(graph)
graph
V(graph)$degree <- degree(graph)
net.data  <- igraph::as_data_frame(graph, what = "both")$edges 
vertices  <- igraph::as_data_frame(graph, what = "both")$vertices 
unique(vertices$group) # ensure group data

mycolor = pal_d3("category20",alpha = 1)(20)
mycolor
cols = mycolor[c(1:7,9:13,14)]
cols <- c("#ed1299", "#09f9f5", "#246b93", "#cc8e12", "#d561dd", "#c93f00", "#ddd53e",
          "#4aef7b", "#e86502", "#9ed84e", "#8249aa", "#99db27", "#e07233", "#ff523f",
          "#ce2523", "#f7aa5d", "#cebb10", "#03827f", "#931635", "#373bbf", "#a1ce4c", "#ef3bb6", "#d66551",
          "#1a918f", "#ff66fc", "#2927c4", "#7149af" ,"#57e559" ,"#8e3af4" ,"#f9a270" ,"#22547f", "#db5e92",
          "#edd05e", "#6f25e8", "#0dbc21", "#280f7a", "#6373ed", "#5b910f" ,"#7b34c1" ,"#0cf29a" ,"#d80fc1",
          "#dd27ce", "#07a301", "#167275", "#391c82", "#2baeb5","#925bea", "#63ff4f")
show_col(cols) 

color = c("positive" ="#D62728FF","negative" ="#2CA02CFF")

bar1 = nrow(vertices[vertices$group == "metabolites",])
bar1 # metabolites
bar2 = nrow(vertices[vertices$group != "metabolites",])
bar2 # otu

## resort id
vertices$id =1
vertices$id[vertices$group == "metabolites"] = seq(1, bar1)
vertices$id[vertices$group != "metabolites"] = seq(1, bar2)
head(vertices)

## calculate angel
vertices$var_angel <- NA
vertices$var_angel[vertices$group == "metabolites"] <- 360 * (vertices$id[vertices$group == "metabolites"]-0.5)/bar1 
vertices$var_angel[vertices$group != "metabolites"] <- 360 * (vertices$id[vertices$group != "metabolites"]-0.5)/bar2

## re-calculate angel
vertices$x <- ifelse(vertices$group == "metabolites",cos(vertices$var_angel)/2,cos(vertices$var_angel))
vertices$y <- ifelse(vertices$group == "metabolites",sin(vertices$var_angel)/2,sin(vertices$var_angel))
head(vertices)

vertices$hjust <- ifelse(vertices$var_angel > 180, 1, 0)
vertices$angle <- ifelse(vertices$var_angel > 180, 90-vertices$var_angel+180, 90-vertices$var_angel)

## inner cycle
vertices$hjust[vertices$group != "metabolites"] = "center" 
vertices$angle[vertices$group != "metabolites"] = 0

net.data <- net.data %>% 
  left_join(vertices,by=c('from'='name')) 

graph <- graph_from_data_frame(net.data, vertices = vertices, directed = FALSE)
graph

#write.csv(vertices,"relation.csv", sep = "\t") #export data
#write.csv(net.data,"netdata0911.csv",quote = FALSE,col.names = NA,row.names = FALSE) #export data

layout1 <- create_layout(
  induced.subgraph(graph,c(V(graph)$name %in% vertices$name[vertices$group != "metabolites"])), 
  layout = 'linear', circular = TRUE)
layout1

## arrangement
layout2 <- create_layout(
  induced.subgraph(graph,c(V(graph)$name %in% vertices$name[vertices$group == "metabolites"])), 
  layout = 'linear', circular = TRUE)
layout2[1:2] <- layout2[1:2]/2 

layout <- create_layout(graph,layout = 'linear', circular = T)
head(layout)

## coordinates
layout$x <- rbind(layout1,layout2)$x
layout$y <- rbind(layout1,layout2)$y
head(layout)

Correlation <- guide_legend(title="Correlation",
                            direction="vertical",
                            order=2,
                            ncol=1,
                            byrow=FALSE,
                            title.theme = element_text(
                              size = 14,
                              face = "bold",
                              colour = "black"),
                            label.theme = element_text(
                              size = 12,colour = "black")
)

width_legend <- guide_legend(title="Pearson |r|",
                             direction="vertical",
                             order=1,
                             ncol=1,
                             byrow=FALSE,
                             title.theme = element_text(
                               size = 14,
                               face = "bold",
                               colour = "black"),
                             label.theme = element_text(
                               size = 12,colour = "black")
)

size_legend <- guide_legend(title="Degree",
                            direction="vertical",
                            order=3,
                            ncol=2, 
                            byrow=FALSE,
                            title.theme = element_text(
                              size = 14,
                              face = "bold",
                              colour = "black"),
                            label.theme = element_text(
                              size = 12,colour = "black")
)

fill_legend <- guide_legend(title="Genus",
                            direction="vertical",
                            order=4,
                            ncol=2,
                            byrow=FALSE,
                            title.theme = element_text(
                              size = 14,
                              face = "bold",
                              colour = "black"),
                            label.theme = element_text(
                              size = 12,
                              colour = "black",face = "italic")
)

set_graph_style(plot_margin = margin(0,0,1,0))

net.cir2 <- ggraph(layout) +
  geom_edge_link(aes(edge_colour = as.factor(linecolor),
                    edge_width = abs(r)
  ),
  edge_alpha=0.6) + 
  scale_edge_colour_manual(values=color,
                           breaks = c("positive","negative"),
                           guide = Correlation) +
  scale_edge_width(
    breaks = seq(0.2,1,0.2),
    label = seq(0.2,1,0.2),
    range = c(0.2,1),
    guide = width_legend )+
  geom_node_point(aes(size=degree, fill=as.factor(type)),
                  shape=21,alpha=1) +
  scale_fill_manual(values=cols,
                    guide = fill_legend) +
  scale_size(
    breaks = seq(3,max(vertices$degree),3),
    label = seq(3,max(vertices$degree),3),
    range = c(5, 15),
    guide = size_legend)+
  geom_node_text(aes(x = x, y=y,
                     label= ifelse(V(graph)$group == "metabolites",as.character(name),""),
                     angle=angle,hjust=hjust,
                     #color = as.factor(group)
  ),
  color = "black",
  size=4.75,# 14pt
  show.legend = FALSE) +
  geom_node_text(aes(x = x*1.4, y=y*1.4,
  label= ifelse(V(graph)$group != "metabolites",as.character(name),""),
  angle=angle,hjust=hjust,
      ),
  color = "black",
  size=4.75,# 14pt
  show.legend = FALSE) +
  scale_color_manual(values = c( "#ff66fc", "#2927c4", "#7149af" ,"#57e559" ,"#8e3af4" ,"#f9a270" ,"#22547f", "#db5e92")) +
  coord_fixed()+
  theme_graph()
net.cir2

cairo_pdf("netcir2.pdf",height=12,width=12,family="Times")
print(net.cir2) 
dev.off()
