library(igraph)
library(networkD3)
library(tidyverse)
library(AnnotationDbi)
library(data.table)
##########################

bplex <- read_tsv("Downloads/BioPlex_293T_Network_10K_Dec_2019 (1).tsv")

edges.47 <- filter(bplex, `SymbolB` == "C6orf47")  %>% select(1, 2, 5, 6, 7, 8, 9)

# MAD2L1, HTR2C AND C16ORF58 HAVE BEEN REMOVED!!!!!!!!!!!!!!!!
edges.else <- filter(bplex, `SymbolA` == "ACP2" |
                       `SymbolB` == "ACP2" |
                       `SymbolA` == "TACR3" |
                       `SymbolB` == "TACR3" |
                       `SymbolA` == "C5AR1" |
                       `SymbolB` == "C5AR1" |
                       `SymbolA` == "SEC62" |
                       `SymbolB` == "SEC62") %>%
  select(1, 2, 5, 6, 7, 8, 9)


  
edges.total <- add_row(edges.47, edges.else) %>% 
  as.data.frame()

edges.d3 <- data.frame(GeneA = edges.total$GeneA - 40,
                       GeneB = edges.total$GeneB - 40,
                       SymbolA = edges.total$SymbolA,
                       SymbolB = edges.total$SymbolB)

gT.df <- rbind(select(edges.total, 1, 3) %>% 
                 rename("Gene" = 'GeneA', "name" = 'SymbolA'),
               select(edges.total, 2, 4) %>%
                 rename("Gene" = 'GeneB', "name" = 'SymbolB')) %>%
  distinct()

nodesT.d3 <- data.frame(Gene = gT.df$Gene, Symbol = gT.df$name)
group.test <- data.frame(group = sample(x = c(1:5), size =263, replace = TRUE))
nodesTest.d3 <- cbind(nodesT.d3, group.test) 


forceNetwork(Links = edges.d3, Nodes = nodesTest.d3, Source="GeneA", Target="GeneB",
             NodeID = "Symbol", linkWidth = 1,linkColour = "#afafaf", 
             fontSize=12, opacity = 0.8, Group = "group")

net.47 <- graph_from_data_frame(edges.47, directed = T)
net.else <- graph_from_data_frame(edges.else, directed = T)
net.total <- graph_from_data_frame(d=edges.total, vertices = nodes.total, directed = T)



##########################

nodes <- read.csv("Downloads/netscix2016 (2)/Dataset1-Media-Example-NODES.csv", 
                  header = T, 
                  as.is = T)
links <- read.csv("Downloads/netscix2016 (2)/Dataset1-Media-Example-EDGES.csv",
                  header = T,
                  as.is = T)

links <- aggregate(links[,3], links[,-3], sum)
links <- links[order(links$from, links$to),]
colnames(links)[4] <- "weight"
rownames(links) <- NULL

net <- graph_from_data_frame(d=links, vertices=nodes, directed=T) 
plot(net, edge.arrow.size = .4, vertex.label = NA)
net <- simplify(net, remove.multiple = F, remove.loops = T)

plot(net, edge.arrow.size = .4, edge.curved = .1)

plot(net, edge.arrow.size = .2, edge.curved = 0, 
     vertex.color = "yellow", vertex.frame.color = "blue",
     vertex.label = V(net)$media, vertex.label.color = "black",
     vertex.label.cex = .7)

colrs <- c("gray50", "tomato", "gold")
V(net)$color <- colrs[V(net)$media.type]
V(net)$size <- V(net)$audience.size*0.7
V(net)$label.color <- "black"
V(net)$label <- NA
E(net)$width <- E(net)$weight/6
E(net)$arrow.size <- 0.2
E(net)$edge.color <- "gray80"
E(net)$width <- 1+E(net)$weight/12

plot(net)

plot(net, edge.color = "orange", vertex.color = "gray50" )

legend(x = -1.5, y=-1.1, c("Newspaper", "Television", "Online News"), pch = 21,
       col ="#777777", pt.bg = colrs, pt.cex = 2, cex = 0.8, bty = "n", ncol = 1)

plot(net, vertex.shape = "none", vertex.label = V(net)$media,
     vertex.label.font = 2, vertex.label.color = "gray40",
     vertex.label.cex = .7, edge.color = "gray85")

edge.start <- ends(net, es=E(net), names = F)[,1]
edge.col <- V(net)$color[edge.start]
plot(net, edge.color = edge.col, edge.curved = .1)

cut.off <- mean(links$weight)

net.sp <- delete_edges(net, E(net)[weight<cut.off])
plot(net.sp)

E(net)$width <- 1.5
plot(net, edge.color = c("dark red", "slategrey")[(E(net)$type == "hyperlink")+1],
     vertex.color = "gray40", layout = layout.circle)

net.m <- net - E(net)[E(net)$type == "hyperlink"]
net.h <- net - E(net)[E(net)$type == "mention"]

nodes.d3 <- cbind(idn=factor(nodes$media, levels=nodes$media), nodes)
#############################

nodes2 <- read.csv("Downloads/netscix2016 (2)/Dataset2-Media-User-Example-NODES.csv", 
                   header = T, 
                   as.is = T)
links2 <- read.csv("Downloads/netscix2016 (2)/Dataset2-Media-User-Example-EDGES.csv", 
                   header = T, 
                   row.names = 1)

links2 <- as.matrix(links2)

net2 <- graph_from_incidence_matrix(links2)

net2.bp <- bipartite.projection(net2)

plot(net2.bp$proj1, 
     vertex.label.colr = "black", 
     vertex.label.dist = 1,
     vertex.size = 7,
     vertex.label = nodes2$media[!is.na(nodes2$media.type)])

plot(net2.bp$proj2, 
     vertex.label.color = "black",
     vertex.label.dist = 1,
     vertex.size = 7,
     vertex.label = nodes2$media[is.na(nodes2$media.type)])


V(net2)$color <- c("steel blue", "orange")[V(net2)$type+1]
V(net2)$shape <- c("square", "circle")[V(net2)$type+1]
V(net2)$label <- ""
V(net2)$label[V(net2)$type==F] <- nodes2$media[V(net2)$type==F] 
V(net2)$label.cex=.4
V(net2)$label.font=2
plot(net2, vertex.label.color = "white", vertex.size = (2 - V(net2)$type)*8)

plot(net2, vertex.label = NA, vertex.size = 7, layout = layout_as_bipartite)

plot(net2, vertex.shape = "none", vertex.label = nodes2$media,
     vertex.label.color= V(net2)$color, vertex.label.font = 2.5, 
     vertex.label.cex = 0.6, edge.color="pink", edge.width = 2)
#############################

net.bg <- sample_pa(80)

V(net.bg)$size <- 8
V(net.bg)$frame.color <- "white"
V(net.bg)$color <- "orange"
V(net.bg)$label <- ""
E(net.bg)$arrow.mode <- 0

plot(net.bg)

plot(net.bg, layout = layout_randomly)

l <- layout_in_circle(net.bg)

l <- layout_on_sphere(net.bg)

l <- layout_with_fr(net.bg)
plot(net.bg, layout = l)

par(mfrow = c(2,2), mar = c(0,0,0,0))
plot(net.bg, layout = layout_with_fr)
plot(net.bg, layout = layout_with_fr)
plot(net.bg, layout = l)
plot(net.bg, layout = l)
dev.off()

l <- layout_with_fr(net.bg)
l <- norm_coords(l, ymin = -1, ymax = 1, xmin = -1, xmax = 1)
plot(net.bg, rescale = F, layout = l)

layouts <- grep("^layout_", ls("package:igraph"), value = TRUE)[-1]
layouts <- layouts[!grepl("bipartite|merge|norm|sugiyama|tree", layouts)]

par(mfrow = c(3,3), mar = c(1,1,1,1))

for(layout in layouts) {
  print(layout)
  l <- do.call(layout, list(net))
  plot(net, edge.arrow.mode = 0, layout = l, main = layout)
}









