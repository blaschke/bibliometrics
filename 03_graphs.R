# Publication authorship
# 
# Import, convert, and tidy bibtex data downloaded from SCOPUS
#
# Copyright (C) 2022 by Dr. Steffen Blaschke (sbl.msc@cbs.dk).
# This work is licensed under a Creative Commons Attribution 4.0 International License (CC BY).
# https://creativecommons.org/licenses/by/4.0/

# load libraries
library(igraph)

# bibliographic coupling matrix
bib_coup <- function(bib) {
  bc <- Matrix::t(cocMatrix(bib, Field = "CR", type = "sparse", n = NULL, sep = ";"))
  bib.coup.m <- Matrix::crossprod(bc, bc)
  bib.coup.m <- normalizeSimilarity(bib.coup.m, type = "salton")
  return(bib.coup.m)
}

# publication authorship matrix
pub_auth <- function(bib) {
  pa <- Matrix::t(cocMatrix(bib, Field = "AU", type = "sparse", n = NULL, sep = ";"))
  pub.auth.m <- Matrix::crossprod(pa, pa)
  pub.auth.m <- normalizeSimilarity(pub.auth.m, type = "salton")
  return(pub.auth.m)
}

## cosine similarity matrix for abstracts
cos_sim <- function(bib) {
  ab <- Matrix::t(cocMatrix(bib, Field = "AB", type = "sparse", n = NULL, sep = " "))
  cos.sim.m <- Matrix::crossprod(ab, ab)
  cos.sim.m <- normalizeSimilarity(cos.sim.m, type = "salton")
  return(cos.sim.m)
}

# compute matrices
# accounting
acc.bib.coup.m <- bib_coup(acc.bib)
acc.pub.auth.m <- pub_auth(acc.bib)
acc.cos.sim.m <- cos_sim(acc.bib)
# astronomy
ast.bib.coup.m <- bib_coup(ast.bib)
ast.pub.auth.m <- pub_auth(ast.bib)
ast.cos.sim.m <- cos_sim(ast.bib)
# gastroenterology
gas.bib.coup.m <- bib_coup(gas.bib)
gas.pub.auth.m <- pub_auth(gas.bib)
gasag.cos.sim.m <- cos_sim(gas.bib)

# convert to igraphs
# accounting
acc.bib.coup.g <- graph_from_adjacency_matrix(acc.bib.coup.m, mode = "upper", weighted = TRUE, diag = FALSE)
acc.pub.auth.g <- graph_from_adjacency_matrix(acc.pub.auth.m, mode = "upper", weighted = TRUE, diag = FALSE)
# astronomy
ast.bib.coup.g <- graph_from_adjacency_matrix(ast.bib.coup.m, mode = "upper", weighted = TRUE, diag = FALSE)
ast.pub.auth.g <- graph_from_adjacency_matrix(ast.pub.auth.m, mode = "upper", weighted = TRUE, diag = FALSE)
# gastroenterology
gas.bib.coup.g <- graph_from_adjacency_matrix(gas.bib.coup.m, mode = "upper", weighted = TRUE, diag = FALSE)
gas.pub.auth.g <- graph_from_adjacency_matrix(gas.pub.auth.m, mode = "upper", weighted = TRUE, diag = FALSE)

# find isolates
# accounting
acc.bib.coup.g.isolates <- length(V(acc.bib.coup.g)[degree(acc.bib.coup.g) == 0])
acc.pub.auth.g.isolates <- length(V(acc.pub.auth.g)[degree(acc.pub.auth.g) == 0])
# astronomy
ast.bib.coup.g.isolates <- length(V(ast.bib.coup.g)[degree(ast.bib.coup.g) == 0])
ast.pub.auth.g.isolates <- length(V(ast.pub.auth.g)[degree(ast.pub.auth.g) == 0])
# gastroenterology
gas.bib.coup.g.isolates <- length(V(gas.bib.coup.g)[degree(gas.bib.coup.g) == 0])
gas.pub.auth.g.isolates <- length(V(gas.pub.auth.g)[degree(gas.pub.auth.g) == 0])

# delete isolates from graph
# accounting
acc.bib.coup.g <- delete.vertices(acc.bib.coup.g, V(acc.bib.coup.g)[degree(acc.bib.coup.g) == 0])
acc.pub.auth.g <- delete.vertices(acc.pub.auth.g, V(acc.pub.auth.g)[degree(acc.pub.auth.g) == 0])
# astronomy
ast.bib.coup.g <- delete.vertices(ast.bib.coup.g, V(ast.bib.coup.g)[degree(ast.bib.coup.g) == 0])
ast.pub.auth.g <- delete.vertices(ast.pub.auth.g, V(ast.pub.auth.g)[degree(ast.pub.auth.g) == 0])
# gastoenterology
gas.bib.coup.g <- delete.vertices(gas.bib.coup.g, V(gas.bib.coup.g)[degree(gas.bib.coup.g) == 0])
gas.pub.auth.g <- delete.vertices(gas.pub.auth.g, V(gas.pub.auth.g)[degree(gas.gpub.auth.g) == 0])

# cluster analysis
# accounting
acc.bib.coup.c <- cluster_fast_greedy(acc.bib.coup.g)
acc.pub.auth.c <- cluster_fast_greedy(acc.pub.auth.g)
# astronomy
ast.bib.coup.c <- cluster_fast_greedy(ast.bib.coup.g)
ast.pub.auth.c <- cluster_fast_greedy(ast.pub.auth.g)
# gastroenterology
gas.bib.coup.c <- cluster_fast_greedy(gas.bib.coup.g)
gas.pub.auth.c <- cluster_fast_greedy(gas.pub.auth.g)

# network statistics
network_statistics <- function(bib.coup.g,
                               pub.auth.g,
                               bib.coup.g.isolates,
                               pub.auth.g.isolates,
                               bib.coup.c,
                               pub.auth.c) {
  tibble(
  vertices = c(vcount(bib.coup.g) + bib.coup.g.isolates, vcount(pub.auth.g) + pub.auth.g.isolates),
  isolates = c(bib.coup.g.isolates, pub.auth.g.isolates),
  edges = c(ecount(bib.coup.g), ecount(pub.auth.g)),
  density = c(graph.density(bib.coup.g), graph.density(pub.auth.g)),
  transitivity = c(transitivity(bib.coup.g), transitivity(pub.auth.g)),
  assortativity.degree = c(assortativity_degree(bib.coup.g), assortativity_degree(pub.auth.g)),
  # components = c(count_components(bib.coup.g), count_components(pub.auth.g)),
  # biggest.component = c(max(components(bib.coup.g)$csize), max(components(pub.auth.g)$csize)),
  clusters = c(length(bib.coup.c), length(pub.auth.c)),
  biggest.cluster = c(max(sizes(bib.coup.c)), max(sizes(pub.auth.c)))
  )
}

acc.bib.network.statistics <- network_statistics(acc.bib.coup.g,
                                                 acc.pub.auth.g,
                                                 acc.bib.coup.g.isolates,
                                                 acc.pub.auth.g.isolates,
                                                 acc.bib.coup.c,
                                                 acc.pub.auth.c
                                                 )

# latex table
xtable::xtable(t(acc.bib.network.statistics), type = "latex")



# cluster distribution ####

# bibliographic coupling in accounting
acc.bib.coup.plot <- tibble(n = as.numeric(sizes(acc.bib.coup.c))) %>%
  # arrange data in descending order
  arrange(desc(n)) %>%
  # add column with cumulative percentages
  mutate(p = cumsum(sort(sizes(acc.bib.coup.c), decreasing = TRUE)) / sum(sizes(acc.bib.coup.c))) %>%
  # add column with cluster number
  mutate(cluster = 1:length(acc.bib.coup.c)) %>%
  # plot
  ggplot(aes(x = cluster, y = n)) +
  geom_bar(stat = "identity", width = 1, fill = "gray") +
  geom_line(mapping = aes(x = cluster, y = p * max(n)), color = "black") +
  scale_x_continuous(breaks = seq(1, length(acc.bib.coup.c), 1)) +
  scale_y_continuous(
    sec.axis = sec_axis(~ ./max(sizes(acc.bib.coup.c)),
    labels = scales::percent,
    name = "Cumulative percentage")
  ) +
  geom_hline(yintercept = 80/(100/max(sizes(acc.bib.coup.c))), color = "black", linetype = "dashed") +
  labs(title = "Bibliographic Coupling in Accounting", x = "Clusters", y = "Articles") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 8),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black")
  )

# article authorship in accounting
acc.pub.auth.plot <- tibble(n = as.numeric(sizes(acc.pub.auth.c))) %>%
  # arrange data in descending order
  arrange(desc(n)) %>%
  # add column with cumulative percentages
  mutate(p = cumsum(sort(sizes(acc.pub.auth.c), decreasing = TRUE)) / sum(sizes(acc.pub.auth.c))) %>%
  # add column with cluster number
  mutate(cluster = 1:length(acc.pub.auth.c)) %>%
  # plot
  ggplot(aes(x = cluster, y = n)) +
  geom_bar(stat = "identity", width = 1, fill = "gray") +
  geom_line(mapping = aes(x = cluster, y = p * max(n)), color = "black") +
  scale_x_continuous(breaks = seq(50, length(acc.pub.auth.c), 50)) +
  scale_y_continuous(
    sec.axis = sec_axis(~ ./max(sizes(acc.pub.auth.c)),
      labels = scales::percent,
      name = "Cumulative percentage")
  ) +
  geom_hline(yintercept = 80/(100/max(sizes(acc.pub.auth.c))), color = "black", linetype = "dashed") +
  labs(title = "Publication Authorship in Accounting", x = "Clusters", y = "Articles") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 8),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black")
  )

# bibliographic coupling in astronomy
ast.bib.coup.plot <- tibble(n = as.numeric(sizes(ast.bib.coup.c))) %>%
  # arrange data in descending order
  arrange(desc(n)) %>%
  # add column with cumulative percentages
  mutate(p = cumsum(sort(sizes(ast.bib.coup.c), decreasing = TRUE)) / sum(sizes(ast.bib.coup.c))) %>%
  # add column with cluster number
  mutate(cluster = 1:length(ast.bib.coup.c)) %>%
  # plot
  ggplot(aes(x = cluster, y = n)) +
  geom_bar(stat = "identity", width = 1, fill = "gray") +
  geom_line(mapping = aes(x = cluster, y = p * max(n)), color = "black") +
  scale_x_continuous(breaks = seq(1, length(ast.bib.coup.c), 1)) +
  scale_y_continuous(
    sec.axis = sec_axis(~ ./max(sizes(ast.bib.coup.c)),
                        labels = scales::percent,
                        name = "Cumulative percentage")
  ) +
  geom_hline(yintercept = 80/(100/max(sizes(ast.bib.coup.c))), color = "black", linetype = "dashed") +
  labs(title = "Bibliographic Coupling in Astronomy", x = "Clusters", y = "Articles") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 8),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black")
  )

# article authorship in astronomy
ast.pub.auth.plot <- tibble(n = as.numeric(sizes(ast.pub.auth.c))) %>%
  # arrange data in descending order
  arrange(desc(n)) %>%
  # add column with cumulative percentages
  mutate(p = cumsum(sort(sizes(ast.pub.auth.c), decreasing = TRUE)) / sum(sizes(ast.pub.auth.c))) %>%
  # add column with cluster number
  mutate(cluster = 1:length(ast.pub.auth.c)) %>%
  # plot
  ggplot(aes(x = cluster, y = n)) +
  geom_bar(stat = "identity", width = 1, fill = "gray") +
  geom_line(mapping = aes(x = cluster, y = p * max(n)), color = "black") +
  scale_x_continuous(breaks = seq(20, length(ast.pub.auth.c), 20)) +
  scale_y_continuous(
    sec.axis = sec_axis(~ ./max(sizes(ast.pub.auth.c)),
                        labels = scales::percent,
                        name = "Cumulative percentage")
  ) +
  geom_hline(yintercept = 80/(100/max(sizes(ast.pub.auth.c))), color = "black", linetype = "dashed") +
  labs(title = "Publication Authorship in Astronomy", x = "Clusters", y = "Articles") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 8),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black")
  )

# bibliographic coupling in gastroenterology
gas.bib.coup.plot <- tibble(n = as.numeric(sizes(gas.bib.coup.c))) %>%
  # arrange data in descending order
  arrange(desc(n)) %>%
  # add column with cumulative percentages
  mutate(p = cumsum(sort(sizes(gas.bib.coup.c), decreasing = TRUE)) / sum(sizes(gas.bib.coup.c))) %>%
  # add column with cluster number
  mutate(cluster = 1:length(gas.bib.coup.c)) %>%
  # plot
  ggplot(aes(x = cluster, y = n)) +
  geom_bar(stat = "identity", width = 1, fill = "gray") +
  geom_line(mapping = aes(x = cluster, y = p * max(n)), color = "black") +
  scale_x_continuous(breaks = seq(5, length(gas.bib.coup.c), 5)) +
  scale_y_continuous(
    sec.axis = sec_axis(~ ./max(sizes(gas.bib.coup.c)),
                        labels = scales::percent,
                        name = "Cumulative percentage")
  ) +
  geom_hline(yintercept = 80/(100/max(sizes(gas.bib.coup.c))), color = "black", linetype = "dashed") +
  labs(title = "Bibliographic Coupling in Gastroenterology", x = "Clusters", y = "Articles") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 8),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black")
  )

# article authorship in gastroenterology
gas.pub.auth.plot <- tibble(n = as.numeric(sizes(gas.pub.auth.c))) %>%
  # arrange data in descending order
  arrange(desc(n)) %>%
  # add column with cumulative percentages
  mutate(p = cumsum(sort(sizes(gas.pub.auth.c), decreasing = TRUE)) / sum(sizes(gas.pub.auth.c))) %>%
  # add column with cluster number
  mutate(cluster = 1:length(gas.pub.auth.c)) %>%
  # plot
  ggplot(aes(x = cluster, y = n)) +
  geom_bar(stat = "identity", width = 1, fill = "gray") +
  geom_line(mapping = aes(x = cluster, y = p * max(n)), color = "black") +
  scale_x_continuous(breaks = seq(10, length(gas.pub.auth.c), 10)) +
  scale_y_continuous(
    sec.axis = sec_axis(~ ./max(sizes(gas.pub.auth.c)),
                        labels = scales::percent,
                        name = "Cumulative percentage")
  ) +
  geom_hline(yintercept = 80/(100/max(sizes(gas.pub.auth.c))), color = "black", linetype = "dashed") +
  labs(title = "Publication Authorship in Gastroenterology", x = "Clusters", y = "Articles") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 8),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black")
  )

# publication-ready plot
ggarrange(
  acc.bib.coup.plot, acc.pub.auth.plot, ast.bib.coup.plot, ast.pub.auth.plot, gas.bib.coup.plot, gas.pub.auth.plot,
  # labels = c("(A)", "(B)", "(C)", "(D)", "(E)", "(F)"),
  # label.y = .25,
  # font.label = list(size = 12, face = "plain", color = "black"),
  ncol = 2, nrow = 3
  )

ggsave("cluster_distributions.pdf", width = 10.5, height = 7.875, units = "in")

# (alternative) plots with tikz
library(tikzDevice)
options(tikzLualatex = '/path/to/lualatex')
tikz(file = "cluster_distributions.tex", width = 5.25, height = 5.25, sanitize = TRUE)
par(
  font.main = 1, # plain font
  cex.main = 1.25, # font size (scale)
  mar = c(0, 1, 2, 1) # margins = c(bottom, left, top, right)
)
layout(matrix(c(1, 2, 3, 4, 5, 6), 3, 2, byrow = TRUE), respect = TRUE)

acc.bib.coup.plot
acc.pub.auth.plot
ast.bib.coup.plot
ast.pub.auth.plot
gas.bib.coup.plot
gas.pub.auth.plot

# default
dev.off()



# inter-cluster cosine similarity ####

# accounting: bibliographic coupling vs cosine similarity of article abstracts
bib.coup.m.temp <- acc.bib.coup.m
bib.coup.m.temp[bib.coup.m.temp > 0] <- 1 # dummy for true value > 0 (existing edge), 0 otherwise
cos.sim.m.temp <- acc.cos.sim.m
cos.sim.m.temp[cos.sim.m.temp == 0] <- 100 # dummy for true value == 0 (existing edge at minimum cosine similarity)
acc.bib.coup.cos.sim.m <- bib.coup.m.temp * cos.sim.m.temp
rm(bib.coup.m.temp, cos.sim.m.temp)

# astronomy: publication authorship vs cosine similarity of article abstracts
pub.auth.m.temp <- acc.pub.auth.m
pub.auth.m.temp[pub.auth.m.temp > 0] <- 1 # dummy for true value > 0 (existing edge), 0 otherwise
cos.sim.m.temp <- acc.cos.sim.m
cos.sim.m.temp[cos.sim.m.temp == 0] <- 100 # dummy for true value == 0 (existing edge at minimum cosine similarity)
acc.pub.auth.cos.sim.m <- pub.auth.m.temp * cos.sim.m.temp
rm(pub.auth.m.temp, cos.sim.m.temp)

# astronomy: bibliographic coupling vs cosine similarity of article abstracts
bib.coup.m.temp <- ast.bib.coup.m
bib.coup.m.temp[bib.coup.m.temp > 0] <- 1 # dummy for true value > 0 (existing edge), 0 otherwise
cos.sim.m.temp <- ast.cos.sim.m
cos.sim.m.temp[cos.sim.m.temp == 0] <- 100 # dummy for true value == 0 (existing edge at minimum cosine similarity)
ast.bib.coup.cos.sim.m <- bib.coup.m.temp * cos.sim.m.temp
rm(bib.coup.m.temp, cos.sim.m.temp)

# astronomy: publication authorship vs cosine similarity of article abstracts
pub.auth.m.temp <- ast.pub.auth.m
pub.auth.m.temp[pub.auth.m.temp > 0] <- 1 # dummy for true value > 0 (existing edge), 0 otherwise
cos.sim.m.temp <- ast.cos.sim.m
cos.sim.m.temp[cos.sim.m.temp == 0] <- 100 # dummy for true value == 0 (existing edge at minimum cosine similarity)
ast.pub.auth.cos.sim.m <- pub.auth.m.temp * cos.sim.m.temp
rm(pub.auth.m.temp, cos.sim.m.temp)

# gastroenterology: bibliographic coupling vs cosine similarity of article abstracts
bib.coup.m.temp <- gas.bib.coup.m
bib.coup.m.temp[bib.coup.m.temp > 0] <- 1 # dummy for true value > 0 (existing edge), 0 otherwise
cos.sim.m.temp <- gas.cos.sim.m
cos.sim.m.temp[cos.sim.m.temp == 0] <- 100 # dummy for true value == 0 (existing edge at minimum cosine similarity)
gas.bib.coup.cos.sim.m <- bib.coup.m.temp * cos.sim.m.temp
rm(bib.coup.m.temp, cos.sim.m.temp)

# gastroenterology: publication authorship vs cosine similarity of article abstracts
pub.auth.m.temp <- gas.pub.auth.m
pub.auth.m.temp[pub.auth.m.temp > 0] <- 1 # dummy for true value > 0 (existing edge), 0 otherwise
cos.sim.m.temp <- gas.cos.sim.m
cos.sim.m.temp[cos.sim.m.temp == 0] <- 100 # dummy for true value == 0 (existing edge at minimum cosine similarity)
gas.pub.auth.cos.sim.m <- pub.auth.m.temp * cos.sim.m.temp
rm(pub.auth.m.temp, cos.sim.m.temp)

# accounting: inter-cluster cosine similarity
acc.bib.coup.cs <- lapply(groups(acc.bib.coup.c), function(x) {
  j <- acc.bib.coup.cos.sim.m[unlist(x), unlist(x)] # filter cluster (groups) from cosine similarity matrix
  j <- j[upper.tri(j)] # upper triangle
  j[j == 0] <- NA # remove all 0s (zeros)
  j[j == 100] <- 0 # replace dummy (100) for true value == 0 (existing edge at minimum cosine similarity)
  return(j)
})

acc.pub.auth.cs <- lapply(groups(acc.pub.auth.c), function(x) {
  j <- acc.pub.auth.cos.sim.m[unlist(x), unlist(x)] # filter cluster (groups) from cosine similarity matrix
  j <- j[upper.tri(j)] # upper triangle
  j[j == 0] <- NA # remove all 0s (zeros)
  j[j == 100] <- 0 # replace dummy (100) for true value == 0 (existing edge at minimum cosine similarity)
  return(j)
})

# astronomy: inter-cluster cosine similarity
ast.bib.coup.cs <- lapply(groups(ast.bib.coup.c), function(x) {
  j <- ast.bib.coup.cos.sim.m[unlist(x), unlist(x)] # filter cluster (groups) from cosine similarity matrix
  j <- j[upper.tri(j)] # upper triangle
  j[j == 0] <- NA # remove all 0s (zeros)
  j[j == 100] <- 0 # replace dummy (100) for true value == 0 (existing edge at minimum cosine similarity)
  return(j)
})

ast.pub.auth.cs <- lapply(groups(ast.pub.auth.c), function(x) {
  j <- ast.pub.auth.cos.sim.m[unlist(x), unlist(x)] # filter cluster (groups) from cosine similarity matrix
  j <- j[upper.tri(j)] # upper triangle
  j[j == 0] <- NA # remove all 0s (zeros)
  j[j == 100] <- 0 # replace dummy (100) for true value == 0 (existing edge at minimum cosine similarity)
  return(j)
})

# gastroenterology: inter-cluster cosine similarity
gas.bib.coup.cs <- lapply(groups(gas.bib.coup.c), function(x) {
  j <- gas.bib.coup.cos.sim.m[unlist(x), unlist(x)] # filter cluster (groups) from cosine similarity matrix
  j <- j[upper.tri(j)] # upper triangle
  j[j == 0] <- NA # remove all 0s (zeros)
  j[j == 100] <- 0 # replace dummy (100) for true value == 0 (existing edge at minimum cosine similarity)
  return(j)
})

gas.pub.auth.cs <- lapply(groups(gas.pub.auth.c), function(x) {
  j <- gas.pub.auth.cos.sim.m[unlist(x), unlist(x)] # filter cluster (groups) from cosine similarity matrix
  j <- j[upper.tri(j)] # upper triangle
  j[j == 0] <- NA # remove all 0s (zeros)
  j[j == 100] <- 0 # replace dummy (100) for true value == 0 (existing edge at minimum cosine similarity)
  return(j)
})

# boxplots (cluster mean)
boxplot <- rbind(
  map(gas.pub.auth.cs, mean, na.rm = TRUE) %>%
    unlist() %>%
    as_tibble() %>%
    mutate(id = factor("Publication Authorship in Gastroenterology"))
  ,
  map(gas.bib.coup.cs, mean, na.rm = TRUE) %>%
    unlist() %>%
    as_tibble() %>%
    mutate(id = factor("Bibliographic Coupling in Gastroenterology"))
  ,
  map(ast.pub.auth.cs, mean, na.rm = TRUE) %>%
    unlist() %>%
    as_tibble() %>%
    mutate(id = factor("Publication Authorship in Astronomy"))
  ,
  map(ast.bib.coup.cs, mean, na.rm = TRUE) %>%
    unlist() %>%
    as_tibble() %>%
    mutate(id = factor("Bibliographic Coupling in Astronomy"))
  ,
  map(acc.pub.auth.cs, mean, na.rm = TRUE) %>%
    unlist() %>%
    as_tibble() %>%
    mutate(id = factor("Publication Authorship in Accounting"))
  ,
  map(acc.bib.coup.cs, mean, na.rm = TRUE) %>%
    unlist() %>%
    as_tibble() %>%
    mutate(id = factor("Bibliographic Coupling in Accounting"))
) %>%
  ggplot(aes(id, value)) +
  geom_boxplot() +
  coord_flip() +
  labs(title = "", x = "", y = "Cosine Similarity") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 8),
    axis.title = element_text(size = 8),
    axis.text.x = element_text(size = 8, angle = 0),
    axis.text.y = element_text(size = 8, angle = 0),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks.y = element_blank()
  ) +
  stat_compare_means(
    method = "wilcox.test",
    method.args = list(alternative = "greater"),
    comparisons = list(
      c("Bibliographic Coupling in Gastroenterology", "Publication Authorship in Gastroenterology"),
      c("Bibliographic Coupling in Astronomy", "Publication Authorship in Astronomy"),
      c("Bibliographic Coupling in Accounting", "Publication Authorship in Accounting")
    )
  )

ggsave("cluster_boxplots.pdf", width = 10.5, height = 5.25, units = "in")

# (alternative) plots with tikz
library(tikzDevice)
options(tikzLualatex = '/path/to/lualatex')
tikz(file = "cluster_boxplots.tex", width = 5.25, height = 2.625, sanitize = TRUE)
par(
  font.main = 1, # plain font
  cex.main = 1.25, # font size (scale)
  mar = c(0, 1, 2, 1) # margins = c(bottom, left, top, right)
)
layout(matrix(c(1, 2, 3, 4, 5, 6), 3, 2, byrow = TRUE), respect = TRUE)

boxplot

# default
dev.off()



# mann-whitney u test
ast.wilcox <- wilcox.test(unlist(map(ast.bib.coup.cs, mean, na.rm = TRUE)), unlist(map(ast.pub.auth.cs, mean, na.rm = TRUE)), alternative = "greater", exact = FALSE, correct = FALSE)
# ast.wilcox <- wilcox.test(sapply(ast.cs.bc, mean), sapply(ast.cs.ac, mean), alternative = "greater", exact = FALSE, correct = FALSE)
# gas.wilcox <- wilcox.test(sapply(gas.cs.bc, mean), sapply(gas.cs.ac, mean), alternative = "greater", exact = FALSE, correct = FALSE)