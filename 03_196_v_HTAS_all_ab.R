library(agricolae)
library(ggplot2)
setwd('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/abudance_comparison/')

data = read.csv2('all_comparisons.csv', sep = ';')
head(data)

gene.colours <- c("#fcc5c0", "#7a0177", "#2b8cbe", "#ccebc5")

# all genes in the same plot 
p = ggplot(data, aes(x=MAUI, y=X196)) + 
  geom_point(aes(color = X)) +
  scale_fill_manual(values = gene.colours) +
  scale_colour_manual(values = gene.colours) +
  scale_y_continuous(expand=c(0.01,0)) +
  scale_x_continuous(expand=c(0.01,0)) +
  theme_bw() +
  labs(x="MAUI-seq", y = "196 genomes") 
p
ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/abudance_comparison/all_genes.pdf', plot = p, width = 10, height = 7, unit = 'cm')

# separate plots 
# rpoB 
rpoB <- subset(data, X == "rpoB")
correlation(rpoB$X196, rpoB$MAUI) 

p <- ggplot(rpoB, aes(x = MAUI, y = X196)) + 
  geom_point(colour = "#2b8cbe") +
  geom_smooth(colour = "#7bccc4", fill = "#ccebc5", method = 'lm') +
  scale_x_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
  labs(x="MAUI-seq", y = "196 genomes") +
  theme_bw() 

p
ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/abudance_comparison/rpoB.pdf', plot = p, width = 7, height = 7, unit = 'cm')

# recA
recA <- subset(data, X == "recA")
correlation(recA$X196, recA$MAUI) 

p <- ggplot(recA, aes(x = MAUI, y = X196)) + 
  geom_point(colour = "#2b8cbe") +
  geom_smooth(colour = "#7bccc4", fill = "#ccebc5", method = 'lm') +
  scale_x_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
  labs(x="MAUI-seq", y = "196 genomes") +
  theme_bw() 

p
ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/abudance_comparison/recA.pdf', plot = p, width = 7, height = 7, unit = 'cm')

# nodA
nodA <- subset(data, X == "nodA")
correlation(nodA$X196, nodA$MAUI) 

p <- ggplot(nodA, aes(x = MAUI, y = X196)) + 
  geom_point(colour = "#7a0177") +
  geom_smooth(colour = "#f768a1", fill = "#fcc5c0", method = 'lm') +
  scale_x_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
  labs(x="MAUI-seq", y = "196 genomes") +
  theme_bw() 

p
ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/abudance_comparison/nodA.pdf', plot = p, width = 7, height = 7, unit = 'cm')

# nodD
nodD <- subset(data, X == "nodD")
correlation(nodD$X196, nodD$MAUI) 
p <- ggplot(nodD, aes(x = MAUI, y = X196)) + 
  geom_point(colour = "#7a0177") +
  geom_smooth(colour = "#f768a1", fill = "#fcc5c0", method = 'lm') +
  scale_x_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
  labs(x="MAUI-seq", y = "196 genomes") +
  theme_bw() 

p
ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/abudance_comparison/nodD.pdf', plot = p, width = 7, height = 7, unit = 'cm')


