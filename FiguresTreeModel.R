# Code for figures from chapter on The Tree Model (Jacques, Pellard & Ryder 2023)
# RJR March 2023

library(phytools)
library(ape)
library(phangorn)
library(ggplot2)
library(cowplot)
library(ggthemes)
library(knitr)

# Import data -------------------------------------------------------------

tt <- read.nexus("SinoTibetanSubset.nex")


# Consensus tree ----------------------------------------------------------

# Consensus topology
cons <- consensus(tt, p = .5, rooted = T)
cons$node.label <- round(cons$node.label * 100, 0)

# Consensus topology + branch lengths. This might take a couple of minutes
ct <- consensus.edges(tt, consensus.tree = cons, rooted = T)
ct$root.edge <- .15

pdf("fig_consensus.pdf", pointsize=10, width = 6, height = 4, family = "URWPalladio")
plot(ct, show.node.label=FALSE, root.edge = TRUE, no.margin = TRUE, font = 1, label.offset = .05)
nodelabels(ct$node.label, frame="none", adj = c(-0.2,.5), cex = .85)
dev.off()
embedFonts("fig_consensus.pdf")
plot_crop("fig_consensus.pdf")


# Densitree ---------------------------------------------------------------

pdf("fig_densitree.pdf", pointsize=10, width = 6.5, height = 7, family = "URWPalladio")
densiTree(tt, consensus = cons, alpha=.005, font = 1, label.offset = .01, cex=.9, scale.bar = TRUE)
title(xlab="time (millenia BP)")
dev.off()
embedFonts("fig_densitree.pdf")
knitr::plot_crop("fig_densitree.pdf")


# Root age plot -----------------------------------------------------------

# Posterior probability of a specific subtree
subtree_support <- function(tips, trees = tt) {
  sapply(trees, is.monophyletic, tips = tips)
}
mean(subtree_support(c("Chepang", "Hayu", "Thulung", "Bokar_Tani", "Yidu", "Tshangla")))
# The 6 leaves listed above form a subtree with posterior probability 0.41

# Creating the plot for the root age depending on the outgroup
# Lists of outgroups
find.outgroup <- function(tree, nl = length(tree$tip.label)) {
  r <- getRoot(tree)
  cr <- tree$edge[tree$edge[, 1] == r, 2]
  s1 <- getDescendants(tree, cr[1])
  s2 <- getDescendants(tree, cr[2])
  s1 <- s1[s1 <= nl]
  s2 <- s2[s2 <= nl]
  # Ensure that s1 is the smallest
  if (length(s1) > length(s2)) s1 <- s2
  return(paste(sort(tree$tip.label[s1]), collapse = " "))
}

outgroup <- sapply(tt, find.outgroup, nl = 22)
head(sort(table(outgroup), dec = T), 3)

outgroup[outgroup == "Beijing_Chinese Guangzhou_Chinese Jieyang_Chinese Xingning_Chinese"] <- "Chinese"
outgroup[outgroup == "Beijing_Chinese Guangzhou_Chinese Jieyang_Chinese Jingpho Rabha Xingning_Chinese"] <- "Chinese-Sal"
outgroup[outgroup == "Bokar_Tani Yidu"] <- "Tani-Yidu"

root.age <- function(tree) {
  node.depth.edgelength(tree)[1]
}
ra <- sapply(tt, root.age)
aged <- data.frame(age = ra * 1000, outgroup = outgroup)
keep <- outgroup %in% c("Chinese", "Tani-Yidu", "Chinese-Sal")

aged2 <- rbind(aged[keep, ], data.frame(age = ra * 1000, outgroup = "all"))

p.age.group <- ggplot(aged[keep, ], aes(x = age)) +
  xlim(0, 15000) +
  geom_density(aes(group = outgroup, colour = outgroup, fill = outgroup), alpha = .3) + 
  scale_fill_few() + scale_color_few() + 
  theme_minimal() +
  theme(legend.position = c(1,.5), legend.justification = c(1,.5))

col <- "purple"
p.age <- ggplot(aged, aes(x = age)) +
  xlim(0, 15000) +
  geom_density(fill = col, colour = col, alpha = .3) + 
  theme_minimal()

plot_grid(p.age, p.age.group, align = "v", axis = "rblt", ncol = 1)

options(scipen = 999)
pdf("fig_ageboth.pdf", pointsize=10, width = 5, height = 5/1.6*2, family = "URWPalladio")
plot_grid(p.age, p.age.group, align="v", axis="rblt", ncol=1)
dev.off()
embedFonts("fig_ageboth.pdf")
knitr::plot_crop("fig_ageboth.pdf")

library(ggridges)
library(tidyverse)
fig_ageoutgroup <- aged2 %>% 
  as_tibble() %>% 
  mutate(outgroup = factor(outgroup, levels = c("Tani-Yidu", "Chinese-Sal", "Chinese", "all"))) %>% 
  ggplot(aes(x = age, y = outgroup, fill = outgroup, height = stat(density))) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2, color = "white") +
  geom_density_ridges(fill = NA, color = "gray40") +
  scale_fill_manual(values = c(rep(few_pal("Light")(2)[1], 3), "grey"), guide = "none") +
  xlim(0, 15000) +
  theme_minimal() +
  xlab("time (years BP)") +
  theme(axis.text.y.left = element_text(size = 11), axis.title = element_text(size = 9))

pdf("fig_ageoutgroup.pdf", pointsize=10, width = 5, height = 5/1.6, family = "URWPalladio")
fig_ageoutgroup
dev.off()
embedFonts("fig_ageoutgroup.pdf")
knitr::plot_crop("fig_ageoutgroup.pdf")

###### MRCA ages

mrca.age <- function(tree, tips) {
  max(dist.nodes(tree)[tips, tips]) / 2
}

mrca.age.names <- function(tree, tipnames) {
  max(cophenetic(tree)[tipnames, tipnames]) / 2
}

mrca.summary <- function(d, l) {
  nt <- length(d)
  res <- rep(NA, nt)
  for (i in 1:nt) {
    if (is.numeric(l)) {
      res[i] <- mrca.age(d[[i]], l)
    } else {
      res[i] <- mrca.age.names(d[[i]], l)
    }
  }

  return(res)
}

# Plot of MRCA age for the following languages: Jingpho, Rabha, 4 Chinese dialects
ll <- tt[[1]]$tip.label
l <- c(7, 12, 13:16) # indices of the languages of interest in vector ll
ages <- mrca.summary(tt, l)
lab <- sapply(tt, is.monophyletic, tips = l)
df <- data.frame(age = ages, monophyletic = lab)

p.age.topo.group <- ggplot(df, aes(x = age)) +
  geom_density(aes(group = monophyletic, colour = monophyletic, fill = monophyletic), alpha = .3)

p.age.topo.all <- ggplot(df, aes(x = age)) +
  geom_density(fill = col, colour = col, alpha = .3)

title <- ggdraw() + draw_label("Age of MRCA", fontface = "bold")

plot_grid(title, p.age.topo.all, p.age.topo.group, align = "v", axis = "rblt", ncol = 1, rel_heights = c(0.2, 1, 1))

library(tidyverse)
fig_agemono <- df %>% 
  mutate(group = ifelse(monophyletic == TRUE, "monophyletic", "paraphyletic")) %>% 
  bind_rows(mutate(df, group = "all")) %>% 
  ggplot(aes(x = age * 1000, y = fct_rev(group), fill = group, height = stat(density))) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2, color = "white") +
  geom_density_ridges(fill = NA, color = "gray40") +
  scale_fill_manual(values = c("grey", rep(few_pal("Light")(2)[1], 2)), guide = "none") +
  xlim(0, 15000) +
  theme_minimal() +
  xlab("time (years BP)") +
  ylab(NULL) +
  theme(axis.text.y.left = element_text(size = 11), axis.title = element_text(size = 9))

# df %>% 
#   mutate(group = ifelse(monophyletic == TRUE, "monophyletic", "paraphyletic")) %>% 
#   bind_rows(mutate(df, group = "all")) %>%
#   ggplot(aes(x = age*1000, fill = group, group = group)) +
#   geom_density(alpha = .5) +
#   scale_fill_few(guide = "none") +
#   xlim(0, 15000) +
#   theme_minimal() +
#   xlab("time (years BP)")

pdf("fig_agemono.pdf", pointsize=10, width = 5, height = 5/1.6, family = "URWPalladio")
fig_agemono
dev.off()
embedFonts("fig_agemono.pdf")
knitr::plot_crop("fig_agemono.pdf")

# Cross-validation of dates -----------------------------------------------
library(HDInterval)

# True age intervals
true.tibetan = c(1000, 1200)
true.burmish = c(700, 900)
true.commonchinese = c(2000, 2200)
true.sinitic = c(2400, 2600)

reconstructed.age = function(file, lang){
  t = read.nexus(file)
  
  #remove 20% for burn-in
  nt = length(t)
  t = t[(nt/2) : nt]
  
  age = rep(NA, length(t))
  for(i in 1:length(t)){
    tt = t[[i]]
    
    if(length(lang) == 1){
      idx = which(tt$tip.label == lang)
    } else {
      idx = getMRCA(tt, lang)
    }
    
    depths = node.depth.edgelength(tt)
    age[i] = max(depths) - depths[idx]
  }
  return(hdi(age))
}


recage.burmish = reconstructed.age("xval/Burmese.nex", 
                                   c("BurmishOldBurmese", "BurmishRangoon"))

recage.commonchinese = reconstructed.age("xval/CommonChinese.nex", 
                                         c("SiniticBeijing", "SiniticGuangzhou", 
                                           "SiniticXingning", "SiniticLonggang",
                                           "SiniticJieyang", "SiniticChaozhou" ))
recage.sinitic = reconstructed.age("xval/Chinese.nex", 
                                   c("SiniticOldChinese" ))
recage.tibetan = reconstructed.age("xval/Tibetan.nex", 
                                   c("TibetanOldTibetan"))


df = as.data.frame(rbind(true.burmish, true.commonchinese, true.sinitic, true.tibetan,
                         recage.burmish, recage.commonchinese, recage.sinitic, recage.tibetan))
df$Type = rep(c("known", "inferred"), each=4)
df$clade = rep(c("Burmish", "Common Chinese", "Sinitic", "Tibetan"), 2)

fig_xages <- ggplot(df, aes(x=clade, y=(lower+upper)/2, ymin=lower, ymax=upper)) +
  geom_linerange(aes(color=Type), position=position_dodge(width=c(0.2)), linewidth=2) +
  ylab("time (years BP)") +
  xlab("node") +
  scale_color_few() +
  ylim(c(0,3000)) +
  theme_minimal() +
  theme(legend.position = "top", axis.text.y.left = element_text(size = 10), axis.title = element_text(size = 9))

pdf("fig_xages.pdf", pointsize=10, width = 5, height = 5/1.6, family = "URWPalladio")
fig_xages
dev.off()
embedFonts("fig_xages.pdf")
plot_crop("fig_xages.pdf")
