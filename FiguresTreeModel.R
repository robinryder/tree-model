# Code for figures from chapter on The Tree Model (Pellard, Ryder & Jacques 2023)
# RJR March 2023

library(phytools)
library(ape)
library(phangorn)
library(tidyverse)
library(ggthemes)
library(ggridges)
library(HDInterval)
library(knitr)

wd <- 5
base_font <- "URW Palladio L"
base_font2 <- "URWPalladio"
base_font_size <- 10
theme_set(
  theme_minimal(base_family = base_font, base_size = base_font_size) +
    theme(
      text = element_text(family = base_font, size = base_font_size, color = "black"),
      axis.text = element_text(size = base_font_size, color = "black"),
      axis.title = element_text(size = base_font_size, color = "black"),
      legend.text = element_text(size = base_font_size, color = "black"),
      plot.margin = margin(0, 0, 0, 3)
    )
)
update_geom_defaults("text", list(family = base_font, size = base_font_size / .pt))

# Import data -------------------------------------------------------------

tt <- read.nexus("SinoTibetanSubset.nex")


# Consensus tree ----------------------------------------------------------

# Consensus topology
cons <- consensus(tt, p = .5, rooted = T)
cons$node.label <- round(cons$node.label * 100, 0)

# Consensus topology + branch lengths. This might take a couple of minutes
ct <- consensus.edges(tt, consensus.tree = cons, rooted = TRUE)
ct$root.edge <- .15

pdf("fig_consensus.pdf", pointsize=10, width = 6.3, height = 4, family = base_font2)
plot(ct, show.node.label=FALSE, root.edge = TRUE, no.margin = TRUE, font = 1, label.offset = .05)
nodelabels(c(NA, ct$node.label[-1]), frame="none", adj = c(1.15,1.25), cex = .9)
dev.off()
embedFonts("fig_consensus.pdf")
plot_crop("fig_consensus.pdf")


# MCC tree --------------------------------------------------------------------------------------------------------

mcct <- ladderize(maxCladeCred(tt))
mcct$node.label <- round(mcct$node.label * 100, 0)
mcct$root.edge <- .15
pdf("fig_mcc.pdf", pointsize=10, width = 6.3, height = 4, family = base_font2)
plot(mcct, show.node.label = FALSE, root.edge = TRUE, no.margin = TRUE, font = 1, label.offset = .05)
nodelabels(c(NA, mcct$node.label[2], NA, mcct$node.label[c(-1:-3)]), frame = "none", adj = c(1.15, 1.25), cex = .9)
nodelabels(c(NA, NA, mcct$node.label[3], rep(NA, length(mcct$node.label) - 3)), frame = "none", adj = c(1.15, -.25), cex = .9)
dev.off()
embedFonts("fig_mcc.pdf")
plot_crop("fig_mcc.pdf")

# Densitree ---------------------------------------------------------------

tt_scaled <- lapply(tt, function(x) {x$edge.length <- x$edge.length * 1000; return(x)})
class(tt_scaled) <- "multiPhylo"
maxBT <- max(phangorn:::getAges(tt_scaled))
label <- rev(pretty(c(maxBT, 0)))
maxBT <- max(label)
pdf("fig_densitree.pdf", pointsize = 10, width = 5.025, height = 6, family = "URWPalladio")
par(mar = c(2, 0, 0, .6), xpd = TRUE)
densiTree(tt_scaled, consensus = cons, alpha = .005, font = 1, label.offset = .01, cex = 1, scale.bar = FALSE)
axis(side = 1, at = seq(0, 1.0, length.out = length(label)), labels = label, line = -1.5, cex = .9)
title(xlab = "years BP", line = 1)
dev.off()
embedFonts("fig_densitree.pdf")
plot_crop("fig_densitree.pdf")


# Root age plot -----------------------------------------------------------

# Posterior probability of a specific subtree
subtree_support <- function(tips, trees = tt) {
  sapply(trees, is.monophyletic, tips = tips)
}
mean(subtree_support(c("Chepang", "Hayu", "Thulung", "Bokar_Tani", "Yidu", "Tshangla")))
# The 6 leaves listed above form a subtree with posterior probability 0.41

# Creating the plot for the root age depending on the outgroup
# Lists of outgroups
find.outgroup <- function(tree) {
  nl <- length(tree$tip.label)
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

outgroup <- sapply(tt, find.outgroup)
outgroup_tb <- table(outgroup) |>
  prop.table() |>
  as_tibble() |>
  arrange(-n) |>
  rename(p = n)

ages_outgroup <- tt |>
  seq_along() |>
  map_df(~
           tibble(
             age = max(node.depth.edgelength(tt[[.x]])),
             outgroup = find.outgroup(tt[[.x]])
           )) |>
  left_join(outgroup_tb) |>
  mutate(outgroup = case_when(
    outgroup == "Beijing_Chinese Guangzhou_Chinese Jieyang_Chinese Xingning_Chinese" ~ "Chinese",
    outgroup == "Beijing_Chinese Guangzhou_Chinese Jieyang_Chinese Jingpho Rabha Xingning_Chinese" ~ "Chinese-Sal",
    outgroup == "Bokar_Tani Yidu" ~ "Tani-Yidu",
    .default = "x"
  )) |>
  mutate(label = paste0(outgroup, "\n(", round(p, 2) * 100, "%)")) %>%
  bind_rows(mutate(., outgroup = "any", label = "any")) |>
  filter(!str_detect(outgroup, "x"))

ages_outgroup %>%
  ggplot(aes(x = age * 1000, y = fct_rev(label), fill = fct_rev(label), height = after_stat(density))) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2, color = "white") +
  geom_density_ridges(fill = NA, color = "gray40") +
  scale_fill_manual(values = c(rep(few_pal("Light")(2)[1], 3), "grey"), guide = "none") +
  scale_x_reverse(limits = c(15000, 0)) +
  scale_y_discrete(expand = expansion(add = c(0.25, 1.4))) +
  xlab("root age (years BP)") +
  ylab("first branch") +
  theme(aspect.ratio = 0.618)
ggsave("fig_ageoutgroup.pdf", width = wd, height = wd * .8, units = "in", device = cairo_pdf)
plot_crop("fig_ageoutgroup.pdf")


###### MRCA ages

getMRCA_age <- function(tree, tips) {
  tips <- if (is.character(tips)) which(tree$tip.label %in% tips) else tips
  mrca <- ifelse(length(tips) > 1, getMRCA(tree, tips), tips)
  root_age <- max(node.depth.edgelength(tree))
  root_age - node.depth.edgelength(tree)[mrca]
}

tips <- tt[[1]]$tip.label |>
  str_subset("Jingpho|Rabha|Chinese")
sinitic_sal <- st_tree |>
  seq_along() |>
  map_df(~
           tibble(
             age = getMRCA_age(tt[[.x]], tips),
             monophyletic = ifelse(is.monophyletic(tt[[.x]], tips), "monophyletic", "paraphyletic")
           )) %>%
  bind_rows(mutate(., monophyletic = "either"))

sinitic_sal |>
  ggplot(aes(x = age * 1000, y = fct_rev(monophyletic), fill = monophyletic, height = after_stat(density))) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2, color = "white") +
  geom_density_ridges(fill = NA, color = "gray40") +
  scale_fill_manual(values = c("grey", rep(few_pal("Light")(2)[1], 2)), guide = "none") +
  scale_x_reverse(limits = c(15000, 0)) +
  scale_y_discrete(expand = expansion(add = c(0.25, 1.25))) +
  xlab("root age (years BP)") +
  ylab("status of Chinese-Sal") +
  theme(aspect.ratio = 0.618)
ggsave("fig_agemono.pdf", width = wd, height = wd * .8, units = "in", device = cairo_pdf)
plot_crop("fig_agemono.pdf")


# Cross-validation of dates -----------------------------------------------
library(HDInterval)

# True age intervals
true.tibetan = c(1000, 1200)
true.burmish = c(700, 900)
true.commonchinese = c(2000, 2200)
true.sinitic = c(2300, 2800)

reconstructed.age = function(file, lang){
  t = read.nexus(file)
  
  #remove 20% for burn-in
  nt = length(t)
  t = t[(nt * .2) : nt]
  
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
  ylab("years BP") +
  xlab("node") +
  scale_color_few() +
  # ylim(c(0,3000)) +
  scale_y_reverse(limits = c(3000, 0)) +
  scale_x_discrete(expand = expansion(add = c(0.25, 1.25))) +
  theme_minimal() +
  theme(legend.position = "top", axis.text.y.left = element_text(size = 10), axis.title = element_text(size = 9))

fig_xages <- df |> 
  mutate(clade = str_replace(clade, " ", "\n")) |> 
  ggplot(aes(y=clade, x=(lower+upper)/2, xmin=lower, xmax=upper)) +
  geom_linerange(aes(color=Type), position=position_dodge(width=c(0.2)), linewidth=2) +
  xlab("years BP") +
  ylab(NULL) +
  scale_color_few(name = NULL) +
  scale_x_reverse(limits = c(3000, 0)) +
  scale_y_discrete(limits=rev, expand = expansion(add = c(0.25, .25))) +
  theme_minimal() +
  theme(legend.position = "top", axis.text.y.left = element_text(size = 10), axis.title = element_text(size = 9))

pdf("fig_xages.pdf", pointsize=10, width = 5, height = 5/1.6, family = "URWPalladio")
fig_xages
dev.off()
embedFonts("fig_xages.pdf")
plot_crop("fig_xages.pdf")
