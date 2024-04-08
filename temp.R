# Settings ----------------------------------------------------------------

library(ape)
library(phangorn)
library(phytools)
library(tidyverse)
library(ggtree)
library(ggthemes)
library(knitr)
library(ggridges)
library(HDInterval)

wd <- 5
base_font <- "URW Palladio L"
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



st_tree <- read.nexus("SinoTibetanSubset.nex")


# Consensus tree ----------------------------------------------------------

st_tree_cs <- consensus(st_tree, p = .5, rooted = TRUE)
st_tree_cs <- consensus.edges(st_tree, consensus.tree = st_tree_cs, rooted = TRUE, method = "mean.edge", if.absent = "ignore")
st_tree_cs$edge.length <- st_tree_cs$edge.length[1:dim(st_tree_cs$edge)[1]]

st_tree_cs |>
  fortify() |>
  mutate(label = ifelse(branch.length == 0 & label == 1, "", label)) |>
  mutate(label = str_replace_all(label, "_", " ")) |>
  ggtree(ladderize = FALSE) +
  geom_tiplab(family = base_font, size = base_font_size / .pt) +
  geom_nodelab(aes(label = round(as.numeric(label) * 100)), family = base_font, size = (base_font_size * .9) / .pt, hjust = 1.2, vjust = -.2) +
  geom_rootedge(.25) +
  theme(plot.margin = margin(t = 0, r = 6.25, b = 0, l = -1.35, unit = "char")) +
  coord_cartesian(clip = "off")

ggsave("st_tree_cs.pdf", width = wd, height = wd, units = "in", device = cairo_pdf)
plot_crop("st_tree_cs.pdf")


cons <- consensus(tt, p = .5, rooted = TRUE) |> ladderize()
cons$node.label <- round(cons$node.label * 100, 0)

# Consensus topology + branch lengths. This might take a couple of minutes
ct <- consensus.edges(tt, consensus.tree = cons, rooted = TRUE)
ct$root.edge <- .15

pdf("fig_consensus.pdf", pointsize = 10, width = 6.3, height = 4, family = "URWPalladio")
plot(ct, show.node.label=FALSE, root.edge = TRUE, no.margin = TRUE, font = 1, label.offset = .05)
nodelabels(c(NA, ct$node.label[-1]), frame="none", adj = c(1.15,1.25), cex = .9)
dev.off()
embedFonts("fig_consensus.pdf")
plot_crop("fig_consensus.pdf")

# MCC tree ----------------------------------------------------------------

st_tree_mcc <- maxCladeCred(st_tree)
st_tree_mcc |>
  fortify() |>
  mutate(label = ifelse(branch.length == 0 & label == 1, "", label)) |>
  mutate(label = str_replace_all(label, "_", " ")) |>
  ggtree(ladderize = TRUE, right = FALSE) +
  geom_tiplab(family = base_font, size = base_font_size / .pt) +
  geom_nodelab(aes(label = round(as.numeric(label) * 100)), family = base_font, size = (base_font_size * .9) / .pt, hjust = 1.2, vjust = -.2) +
  geom_rootedge(.25) +
  theme(plot.margin = margin(t = 0, r = 6.5, b = 0, l = -1.35, unit = "char")) +
  coord_cartesian(clip = "off")
ggsave("st_tree_mcc.pdf", width = wd, height = wd, units = "in", device = cairo_pdf)
plot_crop("st_tree_mcc.pdf")

library(treeio)
# mcct <- read.nexus("mcc_ta.tree") |> ladderize()
mcct_b <- read.beast("mcc_ta.tree")
pst <- tibble(posterior = mcct_b@data$posterior, node = as.integer(mcct_b@data$node)) |> 
  filter(!is.na(posterior)) |>
  arrange(node) |> 
  pull(posterior)
mcct_b@phylo$node.labels <- round(pst,2) * 100
mcct <- ladderize(mcct_b@phylo)
mcct$root.edge <- .15
pdf("fig_mcc.pdf", pointsize=10, width = 6.2, height = 4, family = "URWPalladio")
plot(mcct, show.node.label=FALSE, root.edge = TRUE, no.margin = TRUE, font = 1, label.offset = .05)
nodelabels(c(NA, mcct$node.labels[2:12], NA, mcct$node.labels[14:length(mcct$node.labels)]), frame="none", adj = c(1.15,1.25), cex = .9)
nodelabels(c(rep(NA, 12), mcct$node.label[13], rep(NA, length(mcct$node.label) - 13)), frame="none", adj = c(1.15,-.25), cex = .9)
dev.off()
embedFonts("fig_mcc.pdf")
plot_crop("fig_mcc.pdf")

# mcct <- ladderize(maxCladeCred(tt))
mcct$node.label <- round(mcct$node.label * 100, 0)
mcct$root.edge <- .15
pdf("fig_mcc.pdf", pointsize=10, width = 6.3, height = 4, family = "URWPalladio")
plot(mcct, show.node.label=FALSE, root.edge = TRUE, no.margin = TRUE, font = 1, label.offset = .05)
nodelabels(c(NA, mcct$node.label[2], NA, mcct$node.label[c(-1:-3)]), frame="none", adj = c(1.15,1.25))
nodelabels(c(NA, NA, mcct$node.label[3], rep(NA, length(mcct$node.label) - 3)), frame="none", adj = c(1.15,-.25))
dev.off()
embedFonts("fig_mcc.pdf")
plot_crop("fig_mcc.pdf")

# Densitree ---------------------------------------------------------------


# st_tree_ds <- st_tree[1:1000] |>
#   ggdensitree(alpha = .005) +
#   geom_tiplab(aes(label = str_replace_all(label, "_", " ")), family = base_font, size = base_font_size / .pt) +
#   theme_tree2() +
#   theme(plot.margin = margin(t = 0, r = 6.5, b = 0, l = 0, unit = "char")) +
#   coord_cartesian(clip = "off")
# ggsave("st_tree_ds.pdf", st_tree_ds, width = wd, height = wd, units = "in", device = cairo_pdf)
# plot_crop("st_tree_ds.pdf")

st_tree_scaled <- lapply(st_tree, function(x) {
  x$edge.length <- x$edge.length * 1000
  return(x)
})
class(st_tree_scaled) <- "multiPhylo"
pdf("st_tree_ds.pdf", pointsize = 10, width = wd, height = wd * 1.25, family = "URWPalladio")
par(mar = c(2, 0, 0, .6), oma = c(0, 0, 0, 0), xpd = TRUE)
densiTree(st_tree_scaled, consensus = ladderize(st_tree_cs, right = FALSE), alpha = .005, font = 1, label.offset = .01, cex = 1, scale.bar = TRUE)
title(xlab = "years BP")
dev.off()
embedFonts("st_tree_ds.pdf")
knitr::plot_crop("st_tree_ds.pdf")


tt_scaled <- lapply(tt, function(x) {x$edge.length <- x$edge.length * 1000; return(x)})
class(tt_scaled) <- "multiPhylo"
maxBT <- max(phangorn:::getAges(tt_scaled))
label <- rev(pretty(c(maxBT, 0)))
maxBT <- max(label)
pdf("fig_densitree.pdf", pointsize=10, width = 5.025, height = 6, family = "URWPalladio")
par(mar = c(2, 0, 0, .6), #oma = c(0, 0, 0, 0), 
    xpd = TRUE)
densiTree(tt_scaled, consensus = cons, alpha=.005, font = 1, label.offset = .01, cex=1, scale.bar = FALSE)
axis(side = 1, at = seq(0, 1.0, length.out = length(label)), labels = label, line=-1.5, cex=.9)
title(xlab="years BP", line = 1)
dev.off()
embedFonts("fig_densitree.pdf")
plot_crop("fig_densitree.pdf")

# Outgroup ----------------------------------------------------------------

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

outgroup <- sapply(st_tree, find.outgroup)
outgroup_tb <- table(outgroup) |> 
  prop.table() |> 
  as_tibble() |> 
  arrange(-n) |> 
  rename(p = n)
st_tree |>
  map(function(x) {
    nl <- length(x$tip.label)
    x$tip.label
    # r <- getRoot(x)
    # cr <- x$edge[x$edge[, 1] == r, 2]
    # s1 <- getDescendants(x, cr[1])
    # s2 <- getDescendants(x, cr[2])
    # s1 <- s1[s1 <= nl]
    # s2 <- s2[s2 <= nl]
    # if (length(s1) > length(s2)) s1 <- s2
    # sort(x$tip.label[s1])
  })


sapply(st_tree, function(x) {
  tibble(age = max(node.depth.edgelength(x)), outgroup = find.outgroup(x))
})

ages_outgroup <- st_tree |>
  seq_along() |>
  map_df(~
    tibble(
      age = max(node.depth.edgelength(st_tree[[.x]])),
      outgroup = find.outgroup(st_tree[[.x]])
    )) |>
  left_join(outgroup_tb) |> 
  mutate(outgroup = case_when(
    outgroup == "Beijing_Chinese Guangzhou_Chinese Jieyang_Chinese Xingning_Chinese" ~ "Chinese",
    outgroup == "Beijing_Chinese Guangzhou_Chinese Jieyang_Chinese Jingpho Rabha Xingning_Chinese" ~ "Chinese-Sal",
    outgroup == "Bokar_Tani Yidu" ~ "Tani-Yidu",
    .default = "x"
  )) |> 
  mutate(label = paste0(outgroup, "\n(", round(p,2) * 100, "%)")) %>%
  bind_rows(mutate(., outgroup = "any", label = "any")) |>
  filter(!str_detect(outgroup, "x"))

ages_outgroup %>%
  # mutate(outgroup = factor(outgroup, levels = c("Tani-Yidu", "Chinese-Sal", "Chinese", "any"))) %>%
  ggplot(aes(x = age * 1000, y = fct_rev(label), fill = fct_rev(label), height = after_stat(density))) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2, color = "white") +
  geom_density_ridges(fill = NA, color = "gray40") +
  scale_fill_manual(values = c(rep(few_pal("Light")(2)[1], 3), "grey"), guide = "none") +
  scale_x_reverse(limits = c(15000, 0)) +
  scale_y_discrete(expand = expansion(add = c(0.25, 1.4))) +
  xlab("years BP") +
  ylab("first branch") +
  theme(aspect.ratio = 0.618)
ggsave("fig_ageoutgroup.pdf", width = wd, height = wd * .8, units = "in", device = cairo_pdf)
plot_crop("fig_ageoutgroup.pdf")


# Age and monophyly -------------------------------------------------------

getMRCA_age <- function(tree, tips) {
  tips <- if (is.character(tips)) which(tree$tip.label %in% tips) else tips
  mrca <- ifelse(length(tips) > 1, getMRCA(tree, tips), tips)
  root_age <- max(node.depth.edgelength(tree))
  root_age - node.depth.edgelength(tree)[mrca]
}

tips <- st_tree[[1]]$tip.label |>
  str_subset("Jingpho|Rabha|Chinese")
sinitic_sal <- st_tree |>
  seq_along() |>
  map_df(~
    tibble(
      age = getMRCA_age(st_tree[[.x]], tips),
      monophyletic = ifelse(is.monophyletic(st_tree[[.x]], tips), "monophyletic", "paraphyletic")
    )) %>%
  bind_rows(mutate(., monophyletic = "any"))

sinitic_sal |> 
  count(monophyletic) |> 
  mutate(p = n/2000)

sinitic_sal |>
  ggplot(aes(x = age * 1000, y = fct_rev(monophyletic), fill = monophyletic, height = after_stat(density))) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2, color = "white") +
  geom_density_ridges(fill = NA, color = "gray40") +
  scale_fill_manual(values = c("grey", rep(few_pal("Light")(2)[1], 2)), guide = "none") +
  scale_x_reverse(limits = c(15000, 0)) +
  scale_y_discrete(expand = expansion(add = c(0.25, 1.25))) +
  xlab("years BP") +
  ylab("phyletic status of Chinese-Sal") +
  theme(aspect.ratio = 0.618)
ggsave("fig_agemono.pdf", width = wd, height = wd * .8, units = "in", device = cairo_pdf)
plot_crop("fig_agemono.pdf")


known <- tribble(
  ~language, ~lower, ~upper,
  "Old Tibetan", 1000, 1200,
  "Old Burmese", 700, 900,
  "Common Chinese", 2000, 2200,
  "Old Chinese", 2400, 2600
) |>
  mutate(type = "known")

oldburmese_trees <- read.nexus("xval/Burmese.nex")
commonchinese_trees <- read.nexus("xval/CommonChinese.nex")
oldchinese_trees <- read.nexus("xval/Chinese.nex")
oldtibetan_trees <- read.nexus("xval/Tibetan.nex")


getMRCA_ages <- function(trees, tips, burnin = 0, label = NULL) {
  if (is.null(label)) label <- tips
  trees <- trees[round(length(trees) * burnin):length(trees)]
  seq_along(trees) |>
    map_dbl(~ getMRCA_age(trees[[.x]], tips)) |>
    hdi() |>
    rbind() |>
    as_tibble() |>
    mutate(language = label, type = "inferred")
}

bind_rows(
  known,
  getMRCA_ages(oldchinese_trees, "SiniticOldChinese", .2, "Old Chinese"),
  getMRCA_ages(oldtibetan_trees, "TibetanOldTibetan", .2, "Old Tibetan"),
  getMRCA_ages(oldburmese_trees, "BurmishOldBurmese", .2, "Old Burmese"),
  getMRCA_ages(commonchinese_trees, commonchinese_trees[[1]]$tip.label |> str_subset("Sinitic[^O]"), .2, "Common Chinese")
) |> 
  mutate(language = str_replace(language, " ", "\n")) |> 
  ggplot(aes(y=language, x=(lower+upper)/2, xmin=lower, xmax=upper)) +
  geom_linerange(aes(color=fct_rev(type)), position=position_dodge2(reverse = TRUE, width=c(.35)), linewidth=4) +
  xlab("years BP") +
  ylab(NULL) +
  scale_color_few(name = NULL) +
  scale_x_reverse(limits = c(3000, 0)) +
  scale_y_discrete(limits=rev, expand = expansion(add = c(0.25, .25)))  +
  theme(aspect.ratio = 0.618, legend.position = "top")
ggsave("fig_xages.pdf", width = wd, height = wd * .8, units = "in", device = cairo_pdf)
plot_crop("fig_xages.pdf")


# st_tree |>
#   seq_along() |>
#   map_lgl(~ is.monophyletic(st_tree[[.x]], tips = c("Chepang", "Hayu", "Thulung", "Bokar_Tani", "Yidu", "Tshangla"))) |>
#   mean()
#
# st_tree |>
#   map(~map(.x, print))
#
# sapply(tt, is.monophyletic, tips = c("Chepang", "Hayu", "Thulung", "Bokar_Tani", "Yidu", "Tshangla")) |>
#   mean()

library(patchwork)
mcc_ph <- mcc(st_tree) |> ladderize()
mcc_ta_target <- read.nexus("mcc_ta_target.tree") |> ladderize()
mcc_ta <- read.nexus("mcc_ta.tree") |> ladderize()
mcc_ta_median <- read.nexus("mcc_ta_median.tree") |> ladderize()
mcc_ta_mean <- read.nexus("mcc_ta_median.tree") |> ladderize()

theme_set(theme_tree2() + theme(plot.margin = margin(t = 0, r = 8, b = 0, l = 0, unit = "char"))
)
ph <- ggtree(mcc_ph, ladderize = TRUE) + geom_tiplab() + ggtitle("phangorn\n") + theme_tree2() + coord_fixed(clip = "off")
tatg <- ggtree(mcc_ta_target, ladderize = TRUE) + geom_tiplab() + ggtitle("TreeAnnotator\ntarget heights") + theme_tree2() + coord_fixed(clip = "off")
tach <- ggtree(mcc_ta, ladderize = TRUE) + geom_tiplab() + ggtitle("TreeAnnotator\ncommon ancestors heights") + theme_tree2() + coord_fixed(clip = "off")
tamd <- ggtree(mcc_ta_median, ladderize = TRUE) + geom_tiplab() + ggtitle("TreeAnnotator\nmedian heights") + theme_tree2() + coord_fixed(clip = "off")
tame <- ggtree(mcc_ta_mean, ladderize = TRUE) + geom_tiplab() + ggtitle("TreeAnnotator\nmean heights")  + theme_tree2() + coord_fixed(clip = "off")

cowplot::plot_grid(revts(ph), revts(tamd), revts(tame), revts(tatg), revts(tach))



par(mfrow=c(2,2), xpd=FALSE)
plot(mcc_ph)
add.scale.bar()
legend("bottomright", legend="phangorn", bty="n")
plot(mcc_ta)
legend("bottomright", legend="TreeAnnotator\ncommon ancestors heights", bty="n")
plot(mcc_ta_median)
legend("bottomright", legend="TreeAnnotator\nmedian heights", bty="n")
plot(mcc_ta_mean)
legend("bottomright", legend="TreeAnnotator\nmean heights", bty="n")
dev.off()

