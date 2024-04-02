library(ape)
library(phangorn)
library(phytools)
library(tidyverse)
library(ggtree)
library(ggthemes)
library(knitr)

wd <- 5
base_font <- "URW Palladio L"
base_font_size <- 10
theme_set(
  theme_minimal(base_family = base_font, base_size = base_font_size) +
    theme(text = element_text(family = base_font, size = base_font_size))
)
update_geom_defaults("text", list(family = base_font, size = base_font_size / .pt))


st_tree <- read.nexus("SinoTibetanSubset.nex")


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

st_tree_ds <- st_tree[1:1000] |>
  ggdensitree(alpha = .005) +
  geom_tiplab(aes(label = str_replace_all(label, "_", " ")), family = base_font, size = base_font_size / .pt) +
  theme_tree2() +
  theme(plot.margin = margin(t = 0, r = 6.5, b = 0, l = 0, unit = "char")) +
  coord_cartesian(clip = "off")
ggsave("st_tree_ds.pdf", st_tree_ds, width = wd, height = wd, units = "in", device = cairo_pdf)
plot_crop("st_tree_ds.pdf")

st_tree_scaled <- lapply(st_tree, function(x) {x$edge.length <- x$edge.length * 1000; return(x)})
class(st_tree_scaled) <- "multiPhylo"
pdf("st_tree_ds.pdf", pointsize = 10, width = wd, height = wd * 1.25, family = "URWPalladio")
par(mar = c(2, 0, 0, .6), oma = c(0, 0, 0, 0), xpd=TRUE)
densiTree(st_tree_scaled, consensus = ladderize(st_tree_cs, right = FALSE), alpha = .005, font = 1, label.offset = .01, cex = 1, scale.bar = TRUE)
title(xlab = "years BP")
dev.off()
embedFonts("st_tree_ds.pdf")
knitr::plot_crop("st_tree_ds.pdf")
