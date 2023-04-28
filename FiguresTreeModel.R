# Code for figures from chapter on The Tree Model (Jacques, Pellard & Ryder 2023)
# RJR March 2023

library(phytools)
library(ape)
library(phangorn)
library(ggplot2)
library(cowplot)


tt = read.nexus("SinoTibetanSubset.nex")

# Consensus topology
cons = consensus(tt, p=.5, rooted=T)
cons$node.label = round(cons$node.label, 2)
plot.phylo(cons, show.node.label = T) 

# Consensus topology + branch lengths. This might take a couple of minutes
ct = consensus.edges(tt, consensus.tree = cons, rooted=T)
plot(ct, show.node.label=T) 



plot.phylo(ct, show.node.label = T, main ="Consensus tree")



# Maximum clade credibility tree
mcc = mcc(tt)
mcc$node.label = round(mcc$node.label, 2)
mcc$node.label[getRoot(mcc)] = NA

plot(mcc, show.node.label = T, main="MCC")



densiTree(tt, alpha=.008, consensus = cons, main="Densitree")




# posterior probability of a specific subtree
subtree_support = function(tips, trees = tt){
  sapply(trees, is.monophyletic, tips=tips)
}
mean(subtree_support(c("Chepang","Hayu","Thulung","Bokar_Tani","Yidu","Tshangla")))
# The 6 leaves listed above form a subtree with posterior probability 0.41



# Creating the plot for the root age depending on the outgroup
# Lists of outgroups
find.outgroup = function(tree, nl = length(tree$tip.label)){
  r = getRoot(tree)
  cr = tree$edge[tree$edge[, 1] == r , 2]
  s1 = getDescendants(tree, cr[1])
  s2 = getDescendants(tree, cr[2])
  
  s1 = s1[s1<=nl]
  s2 = s2[s2<=nl]
  
  # Ensure that s1 is the smallest
  if(length(s1) > length(s2)) s1 = s2
  return(paste(sort(tree$tip.label[s1]), collapse = " "))
}

outgroup = sapply(tt, find.outgroup, nl=22)
head(sort(table(outgroup), dec=T), 3)

outgroup[outgroup == "Beijing_Chinese Guangzhou_Chinese Jieyang_Chinese Xingning_Chinese"] = "Chinese"
outgroup[outgroup == "Beijing_Chinese Guangzhou_Chinese Jieyang_Chinese Jingpho Rabha Xingning_Chinese"] = "Chinese + Sal"
outgroup[outgroup == "Bokar_Tani Yidu"] = "Tani-Yidu"

root.age = function(tree){ node.depth.edgelength(tree)[1]}
ra = sapply(tt, root.age)
aged = data.frame(age = ra * 1000, outgroup = outgroup)
keep = outgroup %in% c("Chinese", "Tani-Yidu", "Chinese + Sal")

aged2 = rbind(aged[keep, ], data.frame(age = ra*1000, outgroup="all"))


p.age.group = ggplot(aged[keep,], aes(x=age)) + 
  xlim(0, 15000) +
  geom_density(aes(group=outgroup, colour=outgroup, fill=outgroup), alpha=.3)

col="purple"
p.age = ggplot(aged, aes(x=age)) +
  xlim(0, 15000) +
  geom_density(fill=col, colour=col, alpha=.3)
           
plot_grid(p.age, p.age.group, align="v", axis="rblt", ncol=1)


###### MRCA ages

mrca.age = function(tree, tips){
  max(dist.nodes(tree)[tips, tips])/2
}

mrca.age.names = function(tree, tipnames){
  max(cophenetic(tree)[tipnames, tipnames]) / 2
}

mrca.summary = function(d, l){
  nt = length(d)
  res = rep(NA, nt)
  for(i in 1:nt){
    if(is.numeric(l)){
      res[i] = mrca.age(d[[i]], l)
    }
    else{
      res[i] = mrca.age.names(d[[i]], l)
    }
  }
  
  return(res)
}

# Plot of MRCA age for the following languages: Jingpho, Rabha, 4 Chinese dialects
ll = tt[[1]]$tip.label
l = c(7, 12, 13:16) # indices of the languages of interest in vector ll
ages = mrca.summary(tt, l)
lab = sapply(tt, is.monophyletic, tips=l)
df = data.frame(age=ages, monophyletic=lab)

p.age.topo.group = ggplot(df, aes(x=age)) + 
  geom_density(aes(group=monophyletic, colour=monophyletic, fill=monophyletic), alpha=.3) 

p.age.topo.all = ggplot(df, aes(x=age)) + 
  geom_density(fill=col, colour=col, alpha=.3)

title = ggdraw() + draw_label("Age of MRCA" , fontface='bold')

plot_grid(title, p.age.topo.all, p.age.topo.group, align="v", axis="rblt", ncol=1, rel_heights=c(0.2, 1, 1))



