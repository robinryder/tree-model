# code to reproduce the cross-validation of dates figure
# Robin J. Ryder, April 2023

library(ape)
library(HDInterval)
library(ggplot2)
library(phytools)
library(phangorn)

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


recage.burmish = reconstructed.age("Burmese.nex", 
                      c("BurmishOldBurmese", "BurmishRangoon"))

recage.commonchinese = reconstructed.age("CommonChinese.nex", 
                      c("SiniticBeijing", "SiniticGuangzhou", 
                        "SiniticXingning", "SiniticLonggang",
                        "SiniticJieyang", "SiniticChaozhou" ))
recage.sinitic = reconstructed.age("Chinese.nex", 
                      c("SiniticOldChinese" ))
recage.tibetan = reconstructed.age("Tibetan.nex", 
                      c("TibetanOldTibetan"))


df = as.data.frame(rbind(true.burmish, true.commonchinese, true.sinitic, true.tibetan,
           recage.burmish, recage.commonchinese, recage.sinitic, recage.tibetan))
df$type = rep(c("Constraint", "Reconstructed"), each=4)
df$clade = rep(c("Burmish", "Common Chinese", "Sinitic", "Tibetan"), 2)

ggplot(df, aes(x=clade, y=(lower+upper)/2, ymin=lower, ymax=upper)) +
  geom_linerange(aes(color=type), position=position_dodge(width=c(0.2)), size=1) +
  ylab("Years BP")
