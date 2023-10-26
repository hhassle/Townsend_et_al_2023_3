#load libraries
library(geiger)
library(phytools)
library(tidyverse)
library(viridis)
library(ggridges)
library(Rphylopars)
library(ggplot2)
library(ggpubr)

#set paths
rootPath<-"~/data"
treePath<-"~/trees"

#read in tree
trees<-list.files(treePath,full.name=TRUE)

#get all directories except the root path
directories<-list.dirs(rootPath)[-1]

#get all sites using a crazy regular expression to just isolate the site names
#the \\ is how R escapes regular expressions which makes this look funny
localities<-NULL
for (i in 1:length(directories)){
  localities[i] <- str_extract(directories[i], "(?<=\\/)[^\\/]+$")
}

#function to compute the aic scores, weights, and likelihoods and write into a dataframe for all the files in a directory
compareModels<-function(tree, directory, locality){
  monthPaths<-list.files(directory, pattern = "\\.csv$", full.names = TRUE)
  map_df(monthPaths, ~{
    # Read the file
    file_data <- read.csv(.x,h=T)
    #print to track progress
    print(.x)
    #Fit the models
    star <- phylopars(file_data,tree,model="star",phylo_correlated = TRUE, pheno_correlated = TRUE)
    brown <- phylopars(file_data,tree,model="BM",phylo_correlated = TRUE, pheno_correlated = TRUE)
    lambda <- phylopars(file_data,tree,model="lambda",phylo_correlated = TRUE, pheno_correlated = TRUE)
    tree1 <- force.ultrametric(tree)
    OU <- phylopars(file_data,tree1,model="OU",phylo_correlated = TRUE, pheno_correlated = TRUE)
    # Perform calculations and store the results in a tibble, I'm putting line breaks here for readability
    Summary <- tibble(
      treee = treess,
      locality = locality,
      file = .x,
      x229E_brown = brown$anc_recon[1],
      NL63_brown = brown$anc_recon[2],
      HKU1_brown = brown$anc_recon[3],
      OC43_brown = brown$anc_recon[4],
      CoV2_brown = brown$anc_recon[5],
      CoV2_star = star$anc_recon[5],
      CoV2_lambda = lambda$anc_recon[5],
      CoV2_OU = OU$anc_recon[5]
    )
    
    # Return the calculations
    Summary
  })
}

#run function across directories and summarize results
#empty data frame to store results
result_df <- data.frame()
#I love loops
for (i in 1: length(directories)){
  for (x in 1: length(trees)){
    treess<-trees[x]
    tree<-read.tree(trees[x])
    df<-compareModels(tree, directories[i], localities[i])
    result_df<-rbind(result_df, df)
  }  
}

#write out results to data root directory
setwd(list.dirs(rootPath)[1])
write.csv(result_df, file="SummaryStats.csv")

#plot the AIC weights
data1<- read.csv("~/data/SummaryStats.csv",h=T)
ggplot(subset(data1, treee %in% c("~/trees/IQTREE_time_tree.newick")), aes(x=CoV2_brown, y=CoV2_star)) +
  geom_smooth(method="lm") +
  geom_point(size = 1) + stat_regline_equation(label.x=1, label.y=26) +theme_classic() +scale_x_continuous(breaks=seq(0, 30,1))+scale_y_continuous(breaks=seq(0, 30,10))+
  stat_cor(aes(label=..rr.label..), label.x=1, label.y=24)
ggsave(
  "~/data/brown_star.pdf",
  plot = last_plot(),width = 3, height = 3)
ggplot(data=data1, aes(x=CoV2_brown, y=CoV2_lambda)) +
  geom_smooth(method="lm") +
  geom_point(size = 1) + stat_regline_equation(label.x=1, label.y=26) +theme_classic() +scale_x_continuous(breaks=seq(0, 30,1))+scale_y_continuous(breaks=seq(0, 30,10))+
  stat_cor(aes(label=..rr.label..), label.x=1, label.y=24)
ggsave(
  "~/data/lambda_brown.pdf",
  plot = last_plot(),width = 3, height = 3)
ggplot(subset(data1, treee %in% c("~/trees/IQTREE_time_tree.newick"))
       , aes(x=CoV2_brown, y=CoV2_OU)) +
  geom_smooth(method="lm") +
  geom_point(size = 1) + stat_regline_equation(label.x=1, label.y=26) +theme_classic() +scale_x_continuous(breaks=seq(0, 30,1))+scale_y_continuous(breaks=seq(0, 30,10))+
  stat_cor(aes(label=..rr.label..), label.x=1, label.y=24)
ggsave(
  "~/data/ou_brown.pdf",
  plot = last_plot(),width = 3, height = 3)
ggplot(subset(data1, treee %in% c("~/trees/IQTREE_time_tree.newick")), aes(x=CoV2_star, y=CoV2_lambda)) +
  geom_smooth(method="lm") +
  geom_point(size = 1) + stat_regline_equation(label.x=1, label.y=26) +theme_classic() +scale_x_continuous(breaks=seq(0, 30,1))+scale_y_continuous(breaks=seq(0, 30,10))+
  stat_cor(aes(label=..rr.label..), label.x=1, label.y=24)
ggsave(
  "~/data/lambda_star.pdf",
  plot = last_plot(),width = 3, height = 3)
ggplot(subset(data1, treee %in% c("~/trees/IQTREE_time_tree.newick")), aes(x=CoV2_star, y=CoV2_OU)) +
  geom_smooth(method="lm") +
  geom_point(size = 1) + stat_regline_equation(label.x=1, label.y=26) +theme_classic() +scale_x_continuous(breaks=seq(0, 30,1))+scale_y_continuous(breaks=seq(0, 30,10))+
  stat_cor(aes(label=..rr.label..), label.x=1, label.y=24)
ggsave(
  "~/data/star_ou.pdf",
  plot = last_plot(),width = 3, height = 3)
ggplot(subset(data1, treee %in% c("~/trees/IQTREE_time_tree.newick")), aes(x=CoV2_lambda, y=CoV2_OU)) +
  geom_smooth(method="lm") +
  geom_point( size = 1) + stat_regline_equation(label.x=1, label.y=26) +theme_classic() +scale_x_continuous(breaks=seq(0, 30,1))+scale_y_continuous(breaks=seq(0, 30,10))+
  stat_cor(aes(label=..rr.label..), label.x=1, label.y=24)
ggsave(
  "~/data/lambda_ou.pdf",
  plot = last_plot(),width = 3, height = 3)


