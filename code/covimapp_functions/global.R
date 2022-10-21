
# Load R packages
library(shiny)
library(shinythemes)
library(ggplot2)
library(stringr)
library(esc)
library(meta)
library(tidyverse)
library(ggrepel)
library(plotly)
library(png)
library(grid)
library(mada)
library(ggforce)

##### Load required data for SMD #####
setwd("shiny_covid/data") 
t1 = read.csv("table_s1.csv")
t2 = read.csv("table_s2.csv")

colnames(t1) = gsub("_", " ", colnames(t1))
colnames(t1)[c(8,11)] = c("N COVID-19 cases", "Proteins in Meta-analysis")
colnames(t2) = gsub("_", " ", colnames(t2))

mt = read.csv("meta_results.csv", row.names = 1)
mts = read.csv("meta_sel.csv", row.names = 1)
df = readRDS("study_summaries.RData")
det_prot = read.csv("detected_proteins.csv", row.names = 1)
dorf = read.csv("dorf.csv", row.names = 1)
dorf$k.study = mt$k.study[match(dorf$Gene_name, mt$gene_name)]
drfsig = dorf[which(dorf$lower.Sensitivity>0.5 & dorf$lower.Specificity>0.5),]

##### Load required data for SROC #####
t_summary = readRDS("t_summary.RData")
mtdats = readRDS("mtdats.RData")
fit.reitsmas = readRDS("fitreitsmas.Rdata")
roc_prot = read.csv("roc_prot.csv", row.names = 1)

##### Specify input choices #####
protein_choices = mt$gene_name
ns = length(t2$Author)
studies = t2$Author
studies2 = c("Babacic et al.", "Shen et al.", "Geyer et al.", "Shu et al.", "Di et al.",
             "Messner et al. - Discovery", "Messner et al. - Validation", "Suvarna et al.",
             "Sullivan et al.", "Overmyer et al.", "Zhang et al.")

save.image("~/covimapp/data_files.RData")

