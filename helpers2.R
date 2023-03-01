## ==================================================================================== ##
# START Shiny App for analysis and visualization of transcriptome data.
# Copyright (C) 2016  Jessica Minnier
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# 
# You may contact the author of this code, Jessica Minnier, at <minnier@ohsu.edu>
## ==================================================================================== ##

library(shiny) #devtools::install_github("rstudio/shiny"); devtools::install_github("rstudio/shinyapps")
library(reshape2)
library(ggplot2)
library(ggthemes)
#library(shinyIncubator) #devtools::install_github("rstudio/shinyIncubator")
library(gplots)
#library(rjson)
#library(base64enc)
library(ggvis)
library(dplyr)
library(tidyr)
library(DT) #devtools::install_github('ramnathv/htmlwidgets'); devtools::install_github('rstudio/DT')
library(limma)
#library(DESeq2)
library(edgeR)
library(RColorBrewer)
library(pheatmap)
library(shinyBS)
library(plotly)
library(markdown)
library(NMF)
library(scales)
library(heatmaply)
library(readr)
library(colourpicker)
library(data.table)
library(janitor)

##================================================================================##

source("fun-dotplot2.R")
source("fun-heatmap.R")
source("fun-analyzecounts.R")
source("fun-analysisres.R")
source("fun-groupplots.R")
source("fun-input-analyze-data.R")

#troubleshooting
if(TRUE) {
  seqdata <- read.csv("data/GSE193677_MSCCR_Biopsy_counts3.csv",stringsAsFactors = FALSE)
  seqinfo <- read.csv("data/20230223-GSE193677_phenoData_short.csv",stringsAsFactors = FALSE)
  # rectum.no <- which(seqinfo$regionre.ch1 == "Rectum")
  seqinfo.rectum <- filter(seqinfo, regionre.ch1 == "Rectum")
  sample_column <- data.frame(no = 1:ncol(seqdata), sample.names = colnames(seqdata))
  col.rectum <- which(sample_column$sample.names%in%seqinfo.rectum$sample)
  # seqdata.rectum <- seqdata[, c(1, 2, rectum.no)]
  seqdata.rectum <- seqdata[, c(1, 2, col.rectum)]
  seqdata <- seqdata.rectum
  seqinfo.rectum <- filter(seqinfo, regionre.ch1 == "Rectum")
  seqinfo <- seqinfo.rectum
  seqdata.cpm <- cbind(seqdata.rectum[,1:2], data.frame(cpm(seqdata.rectum[, 3:906], log = FALSE)))
  # seqdata.cpm.sum1 <- sum(seqdata.cpm[,1])
  # seqdata.cpm.sum2 <- sum(seqdata.cpm[,2])
  seqdata.log2cpm <- cbind(seqdata.rectum[,1:2], data.frame(cpm(seqdata.rectum[, 3:906], log = TRUE)))
  # diseases <- unique(seqinfo$ibd_disease.ch1)
  diseases <- factor(c("Control", "UC", "CD", "UC_Pouch", "CD_Pouch"), levels = c("Control", "UC", "CD", "UC_Pouch", "CD_Pouch"))
  Unique.ID <- paste(seqdata$Gene, seqdata$GENEID, sep="_")
  geneids <- data.frame(Unique.ID, Gene.ID = seqdata$Gene, Gene.Symbol = seqdata$GENEID)
  # seqdata <- read.csv("data/mousecounts_example.csv",stringsAsFactors = FALSE)
  # load('data/mousecounts_example_analysis_results.RData')
  # load('data/mousecounts_example_analyzed.RData') #example_data_results
  data_analyzed = list('group_names'=diseases,
                       # 'sampledata'=sampledata,
                       # "results"=results,
                       # "data_long"=data_long, 
                       "seqdata.cpm" = seqdata.cpm,
                       "seqdata.log2cpm" = seqdata.cpm,
                       "geneids"=geneids,
                       "seqinfo"=seqinfo)
                       # "data_results_table"=example_data_results)
  
  # data_results = data_analyzed$results
  
  # test_sel = "CD.Severe.Rectum/Control.NA.rectum"
  # sel_test = test_sel
  # fdrcut = 0.05
  # absFCcut = 1
  group_sel = c("CD.Severe.Rectum","Control.NA.rectum")
  valuename = "log2cpm"
  yname="log2cpm"
  # maxgenes = 200
  # view_group=NULL
  # filter_by_go=FALSE
  # filter_fdr=FALSE
  # filter_maxgene=TRUE
  # filter_cpm=FALSE
  # filter_fc=FALSE
  # fold_change_range=NULL
  # fold_change_groups=NULL
  # group_filter_range =NULL
  # fixed_genes_list=NULL
  # orderby="variance"
  
  # tmpdat = heatmap_subdat(data_analyzed,yname,orderby="variance",maxgenes=100)
  # heatmap_render(data_analyzed,yname,orderby="variance",maxgenes=100)
  # 
  # mydat = heatmap_ggvis_data(
  #   data_analyzed = data_analyzed,
  #   yname = yname,
  #   orderby = "variance",
  #   maxgenes=100)

}
