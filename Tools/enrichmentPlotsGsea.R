# This R script reads lists of GSEA generated enrichment TSV files, 
# and generates a plot with all the results in one panel

############################### Setup
# Turn off warnings
options(warn=-1)
# Load libraries
# If you do not have these packages installed, please install and load them before you proceed
#BiocManager::install('clusterProfiler')
library(foreach)
library(tidyverse)
library(clusterProfiler)
library(multipanelfigure)
library(stringi)

# Set the directory path for input and output
# Place all your GSEA enrichment report TSV files in one folder ('gseareports' in the below example)
#directorypath="/Users/nagarajanv/OneDrive - National Institutes of Health/Ryan_Vijay/NewAnalysis2021/GSEA_Results/tsv/"
#directorypath="/Users/nagarajanv/OneDrive - National Institutes of Health/Ryan_Vijay/NewAnalysisCytotoxic2021/GSEA_Results/tsv/"
directorypath="/Users/nagarajanv/OneDrive - National Institutes of Health/Ryan_Vijay/NewAnalysisCD4CD8_2021/tsv/"


# Initiate list variable to store enrichment lists
datalistentrez=list()

# Assign sample names here
# Sample name should be the prefix of the GSEA enrichment report TSV files
# For example, KO_neg is the sample name, if the file name is KO_neg.tsv
#samples=c("MLN_LAQ","SILP_LAQ","MLN_PBS","SILP_PBS")
#samples=c("IEL_CD4_POS","IEL_DP_POS","IEL_CD8p_CD90p_POS","IEL_CD8p_CD90n_POS","SI_CD4_POS","SI_DP_POS","SI_CD8p_CD90p_POS","MES_CD4_POS","MES_DP_POS","MES_CD8p_CD90p_POS","PLN_CD4_POS","PLN_CD8p_CD90p_POS","IEL_CD4_NEG","IEL_DP_NEG","IEL_CD8p_CD90p_NEG","IEL_CD8p_CD90n_NEG","SI_CD4_NEG","SI_DP_NEG","SI_CD8p_CD90p_NEG","MES_CD4_NEG","MES_DP_NEG","MES_CD8p_CD90p_NEG","PLN_CD4_NEG","PLN_CD8p_CD90p_NEG")
#samples=c("IEL_CD4_POS","IEL_DP_POS","IEL_CD8p_CD90p_POS","IEL_CD8p_CD90n_POS","SI_CD4_POS","SI_DP_POS","SI_CD8p_CD90p_POS","MES_CD4_POS","MES_DP_POS","MES_CD8p_CD90p_POS","PLN_CD4_POS","PLN_CD8p_CD90p_POS")
#samples=c("IEL_CD4","IEL_DP","IEL_CD8p_CD90p","IEL_CD8p_CD90n","SI_CD4","SI_DP","SI_CD8p_CD90p","MES_CD4","MES_DP","MES_CD8p_CD90p","PLN_CD4","PLN_CD8p_CD90p")
#samples=c("PLN_NT_d0_pos","PLN_ST_d0_pos","PLN_LT_d0_pos","MLN_NT_d0_pos","MLN_ST_d0_pos","MLN_LT_d0_pos","SILP_NT_d0_pos","SILP_ST_d0_pos","SILP_LT_d0_pos","IEL_NT_d0_pos","IEL_ST_d0_pos","IEL_LT_d0_pos","PLN_NT_d7_pos","PLN_ST_d7_pos","PLN_LT_d7_pos","MLN_NT_d7_pos","MLN_ST_d7_pos","MLN_LT_d7_pos","SILP_NT_d7_pos","SILP_ST_d7_pos","SILP_LT_d7_pos","IEL_NT_d7_pos","IEL_ST_d7_pos","IEL_LT_d7_pos","PLN_NT_d0_neg","PLN_ST_d0_neg","PLN_LT_d0_neg","MLN_NT_d0_neg","MLN_ST_d0_neg","MLN_LT_d0_neg","SILP_NT_d0_neg","SILP_ST_d0_neg","SILP_LT_d0_neg","IEL_NT_d0_neg","IEL_ST_d0_neg","IEL_LT_d0_neg","PLN_NT_d7_neg","PLN_ST_d7_neg","PLN_LT_d7_neg","MLN_NT_d7_neg","MLN_ST_d7_neg","MLN_LT_d7_neg","SILP_NT_d7_neg","SILP_ST_d7_neg","SILP_LT_d7_neg","IEL_NT_d7_neg","IEL_ST_d7_neg","IEL_LT_d7_neg")
#samples=c("PLN_NT_d0_pos","PLN_ST_d0_pos","PLN_LT_d0_pos","MLN_NT_d0_pos","MLN_ST_d0_pos","MLN_LT_d0_pos","SILP_NT_d0_pos","SILP_ST_d0_pos","SILP_LT_d0_pos","IEL_NT_d0_pos","IEL_ST_d0_pos","IEL_LT_d0_pos","PLN_NT_d7_pos","PLN_ST_d7_pos","PLN_LT_d7_pos","MLN_NT_d7_pos","MLN_ST_d7_pos","MLN_LT_d7_pos","SILP_NT_d7_pos","SILP_ST_d7_pos","SILP_LT_d7_pos","IEL_NT_d7_pos","IEL_ST_d7_pos","IEL_LT_d7_pos")
#samples=c("PLN_NT_d0_neg","PLN_ST_d0_neg","PLN_LT_d0_neg","MLN_NT_d0_neg","MLN_ST_d0_neg","MLN_LT_d0_neg","SILP_NT_d0_neg","SILP_ST_d0_neg","SILP_LT_d0_neg","IEL_NT_d0_neg","IEL_ST_d0_neg","IEL_LT_d0_neg","PLN_NT_d7_neg","PLN_ST_d7_neg","PLN_LT_d7_neg","MLN_NT_d7_neg","MLN_ST_d7_neg","MLN_LT_d7_neg","SILP_NT_d7_neg","SILP_ST_d7_neg","SILP_LT_d7_neg","IEL_NT_d7_neg","IEL_ST_d7_neg","IEL_LT_d7_neg")
#samples=c("PLN_NT_d0","PLN_ST_d0","PLN_LT_d0","MLN_NT_d0","MLN_ST_d0","MLN_LT_d0","SILP_NT_d0","SILP_ST_d0","SILP_LT_d0","IEL_NT_d0","IEL_ST_d0","IEL_LT_d0","PLN_NT_d7","PLN_ST_d7","PLN_LT_d7","MLN_NT_d7","MLN_ST_d7","MLN_LT_d7","SILP_NT_d7","SILP_ST_d7","SILP_LT_d7","IEL_NT_d7","IEL_ST_d7","IEL_LT_d7")

#samples=c("IEL_LT_d0.vs..IEL_NT_d0",	"IEL_LT_d7.vs..IEL_NT_d0",	"IEL_NT_d7.vs..IEL_NT_d0",	"IEL_ST_d0.vs..IEL_NT_d0",	"IEL_ST_d7.vs..IEL_NT_d0",	"MLN_LT_d0.vs..IEL_NT_d0",	"MLN_LT_d7.vs..IEL_NT_d0",	"MLN_NT_d0.vs..IEL_NT_d0",	"MLN_NT_d7.vs..IEL_NT_d0",	"MLN_ST_d0.vs..IEL_NT_d0",	"MLN_ST_d7.vs..IEL_NT_d0",	"PLN_LT_d0.vs..IEL_NT_d0",	"PLN_LT_d7.vs..IEL_NT_d0",	"PLN_NT_d0.vs..IEL_NT_d0",	"PLN_NT_d7.vs..IEL_NT_d0",	"PLN_ST_d0.vs..IEL_NT_d0",	"PLN_ST_d7.vs..IEL_NT_d0",	"SILP_LT_d0.vs..IEL_NT_d0",	"SILP_LT_d7.vs..IEL_NT_d0",	"SILP_NT_d0.vs..IEL_NT_d0",	"SILP_NT_d7.vs..IEL_NT_d0",	"SILP_ST_d0.vs..IEL_NT_d0",	"SILP_ST_d7.vs..IEL_NT_d0")
#samples=c("IEL_LT_d0.vs..SILP_NT_d0",	"IEL_LT_d7.vs..SILP_NT_d0",	"IEL_NT_d0.vs..SILP_NT_d0",	"IEL_NT_d7.vs..SILP_NT_d0",	"IEL_ST_d0.vs..SILP_NT_d0",	"IEL_ST_d7.vs..SILP_NT_d0",	"MLN_LT_d0.vs..SILP_NT_d0",	"MLN_LT_d7.vs..SILP_NT_d0",	"MLN_NT_d0.vs..SILP_NT_d0",	"MLN_NT_d7.vs..SILP_NT_d0",	"MLN_ST_d0.vs..SILP_NT_d0",	"MLN_ST_d7.vs..SILP_NT_d0",	"PLN_LT_d0.vs..SILP_NT_d0",	"PLN_LT_d7.vs..SILP_NT_d0",	"PLN_NT_d0.vs..SILP_NT_d0",	"PLN_NT_d7.vs..SILP_NT_d0",	"PLN_ST_d0.vs..SILP_NT_d0",	"PLN_ST_d7.vs..SILP_NT_d0",	"SILP_LT_d0.vs..SILP_NT_d0",	"SILP_LT_d7.vs..SILP_NT_d0",	"SILP_NT_d7.vs..SILP_NT_d0",	"SILP_ST_d0.vs..SILP_NT_d0",	"SILP_ST_d7.vs..SILP_NT_d0")
#samples=c("IEL_LT_d0.vs..MLN_NT_d0",	"IEL_LT_d7.vs..MLN_NT_d0",	"IEL_NT_d0.vs..MLN_NT_d0",	"IEL_NT_d7.vs..MLN_NT_d0",	"IEL_ST_d0.vs..MLN_NT_d0",	"IEL_ST_d7.vs..MLN_NT_d0",	"MLN_LT_d0.vs..MLN_NT_d0",	"MLN_LT_d7.vs..MLN_NT_d0",	"MLN_NT_d7.vs..MLN_NT_d0",	"MLN_ST_d0.vs..MLN_NT_d0",	"MLN_ST_d7.vs..MLN_NT_d0",	"PLN_LT_d0.vs..MLN_NT_d0",	"PLN_LT_d7.vs..MLN_NT_d0",	"PLN_NT_d0.vs..MLN_NT_d0",	"PLN_NT_d7.vs..MLN_NT_d0",	"PLN_ST_d0.vs..MLN_NT_d0",	"PLN_ST_d7.vs..MLN_NT_d0",	"SILP_LT_d0.vs..MLN_NT_d0",	"SILP_LT_d7.vs..MLN_NT_d0",	"SILP_NT_d0.vs..MLN_NT_d0",	"SILP_NT_d7.vs..MLN_NT_d0",	"SILP_ST_d0.vs..MLN_NT_d0",	"SILP_ST_d7.vs..MLN_NT_d0")
#samples=c("IEL_LT_d0.vs..PLN_NT_d0",	"IEL_LT_d7.vs..PLN_NT_d0",	"IEL_NT_d0.vs..PLN_NT_d0",	"IEL_NT_d7.vs..PLN_NT_d0",	"IEL_ST_d0.vs..PLN_NT_d0",	"IEL_ST_d7.vs..PLN_NT_d0",	"MLN_LT_d0.vs..PLN_NT_d0",	"MLN_LT_d7.vs..PLN_NT_d0",	"MLN_NT_d0.vs..PLN_NT_d0",	"MLN_NT_d7.vs..PLN_NT_d0",	"MLN_ST_d0.vs..PLN_NT_d0",	"MLN_ST_d7.vs..PLN_NT_d0",	"PLN_LT_d0.vs..PLN_NT_d0",	"PLN_LT_d7.vs..PLN_NT_d0",	"PLN_NT_d7.vs..PLN_NT_d0",	"PLN_ST_d0.vs..PLN_NT_d0",	"PLN_ST_d7.vs..PLN_NT_d0",	"SILP_LT_d0.vs..PLN_NT_d0",	"SILP_LT_d7.vs..PLN_NT_d0",	"SILP_NT_d0.vs..PLN_NT_d0",	"SILP_NT_d7.vs..PLN_NT_d0",	"SILP_ST_d0.vs..PLN_NT_d0",	"SILP_ST_d7.vs..PLN_NT_d0")

#samples=c("IEL_ST_d0.vs..IEL_NT_d0","IEL_LT_d0.vs..IEL_NT_d0","IEL_NT_d7.vs..IEL_NT_d0","IEL_ST_d7.vs..IEL_NT_d0","IEL_LT_d7.vs..IEL_NT_d0")
#samples=c("SILP_ST_d0.vs..SILP_NT_d0","SILP_LT_d0.vs..SILP_NT_d0","SILP_NT_d7.vs..SILP_NT_d0","SILP_ST_d7.vs..SILP_NT_d0","SILP_LT_d7.vs..SILP_NT_d0")
#samples=c("MLN_ST_d0.vs..MLN_NT_d0","MLN_LT_d0.vs..MLN_NT_d0","MLN_NT_d7.vs..MLN_NT_d0","MLN_ST_d7.vs..MLN_NT_d0","MLN_LT_d7.vs..MLN_NT_d0")
#samples=c("PLN_ST_d0.vs..PLN_NT_d0","PLN_LT_d0.vs..PLN_NT_d0","PLN_NT_d7.vs..PLN_NT_d0","PLN_ST_d7.vs..PLN_NT_d0","PLN_LT_d7.vs..PLN_NT_d0")

samples=c("IEL_CD4+CD8+_neg","IEL_CD8+CD90-_neg","MES_CD4+CD8+_neg","MES_CD8+_neg","PLN_CD8+_neg","SI_CD4+_neg","IEL_CD4+_neg","IEL_CD8+_neg","MES_CD4+_neg","PLN_CD4+_neg","SI_CD4+CD8+_neg","SI_CD8+_neg","IEL_CD4+CD8+_pos","IEL_CD8+CD90-_pos","MES_CD4+CD8+_pos","MES_CD8+_pos","PLN_CD8+_pos","SI_CD4+_pos","IEL_CD4+_pos","IEL_CD8+_pos","MES_CD4+_pos","PLN_CD4+_pos","SI_CD4+CD8+_pos","SI_CD8+_pos")

############################### Data extraction (for postive files)
# For each sample, open the GSEA report TSV and combine them in one variable 
for(i in samples)
{
  # set filename path
  filenameup=paste(directorypath,i,".tsv",sep="")
  #filenameup=paste(directorypath,i,"_neg.tsv",sep="")
  # set sample lables
  # uplistname=paste(i,"Up",sep="_")
  # read files
  onedataset <- read_delim(filenameup, "\t", escape_double = FALSE, trim_ws = TRUE)
  # sort by pvalues
  onedataset=onedataset[order(onedataset$`NOM p-val`),]
  # keep only the top 10 pathways
  onedataset=as.data.frame(head(onedataset,10))
  # print sample name
  print(i)
  # combine data in datalistentrez variable
  datalistentrez[[i]]=as.data.frame(onedataset)
}
datalistentrez[1]
############################### Extract data for ggplot
# Generate common unique lists
alllists=list()
for (j in seq(1,length(datalistentrez),1))
{
  # extract sample
  dataf=as.data.frame(datalistentrez[j])
  # extract pathway names
  newlist=dataf[,2]
  print(newlist)
  # combine pathway names from all samples
  alllists<-append(alllists,newlist)
}
#alllists<-list("T CELL ACTIVATION")
# generate unique combined pathway names list
alllists<-unlist(alllists)
alllists<-unique(alllists)
#alllists<-alllists[alllists != "T CELL DIFFERENTIATION IN THYMUS"]
#alllists<-alllists[alllists != "POSITIVE THYMIC T CELL SELECTION"]
#alllists<-alllists[alllists != "NEGATIVE THYMIC T CELL SELECTION"]
#alllists<-alllists[alllists != "ANTIGEN PROCESSING AND PRESENTATION"]
#alllists<-alllists[alllists != "HEMOPOIESIS"]

#View(as.data.frame(alllists))
length(alllists)

# Match and extract counts/pvalue for each of the pathways
# Initiate an empty table
updownpathwaystable = data.frame()
samplecount=as.numeric(1)

# For every unique pathway, extract its corresponding data from the enriched pathway list for each of the samples
# Read enrichment data for each list
for (k in seq(1,length(datalistentrez),1))
{
  # extract sample data
  onerowsamplename=names(datalistentrez[k])
  onedataframe=as.data.frame(datalistentrez[k])
  print(onedataframe)
  print(onerowsamplename)
  # Read unique pathway to search/match
  for(d in alllists)
  {
    print(d)
    # extract matching row
    oneline=onedataframe[grep(paste("^",d,"$",sep=""),onedataframe[,2]), ]
    onelinename=unlist(d)
    print(oneline)
    # extract matchin row count
    onelinecount=unlist(oneline[,4])
    print(onelinecount)
    # extract matching row pvalue
    onelinepvalue=unlist(oneline[,7])
    print(onelinepvalue)
    # If a pathway is enriched, extract its data
    if(length(onelinecount) >= 1)
    {
      fordataframe=c(onelinename,onelinecount,onelinepvalue,onerowsamplename)
      updownpathwaystable<-rbind(updownpathwaystable,fordataframe)
      print(fordataframe)
    }
    # If a pathway is not enriched, place 0 and 1 values
    else
    {
      fordataframe=c(onelinename,"0","1",onerowsamplename)
      updownpathwaystable<-rbind(updownpathwaystable,fordataframe)
    }
  }
}
#View(updownpathwaystable)
# Add column headers and remove 'KEGG' from pathway names
colnames(updownpathwaystable)=c("Term","Count","PValue","Sample")
#updownpathwaystable=separate(data=updownpathwaystable, col = Term, into=c(NA,"right"), sep="_", extra = "merge")
head(updownpathwaystable)
# Add column headers again
#colnames(updownpathwaystable)=c("Term","Count","PValue","Sample")

######################## Remove single entries from plot
updownpathwaystablenosingles = data.frame()

for(d in alllists)
{
  print(d)
  # extract matching row
  oneline=updownpathwaystable[grep(paste("^",d,"$",sep=""),updownpathwaystable[,1]), ]
  onelinename=unlist(d)
  #print(oneline)
  # extract matchin row count
  onelinecount=unlist(oneline[,2])
  #print(View(onelinecount))
  onezeros=onelinecount[grep("^0$",onelinecount)]
  onezeros
  print(onezeros)
  print(length(onezeros))
  # extract matching row pvalue
  #onelinepvalue=unlist(oneline[,7])
  #print(onelinepvalue)
  # If a pathway is enriched, extract its data
  # Change to 7 for 4 samples, 9 for 5 samples
  #if(length(onezeros) < 9)
  #{
    #fordataframe=c(onelinename,onelinecount,onelinepvalue,onerowsamplename)
    updownpathwaystablenosingles<-rbind(updownpathwaystablenosingles,oneline)
    print(oneline)
  #}
  # If a pathway is not enriched, place 0 and 1 values
  #else
  #{
  #  fordataframe=c(onelinename,"0","1",onerowsamplename)
  #  updownpathwaystable<-rbind(updownpathwaystable,fordataframe)
  #}
}
#View(updownpathwaystablenosingles)
#Convert to title case
updownpathwaystablenosingles$Term <- stri_trans_totitle(as.character(updownpathwaystablenosingles$Term))
# Order columns based on count
#updownpathwaystablenosingles$Term<-factor(updownpathwaystablenosingles$Term, levels=unique(updownpathwaystablenosingles$Term[order(as.numeric(updownpathwaystablenosingles$Count),decreasing = TRUE)]))
updownpathwaystablenosingles$Term<-factor(updownpathwaystablenosingles$Term, levels=unique(updownpathwaystablenosingles$Term[order(as.numeric(updownpathwaystablenosingles$Count))]))
#View(updownpathwaystablenosingles)


############################Plot
S1=updownpathwaystablenosingles %>% 
  ggplot(aes(x=factor(Sample, levels = samples), y=Term, size=as.numeric(Count), color=as.numeric(PValue), group=Sample)) + 
  geom_point(alpha = 0.6) + 
  theme_classic() +
  # pavlue range for color
  scale_color_gradient(low="red", high="blue", limit=c(0,1), name="PValue")+
  # count range for size
  scale_size(range = c(0, 4), name="Count")+
  #ggtitle("Negative")+
  #ggtitle("Positive")+
  theme(axis.text.x = element_text(size = 8, angle = 90),
    legend.title = element_text(size = 6), 
    legend.text = element_text(size = 6),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_rect(colour = "grey", fill=NA, size=1),
    panel.grid.major = element_line(size = (0.2), colour="grey"))
# Run S1 to plot
S1
# Save the above plot
#ggsave("gsea-IEL_NT_d0_vs_all-neg.png", width = 12, height = 10, units = "cm")

######################
