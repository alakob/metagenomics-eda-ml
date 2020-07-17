library('dplyr')
library(tidyr)
library('readr')

root.detail <- "https://www.ebi.ac.uk/metagenomics/projects/doExportDetails?search=Search&includingChildren=true&biomeLineage=root:"
root        <- "https://www.ebi.ac.uk/metagenomics/projects/doExport?search=Search&includingChildren=true&biomeLineage=root:"
biomes      <- c("Environmental:Terrestrial:Soil", "Engineered","Host-associated:Human",
            "Environmental:Aquatic:Marine", "Host-associated:Plants",
            "Host-associated:Human:Digestive%20system", "Environmental:Aquatic:Freshwater", "Host-associated:Mammals", 
            "Environmental:Terrestrial:Soil:Forest%20soil", "Environmental:Terrestrial:Soil:Grasslands")
##################################################################
# Extract all project metadata .. from MG portal
##################################################################
details <- tbl_df(read.csv(paste(root.detail, biomes[1], sep="")))
simple  <- tbl_df(read.csv(paste(root, biomes[1], sep="")))
simple <- simple %>% rename(Project.Name=Project.name)
biomes.meta <- left_join(simple, details, by=c("Project.Name"="Project.Name")) %>% 
                select(Biome, Project.Name, Study.ID) %>% 
                mutate(Study.ID=as.character(Study.ID))%>% unique() %>% na.omit()

#############################################################################
# READ IN THE PIPELILINES VERSION METAGENOME STUDIES AND METATRANSCRIPTOMICS
#############################################################################

setwd('/home/alakob/Desktop/MGMegaAnalysis')

wgs_2 <- tbl_df(read.delim(file="wgs_pipeline_2.txt", h=T))
wgs_3 <- tbl_df(read.delim(file="wgs_pipeline_3.txt", h=T))
# Are the same studies processed in pipeline 2 and pipeline 3
wgs_3.study<- wgs_3 %>% select(EXT_STUDY_ID) %>% unique()
wgs_2.study<- wgs_2 %>% select(EXT_STUDY_ID) %>% unique()


#################################################################
# Extract MG GO abundance contingency tables ...
#################################################################

getgo <- function(id=i, version=2){
  #id <- 'ERP006678'
  #version <- 2
  #output <- list()
  print(paste("Processing....", id, sep=" ") )
  prefix.down <- "https://www.ebi.ac.uk/metagenomics/projects/" 
  go_start <- "/download/" # append version.. 
  go_end  <- "/export?contentType=text&exportValue=GO_abundances"
  info.close <- "/overview/doExport"
  project <- paste(prefix.down, id , go_start, version,'.0', go_end, sep="")
  project
  project.info <- paste(prefix.down, id, info.close, sep="")
  project <- tbl_df(tryCatch(read.delim(file=project), error=function(e) NULL)) 
  colnames(project) <- gsub("_T", "",gsub("_G","",colnames(project)))
  project.info <- tbl_df(tryCatch(read.csv(file=project.info), error=function(e) NULL))
  if (dim(project.info)[1] <1) {
    print(paste (".... ", id, sep=""))
    output <-id
    return(output)
  }
  if (dim(project)[1] <1) {
    print(paste (".... ", id, sep=""))
    output <- id
    return(output)
  }
  project.info <-  project.info %>% 
                   filter(Release.version==version) %>%
                   select(Run.ID,Experiment.type,Sample.ID, Release.version) 
  project.info <- project.info %>%  
                  mutate(Study.ID=id, Sample.ID=as.character(Sample.ID), Run.ID=as.character(Run.ID))
  project <- project %>% 
             select(-GO) %>%gather(Run.ID, count, -description, -category)
  project.df <- left_join(project, project.info, by="Run.ID")
  project.df <- left_join(project.df, biomes.meta, by="Study.ID")
  project.df <- project.df %>% 
                select(Run.ID, Experiment.type, Release.version, Study.ID ,description, category,count,Biome)
  output <- project.df
  return(output)
}


########################################################
# Extract GO annotation pipeline version 2
########################################################

wgsgov2 <- lapply(as.list(wgs_2.study$EXT_STUDY_ID), function (x) getgo(id=x, version=2))
missingAnnoV2 <- wgsgov2[which(unlist(lapply(wgsgov2, function(x) is.factor(x))))]
project_wo_annoV2 <- as.character(droplevels(unlist(missingAnnoV2)))
wgsgov2_anno <- wgsgov2[-which(unlist(lapply(wgsgov2, function(x) is.factor(x))))]


goV2 <- wgsgov2_anno %>% bind_rows() %>%
  group_by(description,category,Biome) %>% mutate(count=sum(count)) %>%
  ungroup() %>% unique()

system.time(bpV2  <- goV2 %>% filter(category=="biological process") %>% unique() %>% spread(Biome, count, fill=0) %>% select(-category, -Run.ID, -Experiment.type, -Release.version,  -Study.ID) %>% unique())
# user  system elapsed 
# 537.776  10.324 548.327 
system.time(ccpV2 <- goV2 %>% filter(category=="cellular component") %>% unique() %>% spread(Biome, count, fill=0) %>% select(-category, -Run.ID, -Experiment.type, -Release.version,  -Study.ID) %>% unique())
system.time(mfV2  <- goV2 %>% filter(category=="molecular function") %>% unique() %>% spread(Biome, count, fill=0) %>% select(-category, -Run.ID, -Experiment.type, -Release.version,  -Study.ID) %>% unique())


########################################################
# Extract GO annotation pipeline version 3
########################################################

wgsgov3 <- lapply(as.list(wgs_3.study$EXT_STUDY_ID), function (x) getgo(id=x, version=3))
missingAnnoV3 <- wgsgov3[which(unlist(lapply(wgsgov3, function(x) is.factor(x))))]
project_wo_annoV3 <- as.character(droplevels(unlist(missingAnnoV3)))
wgsgov3_anno <- wgsgov3[-which(unlist(lapply(wgsgov3, function(x) is.factor(x))))]


goV3 <- wgsgov3_anno %>% bind_rows() %>%
  group_by(description,category,Biome) %>% mutate(count=sum(count)) %>%
  ungroup() %>% unique()

system.time(bpV3  <- goV3 %>% filter(category=="biological process") %>% unique() %>% spread(Biome, count, fill=0) %>% select(-category, -Run.ID, -Experiment.type, -Release.version,  -Study.ID) %>% unique())
# user  system elapsed 
# 537.776  10.324 548.327 
system.time(ccpV3 <- goV3 %>% filter(category=="cellular component") %>% unique() %>% spread(Biome, count, fill=0) %>% select(-category, -Run.ID, -Experiment.type, -Release.version,  -Study.ID) %>% unique())
system.time(mfV3  <- goV3 %>% filter(category=="molecular function") %>% unique() %>% spread(Biome, count, fill=0) %>% select(-category, -Run.ID, -Experiment.type, -Release.version,  -Study.ID) %>% unique())


goanalysis <- list(wgsgov2=wgsgov2, missingAnnoV2=missingAnnoV2, project_wo_annoV2=project_wo_annoV2, wgsgov2_anno=wgsgov2_anno, goV2 =goV2 , bpV2=bpV2, ccpV2=ccpV2, mfV2=mfV2,
                    wgsgov3=wgsgov3, missingAnnoV3=missingAnnoV3, project_wo_annoV3=project_wo_annoV3, wgsgov3_anno=wgsgov3_anno, goV3 =goV3 , bpV3=bpV3, ccpV3=ccpV3, mfV3=mfV3
                   )
save(goanalysis, file="goAnalysis.Rdata")

####################################################
# Biological processess
#####################################################

library('gplots')
library(reshape2)
library("Hmisc")
library('corrplot')
suppressWarnings(suppressMessages(require(FactoMineR)))   # Multivariate analysis
suppressWarnings(suppressMessages(require(factoextra)))   # Visualize result of
suppressWarnings(suppressMessages(require(XML)))
suppressWarnings(suppressMessages(require(RCurl)))

norm_vec <- function(x) x/sqrt(sum(x^2))

bpV2.mat  <-  bpV2 %>% select(-description) %>% as.matrix()
rownames(bpV2.mat) <- bpV2$description
bpV2.norm <- apply(bpV2.mat, MARGIN = 2, FUN = function(X) (norm_vec(X)))

#######################################################
# Metagenome, Family-rank Vs Biomes 
#######################################################
rc <- rainbow(nrow(bpV2.norm), start=0, end=.3)
cc <- rainbow(ncol(bpV2.norm), start=0, end=.3)
bpV2.map <- heatmap.2(bpV2.norm, col=cm.colors(255), scale="row",
                    RowSideColors=rc, ColSideColors=cc, margin=c(5, 10),
                    xlab="Biomes", ylab= "BP",
                    main="Biological processes vs Biomes",
                    tracecol="steelblue", density="density",
                    srtRow=55, adjRow=c(0, 1), srtCol=10, adjCol=c(1,1)
)


#####################################################
# CA metatrans....
#####################################################
bp.ca <- CA(bp.mat, graph=TRUE, ncp=10)
bp.hcpc <- HCPC(bp.ca, cluster.CA = "columns")

fviz_dend(bp.hcpc, horiz=T, rec_fill=TRUE, rect = TRUE, type="rectangle", k=5, sub=NULL, main="Biological process \nBiomes Clustering") + 
  theme( axis.line.x  = element_line(colour = "grey", size = .3, linetype = "solid"))

bp.scree <- fviz_screeplot(metatrans.ca, barfill='white', addlabels=TRUE) +ggtitle("Family,biomes")+ theme_bw()



###############################################
# BEst contributing .... Metagenomes
###############################################
bp.ctrb  <- fviz_ca_row(bp.ca, alpha.row="contrib", select.row=list(contrib=25))
best.ctrb <- as.character(bp.ctrb$data[,"name"])
bp.ctrb <- bp.mat[as.character(rownames(bp.mat)) %in% best.ctrb, ]
####################################################################################
# Get rid of Biomes that does not does not contain highly contributing families...
####################################################################################
bp.myselect <- apply(bp.ctrb, 2, function(x) all(x==0))
bp.ctrb<- bp.ctrb[,!bp.myselect]
bp.corr <- rcorr((bp.ctrb))

##########################################
# Plot correlation between biomes...
##########################################
corrplot(bp.corr$r, type="lower", order="FPC", method='pie', main="Families correlation",
         p.mat = bp.corr$P, sig.level = 0.05, insig = "blank", tl.col = "black", tl.srt = 45)

###########################################
# CA with the most contributing families
###########################################
bp.best <- CA(bp.ctrb, graph=TRUE, ncp=10)
bp.best.hcpc <- HCPC(bp.best, cluster.CA = "columns")

fviz_dend(bp.best.hcpc, horiz=T, rec_fill=TRUE, rect = TRUE, type="rectangle", k=5, sub=NULL, main="Biomes Clustering on BP") + 
  theme( axis.line.x  = element_line(colour = "grey", size = .3, linetype = "solid"))


##############################################
#  Heatmap of best contributing families. ....
##############################################

bp.ctrb.norm <- apply(bp.ctrb[,-51], MARGIN = 1, FUN = function(X) (norm_vec(X)))

#######################################################
# Biological processess, BP Vs Biomes 
#######################################################
rc <- rainbow(nrow(bp.ctrb.norm), start=0, end=.3)
cc <- rainbow(ncol(bp.ctrb.norm), start=0, end=.3)
bp.map <- heatmap.2(bp.ctrb.norm, col=cm.colors(255), scale="row",
                    RowSideColors=rc, ColSideColors=cc, margin=c(5, 10),
                    ylab="Biomes", xlab= "Biological processes",
                    main="Biological Process vs Biomes",
                    tracecol="steelblue", density="density",
                    srtRow=55, adjRow=c(0, 1), srtCol=10, adjCol=c(1,1)
)



#########################################################
# Molecular functions ...
#########################################################
mf.mat  <-  mf %>% select(-description) %>% as.matrix()
rownames(mf.mat) <- mf$description
mf.norm <- apply(mf.mat[,-51], MARGIN = 1, FUN = function(X) (norm_vec(X)))

#######################################################
# Metagenome, Family-rank Vs Biomes 
#######################################################
rc <- rainbow(nrow(mf.norm), start=0, end=.3)
cc <- rainbow(ncol(mf.norm), start=0, end=.3)
mf.map <- heatmap.2(mf.norm, col=cm.colors(255), scale="row",
                    RowSideColors=rc, ColSideColors=cc, margin=c(5, 10),
                    xlab="Biomes", ylab= "Families",
                    main="Biomes vs Families",
                    tracecol="steelblue", density="density",
                    srtRow=55, adjRow=c(0, 1), srtCol=10, adjCol=c(1,1)
)

########################################################
# Extract GO annotation pipeline version 3
########################################################
wgsgov3 <- lapply(as.list(wgs_3.study$EXT_STUDY_ID), function (x) getgo(id=x, version=3))
missingAnnoV3 <- wgsgov3[which(unlist(lapply(wgsgov3, function(x) is.factor(x))))]
project_wo_annoV3 <- as.character(droplevels(unlist(missingAnnoV3)))
wgsgov3_anno <- wgsgov3[-which(unlist(lapply(wgsgov3, function(x) is.factor(x))))]


go <- go_dat %>% bind_rows() %>%
  group_by(description,category,Biome) %>% mutate(count=sum(count)) %>%
  ungroup() %>% unique()



bp<- go %>% filter(category=="biological process") %>% unique() %>% spread(Biome, count, fill=0) %>% select(-category)
ccp <- go %>% filter(category=="cellular component") %>% unique() %>% spread(Biome, count, fill=0) %>% select(-category)
mf<- go  %>% filter(category=="molecular function") %>% unique() %>% spread(Biome, count, fill=0) %>% select(-category)


table(wgs_3.study$EXT_STUDY_ID %in% wgs_2.study$EXT_STUDY_ID)
table(wgs_2.study$EXT_STUDY_ID %in% wgs_3.study$EXT_STUDY_ID)
wgs_2.study$EXT_STUDY_ID[(wgs_2.study$EXT_STUDY_ID %in% wgs_3.study$EXT_STUDY_ID)]
wgs_3.study$EXT_STUDY_ID[(wgs_3.study$EXT_STUDY_ID %in% wgs_2.study$EXT_STUDY_ID)]

# The following studies are the only one processed by pipeline 2 and pipeline 3
# Out of 91 distinct WGS study_id processes by pipeline 2 and 160 WGS study_id processed by pipeline 3 ,Only 3 study have been processed by both pipeline.
#(ERP001736 ERP015773 ERP016063)
# Regarding your question of comparing GO annotation by both pipeline, we have very sparse data here to work with. 
#I expected to find the annotation of both pipeline 2 and pipeline 3 for a given study, Am i missing something?

############################################
# METATRANSCRIPTOMICS
############################################
mt_2 <- tbl_df(read.delim(file="mt_pipeline_2.txt", h=T))
mt_3 <- tbl_df(read.delim(file="mt_pipeline_3.txt", h=T))
# Are the same studies processed in pipeline 2 and pipeline 3
mt_3.study<- mt_3 %>% select(EXT_STUDY_ID) %>% unique()
mt_2.study<- mt_2 %>% select(EXT_STUDY_ID) %>% unique()

table(mt_3.study$EXT_STUDY_ID %in% mt_2.study$EXT_STUDY_ID)
table(mt_2.study$EXT_STUDY_ID %in% mt_3.study$EXT_STUDY_ID)
mt_2.study$EXT_STUDY_ID[(mt_2.study$EXT_STUDY_ID %in% mt_3.study$EXT_STUDY_ID)]
mt_3.study$EXT_STUDY_ID[(mt_3.study$EXT_STUDY_ID %in% mt_2.study$EXT_STUDY_ID)]



#################################################################
# Extract MG taxonomy abundance contingency tables ...
#################################################################
prefix.down <- "https://www.ebi.ac.uk/metagenomics/projects/" 
suffix.down  <- "/download/2.0/export?contentType=text&exportValue=taxonomy_abundances"
meta.pref <- "https://www.ebi.ac.uk/metagenomics/projects/"
meta.suf <- "/overview/doExport"
#contingencies.bk <- contingencies

mygather <- function(df){
  colnames(df)[1] <- "taxonomy"
  df %>% rename(Taxa=taxonomy) %>%
    gather(Run,count, -Taxa,-Study.ID) %>%
    mutate(Taxa=gsub(".*;f__","", Taxa)) %>% 
    mutate(Taxa=gsub(";.*","", Taxa)) %>% 
    group_by(Taxa, Run) %>% 
    mutate(count=sum(count)) %>% ungroup() %>%
    unique() %>% filter(!Taxa=="") # %>% select(-Study.ID, -Run) %>% group_by(Taxa, Biome) %>% mutate(count=sum(count)) %>% ungroup()
}

contingencies <- list()

system.time(for (i in biomes.meta$Study.ID){
  #i <- "ERP012803"
  print(paste("Processing....", i, sep=" ") )
  project <- paste(prefix.down, i, suffix.down, sep="")
  #project.meta <- paste(meta.pref, i, meta.suf, sep="")
  #project <- tbl_df(tryCatch(read_tsv(file=project), error=function(e) NULL))
  project <- tbl_df(tryCatch(read.delim(file=project), error=function(e) NULL))
  #project.meta <- tbl_df(tryCatch(read.csv(file=project.meta), error=function(e) NULL))
  if(nrow(project)> 1){
    project <- project %>% mutate(Study.ID=i)  %>% mygather()
    project <- left_join(project, biomes.meta, by="Study.ID") %>% 
      select(-Samples, -Run)  %>% 
      group_by(Taxa,Study.ID, Biome) %>% mutate(count=sum(count)) %>% 
      ungroup() %>% unique() %>% filter(!Taxa=="Root")
    contingencies[[i]] <- project 
    
  }
}
)
contingencies
composite <- contingencies %>% bind_rows() %>% group_by(Taxa,Biome,Project.Name,Study.ID) %>% mutate(count=sum(count)) %>% ungroup() %>% unique() %>% spread(Biome, count, fill=0)


contingencies2 <- contingencies %>% bind_rows() %>% group_by(Taxa,Biome,Project.Name,Study.ID) %>% mutate(count=sum(count)) %>% ungroup() %>% unique()


################################################
# Metagenomes projects
################################################
metagenome <-  contingencies2 %>% filter(!grepl("transcript|rna |amplicon",Project.Name)) %>%
              select(-Study.ID,-Project.Name) %>% 
              group_by(Taxa,Biome) %>% mutate(count=sum(count)) %>%
              ungroup() %>% unique() %>% 
              spread(Biome, count, fill=0)
metagenome.mat<- metagenome %>% select(-Taxa) %>% as.matrix()
rownames(metagenome.mat) <- metagenome$Taxa

################################################
# Metatranscriptomic, Amplicon projects
################################################
metatrans <-  contingencies2 %>% filter(grepl("transcript|rna |amplicon",Project.Name)) %>%
  select(-Study.ID,-Project.Name) %>% 
  group_by(Taxa,Biome) %>% mutate(count=sum(count)) %>%
  ungroup() %>% unique() %>% 
  spread(Biome, count, fill=0)
metatrans.mat<- metatrans %>% select(-Taxa) %>% as.matrix()
rownames(metatrans.mat) <- metatrans$Taxa


#################################################
# Calculate norm of vector for normalization....
#################################################

library('gplots')
library(reshape2)
library("Hmisc")
library('corrplot')
suppressWarnings(suppressMessages(require(FactoMineR)))   # Multivariate analysis
suppressWarnings(suppressMessages(require(factoextra)))   # Visualize result of
suppressWarnings(suppressMessages(require(XML)))
suppressWarnings(suppressMessages(require(RCurl)))

norm_vec <- function(x) x/sqrt(sum(x^2))

#######################################################
## HEATMAP normalized pangenome captured pan-genome ..
#                 W/O TARA pan-genome
#######################################################
#comp.mat <- composite %>% select(-Taxa) %>% as.matrix()
#rownames(comp.mat) <- composite$Taxa
metagenome.norm <- apply(metagenome.mat, MARGIN = 2, FUN = function(X) (norm_vec(X)))
metatrans.norm <- apply(metatrans.mat, MARGIN = 2, FUN = function(X) (norm_vec(X)))

#######################################################
# Metagenome, Family-rank Vs Biomes 
#######################################################
rc <- rainbow(nrow(metagenome.norm), start=0, end=.3)
cc <- rainbow(ncol(metagenome.norm), start=0, end=.3)
metagenome.map <- heatmap.2(metagenome.norm, col=cm.colors(255), scale="column",
                       RowSideColors=rc, ColSideColors=cc, margin=c(5, 10),
                       xlab="Biomes", ylab= "Families",
                       main="Biomes vs Families",
                       tracecol="steelblue", density="density",
                       srtRow=55, adjRow=c(0, 1), srtCol=10, adjCol=c(1,1)
)

#######################################################
# Metatrans, Family-rank Vs Biomes 
#######################################################

rc <- rainbow(nrow(metatrans.norm), start=0, end=.3)
cc <- rainbow(ncol(metatrans.norm), start=0, end=.3)
metatrans.map <- heatmap.2(metatrans.norm, col=cm.colors(255), scale="column",
                            RowSideColors=rc, ColSideColors=cc, margin=c(5, 10),
                            xlab="Biomes", ylab= "Families",
                            main="Biomes vs Families",
                            tracecol="steelblue", density="density",
                            srtRow=55, adjRow=c(0, 1), srtCol=10, adjCol=c(1,1)
)

#####################################################
# CA metagenomes....
#####################################################
metagenome.ca <- CA(metagenome.mat, graph=TRUE, ncp=10)
metagenome.hcpc <- HCPC(metagenome.ca, cluster.CA = "columns")

fviz_dend(metagenome.hcpc, horiz=T, rec_fill=TRUE, rect = TRUE, type="rectangle", k=5, sub=NULL, main="Biomes Clustering") + 
  theme( axis.line.x  = element_line(colour = "grey", size = .3, linetype = "solid"))

metagenome.scree<- fviz_screeplot(metagenome.ca, barfill='white', addlabels=TRUE) +ggtitle("Family,biomes")+ theme_bw()

#####################################################
# CA metatrans....
#####################################################
metatrans.ca <- CA(metatrans.mat, graph=TRUE, ncp=10)
metatrans.hcpc <- HCPC(metatrans.ca, cluster.CA = "columns")

fviz_dend(metatrans.hcpc, horiz=T, rec_fill=TRUE, rect = TRUE, type="rectangle", k=5, sub=NULL, main="Biomes Clustering") + 
  theme( axis.line.x  = element_line(colour = "grey", size = .3, linetype = "solid"))

metatrans.scree<- fviz_screeplot(metatrans.ca, barfill='white', addlabels=TRUE) +ggtitle("Family,biomes")+ theme_bw()



###############################################
# BEst contributing .... Metagenomes
###############################################
metagenome.ctrb  <- fviz_ca_row(metagenome.ca, alpha.row="contrib", select.row=list(contrib=50))
best.ctrb <- as.character(metagenome.ctrb$data[,"name"])
metagenome.ctrb <- metagenome.mat[as.character(rownames(metagenome.mat)) %in% best.ctrb, ]
####################################################################################
# Get rid of Biomes that does not does not contain highly contributing families...
####################################################################################
myselect <- apply(metagenome.ctrb, 2, function(x) all(x==0))
metagenome.ctrb<- metagenome.ctrb[,!myselect]
metagenome.corr <- rcorr((metagenome.ctrb))

##########################################
# Plot correlation between biomes...
##########################################
corrplot(metagenome.corr$r, type="lower", order="FPC", method='pie', main="Families correlation",
         p.mat = metagenome.corr$P, sig.level = 0.05, insig = "blank", tl.col = "black", tl.srt = 45)

###########################################
# CA with the most contributing families
###########################################
metagenome.best <- CA(metagenome.ctrb, graph=TRUE, ncp=10)
metagenome.best.hcpc <- HCPC(metagenome.best, cluster.CA = "columns")

fviz_dend(metagenome.best.hcpc, horiz=T, rec_fill=TRUE, rect = TRUE, type="rectangle", k=5, sub=NULL, main="Biomes Clustering") + 
  theme( axis.line.x  = element_line(colour = "grey", size = .3, linetype = "solid"))


##############################################
#  Heatmap of best contributing families. ....
##############################################

metagenome.ctrb.norm <- apply(metagenome.ctrb, MARGIN = 2, FUN = function(X) (norm_vec(X)))

#######################################################
# Metagenome, Family-rank Vs Biomes 
#######################################################
rc <- rainbow(nrow(metagenome.ctrb.norm), start=0, end=.3)
cc <- rainbow(ncol(metagenome.ctrb.norm), start=0, end=.3)
metagenome.map <- heatmap.2(metagenome.ctrb.norm, col=cm.colors(255), scale="column",
                            RowSideColors=rc, ColSideColors=cc, margin=c(5, 10),
                            xlab="Biomes", ylab= "Families",
                            main="Biomes vs Families",
                            tracecol="steelblue", density="density",
                            srtRow=55, adjRow=c(0, 1), srtCol=10, adjCol=c(1,1)
)



#################################################################
# Extract MG GO abundance contingency tables ...
#################################################################
prefix.down <- "https://www.ebi.ac.uk/metagenomics/projects/" 
goslim  <- "/download/2.0/export?contentType=text&exportValue=GO-slim_abundances"
ipr <- "/download/2.0/export?contentType=text&exportValue=IPR_abundances"
go <- "/download/2.0/export?contentType=text&exportValue=GO_abundances"

info.close <- "/overview/doExport"
go_dat <- list()
no_go_no_meta <- list()
no_go_meta <- list()
go_no_meta <- list()


mygather <- function(df){
  colnames(df)[1] <- "taxonomy"
  df %>% rename(Taxa=taxonomy) %>%
    gather(Run,count, -Taxa,-Study.ID) %>%
    mutate(Taxa=gsub(".*;f__","", Taxa)) %>% 
    mutate(Taxa=gsub(";.*","", Taxa)) %>% 
    group_by(Taxa, Run) %>% 
    mutate(count=sum(count)) %>% ungroup() %>%
    unique() %>% filter(!Taxa=="") # %>% select(-Study.ID, -Run) %>% group_by(Taxa, Biome) %>% mutate(count=sum(count)) %>% ungroup()
}

system.time(for (i in biomes.meta$Study.ID){
  i <- 'SRP041174'
  print(paste("Processing....", i, sep=" ") )
  project <- paste(prefix.down, i , go, sep="")
  project.info <- paste(prefix.down, i, info.close, sep="")
  project <- tbl_df(tryCatch(read.delim(file=project), error=function(e) NULL)) # %>% select(-GO,-category)
  colnames(project) <- gsub("_T", "",gsub("_G","",colnames(project)))
  
  project.info <- tbl_df(tryCatch(read.csv(file=project.info), error=function(e) NULL))
  
  if(nrow(project)> 1 & nrow(project.info) > 1){
    print(paste(":) Project....", i, "has GO and meta info....", sep=" "))
    project.info <-  project.info %>% filter(Release.version==2) %>% select(Run.ID,Experiment.type,Sample.ID, Release.version) # %>% select(-GO,-category)
    project.info <- project.info %>% mutate(Study.ID=i, Sample.ID=as.character(Sample.ID), Run.ID=as.character(Run.ID)) #%>% filter(!grepl("transcript|rna |amplicon",Experiment.type))
    project <- project %>% select(-GO) %>%gather(Run.ID, count, -description, -category)
    project.df <- left_join(project, project.info, by="Run.ID")
    project.df <- left_join(project.df, biomes.meta, by="Study.ID")
    project.df <- project.df %>% select(description, category,count,Biome)
    go_dat[[i]] <- project.df 
  }
  else if(nrow(project)> 1 & nrow(project.info) < 1)  {
    print(paste("Project....", i, "Missing meta info....", sep=" "))
    go_no_meta[[i]] <- project
  }
  else if(nrow(project)< 1 & nrow(project.info) < 1)  {
    print(paste("Project....", i, "has no GO and no meta info....", sep=" "))
    no_go_no_meta[[i]] <- "NA"
  }
  else if(nrow(project)< 1 & nrow(project.info) > 1)  {
    print(paste("Project....", i, "has no GO and has meta info....", sep=" "))
    no_go_meta[[i]] <- project.info
  }
}
)


 go <- go_dat %>% bind_rows() %>%
  group_by(description,category,Biome) %>% mutate(count=sum(count)) %>%
  ungroup() %>% unique()



bp<- go %>% filter(category=="biological process") %>% unique() %>% spread(Biome, count, fill=0) %>% select(-category)
ccp <- go %>% filter(category=="cellular component") %>% unique() %>% spread(Biome, count, fill=0) %>% select(-category)
mf<- go  %>% filter(category=="molecular function") %>% unique() %>% spread(Biome, count, fill=0) %>% select(-category)

####################################################
# Biological processess
#####################################################
bp.mat  <-  bp %>% select(-description) %>% as.matrix()
rownames(bp.mat) <- bp$description
bp.norm <- apply(bp.mat, MARGIN = 2, FUN = function(X) (norm_vec(X)))

#######################################################
# Metagenome, Family-rank Vs Biomes 
#######################################################
rc <- rainbow(nrow(bp.norm), start=0, end=.3)
cc <- rainbow(ncol(bp.norm), start=0, end=.3)
bp.map <- heatmap.2(bp.norm, col=cm.colors(255), scale="none",
                            RowSideColors=rc, ColSideColors=cc, margin=c(5, 10),
                            xlab="Biomes", ylab= "BP",
                            main="Biological processes vs Biomes",
                            tracecol="steelblue", density="density",
                            srtRow=55, adjRow=c(0, 1), srtCol=10, adjCol=c(1,1)
)


#####################################################
# CA metatrans....
#####################################################
bp.ca <- CA(bp.mat, graph=TRUE, ncp=10)
bp.hcpc <- HCPC(bp.ca, cluster.CA = "columns")

fviz_dend(bp.hcpc, horiz=T, rec_fill=TRUE, rect = TRUE, type="rectangle", k=5, sub=NULL, main="Biological process \nBiomes Clustering") + 
  theme( axis.line.x  = element_line(colour = "grey", size = .3, linetype = "solid"))

bp.scree <- fviz_screeplot(metatrans.ca, barfill='white', addlabels=TRUE) +ggtitle("Family,biomes")+ theme_bw()



###############################################
# BEst contributing .... Metagenomes
###############################################
bp.ctrb  <- fviz_ca_row(bp.ca, alpha.row="contrib", select.row=list(contrib=25))
best.ctrb <- as.character(bp.ctrb$data[,"name"])
bp.ctrb <- bp.mat[as.character(rownames(bp.mat)) %in% best.ctrb, ]
####################################################################################
# Get rid of Biomes that does not does not contain highly contributing families...
####################################################################################
bp.myselect <- apply(bp.ctrb, 2, function(x) all(x==0))
bp.ctrb<- bp.ctrb[,!bp.myselect]
bp.corr <- rcorr((bp.ctrb))

##########################################
# Plot correlation between biomes...
##########################################
corrplot(bp.corr$r, type="lower", order="FPC", method='pie', main="Families correlation",
         p.mat = bp.corr$P, sig.level = 0.05, insig = "blank", tl.col = "black", tl.srt = 45)

###########################################
# CA with the most contributing families
###########################################
bp.best <- CA(bp.ctrb, graph=TRUE, ncp=10)
bp.best.hcpc <- HCPC(bp.best, cluster.CA = "columns")

fviz_dend(bp.best.hcpc, horiz=T, rec_fill=TRUE, rect = TRUE, type="rectangle", k=5, sub=NULL, main="Biomes Clustering on BP") + 
  theme( axis.line.x  = element_line(colour = "grey", size = .3, linetype = "solid"))


##############################################
#  Heatmap of best contributing families. ....
##############################################

bp.ctrb.norm <- apply(bp.ctrb[,-51], MARGIN = 1, FUN = function(X) (norm_vec(X)))

#######################################################
# Biological processess, BP Vs Biomes 
#######################################################
rc <- rainbow(nrow(bp.ctrb.norm), start=0, end=.3)
cc <- rainbow(ncol(bp.ctrb.norm), start=0, end=.3)
bp.map <- heatmap.2(bp.ctrb.norm, col=cm.colors(255), scale="row",
                            RowSideColors=rc, ColSideColors=cc, margin=c(5, 10),
                            ylab="Biomes", xlab= "Biological processes",
                            main="Biological Process vs Biomes",
                            tracecol="steelblue", density="density",
                            srtRow=55, adjRow=c(0, 1), srtCol=10, adjCol=c(1,1)
)



#########################################################
# Molecular functions ...
#########################################################
mf.mat  <-  mf %>% select(-description) %>% as.matrix()
rownames(mf.mat) <- mf$description
mf.norm <- apply(mf.mat[,-51], MARGIN = 1, FUN = function(X) (norm_vec(X)))

#######################################################
# Metagenome, Family-rank Vs Biomes 
#######################################################
rc <- rainbow(nrow(mf.norm), start=0, end=.3)
cc <- rainbow(ncol(mf.norm), start=0, end=.3)
mf.map <- heatmap.2(mf.norm, col=cm.colors(255), scale="row",
                    RowSideColors=rc, ColSideColors=cc, margin=c(5, 10),
                    xlab="Biomes", ylab= "Families",
                    main="Biomes vs Families",
                    tracecol="steelblue", density="density",
                    srtRow=55, adjRow=c(0, 1), srtCol=10, adjCol=c(1,1)
)


#####################################################
# CA metatrans....
#####################################################
mf.ca <- CA(mf.mat[,-51], graph=TRUE, ncp=10)
mf.hcpc <- HCPC(mf.ca, cluster.CA = "columns")

fviz_dend(mf.hcpc, horiz=T, rec_fill=TRUE, rect = TRUE, type="rectangle", k=5, sub=NULL, main="Biological process \nBiomes Clustering") + 
  theme( axis.line.x  = element_line(colour = "grey", size = .3, linetype = "solid"))

mf.scree <- fviz_screeplot(mf.ca, barfill='white', addlabels=TRUE) +ggtitle("Family,biomes")+ theme_bw()



###############################################
# BEst contributing .... Metagenomes
###############################################
mf.ctrb  <- fviz_ca_row(mf.ca, alpha.row="contrib", select.row=list(contrib=15))
best.ctrb <- as.character(mf.ctrb$data[,"name"])
mf.ctrb <- mf.mat[as.character(rownames(mf.mat)) %in% best.ctrb, ]
####################################################################################
# Get rid of Biomes that does not does not contain highly contributing families...
####################################################################################
mf.myselect <- apply(mf.ctrb, 2, function(x) all(x==0))
mf.ctrb<- mf.ctrb[,!mf.myselect]
mf.corr <- rcorr((mf.ctrb))

##########################################
# Plot correlation between biomes...
##########################################
corrplot(mf.corr$r, type="lower", order="FPC", method='pie', main="Families correlation",
         p.mat = mf.corr$P, sig.level = 0.05, insig = "blank", tl.col = "black", tl.srt = 45)

###########################################
# CA with the most contributing families
###########################################
mf.best <- CA(mf.ctrb, graph=TRUE, ncp=10)
mf.best.hcpc <- HCPC(mf.best, cluster.CA = "columns")

fviz_dend(mf.best.hcpc, horiz=T, rec_fill=TRUE, rect = TRUE, type="rectangle", k=5, sub=NULL, main="Biomes Clustering") + 
  theme( axis.line.x  = element_line(colour = "grey", size = .3, linetype = "solid"))


##############################################
#  Heatmap of best contributing Molecular functions ....
##############################################

mf.ctrb.norm <- apply(mf.ctrb, MARGIN = 2, FUN = function(X) (norm_vec(X)))

#######################################################
# Molecular functions Vs Biomes 
#######################################################
rc <- rainbow(nrow(mf.ctrb.norm), start=0, end=.3)
cc <- rainbow(ncol(mf.ctrb.norm), start=0, end=.3)
mf.map <- heatmap.2(mf.ctrb.norm, col=cm.colors(255), scale="column",
                    RowSideColors=rc, ColSideColors=cc, margin=c(5, 10),
                    xlab="Biomes", ylab= "Families",
                    main="Biomes vs Molecular Functions",
                    tracecol="steelblue", density="density",
                    srtRow=55, adjRow=c(0, 1), srtCol=10, adjCol=c(1,1)
)







####################################################
# Cellular Compartment
#####################################################
cc.mat  <-  ccp %>% select(-description) %>% as.matrix()
rownames(cc.mat) <- ccp$description
cc.norm <- apply(cc.mat, MARGIN = 2, FUN = function(X) (norm_vec(X)))

#######################################################
# Metagenome, Family-rank Vs Biomes 
#######################################################
rc <- rainbow(nrow(cc.norm), start=0, end=.3)
cc <- rainbow(ncol(cc.norm), start=0, end=.3)
cc.map <- heatmap.2(cc.norm, col=cm.colors(255), scale="column",
                    RowSideColors=rc, ColSideColors=cc, margin=c(5, 10),
                    xlab="Biomes", ylab= "Families",
                    main="Biomes vs Families",
                    tracecol="steelblue", density="density",
                    srtRow=55, adjRow=c(0, 1), srtCol=10, adjCol=c(1,1)
)


#####################################################
# CA metatrans....
#####################################################
cc.ca <- CA(cc.mat, graph=TRUE, ncp=10)
cc.hcpc <- HCPC(cc.ca, cluster.CA = "columns")

fviz_dend(cc.hcpc, horiz=T, rec_fill=TRUE, rect = TRUE, type="rectangle", k=5, sub=NULL, main="Cellular Compartment \nBiomes Clustering") + 
  theme( axis.line.x  = element_line(colour = "grey", size = .3, linetype = "solid"))

mf.scree <- fviz_screeplot(cc.ca, barfill='white', addlabels=TRUE) +ggtitle("Family,biomes")+ theme_bw()



###############################################
# BEst contributing .... Metagenomes
###############################################
cc.ctrb  <- fviz_ca_row(cc.ca, alpha.row="contrib", select.row=list(contrib=15))
best.ctrb <- as.character(cc.ctrb$data[,"name"])
cc.ctrb <- cc.mat[as.character(rownames(cc.mat)) %in% best.ctrb, ]
####################################################################################
# Get rid of Biomes that does not does not contain highly contributing families...
####################################################################################
cc.myselect <- apply(cc.ctrb, 2, function(x) all(x==0))
cc.ctrb<- cc.ctrb[,!cc.myselect]
cc.corr <- rcorr((cc.ctrb))

##########################################
# Plot correlation between biomes...
##########################################
corrplot(cc.corr$r, type="lower", order="FPC", method='pie', main="Families correlation",
         p.mat = cc.corr$P, sig.level = 0.05, insig = "blank", tl.col = "black", tl.srt = 45)

###########################################
# CA with the most contributing families
###########################################
cc.best <- CA(cc.ctrb, graph=TRUE, ncp=10)
cc.best.hcpc <- HCPC(cc.best, cluster.CA = "columns")

fviz_dend(cc.best.hcpc, horiz=T, rec_fill=TRUE, rect = TRUE, type="rectangle", k=5, sub=NULL, main="Biomes Clustering") + 
  theme( axis.line.x  = element_line(colour = "grey", size = .3, linetype = "solid"))


########################################################
#  Heatmap of best contributing Molecular functions ....
########################################################

cc.ctrb.norm <- apply(cc.ctrb, MARGIN = 2, FUN = function(X) (norm_vec(X)))

########################################################
# Molecular functions Vs Biomes 
########################################################
rc <- rainbow(nrow(cc.ctrb.norm), start=0, end=.3)
cc <- rainbow(ncol(cc.ctrb.norm), start=0, end=.3)
cc.map <- heatmap.2(cc.ctrb.norm, col=cm.colors(255), scale="column",
                    RowSideColors=rc, ColSideColors=cc, margin=c(5, 10),
                    xlab="Biomes", ylab= "Cellular Compartment",
                    main="Biomes vs CC",
                    tracecol="steelblue", density="density",
                    srtRow=55, adjRow=c(0, 1), srtCol=10, adjCol=c(1,1)
)

############################################################
#                           ____END____
############################################################

mf.mat  <-  mf %>% select(-description) %>% as.matrix()
rownames(mf.mat) <- mf$description


####################################################
# Cellular compartment
#####################################################

cc.mat  <-  cc %>% select(-description) %>% as.matrix()
rownames(cc.mat) <- cc$description



bp<- project2 %>% filter(category=="biological process") %>% unique() 
cp<- project2 %>% filter(category=="cellular component") %>% unique()
mf<- project2 %>% filter(category=="molecular function") %>% unique()


project    
group_by(Taxa,Study.ID,Project.Name, Biome) %>% mutate(count=sum(count)) %>% 
    ungroup() %>% unique() %>% filter(!Taxa=="Root")
  go_dat[[i]] <- project 
  



go_dat <-go_dat %>% bind_rows() %>%
  select(-Study.ID) %>% 
  group_by(Taxa,Biome) %>% 
  mutate(count=sum(count)) %>% 
  ungroup() %>% 
  unique() %>% 
  spread(Biome, count, fill=0)

################################################
#   IPR abundances .....
################################################
ipr_dat  <- list()
system.time(for (i in biomes.meta$Study.ID){
  print(paste("Processing....", i, sep=" ") )
  project <- paste(prefix.down, i, go, sep="")
  project <- tbl_df(tryCatch(read.delim(file=project), error=function(e) NULL)) # %>% select(-GO,-category)
  if(nrow(project)> 1){
    project <- project %>% select(-IPR)
    project <- project %>% mutate(Study.ID=i)  %>% mygather()
    project <- left_join(project, biomes.meta, by="Study.ID") %>% 
      select(-Samples, -Run)  %>% 
      group_by(Taxa,Study.ID,Project.Name, Biome) %>% mutate(count=sum(count)) %>% 
      ungroup() %>% unique() %>% filter(!Taxa=="Root")
    ipr_dat[[i]] <- project 
    
  }
}
)


ipr <-ipr_dat %>% bind_rows() %>%
  select(-Study.ID) %>% 
  group_by(Taxa,Biome) %>% 
  mutate(count=sum(count)) %>% 
  ungroup() %>% 
  unique() %>% 
  spread(Biome, count, fill=0)

load('ipr.Rdata')

ipr.mat <- ipr %>% select(-Taxa) %>% as.matrix()
rownames(ipr.mat) <- ipr$Taxa
ipr.norm <- apply(ipr.mat, MARGIN = 2, FUN = function(X) (norm_vec(X)))


#############################################3
#    DR
##############################################
ipr.ca <- CA(ipr.mat, graph=FALSE, ncp=10)
ipr.ctrb  <- fviz_ca_row(ipr.ca, alpha.row="contrib", select.row=list(contrib=50))
best.ctrb <- as.character(ipr.ctrb$data[,"name"])
ipr.mat.ctrb <- ipr.mat[as.character(rownames(ipr.mat)) %in% best.ctrb, ]
corr.ipr.mat.ctrb <- rcorr(t(ipr.mat.ctrb))

corrplot(corr.ipr.mat.ctrb$r, type="lower", order="FPC", method='pie', main="Families correlation",
         p.mat = corr.ipr.mat.ctrb$P, sig.level = 0.05, insig = "blank", tl.col = "black", tl.srt = 45)




apply(ipr.mat.ctrb, MARGIN = 2, FUN = function(X) table((X)==0))

which(colnames(ipr.mat.ctrb)=="Salt crystallizer pond")
which(colnames(ipr.mat.ctrb)=="Oil-contaminated")
which(colnames(ipr.mat.ctrb)=="Agricultural land")
which(colnames(ipr.mat.ctrb)=="Hot (42-90C)")
ipr.mat.ctrb <- ipr.mat.ctrb[,-c(3,31,42,51)]
ipr.norm <- apply(ipr.mat.ctrb, MARGIN = 2, FUN = function(X) (norm_vec(X)))



##############################################

rc <- rainbow(nrow(ipr.norm), start=0, end=.3)
cc <- rainbow(ncol(ipr.norm), start=0, end=.3)
hv.ipr <- heatmap.2(ipr.norm, col=cm.colors(255), scale="row",
                       RowSideColors=rc, ColSideColors=cc, margin=c(5, 10),
                       xlab="Biomes", ylab= "Families",
                       main="Biomes vs Families",
                       tracecol="steelblue", density="density",
                       srtRow=55, adjRow=c(0, 1), srtCol=10, adjCol=c(1,1)
)

ipr.ca <- CA(ipr.mat, graph=TRUE, ncp=10)
ipr.hcpc <- HCPC(ipr.ca, cluster.CA = "columns")



fviz_dend(ipr.hcpc, horiz=T, rec_fill=TRUE, rect = TRUE, type="rectangle", k=5, sub=NULL, main="Biomes Clustering") + 
  theme( axis.line.x  = element_line(colour = "grey", size = .3, linetype = "solid"))

###############################################
# BEst contributing ....
###############################################
ipr.ctrb  <- fviz_ca_row(ipr.ca, alpha.row="contrib", select.row=list(contrib=25))
best.ctrb <- as.character(ipr.ctrb$data[,"name"])
ipr.mat.ctrb <- ipr.mat[as.character(rownames(ipr.mat)) %in% best.ctrb, ]
corr.ipr.mat.ctrb <- rcorr(t(ipr.mat.ctrb))

corrplot(corr.ipr.mat.ctrb$r, type="lower", order="FPC", method='pie', main="Families correlation",
         p.mat = corr.ipr.mat.ctrb$P, sig.level = 0.05, insig = "blank", tl.col = "black", tl.srt = 45)




apply(ipr.mat.ctrb, MARGIN = 2, FUN = function(X) table((X)==0))

which(colnames(ipr.mat.ctrb)=="Salt crystallizer pond"
which(colnames(ipr.mat.ctrb)=="Oil-contaminated"
which(colnames(ipr.mat.ctrb)=="Agricultural land"
            
X <- ipr.mat.ctrb[,-c(60,51)]

ipr.ca <- CA(X, graph=TRUE, ncp=10)
ipr.hcpc <- HCPC(ipr.ca, cluster.CA = "columns")



fviz_dend(ipr.hcpc, horiz=T, rec_fill=TRUE, rect = TRUE, type="rectangle", k=5, sub=NULL, main="Biomes Clustering") + 
  theme( axis.line.x  = element_line(colour = "grey", size = .3, linetype = "solid"))

comp.norm2 <- apply(comp.mat.ctrb, MARGIN = 2, FUN = function(X) (norm_vec(X)))


rc <- rainbow(nrow(comp.norm2), start=0, end=.3)
cc <- rainbow(ncol(comp.norm2), start=0, end=.3)
hv.notara <- heatmap.2(comp.norm2, col=cm.colors(255), scale="column",
                       RowSideColors=rc, ColSideColors=cc, margin=c(5, 10),
                       xlab="Biomes", ylab= "Families",
                       main="Biomes vs Families",
                       tracecol="steelblue", density="density",
                       srtRow=65, adjRow=c(0, 1), srtCol=45, adjCol=c(1,1)
)

