# last updated: 7/25/2023
# load libraries
require(biomaRt)
library(devtools)
library(edgeR)
library(limma)
library(Glimma)
library(ggplot2)
library(RColorBrewer)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(msigdbr)
library(WriteXLS)
library(scales)
library(wesanderson)
library(gdata)
library(ggpubr)
library(factoextra)
library(tidyr)
library(readxl)
library(pheatmap)
library(openxlsx)
library(dplyr)
library(magrittr)
library(DT)
library(org.Ce.eg.db)
library(gplots)
library(rstudioapi)
library(stringr)
library(data.table)
library(RRHO)
library(gridExtra)
library(dendextend)
library(GO.db)
library(tidyverse)
library(ggrepel)
library(kableExtra)

# create raw counts matrix from STAR output
# use count_result_star folder as directory for download_raw_reads

Download_raw_reads <- function(featurecounts_dir){
  setwd(featurecounts_dir)
  
  temp_data <- read.csv(dir()[!grepl(".summary$",dir())][1],header=T,sep='\t',skip = 1)
  counts_star <- data.frame(ID=temp_data$Geneid)
  rownames(counts_star) <- counts_star$ID
  for (i in dir()[!grepl(".summary$",dir())]){
    temp_name <- strsplit(i,".count")[[1]][1]
    temp_name <- gsub("-","_",temp_name)
    temp_data <- read.csv(i,header=T,sep='\t',skip = 1)
    temp_counts <- temp_data[,7]
    counts_star[temp_name] <- temp_counts
  }
  rm(temp_data)
  counts_star$ID <- rownames(counts_star)
  counts_star <- counts_star[,-1]
  return(counts_star)
}

count_matrix <- Download_raw_reads('/home/wmitchell/greer_RNAseq_celegans/count_result_star')
# write as xls
WriteXLS(as.data.frame(count_matrix), 'genecounts_greer.xls', row.names = T)

########### excluding transcripts with low expression ##########

keep <- rowSums(count_matrix) >= ncol(count_matrix)
counts.keep <- count_matrix[keep,]

# normalize input counts with edgeR

dgeObj <- DGEList(counts.keep)
barplot(dgeObj$samples$lib.size, names=colnames(dgeObj), las=2)

png('boxplot_rawcounts.png', width = 6*400, height= 8*400, pointsize = 8, res = 400)
logcounts <- cpm(dgeObj,log=TRUE)
boxplot(logcounts, xlab="", ylab="Log2 counts per million")
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalized)")
dev.off()

dgeObj <- calcNormFactors(dgeObj, method = 'RLE')
logcounts <- cpm(dgeObj,log=T)

boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (normalized)")

WriteXLS(as.data.frame(logcounts), 'normalized_genecounts_greer.xls', row.names = T)

# for translation efficiency, normalization was performed by dividing raw IP counts by input counts)

############ import counts and perform PCA #####################################################
# set wd
setwd("/home/wmitchell/greer_new/data")
# import filtered normalized genecounts (filtered raw IP counts divided by input counts)
counts = read.xls('genecounts_translation_efficiency.xls', row.names = T)

# replace NA with 0 
counts[is.na(counts)] <- 0

# read metadata
meta <- read.csv("metadata.csv")

# pca
pca <- prcomp(t(counts))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
PCAsummary <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = meta$treatment) 

g <- ggplot(data = PCAsummary, aes_string(x = "PC1", y = "PC2", color = "group")) + 
  geom_point(size = 3) + 
  xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_rect(linetype = "solid", fill = NA, size = 0.75), plot.margin = margin(1,1,1,1,"mm")) +
  theme(axis.ticks.length = unit(1, "mm"), axis.ticks = element_line(size = 0.5), axis.text = element_text(size = rel(0.75))) +
  theme(axis.title = element_text(size = rel(0.75)))+
  theme(legend.text = element_text(size = rel(0.65)), legend.title = element_blank(), legend.key.size = unit(15, "pt"), legend.margin = margin(0,0,0,0,"pt"), legend.box.margin=margin(0,0,0,-10, "pt"))

ggsave(filename="PCA_input.pdf", plot=g, device="pdf", units="in", width=4, height=4)
g

############## differential expression analysis #####################################
design <- model.matrix(~factor(meta$treatment))
colnames(counts[,colnames(counts)%in% meta$sampleID])==meta$sampleID

#### prepare DEGs
fit <- lmFit(counts, design, method = 'robust') #input method = robust
ebayes <- eBayes(fit)
tab <- topTable(ebayes, coef=2, adjust="fdr", n=nrow(counts))
tab$gene = rownames(tab)

####### check the sign of diff expression 
gene1 = tab$gene[tab$logFC>0][1]
counts[grep(gene1, row.names(counts)),]

#### if the sign is reversed, run this
# tab$logFC = tab$logFC * (-1)

write.csv(as.data.frame(tab), paste('DEGs_translation_efficiency','.csv', sep = ''), row.names = T)

################# heatmap of differentially expressed genes (FDR < 0.05) ############################
colors <- colorRampPalette(c("blue","black","yellow"))(256)
plotData <- counts

# keep only genes that are significantly DE by FDR < 0.05
DE_genes <- tab[tab$adj.P.Val <= 0.05,] %>% row.names()

plotData <- plotData[DE_genes,] %>% scale()
distCol <- t(plotData) %>% dist()
hclustCol <- hclust(distCol, method = 'complete') %>% rotate(., order = colnames(plotData)) 


p <- pheatmap(plotData, clustering_distance_rows = 'correlation', clustering_method = "average", scale = 'row', cluster_cols = hclustCol,
              color = colors,  border_color = NA, show_rownames = FALSE, fontsize = 8, silent = TRUE
              
)

png(filename="inputs_heatmap.png", res = 300, height = 4, width = 4, units = 'in')
p
invisible(dev.off())

################# GO enrichment analysis #########################################
tab = read.csv("DEGs_translation_efficiency.csv", header = T, row.names = 1)
tab$rank = ifelse(tab$logFC>0, 1, -1)*(-log10(tab$P.Value))
tab = tab[order(tab$rank, decreasing = T), ]
geneList = tab$rank
names(geneList)= tab$gene
head(geneList)

gse_GO = gseGO(geneList = geneList, OrgDb = org.Ce.eg.db,
               keyType ="WORMBASE",
               pAdjustMethod = "fdr",
               pvalueCutoff  = 0.05)

png('GO_translation_efficiency.png')
dotplot(gse_GO, showCategory = 15, split = '.sign', font.size = 12)+ 
  facet_grid(.~.sign)+ 
  theme(strip.text.x = element_text(size = 16), 
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16))+
  theme_bw()
dev.off()

write.csv(gse_GO@result,  paste('GSEA_GO_translation_efficiency','.csv', sep = ''))

####### convert WORMBASE IDs to SYMBOL #########################
eg = bitr(unique(tab$geneList), fromType="WORMBASE", toType="SYMBOL", OrgDb="org.Ce.eg.db")
tab = cbind(tab, eg$SYMBOL)
write.csv(tab, "DEGs_symbol_translation_efficiency.csv")

############### volcano plot ##################################
# function for outputting tables 
knitr_table <- function(x) {
  x %>% 
    knitr::kable(format = "html", digits = Inf, 
                 format.args = list(big.mark = ",")) %>%
    kableExtra::kable_styling(font_size = 15)
}

# import data
data <- read.csv("DEGs_symbol_translation_efficiency.csv")

head(data) %>% 
  knitr_table()

# simple volcano plot
p1 <- ggplot(data, aes(logFC, -log10(adj.P.Val)))+
  geom_point(size = 2/5) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR"))
p1

# add color to plot
data <- data %>% 
  mutate(
    Change = case_when(logFC >= 1 & -log10(adj.P.Val)	 >= 2 ~ "upregulated",
                       logFC <= -1 & -log10(adj.P.Val)	 >= 2 ~ "downregulated",
                       TRUE ~ "Unchanged")
  )

head(data) %>% 
  knitr_table()

# add FDR labels 
data <- data %>% 
  mutate(
    Significance = case_when(
      abs(logFC) >= 1 & adj.P.Val <= 0.05 & adj.P.Val > 0.01 ~ "FDR 0.05", 
      abs(logFC) >= 1 & adj.P.Val <= 0.01 & adj.P.Val > 0.001 ~ "FDR 0.01",
      abs(logFC) >= 1 & adj.P.Val <= 0.001 ~ "FDR 0.001", 
      TRUE ~ "Unchanged")
  )

# plot
png('translation_efficiency_volcano.png', width = 3.5*400, height= 2.6*400, pointsize = 8, res = 400)
p2 <- ggplot(data, aes(logFC, -log10(adj.P.Val))) +
  geom_point(aes(color = Significance), size = 2/5) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR")) +
  scale_color_viridis_d() +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  theme_bw()
dev.off()

# label top 25 proteins upregulated and downregulated
top <- 25
top_genes <- bind_rows(
  data %>% 
    filter(Change == 'upregulated') %>% 
    arrange(adj.P.Val	, desc(abs(logFC))) %>% 
    head(top),
  data %>% 
    filter(Change == 'downregulated') %>% 
    arrange(adj.P.Val	, desc(abs(logFC))) %>% 
    head(top),
)
top_genes %>% 
  knitr_table()

# plot
p3 <-  p2 +
  geom_label_repel(data = top_genes, max.overlaps = Inf,
                   mapping = aes(logFC, -log10(adj.P.Val)	, label = SYMBOL),
                   size = 2) 

p3

#################### Revigo plot ###########################################
# A plotting R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# --------------------------------------------------------------------------
revigo.names <- c("term_ID","description","frequency","plot_X","plot_Y","log_size","value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0000413","protein peptidyl-prolyl isomerization",0.2195504443282802,-6.48608631751928,-1.2586259452114004,1.3424226808222062,-1.5786845634939757,0.9519227344536538,0.21479701),
                     c("GO:0000902","cell morphogenesis",1.6309461578672242,-0.4941040566494735,7.176873743673375,2.1958996524092336,-1.3404682398737766,0.8591284958134999,0.66739765),
                     c("GO:0000904","cell morphogenesis involved in differentiation",1.3904861474124413,-0.3427683906686398,7.008160237478215,2.1271047983648077,-1.730059439608469,0.7931346485275615,0.39566079),
                     c("GO:0002181","cytoplasmic translation",0.501829587036069,-5.519530011397359,-2.731680712828029,1.6901960800285136,-4.236572006437063,0.8795982491228385,0.12074798),
                     c("GO:0002183","cytoplasmic translational initiation",0.2195504443282802,-5.834836692320991,-2.611711405069071,1.3424226808222062,-1.519294448263424,0.8871294256217143,0.61515595),
                     c("GO:0002376","immune system process",2.6868792472556193,-7.3315540646688095,0.7809737317591259,2.41161970596323,-5.966576244513051,1,-0),
                     c("GO:0003008","system process",6.61787767903816,2.9342430001614184,6.383369398279233,2.802089257881733,-4.236572006437063,0.8650708630919245,0.37048637),
                     c("GO:0003012","muscle system process",0.2927339257710403,3.486109032292473,6.712306143573584,1.462397997898956,-2.827076428020524,0.8620132298048773,0.5177252),
                     c("GO:0006091","generation of precursor metabolites and energy",1.3173026659696812,-3.826938541321201,-5.301988166906822,2.103803720955957,-7.671620396561262,0.9322533509537803,0),
                     c("GO:0006457","protein folding",1.1604809200209096,-7.058724318527822,2.4670272835453524,2.0492180226701815,-7.300162274132754,0.986419055175234,0.01163397),
                     c("GO:0006605","protein targeting",0.9827496079456352,0.05523166760872622,-6.1317626411637125,1.9777236052888478,-1.6992790496484487,0.879897815244999,0.59293207),
                     c("GO:0006805","xenobiotic metabolic process",0.4286461055933089,4.502181084337197,2.0461050190292904,1.6232492903979006,-1.4301278740973578,0.821018937649388,0.46298712),
                     c("GO:0006839","mitochondrial transport",0.5227391531625719,0.21304589682608124,-5.77822357710049,1.7075701760979363,-2.0513575252475325,0.9048185991214867,0.39323715),
                     c("GO:0006888","endoplasmic reticulum to Golgi vesicle-mediated transport",0.6168322007318349,0.49543658772235855,-6.163832477918153,1.7781512503836436,-1.646987171588756,0.8970714042914613,0.52684845),
                     c("GO:0006936","muscle contraction",0.2927339257710403,3.637851606434256,6.492017132390341,1.462397997898956,-2.827076428020524,0.8620132298048773,0.5177252),
                     c("GO:0006937","regulation of muscle contraction",0.5122843700993205,4.403483786685232,-5.689215847370495,1.6989700043360187,-1.333362586166903,0.9562166692363155,0.54147459),
                     c("GO:0006952","defense response",4.087820177731312,5.966038227308893,2.302476415428416,2.593286067020457,-6.373659632624958,0.8737764712592905,0.62801063),
                     c("GO:0006955","immune response",2.6450601150026136,6.855701902079767,2.272807775222166,2.404833716619938,-5.58838029403677,0.8895512170755571,0.25891004),
                     c("GO:0007005","mitochondrion organization",1.704129639309984,-4.254375020390149,5.331218135606481,2.214843848047698,-3.2709024663393107,0.8975129673541473,0.01216601),
                     c("GO:0007155","cell adhesion",1.003659174072138,-5.8194255962742885,1.1768566699065544,1.9867717342662448,-2.37799996382682,0.9866242026050953,0.01144485),
                     c("GO:0007166","cell surface receptor signaling pathway",3.1573444851019343,5.538074556005163,0.8933912536585534,2.481442628502305,-1.5037728586311132,0.861376676306728,0.3723551),
                     c("GO:0007186","G protein-coupled receptor signaling pathway",7.48562467328803,5.736512699497017,0.5674323879798661,2.8555191556678,-3.314830832842499,0.8452847692805744,0.47240221),
                     c("GO:0007218","neuropeptide signaling pathway",1.4741244119184527,6.122241085309467,0.8798359448953806,2.1522883443830563,-3.2709024663393107,0.8732168500147441,0.23852984),
                     c("GO:0007399","nervous system development",3.56508102456874,1.348845462572868,6.494520976576983,2.534026106056135,-3.2709024663393107,0.7636344293969247,0.30719035),
                     c("GO:0007635","chemosensory behavior",2.4464192368008364,4.585598551245444,3.9653216282461905,2.3710678622717363,-2.0950382800746854,0.7743481816585197,0.56809203),
                     c("GO:0008037","cell recognition",0.334553058024046,-0.07137231340272582,-1.193400698328,1.5185139398778875,-1.7029374147517276,0.9879973903978913,0.01019123),
                     c("GO:0008038","neuron recognition",0.2822791427077888,0.957411129591532,6.55095516566233,1.4471580313422192,-1.333362586166903,0.7677042555937347,0.59471112),
                     c("GO:0009141","nucleoside triphosphate metabolic process",0.48092002090956615,-4.894229513839765,-3.921789326953855,1.6720978579357175,-2.1106536375058993,0.8264560864309412,0.67467362),
                     c("GO:0009201","ribonucleoside triphosphate biosynthetic process",0.31364349189754315,-5.163985818576719,-3.6151969419341032,1.4913616938342726,-4.157390760389438,0.778929518412203,0.30305315),
                     c("GO:0009607","response to biotic stimulus",4.056455828541558,6.603458133503103,1.7894145016870024,2.5899496013257077,-6.373659632624958,0.8976022974500109,0.27616896),
                     c("GO:0022618","protein-RNA complex assembly",1.0454783063251438,-3.6675066094069675,5.38063213256391,2.0043213737826426,-1.3404682398737766,0.8997595395104674,0.63157696),
                     c("GO:0022900","electron transport chain",0.6900156821745949,-2.4031586038199353,-2.303482982568032,1.8260748027008264,-7.671620396561262,0.8132050942381421,0.12488009),
                     c("GO:0030030","cell projection organization",2.8541557762676426,-4.133724185998792,4.774773192904503,2.437750562820388,-2.108695529374785,0.9010785762298078,0.37132233),
                     c("GO:0032101","regulation of response to external stimulus",1.2859383167799268,7.084819946607967,-2.763966392136808,2.093421685162235,-1.624221016177851,0.9702029614802237,0.10388624),
                     c("GO:0032543","mitochondrial translation",0.4704652378463147,-5.472470825393914,-2.5313757910419468,1.662757831681574,-1.5851641948251172,0.880217095391878,0.65916203),
                     c("GO:0032989","cellular component morphogenesis",2.007318348144276,-1.9762247450295423,6.018628724429241,2.285557309007774,-1.509736595506118,0.7790035981491742,0.65438983),
                     c("GO:0033108","mitochondrial respiratory chain complex assembly",0.36591740721380034,-3.8897416126069153,5.633474188284792,1.5563025007672873,-2.5684797571711644,0.8889939864504971,0.39557499),
                     c("GO:0040017","positive regulation of locomotion",1.3068478829064296,6.08407885214787,-3.839656258230125,2.100370545117563,-1.4067070530117611,0.9738962801638776,0.10755014),
                     c("GO:0043062","extracellular structure organization",0.35546262415054886,-4.718081616366895,4.592903261890692,1.5440680443502757,-1.5187017876717965,0.9197525772693596,0.27352885),
                     c("GO:0044057","regulation of system process",0.8259278619968635,4.548535030761599,-5.416681908050179,1.9030899869919435,-2.29479823610276,0.9623801045816096,-0),
                     c("GO:0044419","biological process involved in interspecies interaction between organisms",4.056455828541558,-3.1356821276719926,0.47362469733334406,2.5899496013257077,-6.882728704344236,1,-0),
                     c("GO:0045214","sarcomere organization",0.26136957658128596,-2.188129511501771,6.036393119437686,1.414973347970818,-1.4401879968622748,0.7831787928650256,0.57373018),
                     c("GO:0048193","Golgi vesicle transport",1.1918452692106638,1.017671227738087,-6.2151720606353145,2.060697840353612,-1.337472603782109,0.9428220345827734,0.61724123),
                     c("GO:0050804","modulation of chemical synaptic transmission",0.9618400418191323,3.3975394920726423,-2.200689580191268,1.968482948553935,-2.0172992088861,0.9672859719417923,0.09901423),
                     c("GO:0050911","detection of chemical stimulus involved in sensory perception of smell",1.9027705175117617,4.497323144389513,4.300394599682972,2.2624510897304293,-4.570247719997592,0.7169570852389788,0.2470191),
                     c("GO:0051606","detection of stimulus",2.3732357553580763,6.460613153345802,3.0306577337656875,2.357934847000454,-1.880425305129893,0.9044062545692901,0.25487087),
                     c("GO:0061077","chaperone-mediated protein folding",0.2195504443282802,-0.7001754362809086,0.8498684416809856,1.3424226808222062,-1.631115576324238,0.988452214813675,0.00978049),
                     c("GO:0070585","protein localization to mitochondrion",0.3972817564035547,-0.5121030365487556,-6.154717258798767,1.591064607026499,-4.059981844992337,0.8902402283045753,0.65231542),
                     c("GO:0072655","establishment of protein localization to mitochondrion",0.3972817564035547,-0.32609480340689595,-6.0023115315283855,1.591064607026499,-4.059981844992337,0.8702034134751452,0.01036889),
                     c("GO:0090150","establishment of protein localization to membrane",0.6377417668583377,-0.13654887603528224,-6.358507717846144,1.792391689498254,-1.4619114112738243,0.8938514522579448,0.62192981),
                     c("GO:0098542","defense response to other organism",4.035546262415054,5.9713521682798945,2.5250176285175887,2.5877109650189114,-6.431798275933005,0.8221701603220043,-0),
                     c("GO:0099177","regulation of trans-synaptic signaling",0.9618400418191323,3.6744062911377746,-2.3499146840118237,1.968482948553935,-2.0172992088861,0.9672859719417923,0.6357535),
                     c("GO:0120036","plasma membrane bounded cell projection organization",2.7391531625718764,-3.7519341730117994,4.6036638836074815,2.419955748489758,-2.1106536375058993,0.8507596493646798,0.34635388),
                     c("GO:1902884","positive regulation of response to oxidative stress",0.334553058024046,6.731119807157309,-3.4841520823078294,1.5185139398778875,-1.3442142860360293,0.9677755329074159,0.54113732));

one.data <- data.frame(revigo.data);
names(one.data) <- revigo.names;
one.data <- one.data [(one.data$plot_X != "null" & one.data$plot_Y != "null"), ];
one.data$plot_X <- as.numeric( as.character(one.data$plot_X) );
one.data$plot_Y <- as.numeric( as.character(one.data$plot_Y) );
one.data$log_size <- as.numeric( as.character(one.data$log_size) );
one.data$value <- as.numeric( as.character(one.data$value) );
one.data$frequency <- as.numeric( as.character(one.data$frequency) );
one.data$uniqueness <- as.numeric( as.character(one.data$uniqueness) );
one.data$dispensability <- as.numeric( as.character(one.data$dispensability) );

p1 <- ggplot( data = one.data );
p1 <- p1 + geom_point( aes( plot_X, plot_Y, colour = value, size = log_size), alpha = I(0.6) );
p1 <- p1 + scale_colour_gradientn( colours = c("blue", "green", "yellow", "red"), limits = c( min(one.data$value), 0) );
p1 <- p1 + geom_point( aes(plot_X, plot_Y, size = log_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) ));
p1 <- p1 + scale_size( range=c(5, 30)) + theme_bw(); # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
ex <- one.data [ one.data$dispensability < 0.15, ];
p1 <- p1 + geom_text( data = ex, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 3 );
p1 <- p1 + labs (y = "semantic space x", x = "semantic space y");
p1 <- p1 + theme(legend.key = element_blank()) ;
one.x_range = max(one.data$plot_X) - min(one.data$plot_X);
one.y_range = max(one.data$plot_Y) - min(one.data$plot_Y);
p1 <- p1 + xlim(min(one.data$plot_X)-one.x_range/10,max(one.data$plot_X)+one.x_range/10);
p1 <- p1 + ylim(min(one.data$plot_Y)-one.y_range/10,max(one.data$plot_Y)+one.y_range/10);

p1

################################### end of script ####################################################




