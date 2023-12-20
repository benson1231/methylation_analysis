# 0.load package ------------------------------------------------------------
library(sesame)
library(BiocParallell)
library(tidyverse)
library(ggplot2)
library(wheatmap)
library(ggVennDiagram)


# 1.load data and preprocess data using QCDPB -----------------------------
# load sample list
df <- read.csv("/Users/benson/Documents/raw_data/methy1_metal/raw_data/SampleSheet_EPIC V2_N64.csv")
sample_list <- df[8:20, c(1,6,7)] %>% setNames(c("sample","chip","id"))
rownames(sample_list) <- c(1:13)
sample_list <- sample_list %>% mutate(., group = paste0(row.names(sample_list),"_",sample))

# check idat data
idat_dir <- "/Users/benson/Documents/raw_data/methy1_metal/raw_data/ip_L_V_L_p1"
setwd(idat_dir)
list.files(pattern = "*.idat")
searchIDATprefixes(".")

# # get beta value from IDAT file(way 1)
# beta <- do.call(cbind, BiocParallel::bplapply(
#   searchIDATprefixes(idat_dir), function(pfx) {
#     getBetas(prepSesame(readIDATpair(pfx, platform = "EPICv2"), "QCDPB"))
#   }, BPPARAM = BiocParallel::MulticoreParam(2)))
# # 重新命名sample ID
# beta <- beta %>% as.data.frame() %>% 
#   setNames(., sample_list$sample)

# get sigDF from IDAT file and preprocess data use QCDPB preprocessing function code
sdfs <- openSesame(idat_dir, prep = "QCDPB", func = NULL, platform = "EPICv2")
sdfs <- sdfs %>% setNames(.,sample_list$group) # 注意核對sample序號

# saveRDS(sdfs, file = "sdfs.RDS")
sdfs <- "/Users/benson/Documents/project/methy1/sdfs.RDS" %>% readRDS()


# 2.QC ----------------------------------------------------------------------
setwd("/Users/benson/Documents/project/methy1")
# calculate metrics on all IDATs in a specific folder
qcs = openSesame(idat_dir, prep="", func=sesameQC_calcStats)
qcs <- setNames(qcs, sample_list$sample)
sesameQC <- do.call(rbind, lapply(qcs, as.data.frame))
# QC for each parameter
detection_qc <- openSesame(sdfs, prep="", func=sesameQC_calcStats, funs="detection")
numProbes_qc <- openSesame(sdfs, prep="", func=sesameQC_calcStats, funs="numProbes")
intensity_qc <- openSesame(sdfs, prep="", func=sesameQC_calcStats, funs="intensity")
channel_qc <- openSesame(sdfs, prep="", func=sesameQC_calcStats, funs="channel")
dyeBias_qc <- openSesame(sdfs, prep="", func=sesameQC_calcStats, funs="dyeBias")
betas_qc <- openSesame(sdfs, prep="", func=sesameQC_calcStats, funs="betas")
# rename list
detection_qc <- detection_qc %>% setNames(sample_list$sample)
numProbes_qc <- numProbes_qc %>% setNames(sample_list$sample)
intensity_qc <- intensity_qc %>% setNames(sample_list$sample)
channel_qc <- channel_qc %>% setNames(sample_list$sample)
dyeBias_qc <- dyeBias_qc %>% setNames(sample_list$sample)
betas_qc <- betas_qc %>% setNames(sample_list$sample)
# plot QC parameters
sesameQC_plotBar(detection_qc)
ggsave("detection_qc.png")
sesameQC_plotBar(intensity_qc)
ggsave("intensity_qc.png")
sesameQC_plotBar(dyeBias_qc)
ggsave("dyeBias_qc.png")
sesameQC_plotBar(betas_qc)
ggsave("betas_qc.png")

# Dye bias Q-Q plot
pdf("Dye_bias_QQ_plot.pdf")
par(mfrow=c(3,5), mar=c(3,3,2,1))
for (i in 1:13) {
  sesameQC_plotRedGrnQQ(sdfs[[i]], main = names(sdfs)[i])
}
dev.off()

# Intensity-beta plot
pdf("Intensity_beta_plot.pdf")
par(mfrow=c(3,5), mar=c(3,3,2,1))
for (i in 1:13) {
  sesameQC_plotIntensVsBetas(sdfs[[i]], main=names(sdfs)[i]) 
}
dev.off()


# 3.get beta value (way 2)-------------------------------------------------
betas <- openSesame(sdfs, func = getBetas)  # for further use
# saveRDS(betas, file = "betas.RDS")
betas <- "/Users/benson/Documents/project/methy1/betas.RDS" %>% readRDS()

betas_pre <- openSesame(sdfs, func = getBetas, collapseToPfx = TRUE)  # for visualization
# saveRDS(betas_pre, file = "betas_pre.RDS")
betas_pre <- "/Users/benson/Documents/project/methy1/betas_pre.RDS" %>% readRDS()
all <- betas_pre %>% na.omit()

# 4.visualization -----------------------------------------------------------
# create function for chromosome visualization by Gene Name(need betas_pre)
map <- function(x = betas_pre, 
                group = "all",
                gene = gene,
                platform = "EPIC",
                genome = "hg38"){
  if(group == "all"){
    cat(c("-> plotting", gene, "in all group mapping", genome,"\n"))
    visualizeGene(gene, betas = x, platform = platform, genome = genome) +
      WLegendV("betas", TopRightOf("betas"))
  } else if(group == "as"){
    cat(c("-> plotting", gene, "in As group mapping", genome,"\n"))
    x %>% as.data.frame() %>% 
      select(ip_L_V_L_CON, ip_L_V_L_DMS, ip_L_V_L_BAP, ip_L_V_L_AS, ip_L_V_L_AS_BAP) %>%
      as.matrix() %>% 
      visualizeGene(gene, betas = ., platform = platform, genome = genome) +
      WLegendV("betas", TopRightOf("betas")) 
  } else if(group == "co"){
    cat(c("-> plotting", gene, "in Co group mapping", genome,"\n"))
    x %>% as.data.frame() %>% 
      select(ip_L_V_L_CON, ip_L_V_L_DMS, ip_L_V_L_BAP, ip_L_V_L_CO, ip_L_V_L_CO_BAP) %>% 
      as.matrix() %>% 
      visualizeGene(gene, betas = ., platform = platform, genome = genome) +
      WLegendV("betas", TopRightOf("betas"))
  } else if(group == "lcd"){
    cat(c("-> plotting", gene, "in L_Cd ggroup mapping", genome,"\n"))
    x %>% as.data.frame() %>% 
      select(ip_L_V_L_CON, ip_L_V_L_DMS, ip_L_V_L_BAP, ip_L_V_L_LCD, ip_L_V_L_LCD_BAP) %>% 
      as.matrix() %>% 
      visualizeGene(gene, betas = ., platform = platform, genome = genome) +
      WLegendV("betas", TopRightOf("betas")) 
  } else if(group == "hcd"){
    cat(c("-> plotting", gene, "in H_Cd group mapping", genome,"\n"))
    x %>% as.data.frame() %>% 
      select(ip_L_V_L_CON, ip_L_V_L_DMS, ip_L_V_L_BAP, ip_L_V_L_HCD, ip_L_V_L_HCD_BAP) %>%
      as.matrix() %>% 
      visualizeGene(gene, betas = ., platform = platform, genome = genome) +
      WLegendV("betas", TopRightOf("betas"))
  } else{
    cat("error")
  }
}
# use map function to plot gene on chromosome
map(x = betas_pre, group = "all", gene ="EGFR")

# chromosome visualization by Genomic Region 
visualizeRegion(
  "chr7",55019017,55211628, betas_pre, platform='EPIC', genome = "hg38",
  show.probeNames = TRUE) +
  WLegendV("betas", TopRightOf("betas")) 

# chromosome visualization by Probe ID 
visualizeProbes("cg16205058", betas_pre, platform='EPIC', genome = "hg38") +
  WLegendV("betas", TopRightOf("betas")) 


# 5.delta beta value ------------------------------------------------------
delta <- betas %>% as.data.frame(.) %>% 
  mutate(delta_AZA = ip_L_V_L_AZA - ip_L_V_L_CON,
         delta_DAC = ip_L_V_L_DAC - ip_L_V_L_CON,
         delta_AS = ip_L_V_L_AS - ip_L_V_L_CON,
         delta_CO = ip_L_V_L_CO - ip_L_V_L_CON,
         delta_LCD = ip_L_V_L_LCD - ip_L_V_L_CON,
         delta_HCD = ip_L_V_L_HCD - ip_L_V_L_CON,
         delta_BAP = ip_L_V_L_BAP - ip_L_V_L_DMS,
         delta_AS_BAP = ip_L_V_L_AS_BAP - ip_L_V_L_DMS,
         delta_CO_BAP = ip_L_V_L_CO_BAP - ip_L_V_L_DMS,
         delta_LCD_BAP = ip_L_V_L_LCD_BAP - ip_L_V_L_DMS,
         delta_HCD_BAP = ip_L_V_L_HCD_BAP - ip_L_V_L_DMS)

AS_0.1 <- delta %>% filter(abs(delta_AS) > 0.1) %>% 
  select(delta_AS, ip_L_V_L_CON, ip_L_V_L_AS)
CO_0.1 <- delta %>% filter(abs(delta_CO) > 0.1) %>% 
  select(delta_CO, ip_L_V_L_CON, ip_L_V_L_CO)
LCD_0.1 <- delta %>% filter(abs(delta_LCD) > 0.1) %>% 
  select(delta_LCD, ip_L_V_L_CON, ip_L_V_L_LCD)
HCD_0.1 <- delta %>% filter(abs(delta_HCD) > 0.1) %>% 
  select(delta_HCD, ip_L_V_L_CON, ip_L_V_L_HCD)
BAP_0.1 <- delta %>% filter(abs(delta_BAP) > 0.1) %>% 
  select(delta_BAP, ip_L_V_L_DMS, ip_L_V_L_BAP)
AS_BAP_0.1 <- delta %>% filter(abs(delta_AS_BAP) > 0.1) %>% 
  select(delta_AS_BAP, ip_L_V_L_DMS, ip_L_V_L_AS_BAP)
CO_BAP_0.1 <- delta %>% filter(abs(delta_CO_BAP) > 0.1) %>% 
  select(delta_CO_BAP, ip_L_V_L_DMS, ip_L_V_L_CO_BAP)
LCD_BAP_0.1 <- delta %>% filter(abs(delta_LCD_BAP) > 0.1) %>% 
  select(delta_LCD_BAP, ip_L_V_L_DMS, ip_L_V_L_LCD_BAP)
HCD_BAP_0.1 <- delta %>% filter(abs(delta_HCD_BAP) > 0.1) %>% 
  select(delta_HCD_BAP, ip_L_V_L_DMS, ip_L_V_L_HCD_BAP)


# 6.venn diagram of all(up and down probes,>0.1) -----------------------------
as_all <- list(As = rownames(AS_0.1),
               BAP = rownames(BAP_0.1),
               As_BAP = rownames(AS_BAP_0.1))
co_all <- list(Co = rownames(CO_0.1),
               BAP = rownames(BAP_0.1),
               Co_BAP = rownames(CO_BAP_0.1))
lcd_all <- list(L_Cd = rownames(LCD_0.1),
                BAP = rownames(BAP_0.1),
                L_Cd_BAP = rownames(LCD_BAP_0.1))
hcd_all <- list(H_Cd = rownames(HCD_0.1),
                BAP = rownames(BAP_0.1),
                H_Cd_BAP = rownames(HCD_BAP_0.1))
all_probe <- list(As = rownames(AS_0.1),
                  Co = rownames(CO_0.1),
                  L_Cd = rownames(LCD_0.1),
                  H_Cd = rownames(HCD_0.1),
                  BAP = rownames(BAP_0.1))

ggVennDiagram(as_all) + scale_fill_gradient(low="blue",high = "red")
ggVennDiagram(co_all) + scale_fill_gradient(low="blue",high = "red")
ggVennDiagram(lcd_all) + scale_fill_gradient(low="blue",high = "red")
ggVennDiagram(hcd_all) + scale_fill_gradient(low="blue",high = "red")
ggVennDiagram(all_probe) + scale_fill_gradient(low="blue",high = "red")


# 7.抓出Venn Diagram中的gene list ----------------------------------------
as_probe_list <- process_region_data(Venn(as_all))
as_probe_list$item[[7]]

# 8.KnowYourCG Visualization -----------------------------------------------
# test the enrichment over database groups
results <- testEnrichment(as_probe_list$item[[3]], platform = "EPICv2")
head(results)

# plot
KYCG_plotEnrichAll(results)
# dotplot
KYCG_plotDot(results)
# barplot
library(ggplot2)
library(wheatmap)
p1 <- KYCG_plotBar(results, label=TRUE)
p2 <- KYCG_plotBar(results, y="estimate") + ylab("log2(Odds Ratio)") +
  xlab("") + theme(axis.text.y = element_blank())
WGG(p1) + WGG(p2, RightOf(width=0.5, pad=0))
# Waterfall plot
KYCG_plotWaterfall(results)










