# load package -----------------------------------------------------------------
library(sesame)
library(dplyr)
library(illuminaio)
library(ggplot2)

# 工作位置 ---------------------------------------------------------------------
setwd("/Users/benson/Documents/project/methy1")
dir <- "metal_delta_beta"

# 讀檔 -------------------------------------------------------------------------
aza <- readxl::read_xlsx(file.path(dir,"AZA_0.1deltabeta.xlsx"))
dac <- readxl::read_xlsx(file.path(dir,"DAC_0.1deltabeta.xlsx"))
as <- readxl::read_xlsx(file.path(dir,"As_0.1deltabeta.xlsx"))
co <- readxl::read_xlsx(file.path(dir,"Co_0.1deltabeta.xlsx"))
lcd <- readxl::read_xlsx(file.path(dir,"L_Cd_0.1deltabeta.xlsx"))
hcd <- readxl::read_xlsx(file.path(dir,"H_Cd_0.1deltabeta.xlsx"))
bap <- readxl::read_xlsx(file.path(dir,"BAP_0.1deltabeta.xlsx"))
asbap <- readxl::read_xlsx(file.path(dir,"As_BAP_0.1deltabeta.xlsx"))
cobap <- readxl::read_xlsx(file.path(dir,"Co_BAP_0.1deltabeta.xlsx"))
lcdbap <- readxl::read_xlsx(file.path(dir,"L_Cd_BAP_0.1deltabeta.xlsx"))
hcdbap <- readxl::read_xlsx(file.path(dir,"H_Cd_BAP_0.1deltabeta.xlsx"))

# 建立差異表達基因list ---------------------------------------------------------
list <- c(aza$probeID, dac$probeID, as$probeID, co$probeID, lcd$probeID, hcd$probeID,
          bap$probeID, asbap$probeID, cobap$probeID, lcdbap$probeID, hcdbap$probeID) %>% 
  unique() %>% unlist()

# 前處理 -----------------------------------------------------------------------
# 重複probe合併成一個，取平均值
# beta為所有probe與組別的beta value matrix(在raw data資料夾中，手動匯入)
beta <- as.data.frame(beta)
beta$ID <- substr(row.names(beta), 1, 10)  # 去除probe後綴
all <- aggregate(x = beta, by = list(beta$ID),  FUN = mean)
rownames(all) <- all$Group.1 
all <- all[,2:14]
# saveRDS(all,"all.RDS")

# 差異表達部分 
betas <- beta %>% filter(ID %in% list)   # 抓出差異基因表達beta value
# 重複probe合併成一個，取平均值
agg <- aggregate(x = betas, by = list(betas$ID_pre),  FUN = mean)
rownames(agg) <- agg$Group.1 
agg <- agg[,2:14]

# 讀檔(已整理) -----------------------------------------------------------------
# beta為所有probe與組別的beta value matrix(在raw data資料夾中，手動匯入)
beta <- as.data.frame(beta)
agg <- readRDS("/Users/benson/Documents/project/methy1/agg.RDS")  # 差異表達probe
all <- readRDS("/Users/benson/Documents/project/methy1/all.RDS")  # 全部probe

# chromosome visualization by Genomic Region -----------------------------
visualizeRegion(
  'chr19',10260000,10380000, all, platform='EPIC',
  show.probeNames = TRUE)

# chromosome visualization by Gene Name ----------------------------------
visualizeGene("APEX1", betas = all, platform='EPIC')
# 組別分開
all %>% 
  select(ip_L_V_L_CON, ip_L_V_L_DMS, ip_L_V_L_BAP, ip_L_V_L_AS, ip_L_V_L_AS_BAP) %>% 
  visualizeGene("TP53", betas = ., platform='EPIC')
all %>% 
  select(ip_L_V_L_CON, ip_L_V_L_DMS, ip_L_V_L_BAP, ip_L_V_L_CO, ip_L_V_L_CO_BAP) %>% 
  visualizeGene("TP53", betas = ., platform='EPIC')
all %>% 
  select(ip_L_V_L_CON, ip_L_V_L_DMS, ip_L_V_L_BAP, ip_L_V_L_LCD, ip_L_V_L_LCD_BAP) %>% 
  visualizeGene("TP53", betas = ., platform='EPIC')
all %>% 
  select(ip_L_V_L_CON, ip_L_V_L_DMS, ip_L_V_L_BAP, ip_L_V_L_HCD, ip_L_V_L_HCD_BAP) %>% 
  visualizeGene("TP53", betas = ., platform='EPIC')

### gene name output
gene_list <- c("APEX1","TP53","BRCA1")
for(i in gene_list){
  png(paste0(i,"_gene.png"))
  visualizeGene(i, betas = all, platform='EPIC') %>% print()
  dev.off()
}
for(i in gene_list){
  png(paste0(i,"As_gene.png"))
  all %>% 
    select(ip_L_V_L_CON, ip_L_V_L_DMS, ip_L_V_L_BAP, ip_L_V_L_AS, ip_L_V_L_AS_BAP) %>% 
    visualizeGene(i, betas = ., platform='EPIC') %>% print()
  dev.off()
  png(paste0(i,"Co_gene.png"))
  all %>% 
    select(ip_L_V_L_CON, ip_L_V_L_DMS, ip_L_V_L_BAP, ip_L_V_L_CO, ip_L_V_L_CO_BAP) %>% 
    visualizeGene("TP53", betas = ., platform='EPIC') %>% print()
  dev.off()
  png(paste0(i,"L_Cd_gene.png"))
  all %>% 
    select(ip_L_V_L_CON, ip_L_V_L_DMS, ip_L_V_L_BAP, ip_L_V_L_LCD, ip_L_V_L_LCD_BAP) %>% 
    visualizeGene("TP53", betas = ., platform='EPIC') %>% print()
  dev.off()
  png(paste0(i,"H_Cd_gene.png"))
  all %>% 
    select(ip_L_V_L_CON, ip_L_V_L_DMS, ip_L_V_L_BAP, ip_L_V_L_HCD, ip_L_V_L_HCD_BAP) %>% 
    visualizeGene("TP53", betas = ., platform='EPIC') %>% print()
  dev.off()
}

# function for chromosome visualization by Gene Name ---------------------------
save_image <- function(x = all, 
                group = all,
                gene = gene,
                platform = "EPIC"){
  if(group == "all"){
    cat(c("-> plotting", gene, "in all group\n"))
    png(paste(group,gene,".png"))
    visualizeGene(gene, betas = all, platform = platform) %>% print()
    dev.off()
  } else if(group == "as"){
    cat(c("-> plotting", gene, "in As group\n"))
    png(paste(group,gene,".png"))
    all %>% 
      select(ip_L_V_L_CON, ip_L_V_L_DMS, ip_L_V_L_BAP, ip_L_V_L_AS, ip_L_V_L_AS_BAP) %>% 
      visualizeGene(gene, betas = ., platform = platform) %>% print()
    dev.off()
  } else if(group == "co"){
    cat(c("-> plotting", gene, "in Co group\n"))
    png(paste(group,gene,".png"))
    all %>% 
      select(ip_L_V_L_CON, ip_L_V_L_DMS, ip_L_V_L_BAP, ip_L_V_L_CO, ip_L_V_L_CO_BAP) %>% 
      visualizeGene(gene, betas = ., platform = platform) %>% print()
    dev.off()
  } else if(group == "lcd"){
    cat(c("-> plotting", gene, "in L_Cd group\n"))
    png(paste(group,gene,".png"))
    all %>% 
      select(ip_L_V_L_CON, ip_L_V_L_DMS, ip_L_V_L_BAP, ip_L_V_L_LCD, ip_L_V_L_LCD_BAP) %>% 
      visualizeGene(gene, betas = ., platform = platform) %>% print()
    dev.off()
  } else if(group == "hcd"){
    cat(c("-> plotting", gene, "in H_Cd group\n"))
    png(paste(group,gene,".png"))
    all %>% 
      select(ip_L_V_L_CON, ip_L_V_L_DMS, ip_L_V_L_BAP, ip_L_V_L_HCD, ip_L_V_L_HCD_BAP) %>% 
      visualizeGene(gene, betas = ., platform = platform) %>% print()
    dev.off()
  } else{
    cat("error")
  }
}
map <- function(x = all, 
                group = all,
                gene = gene,
                platform = "EPIC"){
  if(group == "all"){
    cat(c("-> plotting", gene, "in all group\n"))
    visualizeGene(gene, betas = all, platform = platform) %>% print()
  } else if(group == "as"){
    cat(c("-> plotting", gene, "in As group\n"))
    all %>% 
      select(ip_L_V_L_CON, ip_L_V_L_DMS, ip_L_V_L_BAP, ip_L_V_L_AS, ip_L_V_L_AS_BAP) %>% 
      visualizeGene(gene, betas = ., platform = platform) %>% print()
  } else if(group == "co"){
    cat(c("-> plotting", gene, "in Co group\n"))
    all %>% 
      select(ip_L_V_L_CON, ip_L_V_L_DMS, ip_L_V_L_BAP, ip_L_V_L_CO, ip_L_V_L_CO_BAP) %>% 
      visualizeGene(gene, betas = ., platform = platform) %>% print()
  } else if(group == "lcd"){
    cat(c("-> plotting", gene, "in L_Cd group\n"))
    all %>% 
      select(ip_L_V_L_CON, ip_L_V_L_DMS, ip_L_V_L_BAP, ip_L_V_L_LCD, ip_L_V_L_LCD_BAP) %>% 
      visualizeGene(gene, betas = ., platform = platform) %>% print()
  } else if(group == "hcd"){
    cat(c("-> plotting", gene, "in H_Cd group\n"))
    all %>% 
      select(ip_L_V_L_CON, ip_L_V_L_DMS, ip_L_V_L_BAP, ip_L_V_L_HCD, ip_L_V_L_HCD_BAP) %>% 
      visualizeGene(gene, betas = ., platform = platform) %>% print()
  } else{
    cat("error")
  }
}
# use map function to plot gene on chromosome
map(x = all, group = "co", gene ="APEX1")

# chromosome visualization by Probe ID -----------------------------------
visualizeProbes(c("cg02382400", "cg03738669"), all, platform='EPIC')


# KnowYourCG Visualization -----------------------------------------------
# test the enrichment over database groups
results <- testEnrichment(rownames(beta)[1:50], platform = "EPIC")
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
# volcano plot
results_2tailed <- testEnrichment(rownames(beta)[1:50], "TFBS", alternative = "two.sided")
KYCG_plotVolcano(results_2tailed)
# Waterfall plot
KYCG_plotWaterfall(results)

# detail of database 
KYCG_listDBGroups("EPIC")
dbs <- KYCG_getDBs("MM285.design")
str(dbs[["PGCMeth"]])


# Gene Enrichment -------------------------------------------------------------
query <- names(sesameData_getProbesByGene("EGFR", "EPIC"))
results <- testEnrichment(query, 
                          KYCG_buildGeneDBs(query, max_distance=100000, platform="EPIC"),
                          platform="EPIC")
results[,c("dbname","estimate","gene_name","FDR", "nQ", "nD", "overlap")]
# lollipop plot
KYCG_plotLollipop(results, label="gene_name")
# volcano plot
KYCG_plotVolcano(results)


# GO/Pathway Enrichment ---------------------------------------------------
genes <- sesameData_getGenesByProbes(query, platform="EPIC")
genes
# perform Gene ontology enrichment analysis
library(gprofiler2)
## use gene name
gostres <- gost(genes$gene_name, organism = "hsapiens")
gostres$result[order(gostres$result$p_value),]
gostplot(gostres)

## use Ensembl gene ID, note we need to remove the version suffix
gene_ids <- sapply(strsplit(names(genes),"\\."), function(x) x[1])
gostres <- gost(gene_ids, organism = "hsapiens")
gostres$result[order(gostres$result$p_value),]
gostplot(gostres)


# Set Enrichment Analysis -------------------------------------------------
query <- KYCG_getDBs("KYCG.MM285.designGroup")[["TSS"]]
res <- testEnrichmentSEA(query, "MM285.seqContextN")
res[, c("dbname", "test", "estimate", "FDR", "nQ", "nD", "overlap")]

query <- KYCG_getDBs("KYCG.MM285.designGroup")[["TSS"]]
db <- KYCG_getDBs("MM285.seqContextN", "distToTSS")
res <- testEnrichmentSEA(aza, db, prepPlot = TRUE)
KYCG_plotSetEnrichment(res[[1]])






