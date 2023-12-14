# load package -----------------------------------------------------------------
library(sesame)
library(dplyr)
library(illuminaio)
library(ggplot2)
sesameDataCache()

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
agg <- readRDS("/Users/benson/Documents/project/methy1/agg.RDS")  # 差異表達probe
all <- readRDS("/Users/benson/Documents/project/methy1/all.RDS")  # 全部probe

# visualization -----------------------------------------------------------
# By Genomic Region
visualizeRegion(
  'chr19',10260000,10380000, all, platform='EPIC',
  show.probeNames = TRUE)

# By Gene Name
visualizeGene("APEX1", betas = all, platform='EPIC')

# By Probe ID
visualizeProbes(c("cg02382400", "cg03738669"), all, platform='EPIC')


### gene name output
gene_list <- c("APEX1","TP53","BRCA1")
for(i in gene_list){
  png(paste0(i,"_gene.png"), height = 800)
  visualizeGene(i, betas = all, platform='EPIC') %>% print()
  dev.off()
}


# GO analysis -------------------------------------------------------------
# GO
results <- testEnrichment(rownames(beta)[1:100], platform = "EPIC")


# dotplot
KYCG_plotDot(results)
# barplot
library(ggplot2)
library(wheatmap)
p1 <- KYCG_plotBar(results, label=TRUE)
p2 <- KYCG_plotBar(results, y="estimate") + ylab("log2(Odds Ratio)") +
  xlab("") + theme(axis.text.y = element_blank())
WGG(p1) + WGG(p2, RightOf(width=0.5, pad=0))

# 
results_pgc <- testEnrichment(betas$ID, platform="EPIC")
head(results_pgc)
KYCG_plotEnrichAll(results_pgc)





