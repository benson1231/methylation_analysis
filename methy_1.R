library(sesame)
library(dplyr)
library(illuminaio)
sesameDataCache()

betas <- sesameDataGet('HM450.10.TCGA.PAAD.normal')
visualizeRegion(
  'chr19',10260000,10380000, betas, platform='HM450',
  show.probeNames = FALSE)

visualizeGene('DNMT1', betas, platform='HM450')


df <- Co_0_1deltabeta$delta_Co %>% as.data.frame()
rownames(df) <- Co_0_1deltabeta$Nameill

visualizeRegion(
  'chr19',10260000,10380000, df, platform='EPIC',
  show.probeNames = FALSE)
visualizeGene('DNMT1', betas, platform='HM450')

sdf = sesameDataGet('EPIC.1.SigDF')
segs <- cnSegmentation(sdf)
visualizeSegments(segs)

library(SummarizedExperiment)
betas = assay(sesameDataGet("MMB.10.SE.tissue"))[,1:2]
compareReference(sesameDataGet("MM285.tissueSignature"), betas)

library(SummarizedExperiment)
df <- rowData(sesameDataGet('MM285.tissueSignature'))
query <- df$Probe_ID[df$branch == "fetal_liver" & df$type == "Hypo"]

results <- testEnrichment(Co_0_1deltabeta$probeID, platform = "EPIC")
KYCG_plotDot(results)
