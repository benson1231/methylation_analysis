# 0.1 ---------------------------------------------------------------------
AZA_0.1 <- delta %>% filter(abs(delta_AZA) > 0.1) %>% 
  select(delta_AZA, ip_L_V_L_CON, ip_L_V_L_AZA) %>% arrange(., desc(delta_AZA))
DAC_0.1 <- delta %>% filter(abs(delta_DAC) > 0.1) %>% 
  select(delta_DAC, ip_L_V_L_CON, ip_L_V_L_DAC) %>% arrange(., desc(delta_DAC))
AS_0.1 <- delta %>% filter(abs(delta_AS) > 0.1) %>% 
  select(delta_AS, ip_L_V_L_CON, ip_L_V_L_AS) %>% arrange(., desc(delta_AS))
CO_0.1 <- delta %>% filter(abs(delta_CO) > 0.1) %>% 
  select(delta_CO, ip_L_V_L_CON, ip_L_V_L_CO) %>% arrange(., desc(delta_CO))
LCD_0.1 <- delta %>% filter(abs(delta_LCD) > 0.1) %>% 
  select(delta_LCD, ip_L_V_L_CON, ip_L_V_L_LCD) %>% arrange(., desc(delta_LCD))
HCD_0.1 <- delta %>% filter(abs(delta_HCD) > 0.1) %>% 
  select(delta_HCD, ip_L_V_L_CON, ip_L_V_L_HCD) %>% arrange(., desc(delta_HCD))
BAP_0.1 <- delta %>% filter(abs(delta_BAP) > 0.1) %>% 
  select(delta_BAP, ip_L_V_L_DMS, ip_L_V_L_BAP) %>% arrange(., desc(delta_BAP))
AS_BAP_0.1 <- delta %>% filter(abs(delta_AS_BAP) > 0.1) %>% 
  select(delta_AS_BAP, ip_L_V_L_DMS, ip_L_V_L_AS_BAP) %>% arrange(., desc(delta_AS_BAP))
CO_BAP_0.1 <- delta %>% filter(abs(delta_CO_BAP) > 0.1) %>% 
  select(delta_CO_BAP, ip_L_V_L_DMS, ip_L_V_L_CO_BAP) %>% arrange(., desc(delta_CO_BAP))
LCD_BAP_0.1 <- delta %>% filter(abs(delta_LCD_BAP) > 0.1) %>% 
  select(delta_LCD_BAP, ip_L_V_L_DMS, ip_L_V_L_LCD_BAP) %>% arrange(., desc(delta_LCD_BAP))
HCD_BAP_0.1 <- delta %>% filter(abs(delta_HCD_BAP) > 0.1) %>% 
  select(delta_HCD_BAP, ip_L_V_L_DMS, ip_L_V_L_HCD_BAP) %>% arrange(., desc(delta_HCD_BAP))

probe_0.1 <- c(rownames(AZA_0.1), rownames(DAC_0.1), 
               rownames(AS_0.1), rownames(CO_0.1),
               rownames(LCD_0.1), rownames(HCD_0.1),
               rownames(BAP_0.1), rownames(AS_BAP_0.1),
               rownames(CO_BAP_0.1), rownames(LCD_BAP_0.1),
               rownames(HCD_BAP_0.1))
bet <- betas %>% na.omit()
probe_0.2_beta <- bet[rownames(bet)%in%probe_0.2,]
ComplexHeatmap::Heatmap(probe_0.2_beta, name = "Beta-value",
                        col = colorRamp2(c(0,0.1,1), c("blue","white", "yellow")),
                        na_col = "black", show_row_names = F,
                        column_title = "Groups", column_title_side = "bottom"
                        )

is_positive <- HCD_BAP_0.1[,1] > 0 
table(is_positive)
