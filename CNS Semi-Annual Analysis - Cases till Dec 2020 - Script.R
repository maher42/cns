## CCHE, Clinical Research Department

## Epidemiology, Histological Classification, Topographical
## Distribution and Clinical Outcome of Pediatric Brain Tumors

## Cases ~ Dec 2020
## FU ~ June 2021

# Import
library(readxl)
CNS <- read_excel("D:/Semiannual Analysis - Casess till end of 2020/Patients till end of 12-2020 - For Analysis.xlsx")

CNS$Extent <- factor(CNS$Extent, levels = c("GTR/NTR", "STR/Debulking", "Biopsy"))


# Slide 2 - Patient Presentation ------------------------------------------
# Average per year
nrow(CNS)/13.5

# Number of patients per year, Table
Years <- table(cut(CNS$date_of_registration, 'year'))

data.frame(Date=format(as.Date(names(Years)), '%Y'),
           Frequency=as.vector(Years))

# Survival
library(ggplot2); library(survival); library(survminer)
summary(survfit(Surv(CNS$OS, CNS$current_status) ~ 1, CNS), times = 60)
summary(survfit(Surv(CNS$EFS, CNS$event) ~ 1, CNS), times = 60)
table(CNS$current_status); table(CNS$event)

# Slide 3 - Diagnosis -----------------------------------------------------

table(CNS$Dx_P_R)
round(prop.table(table(CNS$Dx_P_R))*100, 1)

CNS_R <- subset(CNS, CNS$Dx_P_R == "R")
addmargins(table(CNS_R$pathology_gross))
round(prop.table(with(CNS_R, table(CNS_R$pathology_gross)))*100, 1)

# Slide 4 - Diagnosis by Age Group ----------------------------------------
# Table Pathology by Age Group
addmargins(table(CNS$pathology_gross, CNS$age_group))

library(psych); describe(CNS$age)
round(prop.table(with(CNS, table(CNS$age_group)))*100, 1)
length(which(CNS$age < 1))
(length(which(CNS$age < 1))/nrow(CNS))*100

# Slide 5 - Tumor Site ----------------------------------------------------
addmargins(table(CNS$site_gross))
round(prop.table(with(CNS, table(CNS$site_gross)))*100, 1)
Suprasellar <- subset(CNS, CNS$site_gross == "Suprasellar")
table(Suprasellar$site_systematic)

Spinal <- subset(CNS, CNS$site_systematic == "Spinal")
nrow(Spinal)
length(which(Spinal$pathology_gross == "Astrocytic tumors" | 
               Spinal$pathology_gross == "Ependymomas" |
               Spinal$pathology_gross == "Embryonal"))

# Slide 6 - Low Grade Glioma ----------------------------------------------
LGG <- subset(CNS, pathology_systematic == "Pilocytic, WHO GI" |
                pathology_systematic == "Pilomyxoid" |
                pathology_systematic == "Diffuse Astrocytoma, WHO GII" |
                pathology_systematic == "Ganglioglioma, WHO GI" |
                pathology_systematic == "Ganglioglioma, WHO GII" |
                pathology_systematic == "LGA" |
                pathology_systematic == "SEGA, WHO GI" |
                pathology_systematic == "DNET, WHO GI" |
                pathology_systematic ==
                "Pleomorphic Xantho astrocytoma, WHO GII" |
                pathology_systematic == "Desmoplastic Inf Astrocytoma, WHO GI" |
                pathology_systematic == "Oligoastrocytoma, WHO GII" |
                pathology_systematic == "Other Glioneuronal Tumor" |
                pathology_systematic == "Rosette forming glioneuronal tumor, WHO GI" |
                pathology_systematic == "Chordoid, WHO GII" |
                pathology_systematic == "Angio centric glioma, WHO G I." |
                pathology_systematic == "Gangliocytoma, WHO GI" |
                pathology_systematic == "Papillary Glioneuronal tumor" |
                pathology_systematic == "Ganglioglioma NOS" |
                pathology_systematic == "Oligodendroglioma, WHO GII")

addmargins(table(LGG$pathology_systematic))
round(prop.table(with(LGG, table(pathology_systematic)))*100, 1)

summary(survfit(Surv(LGG$OS, LGG$current_status) ~ 1, LGG), times = 60)
summary(survfit(Surv(LGG$EFS, LGG$event) ~ 1, LGG), times = 60)

# Slide 7 - Pilocytic Astrocytoma -----------------------------------------
Pilocytic <- subset(LGG, pathology_systematic == "Pilocytic, WHO GI")
survdiff(Surv(OS[!Extent=="Unknown cases resected outside CCHE"],
              current_status[!Extent=="Unknown cases resected outside CCHE"]) ~
           Extent[!Extent=="Unknown cases resected outside CCHE"], Pilocytic)
survdiff(Surv(EFS[!Extent=="Unknown cases resected outside CCHE"],
              event[!Extent=="Unknown cases resected outside CCHE"]) ~
           Extent[!Extent=="Unknown cases resected outside CCHE"], Pilocytic)

summary(survfit(Surv(OS, current_status) ~ Extent,
                Pilocytic), times = 60)
summary(survfit(Surv(EFS, event) ~ Extent,
                Pilocytic), times = 60)

library("survminer")
ggsurvplot(survfit(Surv(OS[!Extent=="Unknown cases resected outside CCHE"],
                        current_status[!Extent=="Unknown cases resected outside CCHE"]) ~
                     extent[!extent=="Unknown cases resected outside CCHE"],
                   Pilocytic), pval = TRUE,
           legend.labs = c("GTR / NTR", "STR", "Biopsy"),
           palette = "jco", risk.table = TRUE, xlab = "Time from registration (months)",
           ylab = "OS Probability")

ggsurvplot(survfit(Surv(EFS[!Extent=="Unknown cases resected outside CCHE"],
                        event[!Extent=="Unknown cases resected outside CCHE"]) ~
                     Extent[!Extent=="Unknown cases resected outside CCHE"],
                   Pilocytic), pval = TRUE,
           legend.labs = c("GTR / NTR", "STR", "Biopsy"),
           palette = "jco")


# Slide 8 - High Grade Glioma ---------------------------------------------
HGG <- subset(CNS, pathology_systematic == "Glioblastoma, WHO GIV" |
                pathology_systematic == "AA, WHO GIII" |
                pathology_systematic == "HGA" |
                pathology_systematic == "Gliomatosis cerebri, WHO GIII" |
                pathology_systematic ==
                "Anaplastic Pleomorphic Xantho astrocytoma, WHO GIII" |
                pathology_systematic == "Gliosarcoma, WHO GIV" |
                pathology_systematic == "Anaplastic Ganglioglioma, WHO GIII")

addmargins(table(HGG$pathology_systematic))
round(prop.table(with(HGG, table(pathology_systematic)))*100, 1)
summary(survfit(Surv(OS, current_status) ~ 1,
                HGG), times = 36)
summary(survfit(Surv(EFS, event) ~ 1, HGG), times = 36)
survfit(Surv(OS, current_status) ~ 1, HGG)

# Slide 9 - HGG - Survival and EoR ----------------------------------------
HGG_Extent <- subset(HGG, Extent == "GTR/NTR" | Extent == "STR/Debulking" |
                             Extent == "Biopsy")

survdiff(Surv(OS, current_status) ~ Extent, HGG_Extent)
survdiff(Surv(EFS, event) ~ Extent, HGG_Extent)
summary(survfit(Surv(OS, current_status) ~ Extent,
                HGG_Extent), times = 36)
summary(survfit(Surv(EFS, event) ~ Extent,
                HGG_Extent), times = 36)

ggsurvplot(survfit(Surv(OS, current_status) ~ extent, HGG_extent), pval = TRUE,
           legend.labs = c("GTR / NTR", "STR", "Biopsy"), palette = "jco",
           pval.coord = c(100,0.8), xlab = "Time from registration (months)",
           ylab = "OS Probability")

ggsurvplot(survfit(Surv(EFS, event) ~ extent, HGG_extent), pval = TRUE,
           legend.labs = c("GTR / NTR", "STR", "Biopsy"), palette = "jco",
           pval.coord = c(100,0.8), xlab = "Time from registration (months)",
           ylab = "EFS Probability")

# Slide 10 - Ependymoma ---------------------------------------------------
Ependymoma <- subset(CNS, pathology_gross == "Ependymomas")
addmargins(table(Ependymoma$pathology_systematic))
round(prop.table(with(Ependymoma, table(pathology_systematic)))*100, 1)

summary(survfit(Surv(OS, current_status) ~ 1,
                Ependymoma), times = 60)
summary(survfit(Surv(EFS, event) ~ 1, Ependymoma), times = 60)
Ependymoma_Ext <- subset(Ependymoma, Extent == "GTR/NTR" | Extent == "STR/Debulking")

### Not added bec. p = 0.7 and 0.8; to investigate
survdiff(Surv(OS, current_status) ~ Extent, Ependymoma_Ext)
survdiff(Surv(EFS, event) ~ Extent, Ependymoma_Ext)

# Slide 11 - Anaplastic Ependymoma ----------------------------------------
Anaplastic_Ependymoma <- subset(Ependymoma_Ext, pathology_systematic ==
                                  "Anaplastic Ependymoma, WHO GIII")
summary(survfit(Surv(OS, current_status) ~ Extent, Anaplastic_Ependymoma),
        times = 60)
summary(survfit(Surv(EFS, event) ~ Extent, Anaplastic_Ependymoma),
        times = 60)

survdiff(Surv(OS, current_status) ~ Extent, Anaplastic_Ependymoma)
survdiff(Surv(EFS, event) ~ Extent, Anaplastic_Ependymoma)

ggsurvplot(survfit(Surv(OS, current_status) ~ extent, Anaplastic_Ependymoma),
           pval = TRUE, legend.labs = c("GTR / NTR", "STR"),
           palette = "jco", xlab = "Time from registration (months)",
           ylab = "OS Probability", break.time.by = 24)

ggsurvplot(survfit(Surv(EFS, event) ~ extent, Anaplastic_Ependymoma),
           pval = TRUE, legend.labs = c("GTR / NTR", "STR"),
           palette = "jco", xlab = "Time from registration (months)",
           ylab = "EFS Probability", break.time.by = 24)

# Slide 12 & 13 - Embryonal Tumors: Medulloblastoma > 3 yo ---------------------------------------
MB <- subset(CNS, pathology_systematic ==
                     "Anaplastic/ Large Cell Medulloblastoma, WHO GIV" |
                     pathology_systematic == "Desmoplastic MB, WHO GIV" |
                     pathology_systematic ==
                     "MB with extensive nodularity, WHO GIV" |
                     pathology_systematic == "Classic MB, WHO GIV" |
                     pathology_systematic == "MB with focal nodularity, WHO GIV" |
                     pathology_systematic == "Medulloblastoma NOS, WHO GIV" |
                     pathology_systematic ==
                     "Medulloblastoma with neuronal differentiation, WHO GIV" |
                     pathology_systematic == "Melanotic MB, WHO GIV")
MB_child <- subset(MB, age >= 3)
addmargins(table(MB_child$pathology_systematic))
round(prop.table(with(MB_child, table(pathology_systematic)))*100, 1)

# Slide 14 - Medulloblastoma Risk Stratification --------------------------
round(prop.table(with(MB_child, table(protocol)))*100, 1)
table(MB_child$protocol)

# Slide 15 - Medulloblastoma Overall Survival (5 years) -------------------
MB_child_main3 <- subset(MB_child, pathology_systematic ==
                                        "Anaplastic/ Large Cell Medulloblastoma, WHO GIV" |
                                        pathology_systematic == "Classic MB, WHO GIV" |
                                 pathology_systematic == "Desmoplastic MB, WHO GIV")
summary(survfit(Surv(OS, current_status) ~ pathology_systematic,
                MB_child_main3), times = 60)
ggsurvplot(survfit(Surv(OS, current_status) ~ pathology_systematic, MB_child_main3),
           pval = TRUE, legend.labs = c("Anaplastic MB", "Classic MB", "Desmoplastic MB"),
           palette = "jco", xlab = "Time from registration (months)",
           ylab = "OS Probability", risk.table = TRUE, risk.table.y.text = FALSE)


# Slide 16 - Medulloblastoma EFS (5 years) -------------------------------
summary(survfit(Surv(EFS, event) ~ pathology_systematic,
                MB_child_main3), times = 60)
ggsurvplot(survfit(Surv(EFS, event) ~ pathology_systematic, MB_child_main3),
           pval = TRUE, legend.labs = c("Anaplastic MB", "Classic MB", "Desmoplastic MB"),
           palette = "jco", xlab = "Time from registration (months)",
           ylab = "EFS Probability", risk.table = TRUE, risk.table.y.text = FALSE)

# Slide 17 - Anaplastic MB vs Classic MB ----------------------------------
Anaplastic_vs_Classic <- subset(MB_child, pathology_systematic ==
                                  "Anaplastic/ Large Cell Medulloblastoma, WHO GIV" |
                                  pathology_systematic == "Classic MB, WHO GIV")

survdiff(Surv(OS, current_status) ~ pathology_systematic, Anaplastic_vs_Classic)
survdiff(Surv(EFS, event) ~ pathology_systematic, Anaplastic_vs_Classic)

SRClassic <- subset(MB_child, pathology_systematic ==
                      "Classic MB, WHO GIV" & protocol == "SRMB (CCHE 3-1-2008)")
nrow(SRClassic)

summary(survfit(Surv(OS, current_status) ~ 1, SRClassic), times = 60)
summary(survfit(Surv(EFS, event) ~ 1, SRClassic), times = 60)

# Slide 18 - Infantile Medulloblastoma < 3 yr -----------------------------
Infantile_MB <- subset(MB, age < 3)

addmargins(table(Infantile_MB$pathology_systematic))
round(prop.table(with(Infantile_MB, table(pathology_systematic)))*100, 1)

# Slide 19 - Infantile Medulloblastoma < 3 yr - OS ------------------------
summary(survfit(Surv(OS, current_status) ~ 1, Infantile_MB), times = 60)
summary(survfit(Surv(EFS, event) ~ 1, Infantile_MB), times = 60)
ggsurvplot(survfit(Surv(OS, current_status) ~ 1, Infantile_MB),
           risk.table = TRUE, legend = "none",
           palette = "darkgrey", break.time.by = 12,
           risk.table.pos = "in", ylab = "OS Probability",
           xlab = "Months from registration & N at risk")

# Slide 20 - Infantile Medulloblastoma OS by Pathology --------------------
Infantile_MB_main3 <- subset(Infantile_MB, pathology_systematic ==
                                 "Anaplastic/ Large Cell Medulloblastoma, WHO GIV" |
                                 pathology_systematic == "Classic MB, WHO GIV" |
                                 pathology_systematic == "Desmoplastic MB, WHO GIV")

summary(survfit(Surv(OS, current_status) ~ pathology_systematic,
                Infantile_MB_main3), times = 36)
survdiff(Surv(OS, current_status) ~ pathology_systematic, Infantile_MB_main3)
ggsurvplot(survfit(Surv(OS, current_status) ~ pathology_systematic, Infantile_MB_main3),
           pval = TRUE, legend.labs = c("Anaplastic MB", "Classic MB", "Desmoplastic MB"),
           palette = "jco", xlab = "Time from registration (months)", break.time.by = 24,
           ylab = "OS Probability", risk.table = TRUE, risk.table.y.text = FALSE)

# Slide 21 - Infantile Medulloblastoma EFS --------------------------------
summary(survfit(Surv(EFS, event) ~ pathology_systematic,
                Infantile_MB_main3), times = 36)
survdiff(Surv(EFS, event) ~ pathology_systematic, Infantile_MB_main3)
ggsurvplot(survfit(Surv(EFS, event) ~ pathology_systematic, Infantile_MB_main3),
           pval = TRUE, legend.labs = c("Anaplastic MB", "Classic MB", "Desmoplastic MB"),
           palette = "jco", xlab = "Time from registration (months)", break.time.by = 24,
           ylab = "EFS Probability", risk.table = TRUE, risk.table.y.text = FALSE)

# Slide 22 - Medulloblastoma by Molecular Subgroup ------------------------
addmargins(table(MB$pathology_systematic, MB$MB_subgroup))
Seq_MB <- subset(MB, MB_subgroup == "WNT" | MB_subgroup == "SHH" |
                         MB_subgroup == "Group 3 or 4")
Seq_infMB <- subset(Seq_MB, age < 3)
table(Seq_infMB$MB_subgroup)

# Slide 23 - ATRT ---------------------------------------------------------
ATRT <- subset(CNS, pathology_systematic == "ATRT, WHO GIV")
length(which(ATRT$age >= 3)); length(which(ATRT$age < 3))
length(which(ATRT$age < 1))

# Slide 24 - ATRT Site ----------------------------------------------------
addmargins(table(ATRT$site_gross))

# Slide 25 - ATRT  OS & EFS -----------------------------------------------
summary(survfit(Surv(OS, current_status) ~ 1, ATRT), times = 12)
summary(survfit(Surv(EFS, event) ~ 1, ATRT), times = 12)
ggsurvplot(survfit(Surv(OS, current_status) ~ pathology_systematic, ATRT),
           palette = "jco", xlab = "Months from Registration & N at risk", break.time.by = 12,
           ylab = "OS Probability", risk.table = TRUE, risk.table.pos = "in",
           legend = "none")
ggsurvplot(survfit(Surv(EFS, event) ~ pathology_systematic, ATRT),
           palette = "jco", xlab = "Months from Registration & N at risk", break.time.by = 12,
           ylab = "EFS Probability", risk.table = TRUE, risk.table.pos = "in",
           legend = "none")

# Slide 26 - Other Embryonal Tumors ---------------------------------------
OET <-  subset(CNS, pathology_systematic=="PNET, WHO GIV"|
                       pathology_systematic=="Embryonal tumor, (NOS)"|
                       pathology_systematic=="Medulloepithelioma, WHO GIV"|
                       pathology_systematic=="ETMR, WHO GIV"|
                       pathology_systematic=="Ependymoblastoma, WHO GIV" |
                       pathology_systematic=="ETANTR, WHO GIV")
addmargins(table(OET$pathology_systematic, OET$current_status))
summary(survfit(Surv(OS, current_status) ~ 1, OET), times = 36)
summary(survfit(Surv(EFS, event) ~ 1, OET), times = 36)

ggsurvplot(survfit(Surv(OS, current_status) ~ 1, OET),
           xlab = "OS", ylab = "", legend="none", palette = "darkgrey")

# Slide 27 - Tumors of the Pineal Region ----------------------------------
Pineal <-  subset(CNS, pathology_gross=="Tumors of the Pineal region")

addmargins(table(Pineal$pathology_systematic))
round(prop.table(with(Pineal, table(pathology_systematic)))*100, 1)

Pineoblastoma <- subset(CNS, pathology_systematic == "Pineoblastoma, WHO GIV")
summary(survfit(Surv(OS, current_status) ~ 1, Pineoblastoma), times = 36)
summary(survfit(Surv(EFS, event) ~ 1, Pineoblastoma), times = 36)

# Slide 28 - Intracranial Germ Cell Tumors --------------------------------
Germinoma <- subset(CNS, pathology_systematic == "Germinoma" & protocol=="RTH")
Non_Germinoma <- subset(CNS, pathology_systematic == "Immature Teratoma" |
                                pathology_systematic == "Mixed GCT" |
                                pathology_systematic == "Other NGGCT" |
                                pathology_systematic == "Yolk Sac Tumors" |
                                pathology_systematic == "Embryonal Carcinoma" |
                                pathology_systematic == "Germinoma" & protocol !="RTH")
nrow(Germinoma); nrow(Non_Germinoma)
table(Non_Germinoma$pathology_systematic)

# Slide 29 - Non-Germinoma ------------------------------------------------
addmargins(table(Non_Germinoma$current_status))
addmargins(table(Non_Germinoma$protocol))
summary(survfit(Surv(OS, current_status) ~ 1, Non_Germinoma), times = 60)
summary(survfit(Surv(EFS, event) ~ 1, Non_Germinoma), times = 60)
ggsurvplot(survfit(Surv(OS, current_status) ~ 1, Non_Germinoma),
           xlab = "OS", ylab = "", legend="none", palette = "darkgrey")

# Slide 30 - Germinoma ----------------------------------------------------
addmargins(table(Germinoma$protocol, Germinoma$current_status))
summary(survfit(Surv(OS, current_status) ~ 1, Germinoma), times = 60)
summary(survfit(Surv(EFS, event) ~ 1, Germinoma), times = 60)
ggsurvplot(survfit(Surv(OS, current_status) ~ 1, Germinoma),
           xlab = "OS", ylab = "", legend="none", palette = "darkgrey")
ggsurvplot(survfit(Surv(EFS, event) ~ 1, Germinoma),
           xlab = "EFS", ylab = "", legend="none", palette = "darkgrey")

# Slide 31 - Craniopharyngioma --------------------------------------------
Craniopharyngioma <- subset(CNS, pathology_systematic ==
                                    "Craniopharyngioma, WHO GI")
addmargins(table(Craniopharyngioma$Dx_P_R))
addmargins(table(Craniopharyngioma$Extent))
Craniopharyngioma <- subset(Craniopharyngioma, Dx_P_R == "P")
summary(survfit(Surv(EFS, event) ~ 1, Craniopharyngioma),
        times = 60)
survdiff(Surv(OS, current_status) ~ Extent, Craniopharyngioma)
survdiff(Surv(EFS, event) ~ Extent, Craniopharyngioma)

Craniopharyngioma <- subset(CNS, pathology_systematic ==
                                    "Craniopharyngioma, WHO GI")
table(Craniopharyngioma$rth)
# survdiff(Surv(OS, current_status) ~ rth, Craniopharyngioma)
# survdiff(Surv(EFS, event) ~ rth, Craniopharyngioma)

# Slide 32 - Choroid Plexus Tumors ----------------------------------------
CPT <- subset(CNS, pathology_gross == "Choroid Plexus Tumors")
addmargins(table(CPT$pathology_systematic, CPT$current_status))
addmargins(table(CPT$pathology_systematic, CPT$protocol))

# Slide 33 - Age Distribution of CPT Tumors -------------------------------
length(which(CPT$age < 3))

# Slide 34 - Choroid Plexus Tumors OS & EFS -------------------------------
summary(survfit(Surv(OS, current_status) ~ 1, CPT), times = 60)
summary(survfit(Surv(EFS, event) ~ 1, CPT), times = 60)
ggsurvplot(survfit(Surv(OS, current_status) ~ 1, CPT),
           xlab = "OS", ylab = "", legend="none", palette = "darkgrey")
ggsurvplot(survfit(Surv(EFS, event) ~ 1, Germinoma),
           xlab = "EFS", ylab = "", legend="none", palette = "darkgrey")

# Slide 35 - CPT: Protocol & Current Status -------------------------------
addmargins(table(CPT$protocol, CPT$current_status))

# Slide 36 & 37 - Radiologically Confirmed Diagnoses:  Optic Pathway Glioma -----------------------------------------
nrow(OPG <- subset(CNS, site_gross == "OP" | site_systematic == "Suprasellar - OP"))
table(OPG$genetic_disorder)
table(OPG$current_status);table(OPG$event);table(OPG$protocol)
summary(survfit(Surv(OS, current_status) ~ 1, OPG), times = 60)
summary(survfit(Surv(EFS, event) ~ 1, OPG), times = 60)
ggsurvplot(survfit(Surv(OS, current_status) ~ 1, OPG),
           xlab = "OS", ylab = "", legend="none", palette = "darkgrey")
ggsurvplot(survfit(Surv(EFS, event) ~ 1, OPG),
           xlab = "EFS", ylab = "", legend="none", palette = "darkgrey")

# Slide 38 - Brain Stem Lesions -------------------------------------------
nrow(BSG <- subset(CNS, site_gross == "Brain Stem"))
table(BSG$Dx_P_R); table(BSG$site_systematic)
DIPG <- subset(BSG, site_systematic == "Brain Stem - DIPG")
FBSG <- subset(BSG, site_systematic == "Brain Stem - Focal")
survfit(Surv(OS, current_status) ~ 1, DIPG)
summary(survfit(Surv(OS, current_status) ~ 1, BSG), times = 12)
ggsurvplot(survfit(Surv(OS, current_status) ~ 1, BSG),
           xlab = "OS", ylab = "", legend="none", palette = "darkgrey",
           title = "Brain Stem Glioma: 1-year OS = 46.7%")

summary(survfit(Surv(OS, current_status) ~ 1, BSG), times = 12)

# Slide 39 - Brain Stem Lesions by Pathology ------------------------------
BSG_p <- subset(BSG, Dx_P_R=="P")
nrow(BSG_p)/nrow(BSG)*100
addmargins(table(BSG_p$pathology_systematic),1)

table(BSG_p$site_systematic)

# Slides 40 & 41 - Protocol Survival Rates --------------------------------
lgg_protocol <- subset(CNS, protocol=="LGG (CCHE 3-1-2008)" & site_gross != "Brain Stem")
srmb_protocol <- subset(CNS, protocol=="SRMB (CCHE 3-1-2008)" & age >=3)
infmb_protocol <- subset(CNS, protocol=="InfMB (CCHE 3-1-2008)" & off_protocol=="FALSE")
epnd_protocol <- subset(CNS, protocol=="EPND (CCHE 3-1-2008)")
hrmb_protocol <- subset(CNS, protocol=="HRMB (CCHE 3-1-2008)" & age>=3)
ngcct_protocol <- subset(CNS, protocol=="NGGCT (CCHE 3-1-2008)")
cpt_protocol <- subset(CNS, protocol=="CPT protocol")
hgg_protocol <- subset(CNS, protocol=="HGG (CCHE 3-1-2009)")
atrt_protocol <- subset(CNS, protocol=="ATRT (CCHE 3-1-2010)")

protocols <- rbind(lgg_protocol,srmb_protocol,infmb_protocol,epnd_protocol,
                   hrmb_protocol,ngcct_protocol,
                   cpt_protocol,hgg_protocol,atrt_protocol)

addmargins(table(protocols$protocol, protocols$current_status))
addmargins(table(protocols$protocol, protocols$event))

summary(survfit(Surv(OS, current_status) ~ protocol, protocols),times = 60)
summary(survfit(Surv(EFS, event) ~ protocol,protocols),times = 60)
