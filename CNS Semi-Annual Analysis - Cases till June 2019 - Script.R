## CCHE, Clinical Research Department

## Epidemiology, Histological Classification, Topographical
## Distribution and Clinical Outcome of Pediatric Brain Tumors

## Cases ~ June 2019
## FU ~ Dec 2019

# Import
library(readxl)
CNS <- read_excel("CNS analysis - cases till mid 2019/CNS semi-annual analysis - Cases till June 2019 - Cleaned for analysis.xlsx")

# Average per year
nrow(CNS)/12

# Number of patients per year # Table (1)
Years <- table(cut(CNS$date_of_registration, 'year'))

data.frame(Date=format(as.Date(names(Years)), '%Y'),
           Frequency=as.vector(Years))

# Survival
library(ggplot2); library(survival); library(survminer)
summary(survfit(Surv(CNS$OS, CNS$current_status) ~ 1, CNS), times = 60)
summary(survfit(Surv(CNS$EFS, CNS$event) ~ 1, CNS), times = 60)

# Diagnosis
table(CNS$Dx_P_R)

# Gender
table(CNS$gender)

# Figure (1)
ggplot(CNS, aes(x="M", y="F", fill=gender))+geom_bar(width = 1,
                                                     stat = "identity")+
  coord_polar("y", start=0)+scale_fill_grey()+theme_minimal()+
  xlab(label="")+ylab(label="")

# Table (2)
addmargins(table(CNS$pathology_gross, CNS$gender))

# Age
# Less than 1 year
length(which(CNS$age < 1))
(length(which(CNS$age < 1))/nrow(CNS))*100

library(psych); describe(CNS$age)

# Table (3)
round(prop.table(with(CNS, table(CNS$age_group)))*100, 1)
table(CNS$age_group)

# Geographical Distribution
round(prop.table(with(CNS, table(CNS$geo_area)))*100, 1)
#Figure (2)
CNS <- within(CNS, geo_area <- factor(geo_area, 
                                      levels=names(sort(table(geo_area), 
                                                        decreasing=TRUE))))
ggplot(CNS, aes(x=geo_area, y=frequency(geo_area))) +
  geom_bar(stat="identity")+xlab(label="")+ylab(label="")

# Detailed Pathology
CNS_P <- subset(CNS, Dx_P_R == "P")
# Table (4)
addmargins(table(CNS_P$pathology_systematic, CNS_P$gender))

# Table(5)
addmargins(table(CNS_P$pathology_systematic, CNS_P$age_group))

Markers <- subset(CNS, CNS$Dx_P_R == "R+T")
table(Markers$site_systematic)

# Table (6)
CNS_R <- subset(CNS, CNS$Dx_P_R == "R")
addmargins(table(CNS_R$pathology_systematic, CNS_R$gender))

# Table (7)
addmargins(table(CNS_R$pathology_systematic, CNS_R$age_group))

# Tumor Site
# Table (8)
addmargins(table(CNS$site_systematic))
round(prop.table(with(CNS, table(CNS$site_systematic)))*100, 2)
Suprasellar <- subset(CNS, CNS$site_systematic == "Suprasellar")
table(Suprasellar$site_gross)
round(prop.table(with(Suprasellar, table(Suprasellar$site_gross)))*100, 1)

Spinal <- subset(CNS, CNS$site_systematic == "Spinal")
length(which(Spinal$pathology_gross == "Astrocytic tumors" | 
               Spinal$pathology_gross == "Ependymomas" |
               Spinal$pathology_gross == "Embryonal"))

## Pathology Confirmed Diagnoses
# Glioma
# LGG
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

# Table(9)
addmargins(table(LGG$pathology_systematic))
round(prop.table(with(LGG, table(pathology_systematic)))*100, 1)

summary(survfit(Surv(LGG$OS, LGG$current_status) ~ 1, LGG), times = 60)
summary(survfit(Surv(LGG$EFS, LGG$event) ~ 1, LGG), times = 60)

# Pilocytic Astrocytoma
Pilocytic <- subset(LGG, pathology_systematic == "Pilocytic, WHO GI")
survdiff(Surv(OS[!extent==4], current_status[!extent==4]) ~
           extent[!extent==4], Pilocytic)
survdiff(Surv(EFS[!extent==4], event[!extent==4]) ~
           extent[!extent==4], Pilocytic)

summary(survfit(Surv(OS, current_status) ~ extent,
                Pilocytic), times = 60)
summary(survfit(Surv(EFS, event) ~ extent,
                Pilocytic), times = 60)
# Figure (3)
library("survminer")
ggsurvplot(survfit(Surv(OS[!extent==4], current_status[!extent==4]) ~
                     extent[!extent==4], Pilocytic), pval = TRUE,
           legend.labs = c("GTR / NTR", "STR", "Biopsy"))
# Figure (4)
ggsurvplot(survfit(Surv(EFS[!extent==4], event[!extent==4]) ~
                     extent[!extent==4], Pilocytic), pval = TRUE,
           legend.labs = c("GTR / NTR", "STR", "Biopsy"))

# Diffuse Astrocytoma
Diffuse <- subset(LGG, pathology_systematic == "Diffuse Astrocytoma, WHO GII")
survdiff(Surv(OS[!extent==4], current_status[!extent==4]) ~
           extent[!extent==4], Diffuse)
survdiff(Surv(EFS[!extent==4], event[!extent==4]) ~
           extent[!extent==4], Diffuse)

summary(survfit(Surv(OS, current_status) ~ extent,
                Diffuse), times = 60)
summary(survfit(Surv(EFS, event) ~ extent,
                Diffuse), times = 60)

# HGG
HGG <- subset(CNS, pathology_systematic == "Glioblastoma, WHO GIV" |
                pathology_systematic == "AA, WHO GIII" |
                pathology_systematic == "HGA" |
                pathology_systematic == "Gliomatosis cerebri, WHO GIII" |
                pathology_systematic ==
                "Anaplastic Pleomorphic Xantho astrocytoma, WHO GIII" |
                pathology_systematic == "Gliosarcoma, WHO GIV" |
                pathology_systematic == "Anaplastic Ganglioglioma, WHO GIII")

# Table(10)
addmargins(table(HGG$pathology_systematic))
round(prop.table(with(HGG, table(pathology_systematic)))*100, 1)
HGG <- subset(HGG, extent==1 | extent==2 | extent==3)
summary(survfit(Surv(OS, current_status) ~ 1,
                HGG), times = 36)
summary(survfit(Surv(EFS, event) ~ 1, HGG), times = 36)

survdiff(Surv(OS, current_status) ~ extent, HGG)
survdiff(Surv(EFS, event) ~ extent, HGG)
summary(survfit(Surv(OS, current_status) ~ extent,
                HGG), times = 36)
summary(survfit(Surv(EFS, event) ~ extent,
                HGG), times = 36)
# Figure (5)
ggsurvplot(survfit(Surv(OS, current_status) ~ extent, HGG), pval = TRUE,
           legend.labs = c("GTR / NTR", "STR", "Biopsy"))
# Figure (6)
ggsurvplot(survfit(Surv(EFS, event) ~ extent, HGG), pval = TRUE,
           legend.labs = c("GTR / NTR", "STR", "Biopsy"))

GBM_vs_AA <- subset(HGG, pathology_systematic == "Glioblastoma, WHO GIV" |
                      pathology_systematic == "AA, WHO GIII")
survdiff(Surv(OS, current_status) ~ pathology_systematic, GBM_vs_AA)
survfit(Surv(OS, current_status) ~ pathology_systematic, GBM_vs_AA)

# Ependymoma
Ependymoma <- subset(CNS, pathology_gross == "Ependymomas")
# Table (11)
addmargins(table(Ependymoma$pathology_systematic))
round(prop.table(with(Ependymoma, table(pathology_systematic)))*100, 1)
Ependymoma <- subset(Ependymoma, extent==1 | extent==2 | extent==3)
summary(survfit(Surv(OS, current_status) ~ 1,
                Ependymoma), times = 60)
summary(survfit(Surv(EFS, event) ~ 1, Ependymoma), times = 60)

survdiff(Surv(OS, current_status) ~ extent, Ependymoma)
survdiff(Surv(EFS, event) ~ extent, Ependymoma)

Anaplastic_Ependymoma <- subset(Ependymoma, pathology_systematic ==
                                  "Anaplastic Ependymoma, WHO GIII")
summary(survfit(Surv(OS, current_status) ~ extent, Anaplastic_Ependymoma),
        times = 60)
summary(survfit(Surv(EFS, event) ~ extent, Anaplastic_Ependymoma),
        times = 60)

survdiff(Surv(OS, current_status) ~ extent, Anaplastic_Ependymoma)
survdiff(Surv(EFS, event) ~ extent, Anaplastic_Ependymoma)

# Figure (7)
ggsurvplot(survfit(Surv(OS[!extent==3], current_status[!extent==3]) ~
                     extent[!extent==3], Anaplastic_Ependymoma),
           pval = TRUE, legend.labs = c("GTR / NTR", "STR"))
# Figure (8)
ggsurvplot(survfit(Surv(EFS[!extent==3], event[!extent==3]) ~
                     extent[!extent==3], Anaplastic_Ependymoma),
           pval = TRUE, legend.labs = c("GTR / NTR", "STR"))

# Craniopharyngioma
Craniopharyngioma <- subset(CNS, pathology_systematic ==
                              "Craniopharyngioma, WHO GI")
addmargins(table(Craniopharyngioma$Dx_P_R))
Craniopharyngioma <- subset(Craniopharyngioma, Dx_P_R == "P")
summary(survfit(Surv(OS, current_status) ~ 1, Craniopharyngioma),
        times = 60)
summary(survfit(Surv(EFS, event) ~ 1, Craniopharyngioma),
        times = 60)
addmargins(table(Craniopharyngioma$extent))
survdiff(Surv(OS, current_status) ~ extent, Craniopharyngioma)
survdiff(Surv(EFS, event) ~ extent, Craniopharyngioma)
Craniopharyngioma <- subset(CNS, pathology_systematic ==
                              "Craniopharyngioma, WHO GI")
table(Craniopharyngioma$Radiotherapy)
survdiff(Surv(OS, current_status) ~ Radiotherapy, Craniopharyngioma)
survdiff(Surv(EFS, event) ~ Radiotherapy, Craniopharyngioma)

# Embryonal Tumors
# ATRT
ATRT <- subset(CNS, pathology_systematic == "ATRT, WHO GIV")
addmargins(table(ATRT$protocol, ATRT$current_status))
addmargins(table(ATRT$site_systematic))
length(which(ATRT$age < 1))
length(which(ATRT$age < 5)); length(which(ATRT$age > 5))
summary(survfit(Surv(OS, current_status) ~ 1, ATRT),
        times = 12)
summary(survfit(Surv(EFS, event) ~ 1, ATRT),
        times = 12)

# Medulloblastoma < 3 yo
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
Infantile_MB <- subset(MB, age < 3)
# Table (13)
addmargins(table(Infantile_MB$pathology_systematic))
round(prop.table(with(Infantile_MB, table(pathology_systematic)))*100, 2)
summary(survfit(Surv(OS, current_status) ~ 1, Infantile_MB), times = 60)
# Figure (9)
ggsurvplot(survfit(Surv(OS, current_status) ~ 1, Infantile_MB))

# Table (14)
addmargins(table(Infantile_MB$protocol, Infantile_MB$current_status))

# Medulloblastoma > 3 yo
MB_child <- subset(MB, age >= 3)
# Table (15)
addmargins(table(MB_child$pathology_systematic,
                 MB_child$protocol))

summary(survfit(Surv(OS, current_status) ~ pathology_systematic,
                MB_child), times = 60)
summary(survfit(Surv(EFS, event) ~ pathology_systematic,
                MB_child), times = 60)

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

# Medulloblastoma by Molecular Subgroup (Table not numbered)
addmargins(table(MB$pathology_systematic, MB$MB_subgroup))
Seq_MB <- subset(MB, MB_subgroup == "WNT" | MB_subgroup == "SHH" |
                   MB_subgroup == "Group 3 or 4")
Seq_infMB <- subset(Seq_MB, age < 3)
table(Seq_infMB$MB_subgroup)

# Pineoblastoma
Pineoblastoma <- subset(CNS, pathology_systematic == "Pineoblastoma, WHO GIV")
addmargins(table(Pineoblastoma$protocol, Pineoblastoma$current_status))

summary(survfit(Surv(OS, current_status) ~ 1, Pineoblastoma), times = 36)
summary(survfit(Surv(EFS, event) ~ 1, Pineoblastoma), times = 36)

# PNET
PNET <- subset(CNS, pathology_systematic == "PNET, WHO GIV")
# Table (17)
addmargins(table(PNET$protocol, PNET$current_status))

length(which(PNET$age < 3))
summary(survfit(Surv(OS, current_status) ~ 1, PNET), times = 36)
summary(survfit(Surv(EFS, event) ~ 1, PNET), times = 36)

# Choroid Plexus Tumors
CPT <- subset(CNS, pathology_gross == "Choroid Plexus Tumors")
# Table (18) & (19)
addmargins(table(CPT$pathology_systematic, CPT$current_status))
addmargins(table(CPT$protocol, CPT$current_status))

length(which(CPT$age < 3))
addmargins(table(CPT$protocol, CPT$pathology_systematic))

# Intracranial Germ Cell Tumors
Germinoma <- subset(CNS, pathology_systematic == "Germinoma")
Non_Germinoma <- subset(CNS, pathology_systematic == "Immature Teratoma" |
                          pathology_systematic == "Mixed GCT" |
                          pathology_systematic == "Other NGGCT" |
                          pathology_systematic == "Yolk Sac Tumors" |
                          pathology_systematic == "Embryonal Carcinoma")
GCT <- rbind(Germinoma, Non_Germinoma)
nrow(GCT); nrow(Germinoma); nrow(Non_Germinoma)
table(GCT$Dx_P_R)
table(Non_Germinoma$pathology_systematic)

# Table (20)
addmargins(table(Non_Germinoma$protocol, Non_Germinoma$current_status))

# Table (21)
addmargins(table(Germinoma$protocol, Germinoma$current_status))

table(Germinoma$site_systematic)
table(Germinoma$Dx_P_R)

## Radiologically Confirmed Diagnoses
# Brain Stem Lesions
nrow(BSG <- subset(CNS, site_gross == "Brain Stem"))
table(BSG$Dx_P_R); table(BSG$site_systematic)
DIPG <- subset(BSG, site_systematic == "Brain Stem - DIPG")
FBSG <- subset(BSG, site_systematic == "Brain Stem - Focal")
survfit(Surv(OS, current_status) ~ 1, DIPG)
survfit(Surv(OS, current_status) ~ 1, FBSG)
summary(survfit(Surv(OS, current_status) ~ 1, BSG), times = 12)


# Optic Pathway Glioma
nrow(OPG <- subset(CNS, site_gross == "OP" | site_gross == "Suprasellar - OP"))
table(OPG$genetics); table(OPG$current_status)
table(OPG$event);table(OPG$protocol)

## Protocol Survival Rates
# # # # # # #                            # # # # # # #
# *** The following section will be updated soon ***
# # # # # # #                            # # # # # # #

# 5-year Overall Survival
# Table (22)
summary(survfit(Surv(OS, current_status) ~ protocol=="LGG (CCHE 3-1-2008)" &
                  site_gross != "Brain Stem",
                CNS), times = 60)
survfit(Surv(OS, current_status) ~ protocol=="LGG (CCHE 3-1-2008)" &
          site_gross != "Brain Stem", CNS)
summary(survfit(Surv(OS, current_status) ~ protocol=="SRMB (CCHE 3-1-2008)", CNS),
        times = 60)
summary(survfit(Surv(OS, current_status) ~ protocol=="HRMB (CCHE 3-1-2008)", CNS),
        times = 60)
summary(survfit(Surv(OS, current_status) ~ protocol=="EPND (CCHE 3-1-2008)", CNS),
        times = 60)
summary(survfit(Surv(OS, current_status) ~ protocol=="NGGCT (CCHE 3-1-2008)", CNS),
        times = 60)
summary(survfit(Surv(OS, current_status) ~ protocol=="HGG (CCHE 3-1-2009)", CNS),
        times = 60)
print(survfit(Surv(OS, current_status) ~ protocol, CNS))

# 5-year Event-free Survival
# Table (23)
summary(survfit(Surv(EFS, event) ~ protocol=="LGG (CCHE 3-1-2008)" &
                  site_gross != "Brain Stem",
                CNS), times = 60)
survfit(Surv(EFS, event) ~ protocol=="LGG (CCHE 3-1-2008)" &
          site_gross != "Brain Stem", CNS)
summary(survfit(Surv(EFS, event) ~ protocol=="SRMB (CCHE 3-1-2008)", CNS),
        times = 60)
summary(survfit(Surv(EFS, event) ~ protocol=="HRMB (CCHE 3-1-2008)", CNS),
        times = 60)
summary(survfit(Surv(EFS, event) ~ protocol=="EPND (CCHE 3-1-2008)", CNS),
        times = 60)
summary(survfit(Surv(EFS, event) ~ protocol=="NGGCT (CCHE 3-1-2008)", CNS),
        times = 60)
summary(survfit(Surv(EFS, event) ~ protocol=="HGG (CCHE 3-1-2009)", CNS),
        times = 60)
print(survfit(Surv(EFS, event) ~ protocol, CNS))

# Infantile Medulloblastoma Protocol (OS & EFS); Posterior Fossa & M0
EligibileINFMB <- subset(CNS, inf_protocol_inclusion == 1)

summary(survfit(Surv(OS, current_status) ~ protocol, EligibileINFMB), times = 60)
survfit(Surv(OS, current_status) ~ protocol, EligibileINFMB)

summary(survfit(Surv(EFS, event) ~ protocol, EligibileINFMB), times = 60)
survfit(Surv(EFS, event) ~ protocol, EligibileINFMB)
