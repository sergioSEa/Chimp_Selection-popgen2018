library(qqman)

### Clean data: remove nans, convert negs to 0 and 1 to 0.9999

sch_trogl = read.table('SchweinTrogl.fst',header=TRUE)
# sch_trogl = as.matrix(sch_trogl)

sch_verus = read.table('SchweinVerus.fst',header=TRUE)
# sch_verus = as.matrix(sch_verus)

trogl_verus = read.table('TroglVerus.fst',header=TRUE)
# trogl_verus = as.matrix(trogl_verus)

## FILTER OUT NANS
#Positions to filterout on each dataframe
sch_verus_NaN = sch_verus[is.nan(sch_verus$FST), ]$SNP
sch_trogl_NaN = sch_trogl[is.nan(sch_trogl$FST), ]$SNP
trogl_verus_NaN = trogl_verus[is.nan(trogl_verus$FST), ]$SNP

sch_verus = sch_verus[!(sch_verus$SNP %in% sch_verus_NaN), ]
sch_verus = sch_verus[!(sch_verus$SNP %in% sch_trogl_NaN), ]
sch_verus = sch_verus[!(sch_verus$SNP %in% trogl_verus_NaN), ]

sch_trogl = sch_trogl[!(sch_trogl$SNP %in% sch_verus_NaN), ]
sch_trogl = sch_trogl[!(sch_trogl$SNP %in% sch_trogl_NaN), ]
sch_trogl = sch_trogl[!(sch_trogl$SNP %in% trogl_verus_NaN), ]

trogl_verus = trogl_verus[!(trogl_verus$SNP %in% sch_verus_NaN), ]
trogl_verus = trogl_verus[!(trogl_verus$SNP %in% sch_trogl_NaN), ]
trogl_verus = trogl_verus[!(trogl_verus$SNP %in% trogl_verus_NaN), ]

# Convert negs to 0
sch_trogl$FST[sch_trogl$FST < 0] = 0
sch_verus$FST[sch_verus$FST < 0] = 0
trogl_verus$FST[trogl_verus$FST < 0] = 0

# Convert 1 to 0.9999
sch_trogl$FST[sch_trogl$FST == 1] = 0.99999
sch_verus$FST[sch_verus$FST == 1] = 0.99999
trogl_verus$FST[trogl_verus$FST == 1] = 0.99999

# CALCULATE T VALUES
sch_trogl$T = -log10(1 - sch_trogl$FST)
sch_verus$T = -log10(1 - sch_verus$FST)
trogl_verus$T = -log10(1 - trogl_verus$FST)

PBS = (sch_trogl$T + trogl_verus$T - sch_verus$T)/2
PBSdata = data.frame(PBS=PBS, CHR=sch_trogl$CHR, BP=sch_trogl$SNP)

# PLOT MANHATTAN
manhattan(PBSdata, p = "PBS", logp = FALSE, ylab = "PBS", genomewideline = FALSE, 
          suggestiveline = FALSE, main = "Manhattan plot of PBS")

PBS_relevant = PBSdata[PBSdata$PBS > 1.5, ]



# #### CHECK EVERY SNP IS EQUAL
# resultST = sch_trogl[,2] == sch_verus[,2]
# sum(resultST)/length(resultST)
# 
# resultSV = sch_verus[,2] == sch_verus[,2]
# sum(resultSV)/length(resultSV)
# 
# resultTV = sch_verus[,2] == trogl_verus[,2]
# sum(resultTV)/length(resultTV)



