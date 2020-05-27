
## High frequency mutation ----

files <- list.files(path = "/mnt/raid62/BetaCoV/Person/tangchao/analysis/Mutation/VCFs/20200511/", pattern = "vcf.gz$", full.names = T)
library(vcfR)
library(data.table)
vcfs <- lapply(files, FUN = read.vcfR, verbose = FALSE)
names(vcfs) <- sub(".pass", "", sub(".vcf.gz", "", basename(files)))

vcfs <- vcfs[mapply(vcfR::nrow, vcfs, USE.NAMES = FALSE) > 0]

vcfs_Fix <- lapply(vcfs, function(x) {
  if(nrow(x) > 1) {
    y <- as.data.table(getFIX(x))
  } else {
    y <- as.data.table(t(as.data.frame(getFIX(x))))
  }
  # y[, POS := as.numeric(POS)]
  # y[, SupportFraction := extract.info(x, element = "SupportFraction", as.numeric = T)]
  # y[, BaseCalledFraction := extract.info(x, element = "BaseCalledFraction", as.numeric = T)]
})

for(i in seq_along(vcfs_Fix)) {
  vcfs_Fix[[i]]$Sample <- names(vcfs_Fix)[i]
}

vcfs_Tab <- do.call(rbind, vcfs_Fix)
vcfs_Tab[, ID := paste0(REF, POS, ALT)]
vcfs_Tab[, QUAL := NULL]
vcfs_Tab <- unique(vcfs_Tab)

length(unique(vcfs_Tab$ID))
# [1] 122

table(vcfs_Tab[, .N, by = ID][, N])
#  1  2  3  4  6  7  9 10 12 23 29 
# 87 15  8  5  1  1  1  1  1  1  1 

length(unique(vcfs_Tab[nchar(REF) == 1 & nchar(ALT) == 1, ID]))
# [1] 104

table(vcfs_Tab[nchar(REF) == 1 & nchar(ALT) == 1, .N, by = ID][, N])
#  1  2  3  4  6  7  9 10 29 
# 73 14  8  4  1  1  1  1  1 

length(unique(vcfs_Tab[!(nchar(REF) == 1 & nchar(ALT) == 1), ID]))
# [1] 18
table(vcfs_Tab[!(nchar(REF) == 1 & nchar(ALT) == 1), .N, by = ID][, N])
#  1  2  4 12 23 
# 14  1  1  1  1 

# 1: TGCACCTCATGGTCATGTTATGGTTGAGCTGGTA499T 23
# 2:                             ACC18108AT 12
# 3:                         GAAGATTTTC728G  2
# 4:                              TTG18986T  4


length(unique(vcfs_Tab[ID %in% c("TGCACCTCATGGTCATGTTATGGTTGAGCTGGTA499T", "GAAGATTTTC728G"), Sample]))

vcfs_Tab[ID %in% c('C3717T', 'C8782T', 'A13175G', 'T13243C', 'A16824G', 'C17373T', 'C18060T', 'C18110T', 'C21711T', 'G26144T', 'T28144C', 'G28878A', 'C29095T', 'C29640T', 'G29742A'), .N, by = "ID"]




High_Freq_SNP <- vcfs_Tab[nchar(REF) == 1 & nchar(ALT) == 1, .N, by = "ID"][N >= 1, ID]
High_Freq_Indel <- vcfs_Tab[nchar(REF) > 1 | nchar(ALT) > 1, .N, by = "ID"][N >= 1, ID]

High_Freq <- c(High_Freq_SNP, High_Freq_Indel)

High_Freq <- unique(vcfs_Tab[ID %in% High_Freq, .(POS, ID, REF, ALT)])

High_Freq[, POS := as.numeric(POS)]

High_Freq[ID %in% High_Freq_SNP, start := POS]
High_Freq[ID %in% High_Freq_SNP, end := POS]

High_Freq[ID %in% High_Freq_Indel, start := POS + 1]
High_Freq[ID %in% High_Freq_Indel, end := POS + nchar(REF) - 1]

High_Freq[ID %in% High_Freq_SNP, FinalID := ID]
High_Freq[ID %in% High_Freq_Indel, FinalID := paste0(start, "-", end, ":D")]
High_Freq[nchar(REF) > 1 & nchar(ALT) > 1, FinalID := ID]

High_Freq[ID %in% High_Freq_SNP, Ref := REF]
High_Freq[ID %in% High_Freq_SNP, Alt := ALT]

High_Freq[nchar(ALT) > 1, Alt := substring(ALT, 2)]
High_Freq[nchar(REF) > 1, Ref := substring(REF, 2)]
High_Freq[nchar(REF) > 1 & nchar(ALT) == 1, Alt := "DEL"]

High_Freq <- merge(High_Freq, vcfs_Tab[ID %in% High_Freq$ID, .N, by = "ID"])

write.table(High_Freq, "/mnt/raid62/BetaCoV/Person/tangchao/analysis/Mutation/High_Freq_Mutation/20200511/High_Freq_SNP_and_INDEL.txt", row.names = F, quote = F, sep = "\t")


## Ratio of high frequency mutation ----

### ARTIC ----

vcfs_Fix2 <- lapply(vcfs, function(x) {
  if(nrow(x) > 1) {
    y <- as.data.table(getFIX(x))
  } else {
    y <- as.data.table(t(as.data.frame(getFIX(x))))
  }
  y[, POS := as.numeric(POS)]
  y[, SupportFraction := extract.info(x, element = "SupportFraction", as.numeric = T)]
  y[, BaseCalledFraction := extract.info(x, element = "BaseCalledFraction", as.numeric = T)]
})

for(i in seq_along(vcfs_Fix2)) {
  vcfs_Fix2[[i]]$Sample <- names(vcfs_Fix2)[i]
}

vcfs_Tab2 <- do.call(rbind, vcfs_Fix2)
vcfs_Tab2[, ID := paste0(REF, POS, ALT)]

vcfs_Tab2 <- vcfs_Tab2[ID %in% High_Freq$ID, ]
vcfs_Tab2[, QUAL := NULL]
vcfs_Tab2[, FILTER := NULL]

vcfs_Tab2 <- vcfs_Tab2[, .(SupportFraction = mean(SupportFraction), BaseCalledFraction = mean(BaseCalledFraction)), by = c("Sample", "ID")]

vcfs_Tab2 <- merge(High_Freq, vcfs_Tab2)

write.table(vcfs_Tab2, "/mnt/raid62/BetaCoV/Person/tangchao/analysis/Mutation/High_Freq_Mutation/20200511/High_Freq_SNP_and_INDEL_Ratio_of_ARTIC.txt", row.names = F, quote = F, sep = "\t")
openxlsx::write.xlsx(vcfs_Tab2, "/mnt/raid62/BetaCoV/Person/tangchao/analysis/Mutation/High_Freq_Mutation/20200511/High_Freq_SNP_and_INDEL_Ratio_of_ARTIC.xlsx")


### StrainGV ----

High_Freq[, start := as.numeric(start)]
High_Freq[, end := as.numeric(end)]


setwd("/mnt/raid62/BetaCoV/Person/tangchao/analysis/Assembly/igvtools/fast5_seqtk_v4")
library(data.table)
files <- list.files(".", pattern = ".SummaryStatistic.txt$", recursive = TRUE)

Summa <- lapply(files, function(x) fread(x))
names(Summa) <- sub(".SummaryStatistic.txt", "", basename(files))
Summa <- Summa[sub(".vcf.gz", "", sub(".pass.vcf.gz", "", list.files(path = "/mnt/raid62/BetaCoV/Person/tangchao/analysis/Mutation/VCFs/20200511/", pattern = "vcf.gz$", full.names = F)))]


lapply(seq_len(nrow(High_Freq)), function(i){
  s <- High_Freq[i, start]
  e <- High_Freq[i, end]
  a <- High_Freq[i, Alt]
  
  if(s == e) {
    mapply(Summa, FUN = function(x) {
      as.numeric(x[s, a, with = FALSE])/as.numeric(x[s, Depth])
    }) -> Fraction
    Fraction <- setDT(as.data.frame(Fraction), keep.rownames = "Sample")
    Fraction$ID <- High_Freq[i, ID]
  } else {
    if(a == "DEL") {
      mapply(Summa, FUN = function(x) {
        as.numeric(x[s:e, mean(DEL/Depth)])
      }) -> Fraction
      Fraction <- setDT(as.data.frame(Fraction), keep.rownames = "Sample")
      Fraction$ID <- High_Freq[i, ID]
    } else {
      mapply(Summa, FUN = function(x) {
        (mean(x[s:e, a, with = F][[1]]) + x[s:e, mean(DEL)])/x[s:e, mean(Depth)]
      }) -> Fraction
      Fraction <- setDT(as.data.frame(Fraction), keep.rownames = "Sample")
      Fraction$ID <- High_Freq[i, ID]
    }
  }
  return(Fraction)
}) -> StrainGV_Ratio

StrainGV_Ratio <- do.call(rbind, StrainGV_Ratio)

StrainGV_Ratio <- merge(High_Freq, StrainGV_Ratio)

write.table(StrainGV_Ratio, "/mnt/raid62/BetaCoV/Person/tangchao/analysis/Mutation/High_Freq_Mutation/20200511/High_Freq_SNP_and_INDEL_Ratio_of_StrainGV.txt", row.names = F, quote = F, sep = "\t")
openxlsx::write.xlsx(StrainGV_Ratio, "/mnt/raid62/BetaCoV/Person/tangchao/analysis/Mutation/High_Freq_Mutation/20200511/High_Freq_SNP_and_INDEL_Ratio_of_StrainGV.xlsx")


