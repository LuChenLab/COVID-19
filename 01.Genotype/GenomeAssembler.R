#!/usr/bin/env Rscript
if (!"optparse" %in% installed.packages()){
  stop('There is no package called "optparse"', call.=FALSE)
}
library("optparse")

option_list = list(
  make_option(c("-i", "--igv"), type = "character", default = NULL,
              help = "The path to igvtools executable[default = The path of which igvtools]", metavar = "character"), 
  
  make_option(c("-m", "--minMapQuality"), type = "integer", default = 0,
              help = "The minMapQuality for reads filtering[default = %default]", metavar = "integer"), 
  
  make_option(c("-b", "--BamFile"), type = "character", default = NULL,
              help = "The path of BamFile[default = %default]", metavar = "integer"), 
  
  make_option(c("-r", "--Reference"), type = "character", default = "/mnt/raid62/BetaCoV/Person/tangchao/analysis/Assembly/igvtools/genome/MN908947.fasta",
              help = "The path of reference genome[default = %default]", metavar = "character"),
  
  make_option(c("-T", "--TMPDIR"), type = "character", default = NULL,
              help = "The path of temporary files[default = tempdir()]", metavar = "character"), 
  
  make_option(c("-M", "--major"), type = "double", default = 0.7,
              help = "The threshold of major base sequence, if a base have more than one nucleotide sequence and the there is a major sequence (fraction more than this threshold) will only output the major one[default = %default]", metavar = "double"),
  
  make_option(c("-D", "--degeneracy"), type = "logical", default = FALSE,
              help = "If a base have more than one nucleotide sequence and there are no major sequence, this parameter indicating whether to output degenerate bases except ACTG[default = %default]", metavar = "logical"), 
  
  make_option(c("-d", "--depth"), type = "integer", default = 1000,
              help = "The minimum depth for degeneracy base calculating[default = %default]", metavar = "integer"), 
  
  make_option(c("-t", "--threshold"), type = "double", default = 0.8,
              help = "The threshold of minimum fraction reads indicating a base are insertion and deletion[default = %default]", metavar = "double"), 
  
  make_option(c("-Y", "--homogeneity"), type = "double", default = 0.8,
              help = "The minimum fraction of homogeneity of insertion (If a base have a true homogeneous insertion mean its fraction can't less than this homogeneity value)[default = %default]", metavar = "double"), 
  
  make_option(c("-c", "--MinimumCoverage"), type = "integer", default = 10,
              help = "The minimum coverage of a base, if the effective coverage less than MinimumCoverage, the reads of this base will be invalid reads and genome sequence will use refrence genome at the same coordinate[default = %default]", metavar = "integer"),
  
  make_option(c("-C", "--insertionCalling"), type = "integer", default = 1000,
              help = "The minimum depth for calling a insertion[default = %default]", metavar = "integer"), 

  make_option(c("-q", "--deletionCalling"), type = "integer", default = 100,
              help = "The minimum depth of deletion reads for calling a deletion[default = %default]", metavar = "integer"), 
  
  make_option(c("-R", "--ReplaceAllN"), type = "logical", default = FALSE,
              help = "This parameter indicating whether to replace the all N of the output genome fasta file use refrence file[default = %default]", metavar = "logical"),
  
  make_option(c("-s", "--Trim"), type = "logical", default = FALSE,
              help = "This parameter indicating whether to trim the N at head and tail of the output genome fasta file[default = %default]", metavar = "logical"),
  
  make_option(c("-I", "--ReplaceInternalN"), type = "logical", default = FALSE,
              help = "This parameter indicating whether to replace the N at internal of the output genome fasta file use refrence file[default = %default]", metavar = "logical"),
  
  make_option(c("-H", "--HeaderOfFASTA"), type = "character", default = NULL,
              help = "The title of output genome FASTA file[default = Name of BAM file]", metavar = "character"), 
  
  make_option(c("-O", "--Output"), type = "character", default = NULL,
              help = "The path and prefix of output file[default = %default]", metavar = "character"),
  
  make_option(c("-n", "--homopolymers"), type = "integer", default = 5,
              help = "The number of homopolymers to correct[default = %default]", metavar = "integer"), 
  
  make_option(c("-N", "--cores"), type = "integer", default = 10,
              help = "The number of cores in Parallel job[default = %default]", metavar = "integer")
);


# TMPDIR <- readLines(pipe("echo $TMPDIR"))

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(GenomicAlignments))

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

igv <- opt$igv
minMapQuality <- opt$minMapQuality
BamFile <- opt$BamFile
ref <- opt$Reference
tmpout <- opt$TMPDIR
threshold <- opt$major
degen <- opt$degeneracy
insdel <- opt$threshold
depth <- opt$depth
header <- opt$HeaderOfFASTA
output <- opt$Output
MinimumCoverage <- opt$MinimumCoverage
ReplaceAllN <- opt$ReplaceAllN
Trim <- opt$Trim
ReplaceInternalN <- opt$ReplaceInternalN
homogeneity <- opt$homogeneity
insertionCalling <- opt$insertionCalling
deletionCalling <- opt$deletionCalling
homopolymers <- opt$homopolymers
cores <- opt$cores

## Additional functions: ----
ParseBamRecord <- function(x) {
  chr <- runValue(seqnames(x))
  cgi <- cigar(x)
  seqi <- mcols(x)$seq
  posi <- mcols(x)$pos
  
  # process cigar
  parsed_cgi <- explodeCigarOpLengths(cgi)[[1]]
  names(parsed_cgi) <- explodeCigarOps(cgi)[[1]]
  parsed_cgi <- parsed_cgi[names(parsed_cgi) != "H"]
  
  if(sum(parsed_cgi[names(parsed_cgi) %in% c("S", "M", "I")]) != width(seqi))
    stop(simpleError("qwidth are not equal to S + M + I, maybe P in any CIGAR"))
  
  if(names(parsed_cgi)[1] == "S") {
    seqi <- subseq(seqi, start = as.numeric(parsed_cgi[1]) + 1, end = width(seqi))
    parsed_cgi <- parsed_cgi[-1]
  }
  
  if(is.element("S", names(parsed_cgi))) {
    seqi <- subseq(seqi, start = 1, end = width(seqi) - parsed_cgi["S"])
    parsed_cgi <- parsed_cgi[names(parsed_cgi) != "S"]
  }
  
  if(any(names(parsed_cgi) %in% c("=", "X"))) {
    names(parsed_cgi)[names(parsed_cgi) %in% c("=", "X")] <- "M"
  }
  
  if(!all(names(parsed_cgi) %in% c("M", "I", "D", "N"))) {
    stop(simpleError("Maybe P in your CIGAR"))
  }
  
  GRs <- list()
  for(j in seq_along(parsed_cgi)) {
    if(j == 1) {
      GRs[[j]] <- GRanges(seqnames = chr, ranges = IRanges(start = 1, end = parsed_cgi[j]), cigar = names(parsed_cgi)[j], seq = subseq(seqi, start = 1, end = parsed_cgi[j]))
    } else {
      ci <- names(parsed_cgi)[j]
      tmcg <- parsed_cgi[1:j]
      
      if(ci == "D") {
        start <- sum(tmcg[names(tmcg) != "I"]) - parsed_cgi[j] + 1
        end <- start + parsed_cgi[j] - 1
        GRs[[j]] <- GRanges(seqnames = chr, ranges = IRanges(start = start, end = end), cigar = "D", 
                            seq = seqi)
      }
      
      if(ci == "M") {
        start <- sum(tmcg[names(tmcg) != "I"]) - parsed_cgi[j] + 1
        end <- start + parsed_cgi[j] - 1
        
        start2 <- sum(tmcg[names(tmcg) %in% c("M", "I")]) - parsed_cgi[j] + 1
        end2 <- start2 + parsed_cgi[j] -1
        GRs[[j]] <- GRanges(seqnames = chr, ranges = IRanges(start = start, end = end), cigar = "M", 
                            seq = subseq(seqi, start = start2, end = end2))
      }
      
      if(ci == "I") {
        start <- sum(tmcg[names(tmcg) != "I"]) + 1
        end <- start
        
        start2 <- sum(tmcg[names(tmcg) %in% c("M", "I")]) - parsed_cgi[j] + 1
        end2 <- start2 + parsed_cgi[j] -1
        
        GRs[[j]] <- GRanges(seqnames = chr, ranges = IRanges(start = start, end = end), cigar = "I", 
                            seq = subseq(seqi, start = start2, end = end2))
      }
      
      if(ci == "N") {
        start <- sum(tmcg[names(tmcg) != "I"]) - parsed_cgi[j] + 1
        end <- start + parsed_cgi[j] - 1
        GRs[[j]] <- GRanges(seqnames = chr, ranges = IRanges(start = start, end = end), cigar = "N", 
                            seq = seqi)
      }
    }
  }
  
  GRs <- unlist(GRangesList(GRs))
  GRs <- shift(GRs, shift = posi -1)
  
  GRs$seq <- as.character(GRs$seq)
  if(any(GRs$cigar %in% c("D", "N"))) {
    GRs[GRs$cigar %in% c("D", "N"), ]$seq <- NA
  }
  return(GRs)
}

ParseBamInsertion <- function(bam, cores = 1) {
  chr <- runValue(seqnames(bam))
  ReferenceSpace <- cigarRangesAlongReferenceSpace(cigar(bam), ops = "I", with.ops = TRUE, pos = start(bam))
  seqs <- mcols(bam)$seq
  QuerySpace <- cigarRangesAlongQuerySpace(cigar(bam), ops = "I", with.ops = TRUE)
  parallel::mclapply(seq_along(seqs), FUN = function(x){
    if(length(QuerySpace[[x]]) > 0) {
      subseq(seqs[rep(x, length(QuerySpace[[x]]))], QuerySpace[[x]])
    } else {
      NULL
    }
  }, mc.cores = cores) -> seqi
  
  seqi <- do.call("c", parallel::mclapply(FUN = function(x) as.character(x), seqi, mc.cores = cores))
  return(GRanges(seqnames = chr, ranges = unlist(ReferenceSpace), cigar = "I", seq = seqi))
}

ParseBamDeletion <- function(bam, cores = 1) {
  chr <- runValue(seqnames(bam))
  ReferenceSpace <- cigarRangesAlongReferenceSpace(cigar(bam), ops = "D", with.ops = TRUE, pos = start(bam))
  return(GRanges(seqnames = chr, ranges = unlist(ReferenceSpace)))
}

Fun1 <- function(x){
  if(rowSums(x[, 1:7]) == 0) {
    "N"
  } else {
    if (max(x[, 6:7])/sum(x[, 1:6]) > insdel) {
      if(as.numeric(x[, 6])/sum(x[, 1:6]) > insdel) {
        names(x[, 6])
      } else {
        if(max(x[, 1:4])/sum(x[, 1:4]) > threshold) {
          paste0(names(which.max(x[, 1:4])), "+INS")
        } else {
          if(degen) {
            if (any(x[, 1:4] > depth)) {
              paste0(consensusString(DNAStringSet(rep(names(x[, 1:4])[as.numeric(x[, 1:4]) > depth], as.numeric(x[, 1:4])[as.numeric(x[, 1:4]) > depth]))), "+INS")
            } else {
              if(sum(x[, 1:4] == max(x[, 1:4])) == 1) {
                paste0(names(x[, 1:4])[which(x[, 1:4] == max(x[, 1:4]))], "+INS")
              } else {
                candi <- names(x[, 1:4])[which(x[, 1:4] == max(x[, 1:4]))]
                if(is.element(x[[8]], candi)) {
                  paste0(x[[8]], "+INS")
                } else {
                  paste0(sample(candi, 1), "+INS")
                }
              }
            }
          } else {
            if(sum(x[, 1:4] == max(x[, 1:4])) == 1) {
              paste0(names(x[, 1:4])[which(x[, 1:4] == max(x[, 1:4]))], "+INS")
            } else {
              candi <- names(x[, 1:4])[which(x[, 1:4] == max(x[, 1:4]))]
              if(is.element(x[[8]], candi)) {
                paste0(x[[8]], "+INS")
              } else {
                paste0(sample(candi, 1), "+INS")
              }
            }
          }
        }
      }
    } else {
      if(max(x[,1:4])/sum(x[, 1:4]) > threshold) {
        names(which.max(x[, 1:4]))
      } else {
        if(degen) {
          if (any(x[, 1:4] > depth)) {
            consensusString(DNAStringSet(rep(names(x[, 1:4])[as.numeric(x[, 1:4]) > depth], as.numeric(x[, 1:4])[as.numeric(x[, 1:4]) > depth])))
          } else {
            if(sum(x[, 1:4] == max(x[, 1:4])) == 1) {
              names(x[, 1:4])[which(x[, 1:4] == max(x[, 1:4]))]
            } else {
              candi <- names(x[, 1:4])[which(x[, 1:4] == max(x[, 1:4]))]
              if(is.element(x[[8]], candi)) {
                x[[8]]
              } else {
                sample(candi, 1)
              }
            }
          }
        } else {
          if(sum(x[, 1:4] == max(x[, 1:4])) == 1) {
            names(x[, 1:4])[which(x[, 1:4] == max(x[, 1:4]))]
          } else {
            candi <- names(x[, 1:4])[which(x[, 1:4] == max(x[, 1:4]))]
            if(is.element(x[[8]], candi)) {
              x[[8]]
            } else {
              sample(candi, 1)
            }
          }
        }
      }
    }
  }
}


## Check parameters ----
if(is.null(igv)) 
  igv <- "igvtools"

fil <- pipe(description = paste("which", igv))
if(length(readLines(fil)) == 0)
  stop(simpleError(paste(igv, "are not exist or unexecutable.\n")))
unlink(fil)

if(!file.exists(BamFile))
  stop(simpleError(paste(BamFile, "are not exist.\n")))

if(!file.exists(BamFile))
  stop(simpleError(paste(BamFile, "are not exist.\n")))

if(!file.exists(paste0(BamFile, ".bai")))
  stop(simpleError(paste(paste0(BamFile, ".bai"), "are not exist.\n")))

if(!file.exists(ref))
  stop(simpleError(paste(ref, "are not exist.\n")))

if(is.null(tmpout)) {
  tmpout <- tempdir()
} else {
  if(!dir.exists(tmpout))
    dir.create(tmpout, recursive = TRUE)
}

if(threshold < 0.5 | threshold > 1)
  stop(simpleError("The major base fraction must must be between 0.5 and 1.\n"))

if(insdel < 0.5 | insdel > 1)
  stop(simpleError("The parameter threshold must must be between 0.5 and 1.\n"))

if(insdel < 0.8)
  warning(simpleWarning("It is recommended that parameter threshold is greater than 0.8.\n"))

if(homogeneity < 0.5 | homogeneity > 1)
  stop(simpleError("The parameter homogeneity must must be between 0.5 and 1.\n"))

if(homogeneity < 0.8)
  warning(simpleWarning("It is recommended that parameter homogeneity is greater than 0.8 for a reliable insertion.\n"))

if(is.null(header))
  header <- basename(BamFile)

if(!dir.exists(dirname(output))) 
  dir.create(dirname(output), recursive = TRUE)

tmpf <- paste(tmpout, "igvtools.count.tsv", sep = "/")


logFile <- paste0(output, ".log.final.out")

zz <- file(logFile, open = "wt") # "w" or "wt": Open for writing in text mode.
sink(zz, type = c("message"), append = TRUE)
sink(zz, type = c("output"), append = TRUE)

## Running igvtools ----

sh <- paste(igv, "count --bases --windowSize 1 --minMapQuality", minMapQuality, BamFile, "stdout", ref, ">", tmpf)
system(sh)


## Consensus sequence ----

igvres <- suppressWarnings(data.table::fread(tmpf))
colnames(igvres) <- c("Pos", "A", "C", "G", "T", "N", "DEL", "INS")
igvres[, Pos:=as.character(Pos)]
setkey(igvres, Pos)
geno <- Biostrings::readBStringSet(ref)
leng <- max(c(as.numeric(tail(igvres$Pos, 1)), width(geno)))
igvres <- igvres[as.character(seq_len(leng)), ]
igvres[is.na(igvres)] <- 0

igvres[seq_len(width(geno)), Ref := strsplit(as.character(geno), "")]

igvres[INS < insertionCalling, INS := 0]

igvres[, Seq := Fun1(.SD), by = Pos]

DNA_TAB <- list(A = c("A"), 
                C = c("C"), 
                G = c("G"), 
                `T` = c("T"), 
                M = c("A", "C"), 
                R = c("A", "G"), 
                W = c("A", "T"), 
                S = c("G", "C"), 
                Y = c("C", "T"), 
                K = c("G", "T"), 
                V = c("G", "A", "C"), 
                H = c("A", "T", "C"), 
                D = c("G", "A", "T"), 
                B = c("G", "T", "C"))

igvres[Seq == "DEL" & DEL < deletionCalling, Seq := Ref]

## homopolymers correction ----

rle <- igvres[, S4Vectors::Rle(Ref)]

cg <- paste0(runLength(rle), as.character(factor(runValue(rle), labels = CIGAR_OPS[c(1, 4, 8, 9)])), collapse = "")
cg_gr <- GenomicAlignments::cigarRangesAlongReferenceSpace(cg, with.ops = T)[[1]]
cg_gr <- cg_gr[width(cg_gr) >= homopolymers, ]

for(i in seq_along(cg_gr)) {
  if(length(unique(igvres[start(cg_gr[i]):end(cg_gr[i]), Seq])) != 1) {
    igvres[start(cg_gr[i]):end(cg_gr[i]), Seq := Ref]
  }
}

# correct homopolymers flank base

for(i in seq_along(cg_gr)) {
  ps <- start(cg_gr[i])
  psf <- ifelse(ps == 1, ps, ps - 1)
  psfff <- ifelse(ps == 1, ps, ps - 3)
  pe <- end(cg_gr[i])
  pef <- ifelse(pe == nrow(igvres), pe, pe + 1)
  pefff <- ifelse(pe == nrow(igvres), pe, pe + 3)
  
  if(all(igvres[ps:pe, Seq == Ref])) {
    if(!all(igvres[psfff:pef, Seq == Ref])) {
      igvres[psf, Seq := Ref]
    }
    if(!all(igvres[pef:pefff, Seq == Ref])) {
      igvres[pef, Seq := Ref]
    }
  }
}

## polish of deltion ----

cg <- Rle(igvres[, ifelse(Seq == "DEL", "D", ifelse(nchar(Seq) > 1, "M", ifelse(Seq == Ref, "=", ifelse(Seq == "N", "N", "X"))))])
cg <- paste0(runLength(cg), runValue(cg), collapse = "")

deletion <- cigarRangesAlongReferenceSpace(cg, with.ops = TRUE, ops = c("D"))[[1]]

if(length(deletion) != 0) {
  bh <- scanBamHeader(BamFile, what = "targets")
  chr <- tryCatch(unique(names(bh[[1]]$targets)), error = function(e) "chromosome")
  deletion <- GRanges(seqnames = chr, ranges = deletion)
  
  for(i in seq_along(deletion)) {
    gr <- deletion[i]
    bam <- GenomicAlignments::readGAlignments(file = BamFile, param = ScanBamParam(flag = scanBamFlag(isNotPassingQualityControls = FALSE), 
                                                                                   what = c("qname", "pos", "qwidth", "cigar", "seq"), 
                                                                                   which = gr, 
                                                                                   mapqFilter = minMapQuality))
    bam <- bam[grep("D", cigar(bam)), ]
    
    y <- ParseBamDeletion(bam)
    y <- y[countOverlaps(y, gr, type = "any") > 0]
    y <- y[as.character(y) %in% names(which(prop.table(table(y)) > 0.05))]
    
    if(max(prop.table(table(y))) > 0.5) {
      y <- as(names(which.max(prop.table(table(y)))), "GRanges")
      igvres[start(y):end(y), Seq := "DEL"]
    } else {
      igvres[start(gr):end(gr), Seq := Ref]
    }
  }
}


## Coverage and Efficiency ----

# Output second base 
mclapply(seq_len(nrow(igvres)), FUN = function(i) {
  if(sum(igvres[i, 2:8]) == 0) {
    f = "N"
    s = "N"
  } else {
    f = names(which.max(igvres[i, 2:8]))
    if(sum(igvres[i, 2:8] > 0) == 1) {
      s = "N"
    } else {
      s = names(sort(igvres[i, 2:8], decreasing = T))[2]
    }
  }
  data.frame(First = f, Second = s)
}, mc.cores = cores) -> top_base
top_base <- do.call(rbind, top_base)

igvres <- cbind(igvres, top_base)

# Coverage

igvres[, Depth := rowSums(igvres[, 2:7])]

RightBase <- mcmapply(FUN = function(i){
  n <- gsub("\\+INS", "", igvres[i, Seq])
  if(is.element(n, colnames(igvres))) {
    as.numeric(igvres[i, ..n])
  } else {
    m <- eval(parse(text = paste(paste("DNA_TAB", n, sep = "$"))))
    as.numeric(sum(igvres[i, ..m]))
  }
}, seq_len(nrow(igvres)), mc.cores = cores)

igvres[, Evidence := RightBase]
igvres[, Efficiency := Evidence/Depth]

igvres[Evidence > 0 & Evidence < MinimumCoverage & Ref != Seq, Seq := Ref]
igvres[, Evidence := NULL]

## Fix INS sequence ----

if(any(with(igvres, grepl("INS", Seq)))) {
  pos <- as.numeric(igvres[grepl("INS", Seq), Pos])
  
  bh <- scanBamHeader(BamFile, what = "targets")
  chr <- tryCatch(unique(names(bh[[1]]$targets)), error = function(e) "NC_045512.2")
  
  for(i in seq_along(pos)) {
    Pos <- pos[i]
    gr <- GenomicRanges::GRanges(seqnames = chr, ranges = IRanges(start = Pos, end = Pos), strand = "*")
    bam <- GenomicAlignments::readGAlignments(file = BamFile, param = ScanBamParam(flag = scanBamFlag(isNotPassingQualityControls = FALSE), what = c("qname", "pos", "qwidth", "cigar", "seq"), which = gr, mapqFilter = minMapQuality))
    bam <- bam[grep("I", cigar(bam)), ]
    
    # ParsedReads <- parallel::mclapply(seq_along(bam), FUN = function(x) {
    #   y <- ParseBamRecord(bam[x])
    #   y[y$cigar == "I", ]
    # }, mc.cores = cores)
    # ReadsIns <- unlist(GRangesList(ParsedReads))
    ReadsIns <- ParseBamInsertion(bam = bam, cores = cores)
    PosINS <- ReadsIns[end(ReadsIns) == Pos, ]
    
    if(max(prop.table(table(PosINS$seq))) > homogeneity) {
      is <- names(sort(table(PosINS$seq), decreasing = T)[1])
      igvres[Pos == as.character(pos[i]), Seq := gsub("\\+INS", is, Seq)]
    } else {
      igvres[Pos == as.character(pos[i]), Seq := gsub("\\+INS", "", Seq)]
    }
  }
}

igvres[nchar(Seq) > 1 & !grepl("DEL", Seq), Efficiency := INS/Depth]

cg <- Rle(igvres[, ifelse(Seq == "DEL", "D", ifelse(nchar(Seq) > 1, "M", ifelse(Seq == Ref, "=", ifelse(Seq == "N", "N", "X"))))])
cg <- paste0(runLength(cg), runValue(cg), collapse = "")
deletion <- cigarRangesAlongReferenceSpace(cg, with.ops = TRUE, ops = c("D"))[[1]]

message(paste("Genome coverage (%):", mean(igvres[, Depth] > 0) * 100))
message(paste("deletion:", length(deletion)))
message(paste("insertion:", sum(nchar(igvres$Seq) > 1 & igvres$Seq != "DEL")))
message(paste("SNP:", nrow(igvres[!Seq %in% c("N", "DEL") & stringr::str_extract(Seq, "^[[:alpha:]]") != Ref, ])))

data.frame(igvres[Seq != "N" & Seq != Ref, ])

fwrite(igvres, paste0(output, ".SummaryStatistic.txt"), row.names = F, col.names = T, na = "NA", sep = "\t", quote = FALSE)

if(ReplaceAllN) {
  igvres[Seq == "N", Seq := Ref]
}

if(ReplaceInternalN) {
  igvres[, Pos := as.integer(Pos)]
  st <- igvres[min(which(Seq != "N")), Pos]
  ed <- igvres[max(which(Seq != "N")), Pos]
  igvres[Pos > st & Pos < ed & Seq == "N", Seq := Ref]
}

sequence <- paste0(igvres[Seq != "DEL", Seq], collapse = "")
message(paste("Genome size:", nchar(sequence)))

if(Trim) {
  sequence <- gsub("^N+", "", sequence)
  sequence <- gsub("N+$", "", sequence)
  message(paste("Trimed genome size:", nchar(sequence)))
}

fa <- DNAStringSet(sequence)
names(fa) <- header

writeXStringSet(fa, paste0(output, ".genome.fasta"), format = "fasta")

if(file.remove(tmpf) & file.remove("./igv.log"))
  message("All finished")

sink(type = "message")
close(zz)
## end