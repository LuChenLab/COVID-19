library(data.table)
meta <- as.data.table(readxl::read_xlsx("/mnt/raid62/BetaCoV/Person/tangchao/data/metainfor/final sample information .xlsx"))
meta <- meta[, .(`样本编号`, BamPath, AssemblyPath, Experiment)]
meta[, bamfile := paste0(BamPath, "/", Experiment, ".bam")]

stopifnot(all(file.exists(meta$bamfile)))
mapply(dir.create, unique(meta$AssemblyPath), recursive = TRUE)

meta[, prefix := paste0(AssemblyPath, "/", `样本编号`)]


sh <- meta[, paste("Rscript /mnt/raid62/BetaCoV/Person/tangchao/script/GenomeAssembly/GenomeAssembler.R -b", bamfile, "-O", prefix, "-s TRUE -m 60 -t 0.6 -Y 0.5 -C 100 -H", `样本编号`, 
                   "&& Rscript /mnt/raid62/BetaCoV/Person/tangchao/script/GenomeAssembly/MutationAnnotater.R -S", paste0(prefix, ".SummaryStatistic.txt"))]

p <- unique(meta$AssemblyPath)
for(i in p) {
  write.table(sh[meta[, AssemblyPath] == i], file = paste0(i, "/GenomeAssembler.sh"), row.names = FALSE, col.names = FALSE, quote = F)
}


list.files("/mnt/raid62/BetaCoV/Person/tangchao/analysis/Assembly/igvtools/fast5_seqtk_v4", pattern = "gz$", recursive = T)


cd /mnt/raid62/BetaCoV/Person/tangchao/analysis/Assembly/igvtools/fast5_seqtk_v4
cd 20200223
bash GenomeAssembler.sh
cd ../20200305
bash GenomeAssembler.sh
cd ../20200314
bash GenomeAssembler.sh



cd /mnt/raid62/BetaCoV/Person/tangchao/analysis/Assembly/igvtools/fast5_seqtk_v4
cd 20200320
bash GenomeAssembler.sh
cd ../20200325_0af3696c
bash GenomeAssembler.sh
cd ../20200325_70b09167
bash GenomeAssembler.sh



cd /mnt/raid62/BetaCoV/Person/tangchao/analysis/Assembly/igvtools/fast5_seqtk_v4
cd 20200327_0a3f2e8a
bash GenomeAssembler.sh
cd ../20200327_5c4d299d
bash GenomeAssembler.sh
cd ../20200327_67070962
bash GenomeAssembler.sh



cd /mnt/raid62/BetaCoV/Person/tangchao/analysis/Assembly/igvtools/fast5_seqtk_v4
cd 20200402
bash GenomeAssembler.sh
cd ../20200403_25141bd2
bash GenomeAssembler.sh
cd ../20200403_67512ee6
bash GenomeAssembler.sh



cd /mnt/raid62/BetaCoV/Person/tangchao/analysis/Assembly/igvtools/fast5_seqtk_v4
cd 20200403_ba97463e
bash GenomeAssembler.sh
cd ../20200420
bash GenomeAssembler.sh
cd ../qianmai
bash GenomeAssembler.sh
