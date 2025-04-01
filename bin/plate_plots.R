#!/usr/bin/env Rscript

library("tidyverse")

library("fastqcr")
library("ggplate")
library("RColorBrewer")
library("cowplot")

#Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
id <- args[1]
qc.dir <-"/scratch-cbe/users/elin.axelsson/TagSeq/work/1b/4cf6cc3e3b69ce5d955f374182de4a/"
qc.dir <- "./"
qc <- qc_aggregate(qc.dir)
qc <- qc %>% filter(str_detect(sample, "[A-H][0-9]+")) %>% select(sample, tot.seq, pct.gc, pct.dup) %>% distinct()

# extract part of sample that looks like this "A1", "B2", etc
data_nr_plate = qc %>% mutate(well=str_extract(sample, "[A-H][0-9]+"), `seq.(million)` = as.numeric(tot.seq)/1000000, `gc(%)` = as.numeric(pct.gc), `dup(%)` = as.numeric(pct.dup))


p1 <- plate_plot(
  data = data_nr_plate,
  position = well,
  value = `seq.(million)`,
  plate_size = 96,
  plate_type = "round",
  colour =  c("white",brewer.pal(9,name="Blues")),
  title = "Number of sequences",
  scale = 1.5
)
#dev.off()

#png("plate_dup_mqc.png")
p2 <- plate_plot(
  data = data_nr_plate,
  position = well,
  value = `dup(%)`,
  plate_size = 96,
  plate_type = "round",
  colour =  brewer.pal(9,name="Reds"),
  title = "Percentage of duplicates",
  scale = 1.5)
#dev.off()

#png("plate_gc.png_mqc")
p3 <- plate_plot(
  data = data_nr_plate,
  position = well,
  value = `gc(%)`,
  plate_size = 96,
  plate_type = "round",
  title = "GC content",
  colour = brewer.pal(9,name="YlGnBu"),
  scale = 1.5
)
#dev.off()

sink("sessionInfo.txt")
sessionInfo()
sink()


over_rep_tbl <- tibble()
for (i in dir(qc.dir)) {
  if (str_detect(i, "[A-H][0-9]+")) {
    temp <- qc_read(paste(qc.dir,i,sep=""), modules = "Overrepresented sequences")[[1]]
    if(nrow(temp)>0){
      temp <- temp %>% summarize(tot = sum(Percentage)) %>% pull()
     } else {
      temp <- 0
     }
    well <- str_extract(i, "[A-H][0-9]+")
    over_rep_tbl <- bind_rows(over_rep_tbl, tibble(well = well, proc = temp))
}
}
p4 <- plate_plot(
  data = over_rep_tbl,
  position = well,
  value = proc,
  plate_size = 96,
  plate_type = "round",
  title = "% overrep",
  colour = brewer.pal(9,name="YlOrRd"),
  scale = 1.5
  
)
options(bitmapType='cairo')
dir.create(file.path(id))
png(paste(id,"Plate_QC_mqc.png",sep="/"), width = 1200, height = 1000)
plot_grid(p1, p2, p3, p4, ncol = 2, nrow = 2,greedy = FALSE)
dev.off()
                   