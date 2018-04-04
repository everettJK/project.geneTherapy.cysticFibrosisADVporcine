# Retrieve patient and trial information.
options(stringsAsFactors = FALSE, useFancyQuotes = FALSE)
library(RMySQL)
library(GenomicRanges)
library(BSgenome.Sscrofa.UCSC.susScr3)
library(gtools)
library(reshape2)
library(gt23)
library(ShortRead)
library(parallel)
library(ggseqlogo)
library(tidyverse)

r <- list()


# Characterize the most frequent R2 reads which are desinged to read  
# out of the transposon and into the adjacent genomic DNA.
# .~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~

R2 <- sread(readFastq('/data/sequencing/Illumina-archive/180220_M00281_0323_000000000-BMWYV/Data/Intensities/BaseCalls/Undetermined_S0_L001_R2_001.fastq.gz'))

r$nReads <- length(R2)

o <- vcountPattern('ATTATCTTTCTAGGG', R2, max.mismatch=1)
r$percentReadsWithITR <- round((sum(o) / length(R2))*100, 2)

R2 <- R2[which(o > 0)]

h <- data.frame(head(sort(table(as.character(subseq(R2, 16, 65))), decreasing = TRUE), n = 500))
h$percentReads <- round((h$Freq / length(R2)*100), 2)
h$cumlativePercentage <- cumsum(h$percentReads)

transposonVector <- sread(readFasta('data/transposonVector.ff'))
transposaseVector <- sread(readFasta('data/transposaseVector.ff'))


#   Most abundant read      --------TTAAAAGATCTGGAAGGTGCTGAGGTACGATGAGACCCGCACCAGGTGCA--
#   Transposon vector       --CTAGGGTTAAAAGATCTGGAAGGTGCTGAGGTACGATGAGACCCGCACCAGGTGCAGA
#   Transposase vector      TCCCCGCGGTGGCAGATCTGGAAGGTGCTGAGGTACGATGAGACCCGCACCAGGTGCAGA
#                                    *   *********************************************  


m <- bind_rows(lapply(1:nrow(h), function(x){
  x  <- h[x,]
  n1 <- paste(unlist(start(vmatchPattern(as.character(x$Var1), transposonVector, max.mismatch=2))), collapse = ';')
  n2 <- paste(unlist(start(vmatchPattern(as.character(x$Var1), transposaseVector, max.mismatch=2))), collapse = ';')
  x$n1 <- ifelse(nchar(n1) == 0, NA, n1)
  x$n2 <- ifelse(nchar(n2) == 0, NA, n2)
  x
}))

write.csv(m, file='data/genomicSeqFreqTable.csv', row.names = FALSE)

# A manually update vection of genomicSeqFreqTable.csv (genomicSeqFreqTable_aligned.csv) was created
# in order to highlight the different adenovirus species that were identified.




# Cluster genomic sequences and indentify representative reads.
# .~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~

s <- subseq(R2, 16, 65)
names(s) <- paste0('S', 1:length(s))
writeFasta(s, file='data/genomicSequences_50NT.fasta')
#./software/cdhit/cd-hit-est -i data/genomicSequences_50NT.fasta -o data/genomicSequences_50NT_90 -c 0.90 -d 0 -M 20000 -T 20

parse_cdhitest_output <- function(file){
  clusters <- readChar(file, file.info(file)$size)
  clusters <- unlist(base::strsplit(clusters, '>Cluster'))
  clusters <- clusters[2:length(clusters)]
  lapply(clusters, function(x){
    gsub('\\.\\.\\.', '', unlist(lapply(str_match_all(x, '>([^\\s]+)'), function(y){ y[,2] })))
  })
}

# Identify cluster representatives in the CD-HIT output.
clusterData <- scan('data/genomicSequences_50NT_90.clstr', what = 'character', sep = '\n')
clusterReps <- clusterData[grep('\\*', clusterData)]
clusterReps <- gsub('\\.\\.\\.', '', unlist(lapply(str_match_all(clusterReps, '>([^\\s]+)'), function(y){ y[,2] })))

# Parse the CD-HIT output and use the cluster representatives sequence names for the resulting list names.
clusters <- parse_cdhitest_output('data/genomicSequences_50NT_90.clstr')
names(clusters) <- clusterReps

# Order the sequence read clusters by the number of reads in each from most to lease.
clusters <- clusters[order(unlist(lapply(clusters, length)), decreasing = TRUE)]

# Create a plot showing the number of reads in each cluster.
o <- bind_rows(mapply(function(n, x){ data.frame(Cluster=n, reads=length(x))}, 1:20, clusters[1:20], SIMPLIFY = FALSE))
o$label <- sprintf("%.1f%%", cumsum(unlist(lapply(clusters[1:20], function(x) (length(x) / length(R2))*100))))
readClusters90 <- ggplot(o, aes(Cluster, reads)) +
  theme_bw() +
  geom_bar(stat='identity') +
  labs(x='Sequence cluster', y='Reads in cluster') +
  geom_text(size=3, angle=45, aes(label=label, y=reads+500000))



# Write out the cluster reps.
clusterRepsSeqs <- as.character(s[names(clusters)[1:20]])
write(clusterRepsSeqs, file='data/genomicSequences_50NT_90_clusterReps.txt')



# Retrieve intSites from database.
# .~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~

invisible(sapply(DBI::dbListConnections(MySQL()), dbDisconnect))

# Retrieve subject data from specimen management database.
subjects <- c("p1166", "p1168", "p1169", "p1171", "p9846", "p9850", "p9851")

# Retrieve intSites.
intSites <- getIntSiteData('specimen_management', 'intSites_miseq', patients = subjects)


# Setup parallelization.
CPUs <- 10
cluster <- makeCluster(CPUs)


# Create a splitting vector for parallelization.
intSites$s <- ceiling(seq_along(intSites)/(length(intSites)/CPUs))


# Standardized start positions by patient.
intSites <- unlist(GRangesList(parLapply(cluster, split(intSites, intSites$patient), function(x){
  gt23::standardizeSeqRanges.malani(x, standardizeStart=TRUE, standardizeEnd=FALSE)
})))


# Standardize break points by sample.
intSites <- unlist(GRangesList(parLapply(cluster, split(intSites, intSites$sampleName), function(x){
  gt23::standardizeSeqRanges.malani(x, standardizeStart=FALSE, standardizeEnd=TRUE)
})))


# Merge replicate samples and estimate clonal abundances.
intSites <- unlist(GRangesList(parLapply(cluster, split(intSites, intSites$patient), function(x){
  gt23::calcReplicateAbundances(x)
})))


# Add nearest gene and nearest oncogene annotations.
names(intSites) <- NULL
intSites <- unlist(GRangesList(parLapply(cluster, split(intSites, intSites$s), function(x){
  library(gtools)
  x <- gt23::nearestGenomicFeature(x, genome = 'susScr3')
  x[mixedorder(x$posid),]
})))


# Cleanup
intSites$s <- NULL
names(intSites) <- NULL
stopCluster(cluster)



# Plots and tables.
# .~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~

# Create a list of chromosome lengths for select chromosomes, ie. [['chr1']] <- 248956422
chromosomeLengths <- sapply(rev(paste0("chr", c(seq(1:18), "X", "Y"))),
                            function(x) length(BSgenome.Sscrofa.UCSC.susScr3[[x]]),
                            simplify = FALSE, USE.NAMES = TRUE)

intSiteDistributionPlot(intSites, chromosomeLengths)
intSiteLocPlot <- intSiteDistributionPlot(intSites, chromosomeLengths, alpha = 0.8)


# Create a table of intSites rather than a relative abundance plots because of the low number of intSites.
intSiteTable <- data.frame(intSites) %>% 
                dplyr::group_by(patient, cellType) %>%
                dplyr::summarise(sites = length(unique(posid))) %>% 
                reshape2::dcast(patient~cellType, value.var='sites')

# This table sums unique sites from each replicate. Small differencs from {intSiteTable} arise from the same site being found in multiple replicates.

dbConn  <- dbConnect(MySQL(), group='specimen_management')
sampleData <-  dbGetQuery(dbConn, paste0('select * from gtsp'))
dbDisconnect(dbConn)


runStats <- read.table('data/180220_M00281_0323_000000000-BMWYV.stats', sep = '\t', header = TRUE) %>% 
            dplyr::select(sample, ltredlinkered, vTrimed, readsWithGoodAlgnmts, numUniqueSites) %>%  
            dplyr::mutate(GTSP = gsub('\\-\\d+$', '', sample)) %>%
            dplyr::group_by(GTSP) %>%
            dplyr::select(-sample) %>%
            dplyr::summarise_all(funs(sum)) %>%
            dplyr::ungroup() %>%
            tidyr::drop_na() %>%
            dplyr::mutate(subject  = sampleData[match(GTSP, sampleData$SpecimenAccNum),]$Patient) %>%
            dplyr::mutate(cellType = sampleData[match(GTSP, sampleData$SpecimenAccNum),]$CellType) %>%
            dplyr::select(subject, cellType, ltredlinkered, vTrimed, readsWithGoodAlgnmts, numUniqueSites)

names(runStats) <- c('Subject', 'Sample', 'Filtered reads', 'Non-vector like reads', 'Reads aligning to genome', 'intSites (reps. combined)')


# Heatmap creation.
# .~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~

set.seed(42)
WAS_d0 <- readRDS('data/WAS_d0_intSites.rds')
WAS_d0 <- sample(WAS_d0, length(intSites))

# Write out a sample file for the heatmap maker.
write('sampleName,GTSP,patient',     file = 'output/samples', append = FALSE)
write('piggyBac,piggyBac,xxx', file = 'output/samples', append = TRUE)
write('lentiVirus,lentiVirus,xxx',   file = 'output/samples', append = TRUE)

# Write out site data.
write('seqnames,strand,position,sampleName,refGenome', file = 'output/sites',  append = FALSE)

# Write out piggBac intSites.
df <- data.frame(seqnames(intSites), strand(intSites), start(intSites), 'piggyBac', 'hg38')
names(df) <- c('seqnames', 'strand', 'position', 'sampleName', 'refGenome')
write.table(df, quote = FALSE, row.names = FALSE, col.names = TRUE, file = 'output/sites',  sep=',', append = FALSE)

# Write WAS_d0 intSites.
df <- data.frame(seqnames(WAS_d0), strand(WAS_d0), start(WAS_d0), 'lentiVirus', 'hg38')
write.table(df, quote = FALSE, row.names = FALSE, col.names = FALSE, file = 'output/sites', sep=',', append = TRUE)

# Create output directories.
if(! dir.exists('output/genomicHeatmapMakerOutput'))  dir.create('output/genomicHeatmapMakerOutput')

# Create command line calls to create the heatmaps.
Rscript_path <- '/home/opt/R-3.4.0/bin/Rscript'

comm <- paste(Rscript_path, 'software/genomicHeatmapMaker-from_input/genomic_heatmap_from_file.R', 
              'output/samples',
              '-c software/genomicHeatmapMaker-from_input/INSPIIRED.yml',
              '-o output/genomicHeatmapMakerOutput',
              '-f output/sites',
              '-r hg38')

# Run heat map creation script.
# system(comm)


m <- genomicHeatmap2dataframe('output/genomicHeatmapMakerOutput/main.svg')
m$df$sample <- factor(as.character(m$df$sample), levels=c('lentiVirus', 'piggyBac'))
heatMap1 <- genericHeatmap(m$df, tileColors=m$colorScale$color)




# Characterize the up and downstream intSite sequences.

getSiteFlankingSequences <- function (g, genome, n = 12) 
{
  g.upstream <- flank(g, n, start = TRUE)
  g.downstream <- flank(g, n, start = FALSE)
  i <- which(strand(g.downstream) == "-")
  g.downstream[i] <- shift(g.downstream[i], 1)
  g.downstream[-i] <- shift(g.downstream[-i], -1)
  g.upstream.seq <- BSgenome::getSeq(genome, g.upstream)
  g.downstream.seq <- BSgenome::getSeq(genome, g.downstream)
  list(upstream = g.upstream.seq, downstream = g.downstream.seq)
}

f <- getSiteFlankingSequences(intSites, BSgenome.Sscrofa.UCSC.susScr3)

upStreamDownStreamLogo <- ggseqlogo(paste0(as.character(f$upstream), as.character(f$downstream)),
          seq_type = 'dna', 
          method = 'bit') + 
          theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

f <- getSiteFlankingSequences(intSites, BSgenome.Sscrofa.UCSC.susScr3, n = 4)
r$percentSitesAGGGupstream <- sprintf("%.1f%%", (sum(sapply(f$upstream, function(x){ as.character(x) == 'AGGG' })) / length(intSites))*100)

f <- getSiteFlankingSequences(intSites, BSgenome.Sscrofa.UCSC.susScr3, n = 4)
r$percentSitesTTAAdownstream <- sprintf("%.1f%%", (sum(sapply(f$downstream, function(x){ as.character(x) == 'TTAA' })) / length(intSites))*100)

save.image(file='data/project.RData')
save(r, heatMap1,  runStats, clusterRepsSeqs, intSiteTable, intSiteLocPlot, upStreamDownStreamLogo, readClusters90, file='data/report.RData')


rmarkdown::render('project.Rmd', 
                  params = list('author' = 'John K. Everett, Ph.D. and Frederic Bushman, Ph.D.',
                                'date'  = format(Sys.Date(), format="%B %d, %Y"),
                                'title' = 'Analysis of piggyBac mediated integration in the porcine genome'))



