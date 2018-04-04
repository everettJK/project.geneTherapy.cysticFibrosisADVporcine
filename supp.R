library(ShortRead)
library(stringr)
library(dplyr)
library(parallel)
library(rtracklayer)
library(BSgenome.Sscrofa.UCSC.susScr3)
library(ggplot2)
options(stringsAsFactors = FALSE)
CPU_batches <- 5
CPUs_per_batch <- 20

# Read in the R2 reads, trim them with a sliding window approach and then name them with arbitray ids.
# Focus on only reads that start with the sequencing primer.
# .~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~

R2 <- readFastq('/data/sequencing/Illumina-archive/180220_M00281_0323_000000000-BMWYV/Data/Intensities/BaseCalls/Undetermined_S0_L001_R2_001.fastq.gz')
stat1 <- length(R2)

R2 <- trimTailw(R2, 5, "?", 5)
stat2 <- length(R2)

R2 <- unique(sread(R2))
stat3 <- length(R2)

names(R2) <- paste0('S', 1:length(R2))

R2 <- R2[grep('^ATTAT', R2)]       # 78% of R2 reads start with ATTAT
stat4 <- length(R2)

R2 <- R2[width(R2) >= 50]          # Only consider reads that are at least 70 NT
stat5 <- length(R2)

# extract sequences from NT 25 -> end for clustering since the first ~ 25 NT
# will contain adapter sequences that will not align to genomes or vector experimental vectors.
R2.tails <- subseq(R2, start = 25) 



# Cluster read tails and find representative reads.
# .~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~

writeFasta(R2.tails , file='R2.tails.fasta')

#
# software/cdhit/cd-hit-est -i R2.tails.fasta -o R2.tails -c 0.95 -d 0 -M 20000 -T 30
#

parse_cdhitest_output <- function(file){
  clusters <- readChar(file, file.info(file)$size)
  clusters <- unlist(base::strsplit(clusters, '>Cluster'))
  clusters <- clusters[2:length(clusters)]
  lapply(clusters, function(x){
    gsub('\\.\\.\\.', '', unlist(lapply(str_match_all(x, '>([^\\s]+)'), function(y){ y[,2] })))
  })
}


clusterData <- scan('R2.tails.clstr', what = 'character', sep = '\n')
clusterReps <- clusterData[grep('\\*', clusterData)]
clusterReps <- gsub('\\.\\.\\.', '', unlist(lapply(str_match_all(clusterReps, '>([^\\s]+)'), function(y){ y[,2] })))

clusters    <- parse_cdhitest_output('R2.tails.clstr')
names(clusters) <- clusterReps

R2.tails.clusterReps <- R2.tails[clusterReps]
stat6 <- length(clusterReps)




# Blat cluster representatives against genome.
# .~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~

# Create a batch calculation chunking factor.
elementMetadata(R2.tails.clusterReps)$s1 <- dplyr::ntile(1:length(R2.tails.clusterReps), CPU_batches)
cluster <- makeCluster(CPUs_per_batch)
invisible(lapply(split(R2.tails.clusterReps, elementMetadata(R2.tails.clusterReps)$s1), function(x){
  # Create a chunking vector for this broad group to distribute the job across {CPUs_per_batch} CPUs.
  elementMetadata(x)$s2 <- dplyr::ntile(1:length(x), CPUs_per_batch)
  
  invisible(parLapply(cluster, split(x, elementMetadata(x)$s2), function(x2){
    library(ShortRead)
    f <- paste0('tmpSeq-',  elementMetadata(x2)$s1[1], '-', elementMetadata(x2)$s2[1], '.fasta')
    writeFasta(x2, file=f)
    cmd <- paste0('blat susScr3.2bit ',  f, ' ', f, '.genome.psl -tileSize=11 -stepSize=9 -minIdentity=85 -maxIntron=5 -minScore=27 -dots=1000 -out=psl -noHead')
    system(cmd)
    file.remove(f)
  }))
}))
stopCluster(cluster)



# Collate and clean up.
# .~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~

system('cat *.genome.psl > genome.psl')
system('rm *.genome.psl')



# Determine which read cluster representatives mapped to the genome.
# .~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~

mapReadNames <- function(file){
  b <- read.table(file, header=FALSE, fill=TRUE)
  names(b) <- c('matches', 'misMatches', 'repMatches', 'nCount', 'qNumInsert', 'qBaseInsert', 'tNumInsert', 'tBaseInsert', 'strand',
                'qName', 'qSize', 'qStart', 'qEnd', 'tName', 'tSize', 'tStart', 'tEnd', 'blockCount', 'blockSizes', 'qStarts', 'tStarts')
  b$percentID <- b$matches/b$qSize
  b2 <- b[which(b$percentID >= 0.95 & b$qStart <= 5 & b$qEnd >= (b$qSize - 5)),]
  bind_rows(lapply(split(b2, b2$qName), function(x){x[1,]}))
}

m <- mapReadNames('genome.psl')

stat7 <- length(unique(m$qName))
stat8 <- (length(unique(m$qName)) / length(unique(names(R2.tails.clusters))))*100


clusters.notGenomic <- clusters[! names(clusters) %in% unique(m$qName)]


# !
clusters.notGenomic <- clusters.notGenomic[order(unlist(lapply(clusters.notGenomic, length)), decreasing = TRUE)]



s <- c('CTGGAAGGTGCTGAGGTACGATGAGACCCGCACCAGGTGCAGACCCTGCGAGTGTGGCGGTAAACATATTAGGAACCAGCCTGTGATGCTGGATGTGACCGAGGAGCTGAGGCCCGATCACTTGGTGCTGGCCTGCACCCGCGCTGAGTTT',
       'CTGGAAGGTGCTGAGGTACGATGAGACCCGCACCAGGTGCAGACCCTGCGAGTGTGGCGGTAAACATATTAGGAACCAGCCTGTGATGCTGGATGTGACCGAGGAGCTGAGGCCCGATCACTTGGTGCTGGCCTGCACCCGCGCTGAGTTT',
       unname(as.character(R2.tails[names(clusters.notGenomic)[1:20]])))
names(s) <- c('TransposonVector', 'TransposaseVector', paste0('Cluster', 1:20))

a <- as.character(msa(s, type='dna', order='input', method="Muscle"))
clusterAlignment <- data.frame(seq = names(a), alignment = a)
row.names(clusterAlignment) <- NULL


n <- 0
clusters.notGenomic.hits <- bind_rows(lapply(names(clusters.notGenomic)[1:20], function(x){
  n <<- n+1
  writeFasta(R2.tails[x], file = 'Cluster.fasta')
  system('blat bothVectors.2bit Cluster.fasta Cluster.psl -tileSize=11 -stepSize=9 -minIdentity=85 -maxIntron=5 -minScore=27 -dots=1000 -out=psl -noHead')
  if(file.info("Cluster.psl")$size == 0){
    return(data.frame(cluster = paste('Cluster', n),
                      target    = '-- No hits to vectors --',
                      seqLength = NA,
                      percentID   = NA,
                      percentCoverage = NA,
                      vectorStart = NA))
  }
  b <- read.table('Cluster.psl', header=FALSE, fill=TRUE)
  names(b) <- c('matches', 'misMatches', 'repMatches', 'nCount', 'qNumInsert', 'qBaseInsert', 'tNumInsert', 'tBaseInsert', 'strand',
                'qName', 'qSize', 'qStart', 'qEnd', 'tName', 'tSize', 'tStart', 'tEnd', 'blockCount', 'blockSizes', 'qStarts', 'tStarts')
  file.remove(c('Cluster.fasta', 'Cluster.psl'))
  b$percentID <- round((b$matches/b$qSize)*100, 2)
  b$percentCoverage <- round(((b$qEnd - b$qStart) / b$qSize)*100, 2)
  b$qName <- paste('Cluster', n)
  b$clusterSize <- length(clusters[[x]])
  b <- dplyr::select(b, qName, clusterSize, tName, qSize, percentID, percentCoverage, tStart)
  names(b) <- c('cluster', 'clusterSize', 'target', 'seqLength', 'percentID', 'percentCoverage', 'vectorStart') 
  b
}))


d <- data.frame(cluster = 1:length(clusters.notGenomic), uniqueSeqs = unlist(lapply(clusters.notGenomic, length)))
d <- d[1:50,]
ggplot(d, aes(d$cluster, d$uniqueSeqs)) + 
  theme_bw() +
  geom_bar(stat='identity') +
  labs(x='Cluster', y='Sequences')






























clusters <- scan('R2.tails.notGenomic.clstr', what = 'character', sep = '\n')
clusters <- clusters[grep('\\*', clusters)]
s <- gsub('\\.\\.\\.', '', unlist(lapply(str_match_all(clusters, '>([^\\s]+)'), function(y){ y[,2] })))
R2.tails.notGenomic.clusters <- R2.tails[s]
###file.remove('R2.tails.fasta')













clusters <- readChar('R2.clstr', file.info('R2.clstr')$size)

parse_cdhitest_output <- function(file){
  # cd-hit-est must be ran with -d 0 flag in order to preserve seq names up to first white space.
  clusters <- readChar(file, file.info(file)$size)
  clusters <- unlist(base::strsplit(clusters, '>Cluster'))
  clusters <- clusters[2:length(clusters)]
  lapply(clusters, function(x){
    gsub('\\.\\.\\.', '', unlist(lapply(str_match_all(x, '>([^\\s]+)'), function(y){ y[,2] })))
  })
}

o <- parse_cdhitest_output('m.clstr')
m2 <- subset(m, qName %in% sapply(o, '[', 1))







# Blat each unique sequence against the porcine genome.
# .~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~

setwd('data')









# Identify read ids that map well to the genome.
# This is a bit challenging in that the first ~20 NT are not expected to match.
# We remove 20 NT from the length calculation for percent id hopeful that this will generally catch
# 20 NT not matching at the start of sequences.
# .~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~



writeFasta(R2[unique(m$qName)], file='m.fasta')

# software/cdhit/cd-hit-est -i m.fasta -o m -c 0.95 -d 0




r <- bind_rows(apply(m2, 1, function(x){
       data.frame(a = substr(R2[x['qName']], 1, as.integer(x['qStart'])-1),
                  b = substr(R2[x['qName']], as.integer(x['qStart']),  as.integer(x['qStart'])+25))
}))






# software/cdhit/cd-hit-est -i unique.R2.fasta -o unique.R2.95 -c 0.95 -d 0

parse_cdhitest_output <- function(file){
  # cd-hit-est must be ran with -d 0 flag in order to preserve seq names up to first white space.
  clusters <- readChar(file, file.info(file)$size)
  clusters <- unlist(base::strsplit(clusters, '>Cluster'))
  clusters <- clusters[2:length(clusters)]
  lapply(clusters, function(x){
    gsub('\\.\\.\\.', '', unlist(lapply(str_match_all(x, '>([^\\s]+)'), function(y){ y[,2] })))
  })
}

o <- parse_cdhitest_output('unique.R2.clstr')
readsToBlat <- R2.trimmed.seq.uni[unlist(lapply(o, '[', 1))]
writeFasta(readsToBlat, file = 'clusterSeqs.fasta')

# blat susScr3.2bit clusterSeqs.fasta clusterSeqs.psl -tileSize=11 -stepSize=9 -minIdentity=85 -maxIntron=5 -minScore=27 -dots=1000 -out=psl -noHead


o <- mapReadNames('clusterSeqs.psl')








# Create a chunking vector to break job into {CPU_batches} broad groups.
# We embed the values into the object to be split so that we can retrieve the values during processing.








mappedReads <- mapReadNames('tmp1-5.psl')
writeFasta(unique(R2.trimmed.seq[mappedReads]), file='tmp.ff')









i <- grep('^ATTATC', as.character(R2.trimmed.seq))
R2.trimmed.seq <- R2.trimmed.seq[i]

R2.trimmed.seq <- R2.trimmed.seq[width(R2.trimmed.seq) >= 75]

R2alignments <- bind_rows(lapply(c(18, 25, 50, 75), function(x){
       h <- data.frame(head(sort(table(as.character(subseq(R2.trimmed.seq, 1, x))), decreasing = TRUE), n = 20))
       p <- (h$Freq / length(R2.trimmed.seq))*100
       data.frame(n_Nucleotides  = x,
                  alignment      = as.character(msa(as.character(h$Var1), type='dna', gapOpening = 10, gapExtension = 4, order = 'input')),
                  ltr_TTTCTAGGG  = grepl('TTTCTAGGG', as.character(h$Var1)),
                  numReads       = h$Freq,
                  percentReads   = round(p, 2),
                  cumPercentRead = round(cumsum(p), 2))
}))


vectors <- sread(readFasta('data/bothVectors.ff'))

# Transposon  ATTATCTTTCTAGGGTTAAaaGATCTGGAAGGTGCTGAGGTACGATGAGACCCGCACCAGGTGCAGACCCTGCGAGTGTGGCGGTAAACATATTAGGAACCAGCCTGTGAT
# Transposase ACTCTAGTCCCCGCGGTGGCAGATCTGGAAGGTGCTGAGGTACGATGAGACCCGCACCAGGTGCAGACCCTGCGAGTGTGGCGGTAAACATATTAGGAACCAGCCTGTGAT




# Export the porcine genome to a 2 bit format for Blat.
#export(BSgenome.Sscrofa.UCSC.susScr3, 'data/susScr3.2bit')

# Create a 2bit data file for both vectors.
# faToTwoBit bothVectors.ff bothVectors.2bit

dataPath <- '/data/internal/geneTherapy/processedRuns/180220_M00281_0323_000000000-BMWYV'

R2.reads <- do.call(c, lapply(list.files(dataPath, 
                                         pattern = '^R2\\-\\d+\\.fa$', 
                                         recursive = TRUE,
                                         full.names = TRUE), function(file){
  o <- rev(unlist(base::strsplit(file, '/')))
  r <- sread(readFasta(file))
  names(r) <- paste0(o[2], '.', str_match(o[1], '\\-(\\d+)')[,2], '.', 1:length(r))
  r
}))




setwd('data')
cluster <- makeCluster(CPUs)

# Creat a splitting vector that is associated with the DNA  data object.
elementMetadata(R2.reads)$s <-  ntile(1:length(R2.reads), CPUs)

# Split the sequencing reads by the splitting the vector and blast each chunk against both the genome and the two accompanying vectors.
invisible(parLapply(cluster, split(R2.reads, elementMetadata(R2.reads)$s), function(x){
  library(ShortRead)
  f <- paste0('R2.', elementMetadata(x)$s[1] , '.fasta')
  writeFasta(x, file=f)
  cmd <- paste0('blat susScr3.2bit ',  f, ' ', f, '.genome.psl -tileSize=11 -stepSize=9 -minIdentity=85 -maxIntron=5 -minScore=27 -dots=1000 -out=psl -noHead')
  system(cmd)
  cmd <- paste0('blat bothVectors.2bit ',  f, ' ', f, '.vectors.psl -tileSize=11 -stepSize=9 -minIdentity=85 -maxIntron=5 -minScore=27 -dots=1000 -out=psl -noHead')
  system(cmd)
}))
stopCluster(cluster)


# Concatenate the results and clean up.
system('cat *.genome.psl > R2.genome.psl')
system('cat *.vectors.psl  > R2.vectors.psl ')

setwd('..')

mapReadNames <- function(file){
  b <- read.table(file, header=FALSE)
  names(b) <- c('matches', 'misMatches', 'repMatches', 'nCount', 'qNumInsert', 'qBaseInsert', 'tNumInsert', 'tBaseInsert', 'strand',
                'qName', 'qSize', 'qStart', 'qEnd', 'tName', 'tSize', 'tStart', 'tEnd', 'blockCount', 'blockSizes', 'qStarts', 'tStarts')
  b$percentID <- (b$matches/b$qSize)*100
  b$percentCoverage <- ((b$qEnd - b$qStart)/b$qSize)*100
  i <- which(b$percentID >= 95 & b$tBaseInsert <= 3 & b$qStart <= 5)
  unique(b$qName[i]) # The same sequence can map well to several genomic locations.
}

readsMappingToGenome  <- mapReadNames('data/R2.genome.psl')
readsMappingToVectors <- mapReadNames('data/R2.vectors.psl')

numReadsNotMapping <- n_distinct(names(R2.reads)) - n_distinct(readsMappingToGenome) - n_distinct(readsMappingToVectors)

nonMapped_R2.reads <- R2.reads[! names(R2.reads) %in% c(readsMappingToGenome, readsMappingToVectors)]

writeFasta(unique(nonMapped_R2.reads), file='data/unique_nonMapped_R2_reads.ff')


