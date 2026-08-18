#Perform a quality check of the FASTQ files with FastQC (external, Python)

setwd("D:/ATAC-seq")
#Reading raw chip-seq data
library(ShortRead)
fqSample <- FastqSampler("D:/ATAC-seq/SRR891269_1.fastq.gz", n = 10^6)
fastq <- yield(fqSample)
fastq
readSequences <- sread(fastq)
readQuality <- quality(fastq)
readIDs <- id(fastq)
readSequences

#Quality of data
readQuality <- quality(fastq)
readQualities <- alphabetScore(readQuality)
readQualities[1:10]

library(ggplot2)
toPlot <- data.frame(ReadQ = readQualities)
ggplot(toPlot, aes(x = ReadQ)) + geom_histogram() + theme_minimal()

readSequences <- sread(fastq)
readSequences_AlpFreq <- alphabetFrequency(readSequences)
readSequences_AlpFreq[1:15, ]

summed__AlpFreq <- colSums(readSequences_AlpFreq)
summed__AlpFreq[c("A", "C", "G", "T", "N")]

readSequences_AlpbyCycle <- alphabetByCycle(readSequences)
readSequences_AlpbyCycle[1:4, 1:10]

AFreq <- readSequences_AlpbyCycle["A", ]
CFreq <- readSequences_AlpbyCycle["C", ]
GFreq <- readSequences_AlpbyCycle["G", ]
TFreq <- readSequences_AlpbyCycle["T", ]

lengths <- sapply(list(AFreq, CFreq, GFreq, TFreq), length)
print(lengths)

toPlot <- data.frame(Count = c(AFreq, CFreq, GFreq, TFreq), Cycle = rep(1:32, 4),
                     Base = rep(c("A", "C", "G", "T"), each = 32))


ggplot(toPlot, aes(y = Count, x = Cycle, colour = Base)) + geom_line() + theme_bw()


qualAsMatrix <- as(readQuality, "matrix")
qualAsMatrix[1:2, ]
boxplot(qualAsMatrix[1:1000, ])

#Filtering

fqStreamer <- FastqStreamer("D:/ATAC-seq/SRR891269_1.fastq.gz", n = 1e+05)

TotalReads <- 0
TotalReadsFilt <- 0
while (length(fq <- yield(fqStreamer)) > 0) {
  TotalReads <- TotalReads + length(fq)
  filt1 <- fq[alphabetScore(fq) > 300]
  filt2 <- filt1[alphabetFrequency(sread(filt1))[, "N"] < 10]
  TotalReadsFilt <- TotalReadsFilt + length(filt2)
  writeFastq(filt2, "filtered_SRR891269_1.fastq", mode = "a")
}

TotalReads

TotalReadsFilt

# Creating a reference genome
library(BSgenome.Hsapiens.UCSC.hg38)
BSgenome.Hsapiens.UCSC.hg38

mainChromosomes <- paste0("chr",c(1:21,"X","Y","M"))
mainChrSeq <- lapply(mainChromosomes,
                     function(x)BSgenome.Hsapiens.UCSC.hg38[[x]])
names(mainChrSeq) <- mainChromosomes
mainChrSeqSet <- DNAStringSet(mainChrSeq)
writeXStringSet(mainChrSeqSet,
                "BSgenome.Hsapiens.UCSC.hg38.mainChrs.fa")

#Creating Rsubread index
library(Rsubread)
buildindex("BSgenome.Hsapiens.UCSC.hg38.mainChrs",
           "BSgenome.Hsapiens.UCSC.hg38.mainChrs.fa",
           indexSplit = TRUE,
           memory = 1000)

#Aligning Sequence Reads

read1 <- "D:/ATAC-seq/SRR891269_1.fastq.gz"
read2 <- "D:/ATAC-seq/SRR891269_2.fastq.gz"

#OR

read1 <- readFastq("D:/ATAC-seq/SRR891269_1.fastq.gz")
read2 <- readFastq("D:/ATAC-seq/SRR891269_2.fastq.gz")
id(read1)[1:2]

#Aligning with Rsubread

Rsubread::align(index = "BSgenome.Hsapiens.UCSC.hg38.mainChrs",
                readfile1 = read1,
                readfile2 = read2,
                output_file = "ATAC_50K_2.bam",
                nthreads = 4, type = 1, unique = TRUE, maxFragLength = 2000)

# type = 1 means genomic DNA-seq, whereas type = 0/"rna" is RNA-seq
# Rsubread's default paired-end maximum fragment length is 600 bp
#  To control maximum allowed fragment lengths, I set the maxFragLength parameter to 2000. 
# I also set the unique parameter to TRUE to only include uniquely mapping reads

bam.files <- list.files(pattern = ".BAM$", full.names = TRUE) 
bam.files

library(Rsamtools)
sortBam("ATAC_50K_2.bam", "Sorted_ATAC_50K_2")
indexBam("Sorted_ATAC_50K_2.bam")

# Distribution of mapped reads
# mappedReads <- idxstatsBam(sortedBAM)

# Check the number of mapped reads on every chromosome using the idxstatsBam() function.
mappedReads <- idxstatsBam("Sorted_ATAC_50K_2.bam")
head(mappedReads)

library(ggplot2)
ggplot(mappedReads, aes(seqnames, mapped, fill = seqnames)) + geom_bar(stat = "identity") +
  coord_flip()

# Check:
#   - Number of mapped reads
#   - Chromosomal distribution
#   - Unmapped reads
#
# Before filtering:
#   Examine the overall mapping rate.
#
# A very low mapping rate may indicate:
#   - poor sequence quality
#   - incorrect reference genome
#   - contamination
#   - adapter/technical problems
#   - poor alignment parameters


# ============================================================
# MITOCHONDRIAL READ QC
# ============================================================
# Check the fraction of mapped reads originating from chrM.
#
# High mitochondrial content can indicate poor nuclear library
# Complexity.
#
# First quantify it.
# Then decide whether mitochondrial reads should be removed for downstream nuclear ATAC-seq analysis.

chrM_reads <- mappedReads[
    mappedReads$seqnames == "chrM",
    "mapped"
]

total_mapped <- sum(mappedReads$mapped)

mitochondrial_fraction <- (
    chrM_reads / total_mapped
) * 100

mitochondrial_fraction

# ------------------------------------------------------------
# OPTIONAL: remove chrM reads
# ------------------------------------------------------------
# Use this when the downstream analysis focuses on nuclear chromatin accessibility.
samtools view -b \
  -o Sorted_ATAC_50K_2_noChrM.bam \
  Sorted_ATAC_50K_2.bam \
  chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY
samtools sort \
  -o Sorted_ATAC_50K_2_noChrM_sorted.bam \
  Sorted_ATAC_50K_2_noChrM.bam

samtools index Sorted_ATAC_50K_2_noChrM_sorted.bam

library(Rsamtools)

mappedReads_noChrM <- idxstatsBam(
    "Sorted_ATAC_50K_2_noChrM_sorted.bam"
)

mappedReads_noChrM[
    mappedReads_noChrM$seqnames == "chrM",
]
# For large BAM files, use samtools or another BAM filtering
# approach rather than loading the entire BAM into memory.

#____________Post-alignment processing______________#
## Proper pairs ##
# ATAC-seq fragment analyses generally use properly paired reads.
#
# Check the proportion of properly paired reads first.
#
# Then use properly paired fragments for downstream
# fragment-based analyses.
library(GenomicAlignments)                    
flags = scanBamFlag(isProperPair = TRUE)
myParam = ScanBamParam(flag = flags, what = c("qname", "mapq", "isize"), which = GRanges("chr20",
                                                                                         IRanges(1, 63025520)))
myParam

sortedBAM <- "Sorted_ATAC_50K_2.bam"
atacReads <- readGAlignmentPairs(sortedBAM, param = myParam)
class(atacReads)

# GAlignmentPairs
atacReads[1:2, ]

read1 <- first(atacReads)
read2 <- second(atacReads)
read2[1, ]


# ============================================================
# MAPQ QC + OPTIONAL FILTERING
# ============================================================
# Check MAPQ distribution before choosing a cutoff.
#
# MAPQ represents confidence in read placement.
#
# Common starting thresholds:
#   MAPQ >= 20
#   MAPQ >= 30
#
# The cutoff should depend on the aligner and dataset.
#
# Your Rsubread alignment already requests unique mapping,
# but MAPQ filtering can provide an additional quality filter.

#Retrieve MapQ scores
read1MapQ <- mcols(read1)$mapq
read2MapQ <- mcols(read2)$mapq
read1MapQ[1:2]
#MapQ score frequencies
read1MapQFreqs <- table(read1MapQ)
read2MapQFreqs <- table(read2MapQ)
read1MapQFreqs
read2MapQFreqs
# Plot MapQ scores
library(ggplot2)
toPlot <- data.frame(MapQ = c(names(read1MapQFreqs), names(read2MapQFreqs)), Frequency = c(read1MapQFreqs,
                                                                                           read2MapQFreqs), Read = c(rep("Read1", length(read1MapQFreqs)), rep("Read2",
                                                                                                                                                               length(read2MapQFreqs))))
toPlot$MapQ <- factor(toPlot$MapQ, levels = unique(sort(as.numeric(toPlot$MapQ))))
ggplot(toPlot, aes(x = MapQ, y = Frequency, fill = MapQ)) + geom_bar(stat = "identity") +
  facet_grid(~Read)

# Example:
# Keep fragments where both mates have MAPQ >= 30.
#
# NOTE:
# Apply this only after inspecting the MAPQ distribution.
length(atacReads)
MAPQ_cutoff <- 30

atacReads <- atacReads[
    mcols(first(atacReads))$mapq >= MAPQ_cutoff &
    mcols(second(atacReads))$mapq >= MAPQ_cutoff,
]

length(atacReads)
                     
# ============================================================
#|    MAPQ | Interpretation         | Typical action      |   
#| ------: | ---------------------- | ------------------- |
#|       0 | Ambiguous/poor mapping | Remove              |
#|     1–9 | Very low confidence    | Remove              |
#|   10–19 | Low confidence         | Usually remove      |
#|   20–29 | Moderate               | Depends on analysis |
#| **≥30** | **High confidence**    | **Usually retain**  |
#|     ≥40 | Very high confidence   | Retain              |
# ============================================================

# ============================================================
# FRAGMENT-SIZE QC
# ============================================================
# ATAC-seq should show nucleosome-associated periodicity.
#
# Expected approximate pattern:
#
#   <100 bp        nucleosome-free
#   ~180–250 bp    mononucleosome
#   ~300–500 bp    dinucleosome
#   larger sizes:   higher-order nucleosomes
#
# Check the distribution before filtering.
#Retrieving insert sizes
atacReads_read1 <- first(atacReads)
insertSizes <- abs(elementMetadata(atacReads_read1)$isize)
head(insertSizes)

#Plotting insert sizes
fragLenSizes <- table(insertSizes)
fragLenSizes[1:5]

library(ggplot2)
toPlot <- data.frame(InsertSize = as.numeric(names(fragLenSizes)), Count = as.numeric(fragLenSizes))
fragLenPlot <- ggplot(toPlot, aes(x = InsertSize, y = Count)) + geom_line()
fragLenPlot + theme_bw()

fragLenPlot + scale_y_continuous(trans = "log2") + theme_bw()

fragLenPlot + scale_y_continuous(trans = "log2") + geom_vline(xintercept = c(180,
                                                                             247), colour = "red") + geom_vline(xintercept = c(315, 437), colour = "darkblue") +
  geom_vline(xintercept = c(100), colour = "darkgreen") + theme_bw()

# ============================================================
# FRAGMENT-SIZE FILTERING
# ============================================================
# These filters are NOT general "bad-read" filters.
# They separate biological fragment classes.
#
# Nucleosome-free:
#   <100 bp
#
# Mononucleosome:
#   180–240 bp
#
# Dinucleosome:
#   315–437 bp
#
# Exact boundaries can be adjusted according to the observed
# fragment-size distribution.
# Subsetting ATAC-seq reads by insert sizes

atacReads_NucFree <- atacReads[insertSizes < 100, ]
atacReads_MonoNuc <- atacReads[insertSizes > 180 & insertSizes < 240, ]
atacReads_diNuc <- atacReads[insertSizes > 315 & insertSizes < 437, ]

#Creating BAM files split by insert sizes
nucFreeRegionBam <- gsub("\\.bam", "_nucFreeRegions\\.bam", sortedBAM)
monoNucBam <- gsub("\\.bam", "_monoNuc\\.bam", sortedBAM)
diNucBam <- gsub("\\.bam", "_diNuc\\.bam", sortedBAM)
library(rtracklayer)
export(atacReads_NucFree, nucFreeRegionBam, format = "bam")
export(atacReads_MonoNuc, monoNucBam, format = "bam")
export(atacReads_diNuc, diNucBam, format = "bam")

# ============================================================
# DUPLICATE / LIBRARY COMPLEXITY QC
# ============================================================
# Assess duplicate fragments.
#
# High duplication can indicate:
#   - excessive PCR amplification
#   - low library complexity
#
# However, duplication should not automatically be interpreted
# as a reason to discard the library.

# Creating fragment GRanges
atacReads[1, ]

atacFragments <- granges(atacReads)
atacFragments[1, ]

# We can use the duplicated() function to identify the non-redundant (non-duplicate) fraction of our full length fragments.
duplicatedFragments <- sum(duplicated(atacFragments))
totalFragments <- length(atacFragments)
duplicateRate <- duplicatedFragments/totalFragments
nonRedundantFraction <- 1 - duplicateRate
nonRedundantFraction

# ============================================================
# OPTIONAL DUPLICATE REMOVAL
# ============================================================
# Only remove duplicates after evaluating library complexity.
#
# For quantitative ATAC-seq analyses, duplicate removal is
# commonly considered to prevent PCR-amplified fragments from
# artificially increasing accessibility.
#
# Do NOT remove duplicates simply because the duplicate rate
# is high without checking the overall library quality.

# Creating an open region bigWig
openRegionRPMBigWig <- gsub("\\.bam", "_openRegionRPM\\.bw", sortedBAM)
myCoverage <- coverage(atacFragments, weight = (10^6/length(atacFragments)))
export.bw(myCoverage, openRegionRPMBigWig)


#__________________ ATACseqQC_________________#
# Since this can be fairly memory-heavy, I am just illustrating it here on a BAM file containing just the chromosome 17 reads of the ATACseq data

# Load required package
library(Rsamtools)

# Input BAM file
input_bam <- "Sorted_ATAC_50K_2.bam"
output_bam <- "Sorted_ATAC_50K_2_ch17.bam"

# Define region (chr17)
param <- ScanBamParam(which = GRanges("chr17", IRanges(1, seqlengths(BamFile(input_bam))["chr17"])))

# Filter reads for chr17 and write to new BAM
filterBam(file = input_bam, destination = output_bam, param = param)

# Index the new BAM
indexBam(output_bam)

# ============================================================
#  PCR BOTTLENECK COEFFICIENT
# ============================================================
# Use ATACseqQC to assess library complexity and PCR bottlenecking.
#
# Check:
#   - PCR bottleneck coefficient
#   - Whether multiple QC metrics agree
#
# A poor PBC should be interpreted together with:
#   - duplicate rate
#   - sequencing depth
#   - FRiP
#   - TSS enrichment
#   - fragment-size distribution

#BiocManager::install("ATACseqQC")
library(ATACseqQC)
ATACQC <- bamQC("Sorted_ATAC_50K_2_ch17.bam")
names(ATACQC)

#PCR bottleneck coefficients

# PCR bottleneck coefficients identify PCR bias/overamplification which may have occurred in preparation of ATAC samples.
# The "PCRbottleneckCoefficient_1" is calculated as the number of positions in the genome with exactly 1 read mapped uniquely compared to the number of positions with at least 1 read.
# For example, if we have 20 reads. 16 map uniquely to locations. 4 do not map uniquely; instead, there are 2 locations, both of which have 2 reads. This would lead us to calculate 16/18. We therefore have a PBC1 of 0.889

# ** Values less than 0.7 indicate severe bottlenecking; between 0.7 and 0.9 indicate moderate bottlenecking. Values greater than 0.9 show no bottlenecking. **
ATACQC$PCRbottleneckCoefficient_1

# The "PCRbottleneckCoefficient_2" is our secondary measure of bottlenecking. 
# It is calculated as the number of positions in the genome with exactly 1 read mapped uniquely compared to the number of positions with exactly 2 reads mapping uniquely.
# We can reuse our example. If we have 20 reads, 16 of which map uniquely. 4 do not map uniquely; instead, there are 2 locations, both of which have 2 reads. This would lead us to calculate 16/2. We therefore have a PBC2 of 8.
# ** Values less than 1 indicate severe bottlenecking, between 1 and 3 indicate moderate bottlenecking. Greater than 3 show no bottlenecking.
ATACQC$PCRbottleneckCoefficient_2


#Plotting signal over regions in R
BiocManager::install("soGGi")
library(soGGi)

# ============================================================
# TSS ENRICHMENT QC
# ============================================================
# TSS enrichment is a major ATAC-seq library-level QC metric.
#
# Check:
#   - Aggregate signal around TSSs.
#   - Presence of a clear enrichment peak at the TSS.
#
# A weak/flat profile can indicate poor ATAC-seq enrichment.
#
# IMPORTANT:
#   TSS enrichment is primarily used to evaluate the library.
#   It is not a read-level filtering criterion.

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
TxDb.Hsapiens.UCSC.hg38.knownGene

# Extract gene locations (TSS to TTS) using the genes() function and our TxDb object.
genesLocations <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
genesLocations

tssLocations <- resize(genesLocations, fix = "start", width = 1)
tssLocations

mainChromosomes <- paste0("chr", c(1:22, "X", "Y", "M"))

myindex <- match(seqnames(tssLocations), mainChromosomes)
tssLocations <- tssLocations[!is.na(myindex)]
library(GenomicRanges)
tssLocations <- keepSeqlevels(tssLocations, mainChromosomes, pruning.mode = "coarse")

myindex <- (match(seqnames(tssLocations), mainChromosomes))
tssLocations <- tssLocations[as.numeric(myindex)]
seqlevels(tssLocations) <- mainChromosomes

library(soGGi)
sortedBAM <- "Sorted_ATAC_50K_2.bam"
library(Rsamtools)
# Nucleosome-free
#bamFile parameter and a GRanges to plot over
allSignal <- regionPlot(bamFile = sortedBAM, testRanges = tssLocations)

nucFree <- regionPlot(bamFile = sortedBAM, testRanges = tssLocations, style = "point",
                      format = "bam", paired = TRUE, minFragmentLength = 0, maxFragmentLength = 100,
                      forceFragment = 50)
class(nucFree)
plotRegion(nucFree)

# We can create a plot for our mono-nucleosome signal by adjusting our minFragmentLength and maxFragmentLength parameters those expected for nucleosome length fragments (here 180 to 240)
monoNuc <- regionPlot(bamFile = sortedBAM, testRanges = tssLocations, style = "point",
format = "bam", paired = TRUE, minFragmentLength = 180, maxFragmentLength = 240,
forceFragment = 80)
plotRegion(monoNuc) 

                     
# Peak calling for nucleosome-free regions
# Single-end peak calling (With single-end sequencing from ATAC-seq, we do not know how long the fragments are. 
# To identify open regions, therefore, requires some different parameters for MACS2 peak calling compared to ChIPseq.
# One strategy employed is to shift read 5’ ends by -100 and then extend from this by 200bp. 
# Considering the expected size of our nucleosome-free fragments, this should provide a pile-up over nucleosome regions suitable for MACS2 window size.)

#MACS2 callpeak -t singleEnd.bam --nomodel --shift -100--extsize 200 --format BAM -g MyGenom


#OR (Alternatively, for the nucleosome-occupied data, we can adjust shift and extension to centre the signal on nucleosome centres; nucleosomes wrapped in 147bp of DNA)

#MACS2 callpeak -t singleEnd.bam --nomodel --shift 37--extsize 73 --format BAM -g MyGenome


# Paired-end peak calling (Externally in Python)

macs2 callpeak \
-t Sorted_ATAC_50K_2.bam \
-f BAMPE \
-g hs \
-n Sorted_ATAC_50K_2_Small_Paired \
--outdir PeakDirectory/ATAC_50K_2


#QC
library(ChIPQC)
library(rtracklayer)
library(DT)
library(dplyr)
library(tidyr)
                     
blkList <- import.bed("ENCFF356LFX.bed.gz")
openRegionPeaks <- "D:/ATAC-seq/PeakDirectory/ATAC_50K_2/Sorted_ATAC_50K_2_Small_Paired_peaks.narrowPeak"

# ============================================================
# FRiP QC
# ============================================================
# FRiP = Fraction of Reads in Peaks.
#
# Check:
#   - FRiP for every library.
#   - Relative consistency between biological replicates.
#
# Higher FRiP generally indicates stronger enrichment.
#
# IMPORTANT:
#   FRiP is mainly a library-level QC metric.
#   Do not filter individual reads based on FRiP.
qcRes <- ChIPQCsample("Sorted_ATAC_50K_2.bam",
                      peaks = openRegionPeaks, annotation = "hg38", blacklist = blkList,
                      verboseT = FALSE)

# Summary statistics
summary(qcRes)

# Generate QC report
ChIPQCreport(qcRes, reportName = "ATAC_QC", reportFolder = "QC_Report")

#  Reads in peaks and reads in blacklist from the QCmetrics() function
myMetrics <- QCmetrics(qcRes)
myMetrics[c("RiBL%", "RiP%")]

#  Number of duplicated reads from the flagtagcounts() function
flgCounts <- flagtagcounts(qcRes)
DupRate <- flgCounts["DuplicateByChIPQC"]/flgCounts["Mapped"]
DupRate * 100

#Remove blacklisted peaks
library(GenomicRanges)
library(rtracklayer)

# Extract all MACS2-called peaks as GRanges
MacsCalls <- granges(qcRes)

# Count how many peaks overlap the ENCODE blacklist
blacklist_summary <- data.frame(
  Blacklisted = sum(MacsCalls %over% blkList),
  Not_Blacklisted = sum(!(MacsCalls %over% blkList))
)

# Filter out blacklisted peaks
MacsCalls <- MacsCalls[!(MacsCalls %over% blkList)]

# View the summary
blacklist_summary
length(MacsCalls)  # total number of peaks after blacklist removal


# Define output file path
filtered_narrowpeak <- "PeakDirectory/filtered/Sorted_ATAC_50K_2_Small_Paired_filtered_peaks.narrowPeak"

# Export filtered peaks (GRanges) to narrowPeak format
export(MacsCalls, filtered_narrowpeak, format = "narrowPeak")

                     
#Annotating peaks to genes
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
MacsCalls_Anno <- annotatePeak(MacsCalls, TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb = "org.Hs.eg.db")
MacsCalls_Anno

library(openxlsx)

# Convert annotation object to data frame
MacsCalls_Anno_df <- as.data.frame(MacsCalls_Anno)

# Keep both gene IDs and gene names
# 'SYMBOL' column is automatically added when annoDb = "org.Hs.eg.db" is used
MacsCalls_Anno_df <- MacsCalls_Anno_df %>%
  dplyr::select(seqnames, start, end, width, annotation, geneId, SYMBOL, distanceToTSS, everything())


# Write to Excel
write.xlsx(MacsCalls_Anno_df, file = "Annotated_ATAC_Peaks1.xlsx", rowNames = FALSE)

plotAnnoPie(MacsCalls_Anno)

# Annotating nucleosome free regions
# With this information, we can then subset our peaks/nuc free regions to those only landing in TSS regions (+/- 500)
MacsGR_Anno <- as.GRanges(MacsCalls_Anno)
MacsGR_TSS <- MacsGR_Anno[abs(MacsGR_Anno$distanceToTSS) < 500]
MacsGR_TSS[1, ]

#Functional analysis of peaks
library(rGREAT)
great_Job <- submitGreatJob(MacsCalls, species = "hg38")
availableCategories(great_Job)


great_ResultTable = getEnrichmentTables(great_Job, category = "GO")
names(great_ResultTable)

great_ResultTable[["GO Biological Process"]][1:4, ]


#_________Finding Motifs in ATACseq data_________#

# Cutting sites from ATACseq data
#Finding motifs

library(MotifDb)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
CTCF <- query(MotifDb, c("CTCF"))
CTCF

names(CTCF)

ctcfMotif <- CTCF[[1]]
ctcfMotif[, 1:4]

# Visualising PWMs
library(seqLogo)
seqLogo(ctcfMotif)

#Searching for PWMs in DNAstring
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)

# Initialize list to store hits per chromosome
allChrHits <- list()

# Loop over all chromosomes
for (chrName in seqnames(BSgenome.Hsapiens.UCSC.hg38)) {
  cat("Scanning", chrName, "...\n")
  chrSeq <- BSgenome.Hsapiens.UCSC.hg38[[chrName]]
  
  # Match PWM (adjust min.score as needed)
  hits <- matchPWM(ctcfMotif, chrSeq, min.score = "80%")
  
  #OR
  library(BiocParallel)
  mainChromosomes <- paste0("chr", c(1:22, "X", "Y", "M"))
  bpparam <- MulticoreParam(workers = 4) # or SnowParam for Windows
  
  scanChr <- function(chrName) {
    seq <- BSgenome.Hsapiens.UCSC.hg38[[chrName]]
    hits <- matchPWM(ctcfMotif, seq, min.score = "80%")
    if (length(hits) == 0) return(NULL)
    GRanges(seqnames = chrName, ranges = ranges(hits))
  }
  
  resList <- bplapply(mainChromosomes, scanChr, BPPARAM = bpparam)
  myRes <- do.call(c, resList)
  
  # Convert hits to GRanges if any found
  if (length(hits) > 0) {
    allChrHits[[chrName]] <- GRanges(
      seqnames = chrName,
      ranges = ranges(hits)
    )
  }
}

# Combine all chromosomes into one GRanges object
myRes <- do.call(c, allChrHits)

# Check results
myRes
length(myRes)  # total motif hits genome-wide

library(GenomicRanges)

# If myRes is the GRanges object from genome-wide matchPWM
toCompare <- myRes  # already a GRanges object

# If myRes is a Views/XStringViews object from matchPWM per chromosome
# Example: converting single chromosome hits
# toCompare <- GRanges(seqnames = "chr20", ranges = ranges(myRes))

# Check
toCompare
length(toCompare)


#Shifting reads for cut-sites

#To plot cut-sites we will wish to consider only the 5’ end of reads and will need to adjust for a known o set of 5’ reads to actual T5 cut sites.
#This will involve capturing the 5’end of reads and shifting reads on positive and negative strand by 4bp or -5bp respectively.
#First we read in our nucleosome free region BAM le and extract read pairs.

BAM <- "Sorted_ATAC_50K_2.bam"
atacReads_Open <- readGAlignmentPairs(BAM)
read1 <- first(atacReads_Open)
read2 <- second(atacReads_Open)
read2[1, ]

Firsts <- resize(granges(read1), fix = "start", 1)
First_Pos_toCut <- shift(granges(Firsts[strand(read1) == "+"]), 4)
First_Neg_toCut <- shift(granges(Firsts[strand(read1) == "-"]), -5)
Seconds <- resize(granges(read2), fix = "start", 1)
Second_Pos_toCut <- shift(granges(Seconds[strand(read2) == "+"]), 4)
Second_Neg_toCut <- shift(granges(Seconds[strand(read2) == "-"]), -5)
test_toCut <- c(First_Pos_toCut, First_Neg_toCut, Second_Pos_toCut, Second_Neg_toCut)
test_toCut[1:2, ]

# Coverage for cut-sites

cutsCoverage <- coverage(test_toCut)
cutsCoverage20 <- cutsCoverage["chr20"]
cutsCoverage20[[1]]

# Plotting for cut-sites
CTCF_Cuts_open <- regionPlot(cutsCoverage20, testRanges = toCompare, style = "point",
                             format = "rlelist", distanceAround = 500)
plotRegion(CTCF_Cuts_open, outliers = 0.001) + ggtitle("NucFree Cuts Centred on CTCF") +
  theme_bw()



