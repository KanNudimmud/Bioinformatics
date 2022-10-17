## NGS Analysis Basics
## Package Requirements
source("https://bioconductor.org/biocLite.R")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("Biostrings", "GenomicRanges",
                       "rtracklayer", "systemPipeR", 
                       "seqLogo", "ShortRead"))

### Basic String Matching and Parsing
## String matching
# Generate sample sequence data set
myseq <- c("ATGCAGACATAGTG", "ATGAACATAGATCC", "GTACAGATCAC")

# String searching with regular expression support
myseq[grep("ATG", myseq)] 

# Searches myseq for first match of pattern “AT”
pos1 <- regexpr("AT", myseq) 
as.numeric(pos1); attributes(pos1)$match.length # Returns position information of matches

# Searches myseq for all matches of pattern “AT”
pos2 <- gregexpr("AT", myseq) 
as.numeric(pos2[[1]]); attributes(pos2[[1]])$match.length # Returns positions of matches in first sequence

# String substitution with regular expression support
gsub("^ATG", "atg", myseq) 

## Positional parsing
# Computes length of strings
nchar(myseq) 

# Positional parsing of several fragments from one string
substring(myseq[1], c(1,3), c(2,5)) 

# Positional parsing of many strings
substring(myseq, c(1,4,7), c(2,6,10)) 

### Random Sequence Generation
# Random DNA sequences of any length
rand <- sapply(1:100, function(x) paste(sample(c("A","T","G","C"), 12, replace=TRUE), collapse=""))
rand[1:3]

# Count identical sequences
table(c(rand[1:4], rand[1]))

## Extract reads from reference
library(Biostrings)
ref <- DNAString(paste(sample(c("A","T","G","C"), 100000, replace=T), collapse=""))
randstart <- sample(1:(length(ref)-15), 1000)
randreads <- Views(ref, randstart, width=15)
rand_set <- DNAStringSet(randreads)
unlist(rand_set)
 
### Sequences in Bioconductor
## Sequence Import and Export
# Download the following sequences to your current working directory and then import them into R:
dir.create("data", showWarnings = FALSE)
download.file("https://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Bacteria/Halobacterium_sp_uid217/AE004437.ffn", "data/AE004437.ffn")

# Import FASTA file with readDNAStringSet
myseq <- readDNAStringSet("data/AE004437.ffn")
myseq[1:3]

# Subset sequences with regular expression on sequence name field
sub <- myseq[grep("99.*", names(myseq))]
length(sub)

# Export subsetted sequences to FASTA file
writeXStringSet(sub, file="./data/AE004437sub.ffn", width=80)

## Working with XString Containers
# The XString stores the different types of biosequences in dedicated containers
library(Biostrings)
d <- DNAString("GCATAT-TAC")
d
d[1:4]

# RNA sequences
r <- RNAString("GCAUAU-UAC") 
r <- RNAString(d) # Converts d to RNAString object
r

# Protein sequences
p <- AAString("HCWYHH")
p

# Any type of character strings
b <- BString("I store any set of characters. Other XString objects store only the IUPAC characters.")
b

## Working with XStringSet Containers
# XStringSet containers allow to store many biosequences in one object
dset <- DNAStringSet(c("GCATATTAC", "AATCGATCC", "GCATATTAC")) 
names(dset) <- c("seq1", "seq2", "seq3") # Assigns names
dset[1:2]

# Important utilities for XStringSet containers
width(dset) # Returns the length of each sequences

# The [[ subsetting operator returns a single entry as XString object
d <- dset[[1]] 

# Appends/concatenates two XStringSet objects
dset2 <- c(dset, dset) 

# Converts XStringSet to named vector
dsetchar <- as.character(dset) 

# Collapses many sequences to a single one stored in a DNAString container
dsetone <- unlist(dset) 

# Sequence subsetting by positions
DNAStringSet(dset, start=c(1,2,3), end=c(4,8,5))

## Multiple Alignment Class
# The XMultipleAlignment class stores the different types of multiple sequence alignments
origMAlign <- readDNAMultipleAlignment(filepath = system.file("extdata",
                                                              "msx2_mRNA.aln", package = "Biostrings"), format = "clustal")
origMAlign

## Basic Sequence Manipulations
# Reverse and Complement
randset <- DNAStringSet(rand)
complement(randset[1:2])

reverse(randset[1:2])

reverseComplement(randset[1:2])

## Translate DNA into Protein
translate(randset[1:2])

### Pattern Matching
## Pattern matching with mismatches 
# Find pattern matches in reference
myseq1 <- readDNAStringSet("./data/AE004437.ffn") 
mypos <- matchPattern("ATGGTG", myseq1[[1]], max.mismatch=1) 

# Count only the corresponding matches
countPattern("ATGGCT", myseq1[[1]], max.mismatch=1)

# Count matches in many sequences
vcountPattern("ATGGCT", myseq1, max.mismatch=1)[1:20]

# Results shown in DNAStringSet object
tmp <- c(DNAStringSet("ATGGTG"), DNAStringSet(mypos)) 

# Return a consensus matrix for query and hits
consensusMatrix(tmp)[1:4,] 

# Find all pattern matches in reference
myvpos <- vmatchPattern("ATGGCT", myseq1, max.mismatch=1) 
myvpos # The results are stored as MIndex object.

# Retrieves the result for single entry
Views(myseq1[[1]], start(myvpos[[1]]), end(myvpos[[1]])) 

# Return all matches
sapply(names(myseq1), function(x) 
  as.character(Views(myseq1[[x]], start(myvpos[[x]]), end(myvpos[[x]]))))[1:4] 

## Pattern matching with regular expression support
myseq <- DNAStringSet(c("ATGCAGACATAGTG", "ATGAACATAGATCC", "GTACAGATCAC"))
# String searching with regular expression support
myseq[grep("^ATG", myseq, perl=TRUE)] 

# Searches 'myseq' for first match of pattern "AT"
pos1 <- regexpr("AT", myseq) 
# Returns position information of matches
as.numeric(pos1); attributes(pos1)$match.length 

# Searches 'myseq' for all matches of pattern "AT"
pos2 <- gregexpr("AT", myseq) 
# Match positions in first sequence
as.numeric(pos2[[1]]); attributes(pos2[[1]])$match.length 

# String substitution with regular expression support
DNAStringSet(gsub("^ATG", "NNN", myseq)) 

## PWM Viewing and Searching
# Plot with seqLogo
library(seqLogo) 
pwm <- PWM(DNAStringSet(c("GCT", "GGT", "GCA"))) 
pwm

seqLogo(t(t(pwm) * 1/colSums(pwm)))

# Search sequence for PWM matches with score better than min.score
chr <- DNAString("AAAGCTAAAGGTAAAGCAAAA") 
matchPWM(pwm, chr, min.score=0.9) 

### NGS Sequences
## Sequence and Quality Data: QualityScaleXStringSet
# Phred score interconversion
phred <- 1:9
phreda <- paste(sapply(as.raw((phred)+33), rawToChar), collapse="")
phreda

as.integer(charToRaw(phreda))-33 

# Construct QualityScaledDNAStringSet from scratch
dset <- DNAStringSet(sapply(1:100, function(x) paste(sample(c("A","T","G","C"), 20, replace=T), collapse=""))) # Creates random sample sequence.
myqlist <- lapply(1:100, function(x) sample(1:40, 20, replace=T)) # Creates random Phred score list.
myqual <- sapply(myqlist, function(x) toString(PhredQuality(x))) # Converts integer scores into ASCII characters.
myqual <- PhredQuality(myqual) # Converts to a PhredQuality object.
dsetq1 <- QualityScaledDNAStringSet(dset, myqual) # Combines DNAStringSet and quality data in QualityScaledDNAStringSet object.
dsetq1[1:2]

## Processing FASTQ Files with ShortRead
library(ShortRead)
download.file("http://cluster.hpcc.ucr.edu/~tgirke/HTML_Presentations/Manuals/testdata/samplefastq/data.zip", "data.zip")
unzip("data.zip")

# Important utilities for accessing FASTQ files
fastq <- list.files("data", "*.fastq$"); fastq <- paste("data/", fastq, sep="")
names(fastq) <- paste("flowcell6_lane", 1:length(fastq), sep="_") 
(fq <- readFastq(fastq[1])) # Imports first FASTQ file

# Counts numbers of reads in FASTQ files
countLines(dirPath="./data", pattern=".fastq$")/4 

id(fq)[1] # Returns ID field

sread(fq)[1] # Returns sequence

quality(fq)[1] # Returns Phred scores 

# Coerces Phred scores to numeric matrix
as(quality(fq), "matrix")[1:4,1:12] 

# Constructs a ShortReadQ from components
ShortReadQ(sread=sread(fq), quality=quality(fq), id=id(fq)) 

## FASTQ Quality Reports
# Using systemPipeR
library(systemPipeR)
fqlist <- seeFastq(fastq=fastq, batchsize=800, klength=8) # For real data set batchsize to at least 10^5 
seeFastqPlot(fqlist)

# Using ShortRead
sp <- SolexaPath(system.file('extdata', package='ShortRead'))
fl <- file.path(analysisPath(sp), "s_1_sequence.txt") 
fls <- c(fl, fl) 
coll <- QACollate(QAFastqSource(fls), QAReadQuality(), QAAdapterContamination(), 
                  QANucleotideUse(), QAQualityUse(), QASequenceUse(), QAFrequentSequence(n=10), 
                  QANucleotideByCycle(), QAQualityByCycle())
x <- qa2(coll, verbose=TRUE)
res <- report(x)
if(interactive())
  browseURL(res) 

## Filtering and Trimming FASTQ Files with ShortRead
# Adaptor trimming
fqtrim <- trimLRPatterns(Rpattern="GCCCGGGTAA", subject=fq)
sread(fq)[1:2] # Before trimming

sread(fqtrim)[1:2] # After trimming

# Read counting and duplicate removal
tables(fq)$distribution # Counts read occurences

sum(srduplicated(fq)) # Identifies duplicated reads

fq[!srduplicated(fq)]

# Trimming low quality tails
cutoff <- 30
cutoff <- rawToChar(as.raw(cutoff+33))
sread(trimTails(fq, k=2, a=cutoff, successive=FALSE))[1:2]

# Removal of reads with Phred scores below a threshold value
cutoff <- 30
qcount <- rowSums(as(quality(fq), "matrix") <= 20) 
fq[qcount == 0] # Number of reads where all Phred scores >= 20

# Removal of reads with x Ns and/or low complexity segments
filter1 <- nFilter(threshold=1) # Keeps only reads without Ns
# Removes reads with nucleotide bias, >=20 of any base
filter2 <- polynFilter(threshold=20, nuc=c("A","T","G","C")) 
filter <- compose(filter1, filter2)
fq[filter(fq)]

## Memory Efficient FASTQ Processing
# Streaming through FASTQ files with FastqStreamer and random sampling reads with FastqSampler
fq <- yield(FastqStreamer(fastq[1], 50)) # Imports first 50 reads 
fq <- yield(FastqSampler(fastq[1], 50)) # Random samples 50 reads 

# Streaming through a FASTQ file while applying filtering/trimming functions and writing the results to a new file here SRR038845.fastq_sub in data directory.
f <- FastqStreamer(fastq[1], 50) 
while(length(fq <- yield(f))) {
  fqsub <- fq[grepl("^TT", sread(fq))] 
  writeFastq(fqsub, paste(fastq[1], "sub", sep="_"), mode="a", compress=FALSE)
}
close(f)

### Range Operations
## Range Data Are Stored in IRanges and GRanges Containers
# Construct GRanges Object
library(GenomicRanges); library(rtracklayer)
gr <- GRanges(seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)), ranges = IRanges(1:10, end = 7:16, names = head(letters, 10)), strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)), score = 1:10, GC = seq(1, 0, length = 10)) # Example of creating a GRanges object with its constructor function.

# Import GFF into GRanges Object
gff <- import.gff("http://cluster.hpcc.ucr.edu/~tgirke/Documents/R_BioCond/Samples/gff3.gff") # Imports a simplified GFF3 genome annotation file.
seqlengths(gff) <- end(ranges(gff[which(values(gff)[,"type"]=="chromosome"),])) 
names(gff) <- 1:length(gff) # Assigns names to corresponding slot
gff[1:4,]

seqinfo(gff)

# Coerce GRanges object to data.frame
as.data.frame(gff)[1:4, 1:7]

### Utilities for Range Containers
## Accessor and subsetting methods for GRanges objects
# Subsetting and replacement
gff[1:4]

gff[1:4, c("type", "ID")] 

gff[2] <- gff[3] 

# GRanges objects can be concatenated with the c function
c(gff[1:2], gff[401:402])

# Acessor functions
seqnames(gff)

ranges(gff)

strand(gff)

seqlengths(gff) 

start(gff[1:4])

end(gff[1:4])

width(gff[1:4]) 

# Accessing metadata component
values(gff) # or elementMetadata(gff)

values(gff)[, "type"][1:20] 

gff[values(gff)[ ,"type"] == "gene"] 

## Useful utilities for GRanges objects
# Remove chromosome ranges
gff <- gff[values(gff)$type != "chromosome"]

# Erase the strand information
strand(gff) <- "*" 

# Collapses overlapping ranges to continuous ranges.
reduce(gff) 

# Return uncovered regions
gaps(gff) 

# More intuitive way to get uncovered regions
setdiff(as(seqinfo(gff), "GRanges"), gff) 

# Return disjoint ranges
disjoin(gff)
 
# Returns coverage of ranges
coverage(gff)
 
# Return the index pairings for overlapping ranges
findOverlaps(gff, gff[1:4])

# Counts overlapping ranges
countOverlaps(gff, gff[1:4])[1:40]

# Return only overlapping ranges
subsetByOverlaps(gff, gff[1:4])

## GRangesList Objects
# Stores every range in separate component of a GRangesList object
sp <- split(gff, seq(along=gff)) 

# Stores ranges of each chromosome in separate component
split(gff, seqnames(gff)) 

# Returns data as GRanges object
unlist(sp) 

# Subsetting of GRangesList objects is similar to GRanges objects
sp[1:4, "type"] 

# Looping over GRangesList objects similar to lists
lapply(sp[1:4], length) 

## end