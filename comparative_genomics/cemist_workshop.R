
# SAM format:
# @HD	VN:1.2	SO:coordinate	GO:reference
# @SQ	SN:jgi|Aspell1|394155|fgenesh1_kg.2_#_62_#_Locus3851v1rpkm35.34	LN:2613
# @SQ	SN:jgi|Aspell1|394154|fgenesh1_kg.2_#_61_#_Locus4504v1rpkm27.74	LN:1551
# @SQ	SN:jgi|Aspell1|394153|fgenesh1_kg.2_#_60_#_Locus3620v1rpkm38.19	LN:1248
# @SQ	SN:jgi|Aspell1|394152|fgenesh1_kg.2_#_59_#_Locus4600v2rpkm16.44	LN:828
# @SQ	SN:jgi|Aspell1|313244|e_gw1.2.316.1	LN:1254
# @SQ	SN:jgi|Aspell1|394150|fgenesh1_kg.2_#_57_#_Locus2657v1rpkm58.63	LN:1095
# @SQ	SN:jgi|Aspell1|342026|estExt_Genewise1.C_20127	LN:2100
# @PG	ID:0	VN:2.8.1+	PN:blastn
# XM_025540546.1	0	jgi|Aspell1|394155|fgenesh1_kg.2_#_62_#_Locus3851v1rpkm35.34	124	255	182H611M1I3M1I4M2D1493M1I3M1D1M1D5M2I3M2I249M1I4M2I64M3D43M334H	*	0	0	GGAGGGCTAC...GTCGTTTTGA	*	AS:i:1114	EV:f:0	NM:i:17	PI:f:82.20	BS:f:2058.29
# XM_025614481.1	0	jgi|Aspell1|394155|fgenesh1_kg.2_#_62_#_Locus3851v1rpkm35.34	124	255	188H164M1D1M2D158M2D2M1I2M1I902M1I2M1D160M2D1M1D1016M1D9M1I49M334H	*	0	0	GGAGAG...ACCATAGAAAGAC	*	AS:i:1018	EV:f:0	NM:i:14	PI:f:80.90	BS:f:1881.01


# Note the use of standard SAM/BAM format tags - looking at the first match:
#   AS:i:1808 - integer alignment score, here 1808
#   NM:i:0 - integer edit distance, here 0
# Plus non-standard tags:
#   EV:f:0 - float e-value, here 0
#   PI:f:96.55 - float percent identify, here 96.55
#   BS:f:701.049 - float bit-score, here 701.049
# https://blastedbio.blogspot.com/2015/07/ncbi-working-on-sam-output-from-blast.html
# https://en.wikipedia.org/wiki/SAM_(file_format)

# https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=FAQ#expect

# https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=BlastHelp

nn <- read.table(file = '~/Downloads/Y0B918VA015-Alignment.sam', comment.char = '@')
names(nn) <- c("QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL", 
               "int_aling_score", "evalue", "int_edit_dist", "percent_id", "bitscore"	)

nn$percent_id <- as.numeric( gsub("PI:f:", "", nn$percent_id) )

unique(nn$QNAME)
unique(nn$RNAME)
gg <- aggregate(formula = QNAME ~ RNAME, data = nn, FUN=length) # number of hits per query
gg <- aggregate(x = list(nn$percent_id), by = list(nn$RNAME) , FUN=mean ) # average pident per query

library(ggplot2)

ggplot(data = nn, aes(x = RNAME, y = percent_id)) + geom_boxplot()  + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  

pp <- read.table(file = '~/Downloads/Y0NKASM1015-Alignment (3).txt', comment.char = '#') # hit table (txt)
names(pp) <- c("query" , "subject", "percent_id", "align_len", "mismatches", "gap_opens", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "percent_positives")
gg <- aggregate(formula = subject ~ query, data = pp, FUN=length) # number of hits per query

length( unique(pp$query) )
length( unique(pp$subject) )

ggplot(data = pp, aes(x = reorder(query, percent_id, FUN = median), y = percent_id)) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


nn <- read.table(file = 'sam_Y0B918VA015-Alignment.sam', comment.char = '@')
names(nn) <- c("QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", 
               "PNEXT", "TLEN", "SEQ", "QUAL", 
               "int_aling_score", "evalue", "int_edit_dist", "percent_id", "bitscore"	)
nn$percent_id <- gsub("PI:f:", "", nn$percent_id)
nn$percent_id <- as.numeric( nn$percent_id)
unique(nn$QNAME)
unique(nn$RNAME)

one_2_one_count <- aggregate(formula = QNAME ~ RNAME, data = nn, FUN = length)
one_2_one_mean <- aggregate(formula = percent_id ~ QNAME, data = nn, FUN = mean)
one_2_many_count <- aggregate(formula = QNAME ~ RNAME + POS, data = nn, FUN = length)

tt <- aggregate(formula = RNAME + QNAME ~ , data = nn, FUN = length)


gg <- aggregate(x = list(nn$percent_id), by = list(nn$RNAME) , FUN = mean )



filepath <- "https://bitbucket.org/tcve/cmgbiotech/raw/0e919654c83089c60656a5a4a25b940cfbf8e9f4/data/Aspfo1_GeneCatalog_transcripts_20120615.nt.fasta"
myData <- read.table( file = '~/Dropbox/Drop-WORK/CeMiSt/Workshop/teach_cemist_github/comparative_genomics/Aspfo1_GeneCatalog_transcripts_20120615.nt.fasta.txt', header=TRUE, sep='\t')

source('~/Desktop/cmg_read_fasta.R')

myData <- cmg_read_fasta(filename = filepath, type = 'nuc', orgname = 'test')

boxplot(x = myData$length)

orgA <- cmg_read_fasta(filename = 'Aspalli1_GeneCatalog_proteins_20161018.aa.fasta', type = 'aa', orgname = 'orgA')
orgB <- cmg_read_fasta(filename = 'Aspamb1_GeneCatalog_proteins_20160910.aa.fasta', type = 'aa', orgname = 'orgB')
orgC <- cmg_read_fasta(filename = 'Aspall1_GeneCatalog_proteins_20160910.aa.fasta', type = 'aa', orgname = 'orgC')
orgD <- cmg_read_fasta(filename = 'Aspalbe1_GeneCatalog_proteins_20160720.aa.fasta', type = 'aa', orgname = 'orgD')
allFSA <- rbind(orgA, orgB, orgC, orgD)

library(ggplot2)

ggplot(data = allFSA, aes(x = orgname, y = length, fill = 'blue')) + geom_boxplot()
ggplot(data = allFSA, aes(x = orgname, y = length, fill = 'red')) + geom_boxplot()
ggplot(data = allFSA, aes(x = orgname, y = length)) + geom_boxplot(fill = 'red')
ggplot(data = allFSA, aes(x = orgname, y = length, fill = orgname)) + geom_boxplot()


ggplot(data=allFSA, aes(x = length)) + geom_histogram()

ggplot(data=allFSA, aes(x = length, fill = orgname)) + geom_histogram(position = "dodge")
ggplot(data=allFSA, aes(x = length, fill = orgname)) + geom_histogram() + facet_grid( .~ orgname)
p1 <- ggplot(data=allFSA, aes(x = length, fill = orgname)) + geom_histogram() + facet_grid(orgname ~.)

ggsave(plot = p1, filename = 'proteinLength.pdf', device = 'pdf')

ggplot(data=allFSA, aes(y = length, fill = orgname)) + geom_boxplot() + facet_grid( .~ orgname)
ggplot(data=allFSA, aes(x = length, fill = orgname)) + geom_bar()+ facet_grid( orgname~. )

ggplot(data=allFSA, aes(x = length, fill = orgname)) + geom_density()+ facet_grid( orgname~. )
ggplot(data=allFSA, aes( x = '' , y = length, fill = orgname)) + geom_violin()+ facet_grid( .~orgname )



library(gplots)

organisms <- c('orgA', 'orgB', 'orgC', 'orgD')
measure1 <- c(1,2,1,2)
measure2 <- c(4,3,5,6)
measure3 <- c(2,1,1,1)
measure4 <- c(6,5,6,5)
measureData <- data.frame(organisms, measure1, measure2, measure3, measure4)
matrixData <- measureData[,2:ncol(measureData)]
matrix <- data.matrix(matrixData)
rownames(matrix) <- measureData$organisms

heatmap(matrix)
heatmap(matrix, margins = c(10,10))
heatmap(matrix, margins = c(10,10), col=(cm.colors(5)), xlab='Measures',
        ylab='Organisms', main = 'Example heatmap')


heatmap.2(matrix)
heatmap.2(matrix, margins = c(10,10))

heatmap.2(matrix, margins = c(10,10), col=(cm.colors(5)), xlab='Measures',
          ylab='Organisms', main = 'Example heatmap')

heatmap.2(matrix, margins = c(10,10), col=c('red', 'blue', 'green', 'orange') ,
          Colv = FALSE, dendrogram='row', xlab='Measures', ylab='Organisms',
          main = 'Example heatmap', trace="none")

pcount <- aggregate(allFSA$length, by = list(allFSA$orgname), FUN=summary)
# Agg <- aggregate(df$Result, list(df$Location), summary)

summary(allFSA$length)

pdf(file = 'test.pdf')
matrixData <- pcount[,2:ncol(pcount)]
matrix <- data.matrix(matrixData)
heatmap.2(matrix)
dev.off()

setwd("~/Dropbox/Drop-WORK/CeMiSt/Workshop/teach_cemist_github/comparative_genomics")

fastafiles <- list.files(pattern = 'faa')
fastaEntries <- data.frame(orgname = NA, sequenceid = NA) # cretae data.frame with blank observation

for (file in fastafiles) {
  lines <- read.csv(file = paste0( file), header = FALSE) # get all lines
  headers <- as.data.frame( lines[  grepl(pattern = '>', lines$V1) ,  ]) # select lines that are headers
  names(headers) <- "header"
  
  headers$orgname <- gsub(pattern = '>.+\\s\\[(.+)\\]', replacement = '\\1', x = headers$header) # problem with header format
  headers$orgname <- unique( headers$orgname[ grepl(pattern = 'Aspergillus', x = headers$orgname)] ) # fix problem
  #unique_names <- unique( headers[ !grepl(pattern = '>', x = headers$orgname), ]$orgname ) # fix problem
  #headers$orgname <- unique_names[ grepl(pattern = 'Aspergillus', x = unique_names)] # fix problem
  
  headers$sequenceid <- gsub(pattern = '>(\\S+)\\s.+', replacement = '\\1', x = headers$header) 

  fastaEntries <- rbind(fastaEntries, headers[, c('sequenceid','orgname' )] )
}
fastaEntries <- fastaEntries[!is.na(fastaEntries$orgname), ] # remove the blank observation
unique(fastaEntries$orgname)
fastaEntries_count <- aggregate(formula = sequenceid~orgname, data = fastaEntries, FUN = length )

blastfiles <- list.files(pattern = 'csv')
blastEntries <- data.frame(qorg = NA, qacc = NA, sorg = NA, sacc = NA, percent_id = NA, alignment_length = NA) # cretae data.frame with blank observation

for (file in blastfiles) {
  print(file)
  blastfile <- read.csv(file = file, header = FALSE) 
  names(blastfile) <- c("query_accver" , "subject_accver" , "percent_identity" , "alignment_length" , 
                        "mismatches" , "gap_opens" , "q_start" , "q_end" , "s_start" , "s_end" , 
                        "evalue" , "bit_score" , "percent_positives")
  
  blastfile$subjectOrg <- unique( fastaEntries[fastaEntries$sequenceid %in% blastfile$subject_accver,]$orgname )   
  blastfile$queryOrg <- unique( fastaEntries[fastaEntries$sequenceid %in% blastfile$query_accver,]$orgname )   
  blasthits <- blastfile[, c('queryOrg','query_accver', 'subjectOrg','subject_accver', 'percent_identity', 'alignment_length' ) ]
  names(blasthits) <- c("qorg","qacc","sorg","sacc","percent_id","alignment_length")
  
  blastEntries <- rbind(blastEntries, blasthits)
}

blastEntries <- blastEntries[!is.na(blastEntries$qorg), ] # remove the blank observation
unique(blastEntries$qorg)
unique(blastEntries$sorg)

unique(paste(blastEntries$sorg, blastEntries$qorg))

tt_long <- aggregate(formula = sacc~qorg+sorg, data = blastEntries, FUN = length)
tt_long <- aggregate(formula = qacc~qorg+sorg, data = blastEntries, FUN = length)
tt_long <- aggregate(formula = qacc~qorg+sorg, data = blastEntries, FUN = function(x) length(unique(x)))
tt_long <- merge(x = tt_long, y = fastaEntries_count, by.x = 'qorg', by.y = 'orgname') 
tt_long$percent_qcc <- tt_long$qacc/tt_long$sequenceid*100

library(reshape)
tt_wide <- cast(tt_long, qorg ~ sorg, value = c('qacc') )
tt_wide <- cast(tt_long, qorg ~ sorg, value = c('percent_qcc') )
tt_matrix <- as.matrix( tt_wide[, 2:length(names(tt_wide))] )
colnames(tt_matrix) <- names(tt_wide[, 2:length(names(tt_wide))])
rownames(tt_matrix) <- tt_wide$qorg

heatmap.2(tt_matrix, margins = c(10,10), trace = NULL, col = c('orange', 'red', 'blue'), na.color = 'black')

          