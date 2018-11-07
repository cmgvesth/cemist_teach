
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

