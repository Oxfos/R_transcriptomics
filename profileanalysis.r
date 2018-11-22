# This script is for:
# Creating a list of GO terms for Nasonia vitripennis.
# The used scrip part begins at line 587.


getwd()
setwd("C:/Users/p233711/Dropbox/Others/Gerard_NGS/R")
setwd("C:/Users/Rinaldo C. Bertossa/Dropbox/Others/Gerard_NGS/R")
ls()
getwd()

# Analyzing TestSet_Genes 17-2-2015
set<-read.table('TestSet_Genes.txt', header=TRUE, skip=9, sep="\t")

# inspection of data
head(set)
par(mfrow=c(2,2))
hist(set$Average_Counts, breaks=10000, xlim=c(0,10))
hist(set$RPKM, breaks=10000, xlim=c(0,25))
hist(set$Read_Counts, breaks=10000, xlim=c(0,1000))

# making records subsets
# subset without names
set0<-subset(set, set$Gene=="" & set$RNA.Accession=="" & set$Protein.Accession=="")
# subset of set0 with counts
set01<-subset(set0,set0$Average.Counts!=0)
# subset with names
set1<-subset(set,set$Gene!=""|set$RNA.Accession!=""|set$Protein.Accession!="")

# Exploring the sets
set.l<-nrow(set)
set0.l<-nrow(set0) # the set has 18009 records (14275 records remain)
set01.l<-nrow(set01)
set1.l<-nrow(set1)
set0.l/set.l
set1.l/set.l
set01.l/set.l
set01.l/set0.l


# Which are the biggest values for Average_Counts and RPKM?
max(set0$Average_Counts) # biggest value for average counts is 319.83
median(set0$Average_Counts) # median 0
mean(set0$Average_Counts)  # mean 0.6828491
max(set0$RPKM)           # biggest value for RPKM is 843.2499
median(set0$RPKM)     # median 0
mean(set0$RPKM)       # mean 1.503558
median(set$RPKM)
mean(set$RPKM)
median(set1$RPKM)
mean(set1$RPKM)
median(set01$RPKM)
mean(set01$RPKM)
mean(set$Reference_Length)
mean(set0$Reference_Length)
mean(set01$Reference_Length)
mean(set1$Reference_Length)
median(set$Reference_Length)
median(set0$Reference_Length)
median(set01$Reference_Length)
median(set1$Reference_Length)
max(set1$Reference_Length)
min(set1$Reference_Length)
max(set0$Reference_Length)
min(set0$Reference_Length)
max(set01$Reference_Length)
min(set01$Reference_Length)
dev.off() # resets plot layout
plot(set$Reference_Length, set$Average_Count,xlim=c(0,600000), ylim=c(0,20),cex=0.7)
par(mfrow=c(2,2))
plot(set$Reference_Length, set$RPKM,xlim=c(0,600000), ylim=c(0,20), cex=0.7)
plot(set0$Reference_Length, set0$RPKM,xlim=c(0,600000), ylim=c(0,20), cex=0.7)
plot(set01$Reference_Length, set01$RPKM,xlim=c(0,600000), ylim=c(0,20), cex=0.7)
plot(set1$Reference_Length, set1$RPKM,xlim=c(0,600000), ylim=c(0,20), cex=0.7)
boxplot(set$Reference_Length)
boxplot(set0$Reference_Length)
boxplot(set01$Reference_Length)
boxplot(set0$RPKM)
boxplot(set$RPKM)
boxplot(set)

# statistics on different 'set'
nrow(set)

# subsetting data (only protein accession)
set.prot<-as.data.frame(set$Protein_Accession)
set.pro<-as.data.frame(substr(set.prot, 0, nchar(set.prot)-1))

# Analyzing continuousmRNA 13-2-2015
set<-read.table('TestSet_Genes.txt', header=TRUE, sep="\t")
# subsetting data (only protein accession)
set.prot<-as.data.frame(set$Protein_Accession)
set.pro<-as.data.frame(substr(set.prot, 0, nchar(set.prot)-1))
# inspection of data with histograms
par(mfrow=c(2,2))
hist(set$Average_Counts, breaks=10000, xlim=c(0,50))
hist(set$RPKM, breaks=10000, xlim=c(0,100))
hist(set$Read_Counts, breaks=10000, xlim=c(0,2000))

##############################################################
## Analyzing UniProt and QuickGO Nasonia vitripennis lists  ##
## and importing the unique GO terms into the 'genes' list  ##
##############################################################

# From line 587: I did the whole analysis with a different approach.
# This new approach I now consider to be the official one.


# General outline of procedure:
# There are files containing the GO terms (UniProt and QuickGO) and files containing the list of genes (genes)
# To which the GO terms have to be assigned
#
# Where important information is present:
# UniProt:  GO terms, Gene names, UniProt ID, RefSeq Prot
# Quick GO: GO terms, Gene names, UniProt ID
# genes:    Gene names, RefSeq Prot
# 
# 1. Using the same UniProt ID, QuickGO GO terms are pooled and imported into UniProt
# 2. QuickGO and UniProt GO terms are pooled
# 3. Using the same Gene name, QuickGO GO terms are pooled and imported into UniProt
# 4. QuickGO and UniProt GO are pooled
# 5. Using Gene names, GO terms from UniProt are imported (and pooled) into genes
# 6. Using RefSeq Prot, GO terms from UniProt are imported (and pooled) into genes

#######################################################################################
#                                                                                     #
# 1. GO terms from QuickGO are imported into UniProt using Gene names and UniProt ID  #
#                                                                                     #
#######################################################################################


# Resetting workspace
rm(list=ls())

getwd() # just to check

# Importing the UniProt and QuickGO files
# working at RUG:
Nvu<-read.table('C:/Users/p233711/Dropbox/Others/Gerard_NGS/Data/QuickGO/UniProt_Nvit_20150218.txt', header=TRUE, quote='', stringsAsFactors=FALSE, sep="\t")
Nvq<-read.table('C:/Users/p233711/Dropbox/Others/Gerard_NGS/Data/QuickGO/QuickGO_Nvit_20150215.txt', header=TRUE, quote='', stringsAsFactors=FALSE, sep="\t")
# working at home:
Nvu<-read.table('C:/Users/Rinaldo C. Bertossa/Dropbox/Others/Gerard_NGS/Data/QuickGO/UniProt_Nvit_20150218.txt', header=TRUE, quote='', stringsAsFactors=FALSE, sep="\t")
Nvq<-read.table('C:/Users/Rinaldo C. Bertossa/Dropbox/Others/Gerard_NGS/Data/QuickGO/QuickGO_Nvit_20150215.txt', header=TRUE, quote='', stringsAsFactors=FALSE, sep="\t")

# Importing GO terms from Nvq to Nvu with the same UniProt ID

a<-data.frame(matrix(NA,nrow=nrow(Nvu),ncol=4))
colnames(a)=c("NvqGOs","NvuGOOld","NvuGONew","TEST")
Nvu$Gene_ontology_IDs <- gsub("; ",", ",Nvu$Gene_ontology_IDs) # preparing the df for analysis...
a$NvuGOOld<-Nvu$Gene_ontology_IDs

for (i in 1:nrow(a)) {
  if (Nvu[i,1]=="") {
    next
  }
  else {
    r <- which(Nvq$ID==Nvu$Entry[i])
    r.l <- length(r)
    if (r.l==0) {
      next
    }
    else {
      a$NvqGOs[i] <- toString(unique(Nvq$GO.ID[r]))
    }
  }
}

# 2. Now we combine the newly created GOs list with the old GO in Nvu
# and test whether old GO are the same as new GO

for (i in 1:nrow(Nvu)) {
  a$NvuGONew[i] <- toString(unique(na.omit(unlist(strsplit(c(Nvu$Gene_ontology_IDs[i], a$NvqGOs[i]), ", ")))))
  a$TEST[i] <- a$NvuGONew[i]==a$NvuGOOld[i]
}

b<-subset(a,a$TEST==FALSE)

# There are 1280 records which now have new GO lists in respect of the Nvu GO they had
# We substitute them into Nvu:

Nvu$Gene_ontology_IDs <- a$NvuGONew

# Until now we have imported Nvq GO terms into Nvu according to the same UniProt IDs
# Preparing the lists:

# select only Nasonia vitripennis:
Nvu1<-Nvu[grep("wasp",Nvu$Organism),]

# Checking whether Nvu1 records with the same gene name have the same GO terms:

Nvu1.1 <- data.frame(matrix(NA,nrow=nrow(Nvu1),ncol=4))
colnames(Nvu1.1)<- c("Nvu1$5_GeneNames","DoubleNvuGeneNamesRows","Nvu1$1_Entry","TEST")

Nvu1.1[c(1,3)] <- Nvu1[c(5,1)]

for (i in 1:nrow(Nvu1)) {
  if (Nvu1[i,5]=="") {
    next
  }
  else {
    r <- which (Nvu1[5]==Nvu1[i,5])
    r.l <- length(r)
    if (r.l<2) {
      next
    }
    else {
      Nvu1.1[i,2] <- toString(r)
      t<-NA
      for (ii in 1:(r.l)-1) {
        t[ii] <- Nvu1[r[ii],12]==Nvu1[r[ii+1],12]
        Nvu1.1[i,4] <- toString(t)
      }
    }
  }
}

length(grep("FALSE",Nvu1.1[,4]))
# 78
# > Not all gene names have the same GO terms because they have different UniProt identifiers.
# For instance the fru gene may have different isoforms (and so different uniprot id) but be called always fru
# The same gene names may have different GO terms because they have different UniProt ID.
# I decide here to use only gene names to assign GO terms because this is the division we are using at the end:
# Assigning GO terms to genes, and not to transcripts or 'units less than genes'.
# I will therefore later pool GO terms for the same gene name.
# The same applies actually for records in Nvq (QuickGO)

#########################################################################################
# 3. Loop to import Nvq (QuickGO) GO terms WITH THE SAME GENE NAME as in Nvu1 into Nvu2   #
# and keeping relevant information in one single file                                   #
# (QuickGO GO terms with the same gene name are pooled):                                #
#########################################################################################

# Dataframe for Nvu1 results:
Nvu2<-data.frame(matrix(NA, nrow = nrow(Nvu1), ncol = 7)) 
colnames(Nvu2)<-c("Query:Nvu1Symbol","Nr.IndNvqHits","UniqueNvqHitSymbol","NvqHitRows","NvqGOs","Nvu1GOs","IndHits")

Nvu2[1] <- Nvu1[5] # The symbol(s) searched are recorded into the first column of Nvu2
Nvu2[6] <- Nvu1[12] # Nvu1 GOs are copied into Nvu2, col 6.
Nvu2[7] <- 0

for (i in 1:nrow(Nvu1)) {      # Loop to check every line of Nvq
  if (Nvu1[i,5] == "") {
    next
  }
  else {
    Nvu1i5<-unlist(strsplit(Nvu1[i,5]," "))
    Nvu1i5.l<-length(Nvu1i5)
    for (ii in 1:Nvu1i5.l) { # for every element in Nvu1[i,5]...
      r <- which(Nvq$Symbol==Nvu1i5[ii]) # r are the row numbers of the hits
      l <- length(r) # l stores how many hits are found
      if (l==0) {   # when no hit is found...
        next
      }
      else {
        v <- unique(Nvq$Symbol[r])
        GOn <- unique(Nvq$GO.ID[r]) # pooled(!!) GO terms
        GOp <- unlist(strsplit(Nvu2[i,5],", "))
        Nvu2[i,2] <- toString(na.omit(c(Nvu2[i,2],l))) # number of independent Nvq hits is recorded
        Nvu2[i,3] <- toString(na.omit(c(Nvu2[i,3],v)))  # Nvq hit values added
        Nvu2[i,4] <- toString(na.omit(c(Nvu2[i,4],r))) # Nvq hit rows numbers
        Nvu2[i,5] <- toString(unique(na.omit(c(GOp,GOn))))  # Pooled Nvq GO values added
        Nvu2[i,7] <- Nvu2[i,7]+1
      }
    }
  }
}

which(Nvu2$IndHits>1)
# integer(0) -> no records seem to be hit more than once... which is to be expected,
# because one question (slot) is filled with only one answer.

# I did this before. Now is a bit superflous:
# Check that GO are the same also for duplicates:
Nvu2.1<-data.frame(matrix(NA,nrow=nrow(Nvu2),ncol=8))
colnames(Nvu2.1)<-c("Row#","Genes","NvqGObefore","NvqGOafter","TestNvqGO","Nvu1GObefore","Nvu1GOafter","TestNvu1GO")

Nvu2[,6]<-gsub("; ",", ",Nvu2[,6],perl=TRUE)

Nvu2.1[c(3,4)]<-Nvu2[5] # Nvq GO terms copied
Nvu2.1[c(6,7)]<-Nvu2[6] # Nvu1 GO terms copied

for (c in 1:nrow(Nvu2)) {
  if (Nvu2[c,1]=="") {
    next
  }
  else {
    a<-which(Nvu2[1]==Nvu2[c,1])
    a.l<-length(a)
    if (a.l>1) {
      Nvu2.1[c,1]<-toString(a) # hits row numbers
      Nvu2.1[c,2]<-Nvu2[c,1] # query
      Nvu2.1[c,4]<-toString(unique(unlist(strsplit(Nvu2[a,5], ", ")))) # corresponding Nvq GOs of hits
      Nvu2.1[c,5]<-Nvu2.1[c,3]==Nvu2.1[c,4] # tests whether original Nvu2's Nvq GO list remains the same
      Nvu2.1[c,7]<-toString(unique(unlist(strsplit(Nvu2[a,6], ", ")))) # corresponding Nvu1 GOs of hits
      Nvu2.1[c,8]<-Nvu2.1[c,6]==Nvu2.1[c,7] # tests whether original Nvu2's Nvu GO (unchanged) remain the same
    }
    else {
      next
    }
  }
}

which(Nvu2.1[5]==FALSE)
# integer(0) -> this was expected because Nvq GO's correponding to the same Nvu1 symbols were pooled while importing them
which(Nvu2.1[8]==FALSE)
# Apparently there are instances only for Nvu1 GO terms where the same gene names do not have the same GOs...
# This confirms what found above.


## Check on a previous script and now deleted!: Printing out Nvq Symbols not found in Nvu1
which(Nvq2[,3]==0)   # All Nvq Gene Symbols are found in Nvu1!!! :D
which(Nvq2[,3]==3)   # Check ;)

# with new script above:
# Nvu2 rows with NA hits = 7495
# Nvu2 rows without NA hits = 1687
# Total = 9182
# Total unique gene symbols in Nvq = 1647 > therefore the magnitude is similar to Nvu2 rows with hits.
# In fact, the unique Nvu2 hits found are 1648 (i.e. 1647+NA). We conclude that all Nvq gene symbols are
# found in Nvu1, as found with the previous script.


# Combining Nvu1 (original UniProt) and Nvu2 (UniProt with new findings) in Nvu4; keeping only relevant columns
Nvu3<-cbind(Nvu1[,c(1,4,5,7,12,13)],Nvu2[,c(3,5)])

# Exporting Nvu4 for Gerard
write.table(Nvu3, file="Nvit_GOs_v2.txt", quote=FALSE, sep='\t')
# Check: whether the file is produced correctly:
Nvu3<-read.table('C:/Users/Rinaldo C. Bertossa/Dropbox/Others/Gerard_NGS/R/Nvit_GOs_v2.txt', header=TRUE, quote='', stringsAsFactors=FALSE, sep="\t")
# The file is produced in the correct way

###############################################################
# 4. Check for duplicates GO terms between original Nvu and Nvq; #
# pooling GO terms into one cell and keeping only uniques     #
###############################################################

# have first a look at the class and mode of the data:
sapply(Nvu4,class)
sapply(Nvu4,mode)


# 1. Create a new df where to put go terms
# 2. Copy GO terms from other column if missing
# 3. Chop GO terms within cells into distinct strings
# 4. Combine GO terms into single vector
# 5. Eliminate duplicates and copy to new df

length(which(is.na(Nvu3[,5]) & is.na(Nvu3[,8]))) # There are no rows in which there is NA contemporarily in Nvu4[,5] and Nvu4[,8]
length(which(is.na(Nvu3[,5]) & !is.na(Nvu3[,8]))) # If there is something in Nvu4[,8] then also in Nvu4[,5]
length(which(!is.na(Nvu3[,5]) & is.na(Nvu3[,8]))) # Not always there is something in Nvu4[,8] when in Nvu4[,5]

Nvu4<-data.frame(UniqueGOs=matrix(NA, nrow = nrow(Nvu3), ncol = 1)) # 1. Create a new df with to put go terms

for (i in 1:nrow(Nvu3)) {
  if (is.na(Nvu3[i,8])){  # 2. Copy GO terms from column 5 if missing in column 8
    Nvu4[i,1]<-unique(Nvu3[i,5])
  }
  else if (Nvu3[i,5]=="") {  # 2. Copy GO terms from column 8 if missing in column 5
    Nvu4[i,1]<-unique(Nvu3[i,8]) # this never happens though... see above
  }
  else {
    Nvu3.8<-unlist(strsplit(Nvu3[i,8], ", "))  # 3. Chop GO terms within cells into distinct strings
    Nvu3.5<-unlist(strsplit(Nvu3[i,5], "; "))  # 3. Chop GO terms within cells into distinct strings
    Nvu4[i,1]<-toString(unique(c(Nvu3.5,Nvu3.8)))  # 4. Combine GO terms into single vector + 5. Eliminate duplicates and copy to df
  }
}

# making a unique file and keeping relevant columns
Nvu5<-cbind(Nvu3,Nvu4)
Nvu5<-Nvu5[,c(1,2,4,3,6,9)]

##### START: THIS IS NO MORE NEEDED, ##########################################################
#                                                                                             #
# BECAUSE LATER YOU WILL ANYHOW POOL GO TERMS COMING FROM DIFFERENT Nvu1 Gene records, but    #
# having the same 'gene' name and which had different GO terms (because they had different UniProt ID)

# Now we will pool Nvu GO's of the same gene name together, then
# we test whether records with the same gene name have the same GO's in the old and new version.
# The new version, with pooled GO terms, is copied in the relevant column of the file.

Nvu2.2<-data.frame(matrix(NA,nrow=nrow(Nvu5),ncol=5))
colnames(Nvu2.2)<-c("Row#","Genes","Nvu5GObefore","Nvu5GOafter","TestNvu5GO")
Nvu5[,6]<-gsub("; ",", ",Nvu5[,6],perl=TRUE)

Nvu2.2[c(3,4)]<-Nvu5[6]
for (c in 1:nrow(Nvu5)) {
  if (Nvu5[c,4]=="") {
    next
  }
  else {
    a<-which(Nvu5[4]==Nvu5[c,4])
    a.l<-length(a)
    for (cc in 1:a.l) {
      Nvu2.2[c,1]<-toString(a) # hits row numbers
      Nvu2.2[c,2]<-Nvu5[c,4] # query
      Nvu2.2[c,4]<-toString(unique(unlist(strsplit(Nvu5[a,6], ", ")))) # corresponding Nvu5 GOs of hits
      Nvu2.2[c,5]<-Nvu2.2[c,3]==Nvu2.2[c,4] # tests whether original Nvu5's Nvu1 GO list remains the same
    }
  }
}

a<-which(Nvu2.2[5]==FALSE)
a # Therefore, no different GO terms for the same gene name...
Nvu5[a,6]<-Nvu2.2[a,4] # The pooled GO terms is copied in the old file.
# If you rerun the test, a will be = 0.
#                                                                                            #
#                                                                                            #
##### END: THIS IS NO MORE NEEDED ############################################################


#######################################################################################
#                                                                                     #
# GO terms from UniProt are imported into genes using Gene names and RefSeq Prot   #
#                                                                                     #
#######################################################################################


#############################################################
# Assigning GO terms from Nvu5 to the genes in 'genes' file #
#############################################################

# NOTE: the problem of different GO terms for records having the same gene name falls
# because GO terms from these records will be pooled.
# The 'GO Test' more below will confirm this.

# Import 'genes' file:
getwd()
# Working at work:
genes<-read.table('C:/Users/p233711/Dropbox/Others/Gerard_NGS/Data/NasoniaExpressionReports/Genes/L6_AAAAAT_converted_Output_v20v11_Nvit21_Expression_Report_Genes.txt', skip=9, header=TRUE, quote='', stringsAsFactors=FALSE, sep="\t")
# Working at home:
genes<-read.table('C:/Users/Rinaldo C. Bertossa/Dropbox/Others/Gerard_NGS/Data/NasoniaExpressionReports/Genes/L6_AAAAAT_converted_Output_v20v11_Nvit21_Expression_Report_Genes.txt', skip=9, header=TRUE, quote='', stringsAsFactors=FALSE, sep="\t")

f<-function (x) sub("[.*;].*", "", x, perl=TRUE) # function to 'clean' genes of everything after ';'
genes<-cbind(genes[,1:5],apply(genes[6],2,f),genes[,7:20],stringsAsFactors=FALSE) # genes names prepared for searching

# script to find Nvu5 genes in df 'genes' and assign GO terms:

genes1<-data.frame(matrix(NA, nrow = nrow(genes), ncol = 12)) # Dataframe for Nvu7 results
colnames(genes1)<-c("GeneRowsHit","GeneValuesHit","ValueQuery1","RefSeqRowHit","RefSeqValuesHit","ValueQuery2","Gene'sGOs","RefSeqGOs","#GeneHits","GenesCounter","#RefSeqHits","RefSeqCounter")

genes1[c(1,9,10,11,12)]<-0

nogenes<-0 # initializing counter for 'no genes found'
norefseq<-0 # initializing counter for 'no refseq found'

# 5. First we search dfr 'genes' according to Nvu6's column 4: Genes
for (i in 1:nrow(Nvu5)) {  # for every row of Nvu5...
  if (Nvu5[i,4]=="") {  # if Nvu6[i,4] is empty, go to next cell
    next
  }
  else {
    Nvu5.4<-unlist(strsplit(Nvu5[i,4]," "))  # split the terms in column 4 (Gene_names)
    Nvu5.4l<-length(Nvu5.4)
    for (c in 1:Nvu5.4l) { # For how many terms there are in Nvu6[i,4]...
      genes.r<-which(genes[6]==Nvu5.4[c]) # Look for each of them into column 6 (Gene) of df 'genes'
      genes.l<-length(genes.r)
      if (genes.l==0) { # if nothing is found, go to next value
        nogenes=nogenes+1
      }
      else {
        genes.v<-NA
        genes.v<-genes[genes.r,6]
        for (d in 1:genes.l) {
          genes1[genes.r[d],1]<-paste(genes.r, collapse=", ") # record row numbers of the hit(s) in the 'genes' df
          genes1[genes.r[d],2]<-toString(genes.v) # record value ('gene') of the hit(s)
          genes1[genes.r[d],3]<-Nvu5[i,4] # record query
          genes1[genes.r[d],7]<-toString(unique(na.omit(unlist(strsplit(c(genes1[genes.r[d],7],Nvu5[i,6]), ", "))))) # record corresponding GO terms !!! It overrides even with empty content!!!
          genes1[genes.r[d],9]<-genes.l # record number of hits
          genes1[genes.r[d],10]<-genes1[genes.r[d],10]+1 # record if hit more than once independently
        }
      }  
    }
  }  
}
# 6. Then we search dfr 'genes' according to Nvu6's column 5: RefSeq terms
for (i in 1:nrow(Nvu5)) {  # for every row of Nvu6...
  if (Nvu5[i,5]=="") {
    next
  }
  else {
    Nvu5.5<-unlist(strsplit(Nvu5[i,5],";"))  # split the terms in column 5 (RefSeq)
    q<-sub("[.*].*", "", Nvu5.5, perl=TRUE) # removes '.' and anything after it
    Nvu5.5<-q
    Nvu5.5l<-length(Nvu5.5)
    for (c in 1:Nvu5.5l) { # For how many terms there are...
      genes.r<-grep(Nvu5.5[c], genes[,10], fixed=TRUE) # Look for each of them into column 6 (Gene) of df 'genes'
      genes.r.l<-length(genes.r)
      if (genes.r.l==0) {
        norefseq<-norefseq+1
      }
      else {
        genes.v<-grep(Nvu5.5[c], genes[,10], fixed=TRUE, value=TRUE)
        for (d in 1:genes.r.l) {
          genes1[genes.r[d],4]<-toString(genes.r) # row hit recorded
          genes1[genes.r[d],5]<-toString(genes.v) # value hit recorded
          genes1[genes.r[d],6]<-toString(q) # value query
          genes1[genes.r[d],8]<-Nvu5[i,6] # GOs: same problem as above if GOs different
          genes1[genes.r[d],11]<-genes.r.l # record number of hits
          genes1[genes.r[d],12]<-genes1[genes.r[d],12]+1 # record if hit more than once
        }
      }  
    }
  } 
}

nogenes
norefseq

# We check now queries with multiple hits and hits hitted more than once with distinct queries:
# First we check for 'RefSeq' then for 'genes':

# RefSeq:

# Fortunately there are no multiple hits for RefSeq queries:
which(genes1[11]>1)
# integer(0)

# Fortunately there are no hits hit more than once independently (with distinct queries) for 'RefSeq':
which(genes1[12]>1)
#integer(0)

# Genes:

# Strangely there are multiple hits for 'genes' queries:
a<-which(genes1[9]>1)
a
#[1]  5441  5524  5734  5847  5882  6052  6160  6176  7678 15487 15602 18989 19048 19231 25036 25077 25115 25208 25290 25446 26590 26898 27597
#[24] 27702 28397 29339 30524
# Which implies that there are duplicates also in the 'genes' list.
# Now we have to be sure that they all got the same GO terms.

# 'GO Test': Checking that multiple hits have the same GOs:

morehits<-which(genes1[9]>1)
morehits.l<-length(morehits)
more<-data.frame(matrix(NA, nrow = morehits.l, ncol = 6))
colnames(more)<-c("RowsMorehits","#TotalHits","genes1$ValueQuery1","Row#ofNvu5Query","QueryInNvu5$Gene_names","GOTests")

more[1]<-morehits # row numbers of records hit more than once are recorded
more[2]<-genes1[(morehits),9] # the number of ultiple hits are recorded
more[3]<-genes1[(morehits),3] # original Nvu5$4 query to find them
more<-more[order(more[2],more[3]),] # records are ordered according to # of hits

for (c in 1:morehits.l) {
  more.t<-which(genes1[3]==more[c,3]) # all genes1 records with the same original Nvu5$4 query1
  more[c,4]<-toString(more.t) # rows of genes1 hits with the same original Nvu5$4 query1
  more[c,5]<-toString(genes1[more.t,3]) # values of the original query1 value in Nvu5$4: wrong because row names used are of genes1!!
  more.t.l<-length(more.t)
  t<-NA
  for (cc in 1:((more.t.l)-1)) {
    t[cc]<-genes1[more.t[cc],7]==genes1[more.t[cc+1],7] # we test whether the GO terms of genes1 records found with the same query1 are the same
    more[c,6]<-toString(t) # the result (logical) of the test is stored in col 5 of df more
  }
}

# All multiple hits have the same GO terms.

# All GO terms for the same gene record are pooled:

for (i in 1:nrow(genes1)) {
  GO <- c(genes1[i,7],genes1[i,8])
  genes1[i,12] <- toString(unique(na.omit(unlist(strsplit(GO, ", "))))) # we put the unique GO into column 12 which is not used anymore
}

# Now we produce a unique file 'genes2' with relevant info:

genes2 <- cbind(genes[c(4,5,6,8,9,10,13:20)],genes1[12])

colnames(genes2)[15]<-c("UniqueGOs")

# We check that duplicates do have the same GO terms:

e<-NA
for (i in 1:nrow(genes2)) {
  if (genes2[i,3]=="") {
    next
  }
  else {
    e[i] <- length(which(genes2[3]==genes2[i,3]))
  }
}

a<-which(e>1)
b<-genes2[a,]
b0<-subset(b,b[15]=="")
b1<-subset(b,b[15]!="")

b0<-b0[order(b0[3]),c(3,15)]
b1<-b1[order(b1[3]),c(3,15)]

# > Duplicates do have the same GO terms

# Everything once again:
#
# 1. Pool GO terms in Nvq (=QuickGO GO terms) according to UniProt ID and add this column to Nvq: NvqGO.UP
# 2. Pool GO terms in Nvq according to Symbol and add this column to Nvq: NvqGO.Symbol
# 3. Pool Nvq GO terms
# 4. Import pooled GO terms into Nvu (=UniProt GO terms) (into new column) according to common UniProt
# 5. Import pooled GO terms into Nvu (into new column) according to common Gene names
# 6. Pool all Nvu GO terms: 1st time
# 7. Pool Nvu GO terms in Nvu according to UniProt ID and add this column to Nvu: no need: Nvu records have unique UniProt
# 8. Pool Nvu GO terms in Nvu according to Gene_name and add this column to Nvu: NvuGO.Gene
# 9. Pool all Nvu GO terms: 2nd time
# 10. Complete 'genes' data with 'Nvit_2.1_NCBI_ProteinTable449_202468_20150331.txt' (>called 'genes1')
# 11. Import Nvu GOs into genes1 (into new column) according to common Gene name
# 12. Import Nvu GOs into genes1 (into new column) according to common RefSeq
# 13. Pool all genes1 GOs columns
# 14. Pooling GO terms in genes1 according to Gene
# 15. Pooling GO terms in genes1 according to GeneID
# 16. Pooling GO terms in genes1 according to RefSeq (Protein.Accession)
# 17. We give a new name to the df: genes2
# 18. Preparing a file to send to Gerard
# 19. List of genes with associated GO terms

# We import the files first:

# Resetting workspace
rm(list=ls())

getwd() # just to check

# Importing the UniProt and QuickGO files
# working at RUG:
Nvu<-read.table('C:/Users/p233711/Dropbox/Others/Gerard_NGS/Data/GO/UniProt_Nvit_20150218.txt', header=TRUE, quote='', stringsAsFactors=FALSE, sep="\t")
Nvq<-read.table('C:/Users/p233711/Dropbox/Others/Gerard_NGS/Data/GO/QuickGO_Nvit_20150215.txt', header=TRUE, quote='', stringsAsFactors=FALSE, sep="\t")
# working at home:
Nvu<-read.table('C:/Users/Rinaldo C. Bertossa/Dropbox/Others/Gerard_NGS/Data/GO/UniProt_20150331.txt', header=TRUE, quote='', stringsAsFactors=FALSE, sep="\t")
Nvq<-read.table('C:/Users/Rinaldo C. Bertossa/Dropbox/Others/Gerard_NGS/Data/GO/QuickGO_UniProt_associations_20150331.txt', header=TRUE, quote='', stringsAsFactors=FALSE, sep="\t")

Nvu<-subset(Nvu,Nvu$Organism=="Nasonia vitripennis (Parasitic wasp)")


# 1. Pool GO terms in Nvq according to UniProt ID and add this column to Nvq: NvqGO.UP

NvqGOUP <- data.frame(matrix(NA,nrow=(nrow(Nvq)),ncol=1))
colnames(NvqGOUP)<-"NvqGO.UP"

for (i in 1:nrow(Nvq)) {
  if (Nvq$ID[i]=="") {
    next
  }
  else {
    r <- which(Nvq$ID==Nvq$ID[i])
    NvqGOUP$NvqGO.UP[i] <- toString(unique(Nvq$GO.ID[r]))
  }  
}

Nvq <- cbind(Nvq,NvqGOUP)


# 2. Pool GO terms in Nvq according to Symbol and add this column to Nvq: NvqGO.Symbol

NvqGOSymbol <- data.frame(matrix(NA,nrow=(nrow(Nvq)),ncol=1))
colnames(NvqGOSymbol)<-"NvqGO.Symbol"

for (i in 1:nrow(Nvq)) {
  if (Nvq$Symbol[i]=="-") {
    next
  }
  else {
    r <- which(Nvq$Symbol==Nvq$Symbol[i])
    NvqGOSymbol$NvqGO.Symbol[i] <- toString(unique(Nvq$GO.ID[r]))
  }  
}

Nvq <- cbind(Nvq,NvqGOSymbol)

# 3. Pool all Nvq GO terms together

NvqGOPool <- data.frame(matrix(NA,nrow=nrow(Nvq),ncol=1))
colnames(NvqGOPool) <- "NvqGOPool"

for (i in 1:nrow(Nvq)) {
  a <- Nvq$GO.ID[i]
  b <- unlist(strsplit(Nvq$NvqGO.UP[i],", "))
  c <- unlist(strsplit(Nvq$NvqGO.Symbol[i],", "))
  NvqGOPool$NvqGOPool[i] <- toString(na.omit(unique(c(a,b,c))))
}

Nvq <- cbind(Nvq,NvqGOPool)


# 4. Import pooled QuickGO GO terms into Nvu (into new column) according to common UniProt

NvuqGO <- data.frame(matrix(NA,nrow=(nrow(Nvu)),ncol=2))
colnames(NvuqGO)<-c("NvqGOPool.UP","NvqGOPool.Gene")

for (i in 1:nrow(Nvu)) {
  if (Nvu$Entry[i]=="") { # actually all records have a UniProt ID and therefore this line is superflous
    next
  }
  else {
    r <- which(Nvq$ID==Nvu$Entry[i])
    l <- length(r)
    if (l==0) {
      next
    }
    else {
      NvuqGO$NvqGOPool.UP[i] <- Nvq$NvqGOPool[r[1]] # No need for pooling because already pooled at origin (Nvq)
    }
  }  
}


# 5. Import pooled QuickGO GO terms into Nvu (into new column) according to common Gene names

for (i in 1:nrow(Nvu)) {
  if (Nvu$Gene.names[i]=="") { # actually all records have a UniProt ID and therefore this line is superflous
    next
  }
  else {
    NvuGene<-unlist(strsplit(Nvu$Gene.names[i]," "))  # split the terms in column 4 (Gene_names)
    NvuGene.l<-length(NvuGene)
    for (c in 1:NvuGene.l) {
      r <- which(Nvq$Symbol==NvuGene[c])
      l <- length(r)
      if (l==0) {
        next
      }
      else {
        NvuqGO$NvqGOPool.Gene[i] <- Nvq$NvqGOPool[r[1]]
      }
    }
  }  
}

Nvu <- cbind(Nvu, NvuqGO)


# 6. Pool Nvu GO columns: 1st time

NvuGOPool1 <- data.frame(matrix(NA,nrow=nrow(Nvu),ncol=1))
colnames(NvuGOPool1) <- "NvuGOPool1"

for (i in 1:nrow(Nvu)) {
  a <- unlist(strsplit(Nvu$Gene.ontology.IDs[i],"; "))
  b <- unlist(strsplit(Nvu$NvqGOPool.UP[i],", "))
  c <- unlist(strsplit(Nvu$NvqGOPool.Gene[i],", "))
  NvuGOPool1$NvuGOPool1[i] <- toString(na.omit(unique(c(a,b,c))))
}

Nvu <- cbind(Nvu,NvuGOPool1)


# 7. Pool Nvu GO terms according to UniProt ID and add this column to Nvu

# No need: Nvu records have unique UniProt


# 8. Pool Nvu GO terms in Nvu according to Gene_name and add this column to Nvu: NvuGO.Gene

NvuGOGene <- data.frame(matrix(NA,nrow=(nrow(Nvu)),ncol=1))
colnames(NvuGOGene)<-"NvuGO.Gene"

for (i in 1:nrow(Nvu)) {
  if (Nvu$Gene.names[i]=="") {
    next
  }
  else {
    NvuGene<-unlist(strsplit(Nvu$Gene.names[i]," "))  # split the terms in column 'Gene.names'
    NvuGene.l<-length(NvuGene)
    for (c in 1:NvuGene.l) {
      r <- which(Nvu$Gene.names==NvuGene[c])
      l <- length(r)
      if (l<2) {
        next
      }
      else {
        NvuGOGene$NvuGO.Gene[i] <- toString(unique(unlist(strsplit(Nvu$NvuGOPool1[r],", "))))
      }
    }  
  }
}

Nvu <- cbind(Nvu,NvuGOGene)


# 9. Pool Nvu GO columns: 2nd time

NvuGOPool2 <- data.frame(matrix(NA,nrow=nrow(Nvu),ncol=1))
colnames(NvuGOPool2) <- "NvuGOPool2"

for (i in 1:nrow(Nvu)) {
  a <- unlist(strsplit(Nvu$NvuGOPool1[i],", "))
  b <- unlist(strsplit(Nvu$NvuGO.Gene[i],", "))
  NvuGOPool2$NvuGOPool2[i] <- toString(na.omit(unique(c(a,b))))
}

Nvu <- cbind(Nvu,NvuGOPool2)


# Creating and saving a final Nvu version (containing the Nasonia GO terms so far)

colnames(Nvu)
colnames <- c("Entry","Entry.name","Status","Protein.names","Gene.names","Length","Alternative.products","Cross.reference..RefSeq.","NvuGOPool2")
NvuFinal <- Nvu[,colnames]
write.table(NvuFinal, file="NvuFinal.txt", quote=FALSE, sep='\t')
# Test it:
ts<-read.table('C:/Users/Rinaldo C. Bertossa/Dropbox/Others/Gerard_NGS/R/NvuFinal.txt', header=TRUE, quote='', stringsAsFactors=FALSE, sep="\t")
# The file is produced in the correct way


# Import 'genes' file (=the list of genes as produced by NextGENe's 'Expression Report' for 'genes'):
getwd()
# Working at work:
genes<-read.table('C:/Users/p233711/Dropbox/Others/Gerard_NGS/Data/NasoniaExpressionReports/Genes/L6_AAAAAT_converted_Output_v20v11_Nvit21_Expression_Report_Genes.txt', skip=9, header=TRUE, quote='', stringsAsFactors=FALSE, sep="\t")
# Working at home:
genes<-read.table('C:/Users/Rinaldo C. Bertossa/Dropbox/Others/Gerard_NGS/Data/NasoniaExpressionReports/Genes/L6_AAAAAT_converted_Output_v20v11_Nvit21_Expression_Report_Genes.txt', skip=9, header=TRUE, quote='', stringsAsFactors=FALSE, sep="\t")

# function to remove '; +/-' from genes$Gene:
f1<-function (x) sub("[.*;].*", "", x, perl=TRUE)
# genes names prepared for searching:
genes<-cbind(genes[1:5],apply(genes[6],2,f1),genes[7:20],stringsAsFactors=FALSE) 


# 10. Complete 'genes' data with 'Nvit_2.1_NCBI_ProteinTable449_202468_20150331.txt'

NvitNCBI<-read.table('C:/Users/Rinaldo C. Bertossa/Dropbox/Others/Gerard_NGS/Data/Nvit_2.1_NCBI_ProteinTable449_202468_20150331.txt', header=TRUE, stringsAsFactors=FALSE, sep="\t")


# 10.1 searching in NvitNCBI$Protein.product with genes$Protein.Accession for extra Gene names and GeneID

com <- data.frame(matrix(nrow=nrow(genes),ncol=(6)))
colnames(com) <- c("NrHits","ExtraLocus","GeneID","genes.Proteins","ProtAcc.Gene","ProtAcc.GeneID")

com$genes.Proteins <- genes$Protein.Accession

for (i in 1:nrow(genes)) {
  if (genes$Protein.Accession[i]=="") {
    next
  }
  else {
    q <- sub("[.*].*", "", genes$Protein.Accession[i], perl=TRUE)
    r <- grep(q, NvitNCBI$Protein.product, fixed=TRUE)
    vLocus <- NvitNCBI$Locus[r]
    vLocus <- vLocus[vLocus[]!="-"]
    com$NrHits[i] <- length(r)
    com$ExtraLocus[i] <- toString(unique(vLocus))
    vID <- NvitNCBI$GeneID[r]
    com$GeneID[i] <- toString(unique(vID))
  }
}

length(which(com$NrHits==1)) # 7451 genes$Protein.Accession records matched in NvitNCBI
length(which(com$NrHits==2)) # At least 24 NvitNCBI$Protein.product records are repeated
grep(",",com$ExtraLocus,fixed=TRUE) # 0 > there are no extra 'Gene' descriptors in NvitNCBI corresponding to the same proteins.


# 10.2 searching NvitNCBI$Locus with genes$Gene for extra Protein names and GeneID

for (i in 1:nrow(genes)) {
  if (genes$Gene[i]=="") {
    next
  }
  else {
    q <- sub("[.*;].*", "", genes$Gene[i], perl=TRUE)
    r <- which(NvitNCBI$Locus==q)
    v <- NvitNCBI$Protein.product[r]
    com$NrHits[i] <- length(r)
    com$GeneID[i] <- toString(unique(na.omit(unlist(strsplit(c(com$GeneID[i],NvitNCBI$GeneID[r]),", ")))))
    com$ProtAcc.Gene[i] <- toString(unique(v))
  }
}

r <- head(grep(",", com$ProtAcc.Gene, fixed=TRUE))
length(r)
com[r,] # THere are 6 cases where extra protein accessions are added...
grep(",", com$GeneID, fixed=TRUE) # no multiple GeneID are found...


# 10.3 searching NvitNCBI$GeneID with com$GeneID for extra Protein names

com$GeneID[is.na(com$GeneID)]<-""
com$NrHits <- NA

for (i in 1:nrow(com)) {
  if (com$GeneID[i]=="") {
    next
  }
  else {
    r <- which(NvitNCBI$GeneID==com$GeneID[i])
    v <- NvitNCBI$Protein.product[r]
    com$NrHits[i] <- length(r)
    com$ProtAcc.GeneID[i] <- toString(unique(v))
  }
}

length(which(com$NrHits>1))
# Evidently new Protein accessions are found. THey are added to genes1:

# 10.4 'ExtraLocus', 'GeneID', 'ProtAcc.Gene' and 'ProtAcc.GeneID' are added to 'genes' and where needed pooled

com$ExtraLocus[is.na(com$ExtraLocus)]<-""
com$GeneID[is.na(com$GeneID)]<-""
com$ProtAcc.Gene[is.na(com$ProtAcc.Gene)]<-""
com$ProtAcc.GeneID[is.na(com$ProtAcc.GeneID)]<-""

length(which(com$ExtraLocus!=""))

sum <- data.frame(matrix(nrow=nrow(genes),ncol=2))
colnames(sum) <- c("Gene1","Protein.Accession1")

for (i in 1:nrow(genes)) {
   sum$Gene1[i] <- toString(unique(na.omit(unlist(strsplit(c(genes$Gene[i],com$ExtraLocus[i]), ", ")))))
   sum$Protein.Accession1[i] <- toString(unique(na.omit(unlist(strsplit(c(genes$Protein.Accession[i],com$ProtAcc.Gene[i],com$ProtAcc.GeneID[i]), ", ")))))
}

genes1 <- cbind(genes[1:5],com[3],sum[1],genes[7:9],sum[2],genes[11:20])



# 11. Import Nvu GOs into genes1 (into new column) according to common Gene name

GOgenes <- data.frame(matrix(nrow=nrow(genes1),ncol=1))
colnames(GOgenes) <- "GO.genes"

for (i in 1:nrow(Nvu)) {
  if (Nvu$Gene.names[i]=="") {
    next
  }
  else {
    gen<-unlist(strsplit(Nvu$Gene.names[i]," "))  # split the terms in column 4 (Gene_names)
    gen.l<-length(gen)
    for (c in 1:gen.l) {
      r <- which(genes1$Gene1==gen[c])
      l <- length(r)
      if (l==0) {
        next
      }
      else {
        GOgenes$GO.genes[r] <- Nvu$NvuGOPool2[i]
      }
    }
  }  
}

genes1 <- cbind(genes1,GOgenes)


# 12. Import Nvu GOs into genes1 (into new column) according to common RefSeq

GORef <- data.frame(matrix(nrow=nrow(genes1),ncol=1))
colnames(GORef) <- "GO.Ref"

for (i in 1:nrow(Nvu)) {
  if (Nvu$Cross.reference..RefSeq.[i]=="") {
    next
  }
  else {
    gen <- unlist(strsplit(Nvu$Cross.reference..RefSeq.[i],";"))  # split the RefSeq terms
    gen <- gsub("[.*].*", "", gen, perl=TRUE)
    gen.l<-length(gen)
    for (c in 1:gen.l) {
      r <- grep(gen[c], genes1$Protein.Accession1, fixed=TRUE)
      l <- length(r)
      if (l==0) {
        next
      }
      else {
        GORef$GO.Ref[r] <- toString(unique(na.omit(unlist(strsplit(c(GORef$GO.Ref[r],Nvu$NvuGOPool2[i]), ", ")))))
      }
    }
  }  
}

genes1 <- cbind(genes1,GORef)


# 13. Pool all genes1 GOs columns

genesGOPool <- data.frame(matrix(nrow=nrow(genes1),ncol=1))
colnames(genesGOPool) <- "GOPool"

for (i in 1:nrow(genes1)) {
  a <- unlist(strsplit(genes1$GO.genes[i],", "))
  b <- unlist(strsplit(genes1$GO.Ref[i],", "))
  genesGOPool$GOPool[i] <- toString(na.omit(unique(c(a,b))))
}

genes1 <- cbind(genes1,genesGOPool)


# 14. Pooling GO terms in genes1 according to Gene

GOgenes <- data.frame(matrix(nrow=nrow(genes1),ncol=3))
colnames(GOgenes) <- c("GO.Gene","GO.GeneID","GO.ProtAcc")

for (i in 1:nrow(genes1)) {
  if (genes1$Gene1[i]=="") {
    next
  }
  else {
    r <- which(genes1$Gene1==genes1$Gene1[i])
    GOgenes$GO.Gene[i] <- toString(na.omit(unique(unlist(strsplit(genes1$GOPool[r],", ")))))
  }
}

# We pool the result with the original GO in genes1$GOPool:
for (i in 1:nrow(genes1)) {
  genes1$GOPool[i] <- toString(unique(na.omit(unlist(strsplit(c(genes1$GOPool[i], GOgenes$GO.gene[i]), ", ")))))
}


# 15. Pooling GO terms in genes1 according to GeneID

for (i in 1:nrow(genes1)) {
  if (genes1$GeneID[i]=="") {
    next
  }
  else {
    r <- which(genes1$GeneID==genes1$GeneID[i])
    GOgenes$GO.GeneID[i] <- toString(na.omit(unique(unlist(strsplit(genes1$GOPool[r],", ")))))
  }
}

# We pool the result with the original GO in genes1$GOPool:
for (i in 1:nrow(genes1)) {
  genes1$GOPool[i] <- toString(unique(na.omit(unlist(strsplit(c(genes1$GOPool[i], GOgenes$GO.GeneID[i]), ", ")))))
}


# 16. Pooling GO terms in genes1 according to RefSeq (Protein.Accession)

for (i in 1:nrow(genes1)) {
  if (genes1$Protein.Accession1[i]=="") {
    next
  }
  else {
    v <- unique(unlist(strsplit(genes1$Protein.Accession1[i], ", ")))
    v.l <- length(v)
    for (ii in 1:v.l) {
      r <- grep(v[ii], genes1$Protein.Accession1, fixed=TRUE)
      GOgenes$GO.ProtAcc[i] <- toString(na.omit(unique(unlist(strsplit(genes1$GOPool[r],", ")))))
    }
  }
}

# We pool the result with the original GO in genes1$GOPool:
for (i in 1:nrow(genes1)) {
  genes1$GOPool[i] <- toString(unique(na.omit(unlist(strsplit(c(genes1$GOPool[i], GOgenes$GO.ProtAcc[i]), ", ")))))
}


# 17. We give a new name to the df: genes2

genes2 <- cbind(genes1[1:21],genes1[24])

# We save this file for later

write.table(genes2, file="genes2.txt", quote=FALSE, sep='\t')
# Test it:
ts<-read.table('C:/Users/Rinaldo C. Bertossa/Dropbox/Others/Gerard_NGS/R/genes2.txt', header=TRUE, quote='', stringsAsFactors=FALSE, sep="\t")
# The file is produced in the correct way
# genes2 contains the genes list as produced by NextGENe with associated IDs and GOs found so far

# 18. Preparing a file to send to Gerard

write.table(genes2, file="Nvit_GOs_v4.txt", quote=FALSE, sep='\t')
# Check: whether the file is produced correctly:
Nvit_GOs_v4<-read.table('C:/Users/Rinaldo C. Bertossa/Dropbox/Others/Gerard_NGS/R/Nvit_GOs_v4.txt', header=TRUE, quote='', stringsAsFactors=FALSE, sep="\t")
# The file is produced in the correct way

# Gerard wants one gene per line...
# 19. List of single genes with associated single GO terms

# 19.1 We first remove duplicates

which(genes2$GOPool!=""&genes2$Gene1=="")
# If there is a GO term there is also a Gene name

GenesSH<-genes2[genes2$Gene1!="",] # 'SH' = Shortlist

du <- data.frame(matrix(nrow=nrow(GenesSH),ncol=1))
colnames(du) <- c("duplicates")

for (i in 1:nrow(GenesSH)) {
  v <- which(GenesSH$Gene1==GenesSH$Gene1[i])
  du[v,1] <- rank(v)
}

length(which(du$duplicates==2))
# 75 genes have at least a duplicate

GenesSHu <- cbind(GenesSH,du) # 'u' = unique 
GenesSHu <- GenesSHu[GenesSHu$du==1,]
which(GenesSHu$du>1)


# 19.2 We create the dataframe

Df <- data.frame(Gene=NA,GO=NA)

for (i in 1:nrow(GenesSHu)) {
  if (GenesSHu$GOPool[i]=="") {
    next
  }
  else {
    GO <- unique(unlist(strsplit(GenesSHu$GOPool[i], ", ")))
    df <- data.frame(Gene=GenesSHu$Gene1[i],GO)
    Df <- rbind(Df,df)
  }
}

Df <- Df[-1,]
rownames(Df) <- NULL

write.table(Df, file="Nvit_GOs_v5.txt", quote=FALSE, sep='\t')
# Check: whether the file is produced correctly:
Nvit_GOs_v5<-read.table('C:/Users/Rinaldo C. Bertossa/Dropbox/Others/Gerard_NGS/R/Nvit_GOs_v5.txt', header=TRUE, quote='', stringsAsFactors=FALSE, sep="\t")
# The file is produced in the correct way

# Nr. of genes in the list: 5205
# Nr. of rows: 21076
# ~4.05 GO terms per gene

####################################
# 19. OLDER version: List of genes with associated GO terms

Nvit_GOs_v3<-read.table('C:/Users/Rinaldo C. Bertossa/Dropbox/Others/Gerard_NGS/R/Nvit_GOs_v3.txt', header=TRUE, quote='', stringsAsFactors=FALSE, sep="\t")
# The file is produced in the correct way

# 19.1 We first remove duplicates

which(Nvit_GOs_v3$GO.ID!=""&Nvit_GOs_v3$Gene=="")
# If there is a GO term there is also a Gene name

GenesSH<-Nvit_GOs_v3[Nvit_GOs_v3$Gene!="",] # 'SH' = Shortlist

du <- data.frame(matrix(nrow=nrow(GenesSH),ncol=1))
colnames(du) <- c("duplicates")

for (i in 1:nrow(GenesSH)) {
  v <- which(GenesSH$Gene==GenesSH$Gene[i])
  du[v,1] <- rank(v)
}

length(which(du$duplicates==2))
# 75 genes have at least a duplicate

GenesSHu <- cbind(GenesSH,du) # 'u' = unique 
GenesSHu <- GenesSHu[GenesSHu$du==1,]
which(GenesSHu$du>1)


# 19.2 We create the dataframe

Df <- data.frame(Gene=NA,GO=NA)

for (i in 1:nrow(GenesSHu)) {
  if (GenesSHu$GO.ID[i]=="") {
    next
  }
  else {
    GO <- unique(unlist(strsplit(GenesSHu$GO.ID[i], ", ")))
    df <- data.frame(Gene=GenesSHu$Gene[i],GO)
    Df <- rbind(Df,df)
  }
}

Df <- Df[-1,]
rownames(Df) <- NULL
