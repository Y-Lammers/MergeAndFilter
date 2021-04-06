#!/usr/bin/env Rscript

# R script for merging two sets of OBITools output TSV files. The script appies
# a number of basic filtering criteria and outputs three tables. One table
# with the filtered results, a summary table for the filterd results and a 
# table that contains the number and proportion of reads for each sample
# in the datasets during the filtering.

# USAGE: MergeAndFilter.R [number of input files (1)] [obitools output x]
# [obitools table x name] [count file] [synthetic blacklist file] 
# [region blacklist file] [base output name] [min iden (1)] [min reads (3)] 
# [min total reads (10)] [min total repeats (3)]

# note: at least one obitools tables, count, blacklist files and output name
# are required. The remaining settings will default to the values below.

# Contact: youri.lammers@gmail.com
# Version: 2.1.4

# set arguments

# get the number of input files (default 1)
if (is.na(commandArgs(trailingOnly = TRUE)[1])) {
	datasets=1
} else {
	datasets=as.numeric(commandArgs(trailingOnly = TRUE)[1])
}

# get the obitools output files and names in a separate list
obifiles <- commandArgs(trailingOnly = TRUE)[2:((datasets*2)+1)]

# get the count table, blacklists and output name
count_table=commandArgs(trailingOnly = TRUE)[(2*datasets)+2]
synthetic_blacklist_name=commandArgs(trailingOnly = TRUE)[(2*datasets)+3]
region_blacklist_name=commandArgs(trailingOnly = TRUE)[(2*datasets)+4]
output_name=commandArgs(trailingOnly = TRUE)[(2*datasets)+5]

# minimum identity
if (is.na(commandArgs(trailingOnly = TRUE)[(2*datasets)+6])) {
	min_iden=1
} else {
	min_iden=as.numeric(commandArgs(trailingOnly = TRUE)[(2*datasets)+6])
}

# minimum reads
if (is.na(commandArgs(trailingOnly = TRUE)[(2*datasets)+7])) {
	min_reads=3
} else {
	min_reads=as.numeric(commandArgs(trailingOnly = TRUE)[(2*datasets)+7])
}

# minimum total reads
if (is.na(commandArgs(trailingOnly = TRUE)[(2*datasets)+8])) {
	min_total_reads=10
} else {
	min_total_reads=as.numeric(commandArgs(trailingOnly = TRUE
	)[(2*datasets)+8])
}

# minimum total repeats
if (is.na(commandArgs(trailingOnly = TRUE)[(2*datasets)+9])) {
	min_total_rep=3
} else {
	min_total_rep=as.numeric(commandArgs(trailingOnly = TRUE
	)[(2*datasets)+9])
}


#################################
# Values for testing or running #
# the R script manually         #
#################################

#obifiles=c("obifile 1","name 1","obifile 2","name 2")
#datasets=length(obifiles)/2
#count_table="count file"
#synthetic_blacklist_name="synthetic blacklist"
#region_blacklist_name="regional blacklist"
#output_name="output name"
#min_iden=1
#min_reads=3
#min_total_reads=10
#min_total_rep=3



############################
# Read the obitools tables #
############################

# lists used for renaming and ordering of the info columns

# list of columns to rename
rename <- c("best_identity.*","family","family_name","genus","genus_name",
	"rank","scientific_name","taxid","species_list.*")

# list of columns and order for the id, count and replicate information
first_set <- c("id","pre_read","pre_rep","post_read","post_rep",
		"read_ratio","rep_ratio")

# loop through the number of files provided
for (i in 1:datasets) { 

	# Load the OBITools table
	obi = read.table(obifiles[((i-1)*2)+1],header=TRUE,sep="\t",quote="",comment.char="")

	# get the table name and clean it for odd characters
	name = make.names(obifiles[((i-1)*2)+2])

	# Rename the informative columns, use grep to locate in case
	# the columns have been ordered. Store the location so that the
	# columns can be reordered
	col_info <- c()
	for (col in rename) {

		# store the column location
		col_info <- c(col_info,grep(paste("^",col,"$",sep=""),
			colnames(obi)))

		# get the column name and rename it
		colnames(obi)[grep(paste("^",col,"$",sep=""),colnames(obi))] <-
			sub("best_identity","iden",sub(".*","",paste(name,"_",
			col,sep=""),fixed=TRUE))
	}

	# order the rows by the barcode ID (incase the input files have 
	# been order differently)
	obi <- obi[order(obi$id),]

	# Extract the obitools information, if some is already stored,
	# add the information, otherwise, create a new dataframe.
	if (exists("obi_info")) {

		# add the obi info to the table
		obi_info <- cbind(obi_info,obi[,col_info])

	} else {

		# add blank columns for the replicate info, copy the count info
		# and add a blank column for the post filter count data
		obi$pre_read <- obi$count
		obi[,c("pre_rep","post_read","post_rep","read_ratio",
			"rep_ratio")] <- NA

		# if its the first obitools dataset, get both the count info
		# identification info and the sequences
		obi_info <- obi[,first_set]
		obi_info <- cbind(obi_info,obi[,col_info])
		sequences <- obi[,"sequence",drop=FALSE]
	}

	# store the actual occurence and obiclean data, only do this for
	# the last input file in order to reduce the max memory usage
	if (i == datasets) {

		# extract the sample and obiclean data
		sample <- obi[,which(grepl("sample.",colnames(obi),fixed=TRUE))]
		clean <- obi[,which(grepl("obiclean_status.",
			colnames(obi),fixed=TRUE))]

		# combine the raw sample and obiclean data
		rcombi=cbind(sample,clean)

		# remove the separate variables
		rm(sample, clean)
	}

	# clear the memory
	rm(obi)
}



###################################################
# Read the counts data and generate a blank table #
###################################################

# read the raw count data
counts = read.table(count_table,header=TRUE,sep="\t")

# reformat the sample names in the counts table so that they
# match with the OBITools input
counts[,1] <- gsub("[][():+-]",".",counts[,1])

# extract the sample names
csample <- unique(sub("(.*[^0-9])([0-9]+$)","\\1",counts[,1]))

# get the maximum number of repeats in the library
repeats <- sort(unique(sub("(.*[^0-9])([0-9]+$)","\\2",counts[,1]),na.rm=TRUE))

# create an empty dataframe for the sample stat information
samplestat <- data.frame(matrix(NA,nrow=length(csample),
	ncol=((7*(length(repeats)+2))+10)))

# fix the row and column names
# add the rownames based on the sample names
rownames(samplestat) <- csample

# create the column names for the repeats based on 
# the number of samples
tcolnames <- c()
for (cat in c("raw","prop_raw","prop_filt","prop_noniden",
	"prop_synthetic_blacklist","prop_region_blacklist","prop_retained")) {
	for (num in repeats) {
		tcolnames <- c(tcolnames, paste(cat,num,sep=""))
	}
	tcolnames <- c(tcolnames,paste(cat,"_avg",sep=""))
	tcolnames <- c(tcolnames,paste(cat,"_sd",sep=""))
}
				
# add the additional columns for the various sample averages
tcolnames <- c(tcolnames,c("avg_rep","sd_rep","avg_filt_rep","sd_filt_rep",
	"overlap","avg_seq_length","sd_seq_length","avg_single_seq_length",
	"sd_single_seq_length","single_seq_count"))

# add the column names
colnames(samplestat) <- tcolnames



#########################################
# Add blank columns for missing repeats #
#########################################

# for sample in the count file
for (sample in counts[,1]) {

	# get the column position
	pos <- grep(paste("^sample.",sample,sep=""),colnames(rcombi))

	# if the position does not exist, add the sample and obiclean sample
	# to the raw combined dataframe
	if(length(pos) == 0) {

		rcombi[,paste("sample.",sample,sep="")] <- 0
		rcombi[,paste("obiclean_status.",sample,
			sep="")] <- as.factor(NA)
	}
}

#sort the new dataframe, so the the columns are grouped together
rcombi <- rcombi[,order(names(rcombi))]



#####################################################
# Add the raw counts data to the sample stat table, #
# both as raw counts and proportion of the sample.  #
#####################################################

# Parse through the sample names
for (sample in csample) {

	# get the rows with data
	pos <- grep(paste("^",sample,"[0-9]+$",sep=""),counts$Sample.name)

	# get the row for the sample in the samplestat table
	samplepos <- grep(paste("^",sample,"$",sep=""),rownames(samplestat))

	# get the total raw read count for the sample
	totraw <- sum(counts$Read.count[pos])

	# create an empty vector that stores the repeat
	# proportions in order to calc the standard deviation
	sampleprop <- c()	

	# calculate the proportion of raw reads per repeat
	# and add it to the sample stat table, as well as
	# the raw number itself.
	for (rep in pos) {

		# get the repeat number so it can be stored in
		# the proper column
		repname <- as.character(counts$Sample.name[rep])
		repchar <- sub("(.*[^0-9])([0-9]+$)","\\2",repname)

		# calc the read proportion for the repeat
		prep <-	counts[grep(paste("^",sample,repchar,"$",sep=""),
				counts[,1]),2]/totraw

		# add the proportion to the sampleprop vector
		sampleprop <- c(sampleprop, prep)

		# add the raw repeat counts and proportional repeat
		# counts to the sample stat table
		samplestat[samplepos,grep(paste("raw",repchar,"$",sep=""),
			colnames(samplestat))] <-
			counts[grep(paste("^",sample,repchar,"$",sep=""),
				counts[,1]),2]
		samplestat[samplepos,grep(paste("prop_raw",repchar,"$",sep=""),
			colnames(samplestat))] <- prep
	}

	# calculate the average and standard deviation for the
	# read counts and proportion of reads and add them
	# to the samplestat table
	samplestat$raw_avg[samplepos] <- 
		mean(counts$Read.count[pos],na.rm=TRUE)
	samplestat$raw_sd[samplepos] <- 
		sd(counts$Read.count[pos],na.rm=TRUE)
	samplestat$prop_raw_avg[samplepos] <- mean(sampleprop,na.rm=TRUE)
	samplestat$prop_raw_sd[samplepos] <- sd(sampleprop,na.rm=TRUE)
}



#######################
# Detect homopolymers #
#######################

# create two new rows for the homopolymer info
sequences$homopolymer <- NA
sequences$homopolymer_type <- NA

# go through the sequences
for (seq in 1:nrow(sequences)) {

	# search for homopolymers that are 5bp or longer
	rep <- regexpr("(.)\\1{4,}",sequences[seq,"sequence"])
	if(rep!=-1) {
		# if present: add the info to the columns
		sequences[seq,"homopolymer"] <- attr(rep,"match.length")
		sequences[seq,"homopolymer_type"] <- substr(
			sequences[seq,"sequence"],rep,rep)
	}
}



########################################
# remove both internal and sequence    #
# below the minimum sequence threshold #
########################################

# copy the raw combined sequence table
tcombi <- rcombi

# detect the sample columns
samples <- grep("sample.",colnames(tcombi),fixed=TRUE)

# for sample in the sample list
for (s in 1:length(samples)) {

	# find the matching obiclean column
	obiclean <- grep(paste("^",sub("sample","obiclean_status",
		colnames(tcombi)[samples[s]]),"$",sep=""), colnames(tcombi))

	# remove low sequence occurences and internal sequences
	tcombi[which(tcombi[,samples[s]]<min_reads | 
		tcombi[,obiclean]=="i"), samples[s]] <- 0
}	
	
# recalculate the total count
for (seq in 1:nrow(tcombi)) {
	obi_info[seq,"post_read"] <- sum(tcombi[seq,samples])
}



################################################
# Remove sequences below the minimum total     #
# count or below the minimum number of repeats #
################################################

# remove the low total read count sequences
subset = which(obi_info[,"post_read"]>=min_total_reads)
rcombi=rcombi[subset,samples]
tcombi=tcombi[subset,samples]
obi_info=obi_info[subset,]
sequences=sequences[subset,]

# count the total number of repeats
obi_info[,"post_rep"] <- rowSums(tcombi>0)
obi_info[,"pre_rep"] <- rowSums(rcombi>0)

# calculate the pre-post filtering ratios
obi_info[,"read_ratio"] <- obi_info$post_read / obi_info$pre_read
obi_info[,"rep_ratio"] <- obi_info$post_rep / obi_info$pre_rep

# remove the low total repeat sequences
subset = which(obi_info[,"post_rep"]>=min_total_rep)
tcombi=tcombi[subset,]
obi_info=obi_info[subset,]
sequences=sequences[subset,]


##################################
# remove sequences with both     #
# identities below the threshold #
##################################

# copy the technical filtered combined sequence table
icombi <- tcombi

# find the rows that match the minimum identity
subset = which(rowSums(obi_info[,grep("_iden",
	colnames(obi_info)),drop=FALSE]>=min_iden)>=1)

# apply the filtering to the data
icombi=icombi[subset,]
obi_info=obi_info[subset,]
sequences=sequences[subset,,drop=FALSE]



###########################################
# Apply the synthetic blacklist filtering #
###########################################

# Read the blacklist
synthetic_blacklist = read.table(synthetic_blacklist_name,sep="\t",fill=TRUE)

# copy the identity filtered combined sequence table
sbcombi <- icombi

# find the rows that match sequences in the blacklist
# and select the non matching rows
subset = which(!sequences[,"sequence"] %in% synthetic_blacklist[,1])

# apply the filtering to the data
sbcombi=sbcombi[subset,]
obi_info=obi_info[subset,]
sequences=sequences[subset,,drop=FALSE]



########################################
# Apply the region blacklist filtering #
########################################

# Read the blacklist #
region_blacklist = read.table(region_blacklist_name,sep="\t",fill=TRUE)

# copy the identity filtered combined sequence table
rbcombi <- sbcombi

# find the rows that match sequences in the blacklist
# and select the non matching rows
subset = which(!sequences[,"sequence"] %in% region_blacklist[,1])

# apply the filtering to the data
rbcombi=rbcombi[subset,]
obi_info=obi_info[subset,]
sequences=sequences[subset,,drop=FALSE]



##########################
# write the cleaned data #
##########################

# get the final table
finaltable=cbind(obi_info,sequences,rbcombi)

# resort the final table based on the total read count
finaltable <- finaltable[order(-finaltable[,"post_read"]),]

write.table(finaltable, file=paste(output_name,"_filtered.tsv",sep=""),
	quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)



############################################
# get the sample repeat count, read sums   #
# and the (weighted) proportion of repeats #
############################################

# copy the final filtered table for summarization
scombi <- rbcombi

# get the unique sample names (by removing the last character
# from the sample name, which is presumed to the repeat number.)
usamples <- unique(sub("(.*[^0-9])([0-9]+$)","\\1",colnames(scombi)))

# loop through the samples
for (us in usamples) {

	# get the relevant columns and column count
	pos <- grep(paste("^",us,"[0-9]+$",sep=""),colnames(scombi))
	lpos <- length(pos)

	# Calculate the proportion of reads per repeat.
	pprop <- colSums(scombi[,pos,drop=FALSE])/sum(scombi[,pos])

	# get the read proportions within a sample
	temp <- as.data.frame(lapply(scombi[,pos,drop=FALSE],
		function(x) prop.table(x)))
	
	# replace empty NaN columns with zeros, so the that the average
	# calculation doesn't break
	temp[is.na(temp)] <- 0

	# get the sample sums, replicates and the proportion of replicates
	totsumv <- rowSums(scombi[,pos,drop=FALSE])
	totrepv <- rowSums(scombi[,pos,drop=FALSE]>0)
	proprepv <- totrepv/lpos

	# and the mean plus standard deviation		
	avgpreadv <- rowMeans(temp,na.rm=TRUE)
	sdpreadv <- apply(temp,1,sd,na.rm=TRUE)

	# remove NAs from single rep samples.		
	sdpreadv[is.na(sdpreadv)] <- 0

	# get the number of barcodes
	barcodecount <- sum(rowSums(scombi[,pos,drop=FALSE]>0)>0)

	# store the results
	scombi[,paste("totread_",us,sep="")] <- totsumv
	scombi[,paste("avgpread_",us,sep="")] <- avgpreadv
	scombi[,paste("sdpread_",us,sep="")] <- sdpreadv
	scombi[,paste("totrep_",us,sep="")] <- totrepv
	scombi[,paste("proprep_",us,sep="")] <- proprepv

	# calculate the mean number of repeats for the sample
	# (ignoring 0 repeats) as well as the proportion of 
	# mean repeats
	mean_rep <- mean(totrepv[totrepv>0])
	rmean_rep <- mean_rep / lpos

	# if the proportion of repeats is higher than 0.33 AND there are
	# atleast 10 different barcodes in the sample:
	# calculate the weighted proportion of repeats per sequence.
	# if proportion is lower, use the regular proportion.
	if ((rmean_rep >= 0.33) & (barcodecount>=10)) {
	
		# calculate the weighted proportion of repeats
		wproprepv <- rowSums(t(t(scombi[,pos,drop=FALSE]>0)*pprop))

		# add the weighted proportion to the dataframe
		scombi[,paste("weightrep_",us,sep="")] <- wproprepv

	} else {

		# reuse the regular proprep if the ratio is too low
		scombi[,paste("weightrep_",us,sep="")] <- 
			scombi[,paste("proprep_",us,sep="")]
	}
}



################################
# Write the summarized results #
################################

# get the positions for the summarized columns
pos <- grep("totread_|avgpread_|sdpread_|totrep_|proprep_|weightrep_",
	colnames(scombi))

# get the summary subset
summary <- scombi[,pos]

# get the final table
finaltable=cbind(obi_info,sequences,summary)

# re-sort the final table based on the total read count
finaltable <- finaltable[order(-finaltable[,"post_read"]),]

# write the table
write.table(finaltable, file=paste(output_name,"_summary.tsv",sep=""),
	quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)



#############################################
# Calculate the proportion of reads removed #
# and retained for the filtering steps      #
#############################################

# create a new data frame that will hold the raw sequence lengths
lengthpile <- data.frame(matrix(NA,nrow=length(rownames(samplestat)),ncol=2))
colnames(lengthpile) <- c("reps","single")
rownames(lengthpile) <- rownames(samplestat)

# loop through the samples based on the counts table
for (us in csample) {

	# get the row in the samplestat table
	samplepos <- grep(paste("^",us,"$",sep=""),rownames(samplestat))

	# get the sample names
	repnames <- counts[grep(paste("^",us,"[0-9]+$",sep=""),counts[,1]),1]

	# create empty vectors for storing the proportional data
	# in order to calculate the mean and standard deviation
	tfilt <- c()
	ifilt <- c()
	sbfilt <- c()
	rbfilt <- c()
	retain <- c()

	# calculate the proportion of raw reads per repeat
	# and add it to the sample stat table, as well as
	# the raw number itself.
	for (rep in repnames) {

		# get the repeat number so it can be stored in
		# the proper column
		repname <- paste("sample.",rep,sep="")
		repchar <- sub("(.*[^0-9])([0-9]+$)","\\2",repname)

		# get the repeat location in the rcombi table
		repc <- grep(repname,colnames(rcombi))

		# Check if the repeat still exists in the OBITools
		# output, if yes, compute the proportions, if not
		# all data is filtered during OBITools itself and the
		# proportions default to a 100% removed.
		if (sum(rcombi[,repc])>0) {

			# get the raw number of reads for the repeat
			raw <- samplestat[samplepos,paste("raw",
				repchar,sep="")]

			# get the proportion of technical filtered reads
			pt <- (raw-sum(tcombi[,repc]))/raw

			# get the proportion of identity filtered reads
			pi <- (sum(tcombi[,repc])-sum(icombi[,repc]))/raw

			# get the proportion of synthetic blacklist
			# filtered reads
			psb <- (sum(icombi[,repc])-sum(sbcombi[,repc]))/raw

			# get the proportion of region blacklist 
			# filtered reads
			prb <- (sum(sbcombi[,repc])-sum(rbcombi[,repc]))/raw

			# get the proportion of retained reads
			pr <- sum(rbcombi[,repc])/raw

		} else {

			# set the proportions to the default when the
			# sample is not in the OBITools output
			pt <- 1
			pi <- 0
			psb <- 0
			prb <- 0
			pr <- 0
		}

		# add the proprotion of technical filtered reads
		# to the samplestat table and vector
		samplestat[samplepos,grep(paste("prop_filt",repchar,"$",
			sep=""),colnames(samplestat))] <- pt
		tfilt <- c(tfilt,pt)

		# add the proprotion of identity filtered reads
		# to the samplestat table and vector
		samplestat[samplepos,grep(paste("prop_noniden",repchar,"$",
			sep=""),colnames(samplestat))] <- pi
		ifilt <- c(ifilt,pi)

		# add the proprotion of synthetic blacklist filtered reads
		# to the samplestat table and vector
		samplestat[samplepos,grep(paste("prop_synthetic_blacklist",
			repchar,"$",sep=""),colnames(samplestat))] <- psb
		sbfilt <- c(sbfilt,psb)

		# add the proprotion of region blacklist filtered reads
		# to the samplestat table and vector
		samplestat[samplepos,grep(paste("prop_region_blacklist",
			repchar,"$",sep=""),colnames(samplestat))] <- prb
		rbfilt <- c(rbfilt,prb)

		# add the proprotion of retained reads to the
		# samplestat table and vector
		samplestat[samplepos,grep(paste("prop_retained",repchar,"$",
			sep=""),colnames(samplestat))] <- pr
		retain <- c(retain,pr)
	}

	# get the relevant columns
	pos <- grep(paste("^sample.",us,"[0-9]+$",sep=""),colnames(rbcombi))

	# create the vectors for the sequence length information
	replength <- c()
	singlelength <- c()

	# loop through the sequences, and get the sequence lengths
	# for each barcode
	for (seq in 1:nrow(rbcombi)) {

		# get the number of repeats for the sample
		repcount <- sum(rbcombi[seq,pos]>=1)

		# get the seq info
		if (repcount >= 1) {

			# add multiple seq lengths, based on reps
			replength <- c(replength, replicate(repcount,nchar(
				as.character(sequences[seq,"sequence"]))))

			# add a single seq length per sample
			singlelength <- c(singlelength, nchar(as.character(
				sequences[seq,"sequence"])))
		} 
	}

	# save the average sequence lenghts
	lengthpile[samplepos,"reps"] <- paste(replength,collapse=",")
	lengthpile[samplepos,"single"] <- paste(singlelength,collapse=",")

	# calculate the mean and standard deviation for the
	# proportion of technical filtered reads and add them
	# to the samplestat table
	samplestat[samplepos,grep("prop_filt_avg",colnames(samplestat),
		fixed=TRUE)] <- mean(tfilt)
	samplestat[samplepos,grep("prop_filt_sd",colnames(samplestat),
		fixed=TRUE)] <- sd(tfilt)

	# the same for the identity filtered reads
	samplestat[samplepos,grep("prop_noniden_avg",colnames(samplestat),
		fixed=TRUE)] <- mean(ifilt)
	samplestat[samplepos,grep("prop_noniden_sd",colnames(samplestat),
		fixed=TRUE)] <- sd(ifilt)

	# and the synthetic blacklist filtered reads
	samplestat[samplepos,grep("prop_synthetic_blacklist_avg",
		colnames(samplestat),fixed=TRUE)] <- mean(sbfilt)
	samplestat[samplepos,grep("prop_synthetic_blacklist_sd",
		colnames(samplestat),fixed=TRUE)] <- sd(sbfilt)

	# and the region blacklist filtered reads
	samplestat[samplepos,grep("prop_region_blacklist_avg",
		colnames(samplestat),fixed=TRUE)] <- mean(rbfilt)
	samplestat[samplepos,grep("prop_region_blacklist_sd",
		colnames(samplestat),fixed=TRUE)] <- sd(rbfilt)

	# and the proportion of retained reads
	samplestat[samplepos,grep("prop_retained_avg",colnames(samplestat),
		fixed=TRUE)] <- mean(retain)
	samplestat[samplepos,grep("prop_retained_sd",colnames(samplestat),
		fixed=TRUE)] <- sd(retain)

	# get the number of unique barcodes in the sample
	barlength <- length(singlelength)

	# add the number of unique barcodes in the sample
	samplestat[samplepos,grep("single_seq_count",colnames(samplestat),
		fixed=TRUE)] <- barlength

	# if there are barcodes, calculate the mean and sd, otherwise, skip
	if (barlength > 0) {

		# and the sequence length across all replicates
		samplestat[samplepos,grep("avg_seq_length",
			colnames(samplestat),fixed=TRUE)] <- mean(replength)
		samplestat[samplepos,grep("sd_seq_length",
			colnames(samplestat),fixed=TRUE)] <- sd(replength)

		# and the sequence lenght per sample (no duplicates)
		samplestat[samplepos,grep("avg_single_seq_length",
			colnames(samplestat),fixed=TRUE)] <- mean(singlelength)
		samplestat[samplepos,grep("sd_single_seq_length",
			colnames(samplestat),fixed=TRUE)] <- sd(singlelength)
	}
}

# preparte the sequence length output
lengthpile_mod <- lengthpile
lengthpile_mod <- rbind(colnames(lengthpile_mod),lengthpile_mod)
rownames(lengthpile_mod)[1] <- ""

# write the sequence length table
write.table(lengthpile_mod, file=paste(output_name,"_lengths.tsv",sep=""),
        quote=FALSE, sep="\t", col.names=FALSE, row.names=TRUE)



#######################################
# get the average repeats for the top #
# 10 filtered and unfiltered taxa     #
#######################################

# loop through the samples
for (us in usamples) {

	# get the relevant columns and column count
	pos <- grep(paste(us,"[0-9]+$",sep=""),colnames(rcombi))
	lpos <- length(pos)

	# create two subsets for the raw and filtered data
	rsubset <- rcombi[,pos,drop=FALSE]
	fsubset <- rbcombi[,pos,drop=FALSE]

	# sort the subsets and select the top 10 most read abundant sequences
	rsubset <- rsubset[order(rowSums(-rsubset)),,
		drop=FALSE][1:10,,drop=FALSE]
	fsubset <- fsubset[order(rowSums(-fsubset)),,
		drop=FALSE][1:10,,drop=FALSE]
	
	# remove NAs
	rsubset[is.na(rsubset)] <- 0
	fsubset[is.na(fsubset)] <- 0

	# calculate the average repeats for the
	# top 10 most abundant sequences
	ravgreps <- mean(rowSums(rsubset>0))/lpos
	rsdreps <- sd(rowSums(rsubset>0)/lpos)
	favgreps <- mean(rowSums(fsubset>0))/lpos
	fsdreps <- sd(rowSums(fsubset>0)/lpos)

	# Remove the rows without reads to avoid excidental overlap
	rsubset <- rsubset[which(rowSums(rsubset)>0),,drop=FALSE]
	fsubset <- fsubset[which(rowSums(fsubset)>0),,drop=FALSE]

	# calculate the overlap between the two sets
	overlap <- length(intersect(rownames(rsubset)[!grepl("NA",
	rownames(rsubset))],rownames(fsubset)[!grepl("NA",rownames(fsubset))]))

	# clean the sample name so it matches the samplestat names
	cus <- sub("sample.","^",us)

	# add the average number of repeats and the overlap 
	# between them to the samplestat table
	samplestat$avg_rep[grep(cus,rownames(samplestat))] <- ravgreps
	samplestat$sd_rep[grep(cus,rownames(samplestat))] <- rsdreps
	samplestat$avg_filt_rep[grep(cus,rownames(samplestat))] <- favgreps
	samplestat$sd_filt_rep[grep(cus,rownames(samplestat))] <- fsdreps
	samplestat$overlap[grep(cus,rownames(samplestat))] <- overlap 
}



###############################
# Write the sample stat table #
###############################

# copy the sample stat table
samplestat_mod <- samplestat

# move the column names to their own row and add them to the
# samplestat_mod table. This is to avoid issues with the 
# write.table function and it not correctly outputting the header.
samplestat_mod <- rbind(colnames(samplestat_mod),samplestat_mod)
rownames(samplestat_mod)[1] <- ""

# write the table
write.table(samplestat_mod, file=paste(output_name,"_samplestat.tsv",sep=""),
	quote=FALSE, sep="\t", col.names=FALSE, row.names=TRUE)

