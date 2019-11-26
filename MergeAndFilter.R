#!/usr/bin/env Rscript

# R script for merging two sets of OBITools output TSV files. The script appies
# a number of basic filtering criteria and outputs three tables. One table
# with the filtered results, a summary table for the filterd results and a 
# table that contains the number and proportion of reads for each sample
# in the datasets during the filtering.

# USAGE: MergeAndFilter.R [obitools output 1] [obitools output 2] 
# [count file] [synthetic blacklist file] [region blacklist file] 
# [base output name] [obitools table 1 name] [obitools table 2 name]
# [min iden] [min reads] [min total reads] [min total repeats]

# note: The two obitools tables, count, blacklist files and output name
# are required. The remaining settings will default to the values below.

# Contact: youri.lammers@gmail.com
# Version: 1.5.5

# set arguments
obi1file=commandArgs(trailingOnly = TRUE)[1]
obi2file=commandArgs(trailingOnly = TRUE)[2]
count_table=commandArgs(trailingOnly = TRUE)[3]
synthetic_blacklist_name=commandArgs(trailingOnly = TRUE)[4]
region_blacklist_name=commandArgs(trailingOnly = TRUE)[5]
output_name=commandArgs(trailingOnly = TRUE)[6]

# obitools table 1 name
if (is.na(commandArgs(trailingOnly = TRUE)[7])) {
	obi1name="obi1"
} else
	obi1name=commandArgs(trailingOnly = TRUE)[7]

# obitools table 2 name
if (is.na(commandArgs(trailingOnly = TRUE)[8])) {
	obi2name="obi2"
} else
	obi2name=commandArgs(trailingOnly = TRUE)[8]

# minimum identity
if (is.na(commandArgs(trailingOnly = TRUE)[9])) {
	min_iden=1
} else
	min_iden=as.numeric(commandArgs(trailingOnly = TRUE)[9])

# minimum reads
if (is.na(commandArgs(trailingOnly = TRUE)[10])) {
	min_reads=3
} else
	min_reads=as.numeric(commandArgs(trailingOnly = TRUE)[10])

# minimum total reads
if (is.na(commandArgs(trailingOnly = TRUE)[11])) {
	min_total_reads=10
} else
	min_total_reads=as.numeric(commandArgs(trailingOnly = TRUE)[11])

# minimum total repeats
if (is.na(commandArgs(trailingOnly = TRUE)[12])) {
	min_total_rep=3
} else
	min_total_rep=as.numeric(commandArgs(trailingOnly = TRUE)[12])



#################################
# Values for testing or running #
# the R script manually         #
#################################

#obi1file="path to the first obitools file"
#obi2file="path to the second obitools file"
#count_table="path to the sample counts table"
#output_name="base output name"
#synthetic_blacklist_name="path to the synthetic blacklist"
#region_blacklist_name="path to the regional blacklist"
#obi1name="dataset 1"
#obi2name="dataset 2"
#min_iden=1
#min_reads=3
#min_total_reads=10
#min_total_rep=3


############################
# Read the obitools tables #
############################

# Load the first OBITools table
obi1 = read.table(obi1file,header=TRUE,sep="\t")

# rename the columns for the first obitools file
colnames(obi1)[3] <- paste(obi1name,"_iden",sep="")
colnames(obi1)[6:9] <- c(paste(obi1name,"_family",sep=""),
	paste(obi1name,"_family_name",sep=""),
	paste(obi1name,"_genus",sep=""),
	paste(obi1name,"_genus_name",sep=""))
colnames(obi1)[(dim(obi1)[2]-6):(dim(obi1)[2]-1)] <- c(
	paste(obi1name,"_rank",sep=""),
	paste(obi1name,"_scientific_name",sep=""),
	paste(obi1name,"_species",sep=""),
	paste(obi1name,"_species_list",sep=""),
	paste(obi1name,"_species_name",sep=""),
	paste(obi1name,"_taxid",sep=""))

# order by ID
obi1 <- obi1[order(obi1$id),]

# get the information columns and the sequences
obi1_info <- obi1[,c(1,3,5,7,9,(dim(obi1)[2]-6):(dim(obi1)[2]-1))]
sequences <- obi1[,dim(obi1)[2],drop=FALSE]

# clear the memory
rm(obi1)

# read the second obitools table
obi2 = read.table(obi2file,header=TRUE,sep="\t")

# rename the columns for the second obitools file
colnames(obi2)[3] <- paste(obi2name,"_iden",sep="")
colnames(obi2)[6:9] <- c(paste(obi2name,"_family",sep=""),
	paste(obi2name,"_family_name",sep=""),
	paste(obi2name,"_genus",sep=""),
	paste(obi2name,"_genus_name",sep=""))
colnames(obi2)[(dim(obi2)[2]-6):(dim(obi2)[2]-1)] <- c(
	paste(obi2name,"_rank",sep=""),
	paste(obi2name,"_scientific_name",sep=""),
	paste(obi2name,"_species",sep=""),
	paste(obi2name,"_species_list",sep=""),
	paste(obi2name,"_species_name",sep=""),
	paste(obi2name,"_taxid",sep=""))

# sort the second dataset by the ID
obi2 <- obi2[order(obi2$id),]

# get the info
obi2_info <- obi2[,c(3,7,9,(dim(obi2)[2]-6):(dim(obi2)[2]-1))]

# get the sample and obiclean info
sample <- obi2[,which(grepl("sample.",colnames(obi2),fixed=TRUE))]
clean <- obi2[,which(grepl("obiclean_status.",colnames(obi2),fixed=TRUE))]

# combine the raw sample and obiclean data
rcombi=cbind(sample,clean)

# clear the memory
rm(obi2)


###################################################
# Read the counts data and generate a blank table #
###################################################

# read the raw count data
counts = read.table(count_table,header=TRUE,sep="\t")

# reformat the sample names in the counts table so that they
# match with the OBITools input
counts[,1] <- gsub("[][():+-]",".",counts[,1])

# extract the sample names
csample <- unique(gsub(".{1}$",'',counts[,1]))

# get the maximum number of repeats in the library
repeats <- max(as.numeric(gsub("^.*rpt","",counts[,1])))

# create an empty dataframe for the sample stat information
samplestat <- data.frame(matrix(NA,nrow=length(csample),
	ncol=(7*(repeats+2)+10)))

# fix the row and column names
# add the rownames based on the sample names
rownames(samplestat) <- csample

# create the column names for the repeats based on 
# the number of samples
tcolnames <- c()
for (cat in c('raw','prop_raw','prop_filt','prop_noniden',
	'prop_synthetic_blacklist','prop_region_blacklist','prop_retained')){
	for (num in 1:repeats){
		tcolnames <- c(tcolnames, paste(cat,num,sep=''))
	}
	tcolnames <- c(tcolnames,paste(cat,'_avg',sep=''))
	tcolnames <- c(tcolnames,paste(cat,'_sd',sep=''))
}

# add the additional columns for the various sample averages
tcolnames <- c(tcolnames,'avg_rep')
tcolnames <- c(tcolnames,'sd_rep')
tcolnames <- c(tcolnames,'avg_filt_rep')
tcolnames <- c(tcolnames,'sd_filt_rep')
tcolnames <- c(tcolnames,'overlap')
tcolnames <- c(tcolnames,'avg_seq_length')
tcolnames <- c(tcolnames,'sd_seq_length')
tcolnames <- c(tcolnames,'avg_single_seq_length')
tcolnames <- c(tcolnames,'sd_single_seq_length')
tcolnames <- c(tcolnames,'single_seq_count')


# add the column names
colnames(samplestat) <- tcolnames



#####################################################
# Add the raw counts data to the sample stat table, #
# both as raw counts and proportion of the sample.  #
#####################################################

# Parse through the sample names
for (sample in csample){

	# get the rows with data
	pos <- grep(paste("^",sample,sep=""),counts$Sample.name)

	# get the row for the sample in the samplestat table
	samplepos <- grep(paste("^",sample,sep=""),rownames(samplestat))

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
		repchar <- substr(repname, nchar(repname), nchar(repname))

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
for (seq in 1:nrow(sequences)){

	# search for homopolymers that are 5bp or longer
	rep <- regexpr("(.)\\1{4,}",sequences[seq,"sequence"])
	if(rep!=-1){
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
for (s in 1:length(samples)){

	# find the matching obiclean column
	obiclean <- grep(paste("^",sub("sample","obiclean_status",
		colnames(tcombi)[samples[s]]),"$",sep=""), colnames(tcombi))

	# remove low sequence occurences and internal sequences
	tcombi[which(tcombi[,samples[s]]<min_reads | 
		tcombi[,obiclean]=="i"), samples[s]] <- 0

}	
	
# recalculate the total count
for (seq in 1:nrow(tcombi)){
	obi1_info[seq,"count"] <- sum(tcombi[seq,samples])
}



################################################
# Remove sequences below the minimum total     #
# count or below the minimum number of repeats #
################################################

# remove the low total read count sequences
subset = which(obi1_info[,"count"]>=min_total_reads)
rcombi=rcombi[subset,samples]
tcombi=tcombi[subset,samples]
obi1_info=obi1_info[subset,]
obi2_info=obi2_info[subset,]
sequences=sequences[subset,]

# count the total number of repeats
obi1_info$total_rep <- NA
obi1_info[,"total_rep"] <- rowSums(tcombi[,samples]>0)

# remove the low total repeat sequences
subset = which(obi1_info[,"total_rep"]>=min_total_rep)
tcombi=tcombi[subset,]
obi1_info=obi1_info[subset,]
obi2_info=obi2_info[subset,]
sequences=sequences[subset,]



##################################
# remove sequences with both     #
# identities below the threshold #
##################################

# copy the technical filtered combined sequence table
icombi <- tcombi

# find the rows that match the minimum identity
subset = which(obi1_info[,paste(obi1name,"_iden",sep="")]>=min_iden|
	obi2_info[,paste(obi2name,"_iden",sep="")]>=min_iden)

# apply the filtering to the data
icombi=icombi[subset,]
obi1_info=obi1_info[subset,]
obi2_info=obi2_info[subset,]
sequences=sequences[subset,,drop=FALSE]



###########################################
# Apply the synthetic blacklist filtering #
###########################################

# Read the blacklist #
synthetic_blacklist = read.table(synthetic_blacklist_name,sep="\t",fill=TRUE)

# copy the identity filtered combined sequence table
sbcombi <- icombi

# find the rows that match sequences in the blacklist
# and select the non matching rows
subset = which(!sequences[,"sequence"] %in% synthetic_blacklist[,1])

# apply the filtering to the data
sbcombi=sbcombi[subset,]
obi1_info=obi1_info[subset,]
obi2_info=obi2_info[subset,]
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
obi1_info=obi1_info[subset,]
obi2_info=obi2_info[subset,]
sequences=sequences[subset,,drop=FALSE]



##########################
# write the cleaned data #
##########################

# get the samples from the final filtered combi table, i.e. remove the
# obiclean information since this is no longer needed

# get the final table
finaltable=cbind(obi1_info,obi2_info,sequences,rbcombi)

# resort the final table based on the total read count
finaltable <- finaltable[order(-finaltable[,"count"]),]

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
usamples <- unique(gsub(".{1}$",'',colnames(scombi[samples])))

# loop through the samples
for (us in usamples){

	# get the relevant columns and column count
	pos <- grep(paste("^",us,sep=""),colnames(scombi))
	lpos <- length(pos)

	# Calculate the sum of reads for each repeat, as well as the proportion
	# of reads for each repeat. If the number of repeats is 1, set the
	# proportion to 1.
	if (lpos > 1){
		pprop <- colSums(scombi[,pos])/sum(scombi[,pos])
	} else {
		pprop <- 1
	}

	# create four new columns for total read sum, repeat count, 
	# proportion of repeats and weighted proportion of repeats
	scombi[[paste("totread_",us,sep="")]] <- NA
	scombi[[paste("avgpread_",us,sep="")]] <- NA
	scombi[[paste("sdpread_",us,sep="")]] <- NA
	scombi[[paste("totrep_",us,sep="")]] <- NA
	scombi[[paste("proprep_",us,sep="")]] <- NA
	scombi[[paste("weightrep_",us,sep="")]] <- NA


	# get the read proportions withing a sample
	temp <- as.data.frame(lapply(scombi[,pos],function(x) prop.table(x)))

	# calculate the totals, means and sd of the proportions, 
	# ignore one rep samples
	if (lpos == 1){

		# get the sample sums, reps are just straight copies
		totsumv <- scombi[,pos]
		totrepv <- scombi[,pos]
		totrepv[which(totrepv>0)] <- 1
		proprepv <- totrepv

		# and the mean plus standard deviation		
		avgpreadv <- as.double(temp)
		sdpreadv <- replicate(length(avgpreadv), 0)

		# get the number of barcodes
		barcodecount <- sum(scombi[,pos]>0)

	} else {

		# get the sample sums, replicates and the proportion of replicates
		totsumv <- rowSums(scombi[,pos])
		totrepv <- rowSums(scombi[,pos]>0)
		proprepv <- totrepv/lpos

		# and the mean plus standard deviation		
		avgpreadv <- rowMeans(temp,na.rm=TRUE)
		sdpreadv <- apply(temp,1,sd,na.rm=TRUE)

		# get the number of barcodes
		barcodecount <- sum(rowSums(scombi[,pos]>0)>0)

	}

	# store the results
	scombi[,paste("totread_",us,sep="")] <- totsumv
	scombi[,paste("avgpread_",us,sep="")] <- avgpreadv
	scombi[,paste("sdpread_",us,sep="")] <- sdpreadv
	scombi[,paste("totrep_",us,sep="")] <- totrepv
	scombi[,paste("proprep_",us,sep="")] <- proprepv


	# calculate the mean number of repeats for the sample
	# (ignoring 0 repeats) as well as the proportion of 
	# mean repeats
	if (sum(scombi[,paste("proprep_",us,sep="")]) > 0){
		mean_rep <- mean(scombi[scombi[,paste("totrep_",us,sep="")]>0
			,paste("totrep_",us,sep="")], na.rm=TRUE)
		rmean_rep <- mean_rep / lpos
	} else {
		mean_rep <- 0
		rmean_rep <- 0
		
	}

	# if the proportion of repeats is higher than 0.33 AND there are
	# atleast 10 different barcodes in the sample:
	# calculate the weighted proportion of repeats per sequence.
	# if proportion is lower, use the regular proportion.
	if ((rmean_rep >= 0.33) & (barcodecount>=10)){

		# calculate the weighted proportion of repeats
		wproprepv <- rowSums(t(t(scombi[,pos]>0)*pprop))

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
finaltable=cbind(obi1_info,obi2_info,sequences,summary)

# re-sort the final table based on the total read count
finaltable <- finaltable[order(-finaltable[,"count"]),]

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
for (us in csample){

	# get the row in the samplestat table
	samplepos <- grep(paste("^",us,"$",sep=""),rownames(samplestat))

	# get the sample names
	repnames <- counts[grep(paste("^",us,sep=""),counts[,1]),1]

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
	for (rep in repnames){

		# get the repeat number so it can be stored in
		# the proper column
		repname <- paste("sample.",rep,sep="")
		repchar <- substr(repname, nchar(repname), nchar(repname))
		repnum <- as.numeric(repchar)

		# check if the repeat still exists in the OBITools
		# output, if yes, compute the proportions, if not
		# all data is filtered during OBITools itself and the
		# proportions default to a 100% removed.
		if (repname %in% colnames(rcombi)) {

			# get the repeat location in the rcombi table
			repc <- grep(repname,colnames(rcombi))

			# get the raw number of reads for the repeat
			raw <- samplestat[samplepos,repnum]

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
	pos <- grep(paste("^sample.",us,sep=""),colnames(rbcombi))

	# create the vectors for the sequence length information
	replength <- c()
	singlelength <- c()

	# loop through the sequences, and get the sequence lengths
	# for each barcode
	for (seq in 1:nrow(rbcombi)){

		# get the number of repeats for the sample
		repcount <- sum(rbcombi[seq,pos]>=1)

		# get the seq info
		if (repcount >= 1){

			# add multiple seq lengths, based on reps
			replength <- c(replength, replicate(repcount,nchar(
				as.character(sequences[seq,"sequence"]))))

			# add a single seq length per sample
			singlelength <- c(singlelength, nchar(as.character(
				sequences[seq,"sequence"])))

		} 
	}

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

	# and the number of unique barcodes in the sample
	samplestat[samplepos,grep("single_seq_count",colnames(samplestat),
		fixed=TRUE)] <- length(singlelength)

}

# preparte the sequence length output
lengthpile_mod <- lengthpile
lengthpile_mod <- rbind(colnames(lengthpile_mod),lengthpile_mod)
rownames(lengthpile_mod)[1] <- ''

# write the sequence length table
write.table(lengthpile_mod, file=paste(output_name,"_lengths.tsv",sep=""),
        quote=FALSE, sep="\t", col.names=FALSE, row.names=TRUE)



#######################################
# get the average repeats for the top #
# 10 filtered and unfiltered taxa     #
#######################################

# loop through the samples
for (us in usamples){

	# get the relevant columns and column count
	pos <- grep(us,colnames(rcombi),fixed=TRUE)
	lpos <- length(pos)

	# sort the subset columns based on the readcount
	# across the rows (note: just order by count
	# if there is only one column).
	if (lpos == 1){

		# create two subsets for the raw and filtered data
		# in order to preserve the dimensions, add an
		# additional obiclean table that will not be used
		rsubset <- rcombi[,c(pos,(dim(rcombi)[2]-1))]
		fsubset <- rbcombi[,c(pos,(dim(rbcombi)[2]-1))]

		# sort and select the top 10
		rsubset <- rsubset[order(rsubset[,1], decreasing=TRUE),][1:10,]
		fsubset <- fsubset[order(fsubset[,1], decreasing=TRUE),][1:10,]

		# remove NAs
		rsubset[is.na(rsubset)] <- 0
		fsubset[is.na(fsubset)] <- 0

		# calculate the average repeats for the
		# top 10 most abundant sequences
		ravgreps <- mean(rsubset[,1]>0)
		rsdreps <- sd(rsubset[,1]>0)
		favgreps <- mean(fsubset[,1]>0)
		fsdreps <- sd(fsubset[,1]>0)


	} else {

		# create two subsets for the raw and filtered data
		rsubset <- rcombi[,pos]
		fsubset <- rbcombi[,pos]

		# sort and select the top 10
		rsubset <- rsubset[order(rowSums(-rsubset)),][1:10,]
		fsubset <- fsubset[order(rowSums(-fsubset)),][1:10,]
	
		# remove NAs
		rsubset[is.na(rsubset)] <- 0
		fsubset[is.na(fsubset)] <- 0

		# calculate the average repeats for the
		# top 10 most abundant sequences
		ravgreps <- mean(rowSums(rsubset>0))/lpos
		rsdreps <- sd(rowSums(rsubset>0)/lpos)
		favgreps <- mean(rowSums(fsubset>0))/lpos
		fsdreps <- sd(rowSums(fsubset>0)/lpos)

	}


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
rownames(samplestat_mod)[1] <- ''

# write the table
write.table(samplestat_mod, file=paste(output_name,"_samplestat.tsv",sep=""),
	quote=FALSE, sep="\t", col.names=FALSE, row.names=TRUE)

