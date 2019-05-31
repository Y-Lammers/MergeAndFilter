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
# Version: 1.3.4

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

#obi1file="ECOGEN-ECOG-1/ECOG-1-2_tag.ali.frm.uniq.index_swap.c2.cl.arctborbryo-iden.ann.sort.tsv"
#obi2file="ECOGEN-ECOG-1/ECOG-1-2_tag.ali.frm.uniq.index_swap.c2.cl.NCBI-iden.ann.sort.tsv"
#count_table="ECOGEN-ECOG-1/ECOG-1-2_tag.ali.frm.uniq.counts.tsv"
#synthetic_blacklist_name="MergeAndFilter/synthetic_blacklist.tsv"
#region_blacklist_name="MergeAndFilter/N-Norway_blacklist.tsv"
#output_name="test"
#obi1name="arctborbryo"
#obi2name="ncbi"
#min_iden=1
#min_reads=3
#min_total_reads=10
#min_total_rep=3



############################
# Read the obitools tables #
############################

obi1 = read.table(obi1file,header=TRUE,sep="\t")
obi2 = read.table(obi2file,header=TRUE,sep="\t")

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

# sort both datasets by the ID
obi1 <- obi1[order(obi1$id),]
obi2 <- obi2[order(obi2$id),]

# get the information columns
obi1_info <- obi1[,c(1,3,5,7,9,(dim(obi1)[2]-6):(dim(obi1)[2]-1))]
obi2_info <- obi2[,c(3,7,9,(dim(obi2)[2]-6):(dim(obi2)[2]-1))]
sequences <- obi1[,dim(obi1)[2],drop=FALSE]

# get the sample and obiclean info
sample <- obi1[,which(grepl("sample.",colnames(obi1),fixed=TRUE))]
clean <- obi1[,which(grepl("obiclean_status.",colnames(obi1),fixed=TRUE))]

# combine the raw sample and obiclean data
rcombi=cbind(sample,clean)



###################################################
# Read the counts data and generate a blank table #
###################################################

# read the raw count data
counts = read.table(count_table,header=TRUE,sep="\t")

# reformat the sample names in the counts table so that they
# match with the OBITools input
counts[,1] <- gsub("[:-]",".",counts[,1])

# extract the sample names
csample <- unique(gsub(".{1}$",'',counts[,1]))

# get the maximum number of repeats in the library
repeats <- c()
for (sample in csample){
	rep <- length(grep(paste("^",sample,sep=""),counts[,1]))
	repeats <- c(repeats,rep)
}
repeats <- max(repeats)

# create an empty dataframe for the sample stat information
samplestat <- data.frame(matrix(NA,nrow=length(csample),
	ncol=(7*(repeats+2)+3)))

# fix the row and column names
# add the rownames based on the sample names
rownames(samplestat) <- csample

# create the column names based on the number of samples
tcolnames <- c()
for (cat in c('raw','prop_raw','prop_filt','prop_noniden',
	'prop_synthetic_blacklist','prop_region_blacklist','prop_retained')){
	for (num in 1:repeats){
		tcolnames <- c(tcolnames, paste(cat,num,sep=''))
	}
	tcolnames <- c(tcolnames,paste(cat,'_avg',sep=''))
	tcolnames <- c(tcolnames,paste(cat,'_sd',sep=''))
}
tcolnames <- c(tcolnames,'avg_rep')
tcolnames <- c(tcolnames,'avg_filt_rep')
tcolnames <- c(tcolnames,'overlap')

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

	# for sequence in the sequence list
	for (seq in 1:nrow(tcombi)){

		# if the read number is below the
		# threshold, set it to zero
		if(tcombi[seq,samples[s]]<min_reads){
			tcombi[seq,samples[s]] <- 0
		} else
			# if the read number is above the
			# threshold, check if it is an internal
			# sequence, if so, set it to zero
			if(tcombi[seq,obiclean]=="i"){
				tcombi[seq,samples[s]] <- 0
			}
	}
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
tcombi=tcombi[subset,]
obi1_info=obi1_info[subset,]
obi2_info=obi2_info[subset,]
sequences=sequences[subset,]

# count the total number of repeats
obi1_info$total_rep <- NA
for (seq in 1:nrow(tcombi)){
	obi1_info[seq,"total_rep"] <- sum(tcombi[seq,samples]>0)
}

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
synthetic_blacklist = read.table(synthetic_blacklist_name,sep="\t")

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
region_blacklist = read.table(region_blacklist_name,sep="\t")

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
clean_sample <- rbcombi[,which(grepl("sample.",colnames(rbcombi),fixed=TRUE))]

# get the final table
finaltable=cbind(obi1_info,obi2_info,sequences,clean_sample)

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
		ssum <- colSums(scombi[,pos])
		pprop <- colSums(scombi[,pos])/sum(scombi[,pos])
	} else {
		ssum <- sum(scombi[,pos])
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


	# loop through the sequences in order to calculate
	# the total sum, count and proportional repeats per sequence
	for (seq in 1:nrow(scombi)){

		# get the total sum, count and prop count
		totsum <- sum(scombi[seq,pos])
		totrep <- sum(scombi[seq,pos]>0)
		proprep <- totrep/lpos

		# calculate the proportion of reads assigned to a taxa
		# per repeat. Use the proportion to work out the average
		# and standard deviation for each sample.
		if (length(scombi[seq,pos]) == 1){
			avgpread <- scombi[seq,pos]/ssum
			sdpread <- 0
		} else {
			avgpread <- rowMeans(scombi[seq,pos]/ssum,na.rm=TRUE)
			sdpread <- sd(scombi[seq,pos]/ssum,na.rm=TRUE)
		}

		# add the values to the new columns
		scombi[seq,paste("totread_",us,sep="")] <- totsum
		scombi[seq,paste("avgpread_",us,sep="")] <- avgpread
		scombi[seq,paste("sdpread_",us,sep="")] <- sdpread
		scombi[seq,paste("totrep_",us,sep="")] <- totrep
		scombi[seq,paste("proprep_",us,sep="")] <- proprep

	}

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

	# if the proportion of repeats is higher than 0.33, 
	# calculate the weighted proportion of repeats per sequence.
	# if proportion is lower, use the regular proportion.
	if (rmean_rep >= 0.33){

		# loop through the sequences again, this time calculating
		# the weighted proportion of repeats if needed
		for (seq in 1:nrow(scombi)){

			# recalc the proportion based on the number of
			# reads in each repeat across the PCR replicates		
			wproprep <- sum((scombi[seq,pos]>0)*pprop)

			# add the weighted proportion to the dataframe
			scombi[seq,paste("weightrep_",us,sep="")] <- wproprep						
		}

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

# loop through the samples based on the counts table
for (us in csample){

	# get the row in the samplestat table
	samplepos <- grep(paste("^",us,"$",sep=""),rownames(samplestat))
	
	# get the relevant columns and column count
	#rep <- grep(us,colnames(rcombi),fixed=TRUE)

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
	for (rep in repnames) {

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

}



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

		rsubset <- rsubset[order(rsubset[,1], decreasing=TRUE),]
		fsubset <- fsubset[order(fsubset[,1], decreasing=TRUE),]

		# calculate the average repeats for the
		# top 10 most abundant sequences
		ravgreps <- mean(rsubset[1:10,1]>0)
		favgreps <- mean(fsubset[1:10,1]>0)

		# calculate the overlap between the two sets
		overlap <- length(intersect(rownames(rsubset[1:10,]),
		rownames(fsubset[1:10,])))


	} else {

		# create two subsets for the raw and filtered data
		rsubset <- rcombi[,pos]
		fsubset <- rbcombi[,pos]

		rsubset <- rsubset[order(rowSums(-rsubset)),]
		fsubset <- fsubset[order(rowSums(-fsubset)),]

		# calculate the average repeats for the
		# top 10 most abundant sequences
		ravgreps <- mean(rowSums(rsubset[1:10,]>0))
		favgreps <- mean(rowSums(fsubset[1:10,]>0))

		# calculate the overlap between the two sets
		overlap <- length(intersect(rownames(rsubset[1:10,]),
		rownames(fsubset[1:10,])))

	}

	# clean the sample name so it matches the samplestat names
	cus <- sub("sample.","^",us)

	# add the average number of repeats and the overlap 
	# between them to the samplestat table
	samplestat$avg_rep[grep(cus,rownames(samplestat))] <- ravgreps
	samplestat$avg_filt_rep[grep(cus,rownames(samplestat))] <- favgreps
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

