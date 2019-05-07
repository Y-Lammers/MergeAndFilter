#!/usr/bin/env Rscript

# R script for merging two sets of OBITools output. Applying some basic
# filters and producing a file that summarizes the PCR repeats.

# USAGE: MergeAndFilter.R [arctborbryo file] [ncbi file] [count file]
# [blacklist file] [base output name] [min iden] [min reads] 
# [min total reads] [min total repeats]

# note: The arctborbryo, ncbi, count and blacklist files and output name
# are required. The remaining settings will default to the values below.

# Contact: youri.lammers@gmail.com
# Version: 1.3.1

# set arguments
arct_name=commandArgs(trailingOnly = TRUE)[1]
ncbi_name=commandArgs(trailingOnly = TRUE)[2]
count_table=commandArgs(trailingOnly = TRUE)[3]
blacklist_name=commandArgs(trailingOnly = TRUE)[4]
output_name=commandArgs(trailingOnly = TRUE)[5]

# minimum identity
if (is.na(commandArgs(trailingOnly = TRUE)[6])) {
	min_iden=1
} else
	min_iden=as.numeric(commandArgs(trailingOnly = TRUE)[6])

# minimum reads
if (is.na(commandArgs(trailingOnly = TRUE)[7])) {
	min_reads=3
} else
	min_reads=as.numeric(commandArgs(trailingOnly = TRUE)[7])

# minimum total reads
if (is.na(commandArgs(trailingOnly = TRUE)[8])) {
	min_total_reads=10
} else
	min_total_reads=as.numeric(commandArgs(trailingOnly = TRUE)[8])

# minimum total repeats
if (is.na(commandArgs(trailingOnly = TRUE)[9])) {
	min_total_rep=3
} else
	min_total_rep=as.numeric(commandArgs(trailingOnly = TRUE)[9])



#################################
# Values for testing or running #
# the R script manually         #
#################################

#arct_name="AOHL3_8.ali.frm.uniq.tagswap.c2.cl.arctborbryo-iden.ann.sort.tsv"
#ncbi_name="AOHL3_8.ali.frm.uniq.tagswap.c2.cl.NCBI-iden.ann.sort.tsv"
#count_table="AOHL3_8.ali.frm.uniq.counts.tsv"
#blacklist_name="MergeAndFilter/blacklist.tsv"
#output_name="test"
#min_iden=1
#min_reads=3
#min_total_reads=10
#min_total_rep=3



############################
# Read the obitools tables #
############################

arct = read.table(arct_name,header=TRUE,sep="\t")
ncbi = read.table(ncbi_name,header=TRUE,sep="\t")

# rename the arct iden specific columns
colnames(arct)[3] <- "arct_iden"
colnames(arct)[6:9] <- c("arct_family","arct_family_name",
	"arct_genus","arct_genus_name")
colnames(arct)[(dim(arct)[2]-6):(dim(arct)[2]-1)] <- c(
	"arct_rank","arct_scientific_name","arct_species",
	"arct_species_list","arct_species_name","arct_taxid")

# rename the ncbi iden specific columns
colnames(ncbi)[3] <- "ncbi_iden"
colnames(ncbi)[6:9] <- c("ncbi_family","ncbi_family_name",
	"ncbi_genus","ncbi_genus_name")
colnames(ncbi)[(dim(ncbi)[2]-6):(dim(ncbi)[2]-1)] <- c(
	"ncbi_rank","ncbi_scientific_name","ncbi_species",
	"ncbi_species_list","ncbi_species_name","ncbi_taxid")

# sort both datasets by the ID
arct <- arct[order(arct$id),]
ncbi <- ncbi[order(ncbi$id),]

# get the arct and ncbi info columns
arct_info <- arct[,c(1,3,5,7,9,(dim(arct)[2]-6):(dim(arct)[2]-1))]
ncbi_info <- ncbi[,c(3,7,9,(dim(arct)[2]-6):(dim(arct)[2]-1))]
sequences <- arct[,dim(arct)[2],drop=FALSE]

# get the sample and obiclean info
sample <- arct[,which(grepl("sample.",colnames(arct),fixed=TRUE))]
clean <- arct[,which(grepl("obiclean_status.",colnames(arct),fixed=TRUE))]

# combine the raw sample and obiclean data
rcombi=cbind(sample,clean)



######################
# Read the blacklist #
######################

blacklist = read.table(blacklist_name,sep="\t")



###################################################
# Read the counts data and generate a blank table #
###################################################

# read the raw count data
counts = read.table(count_table,header=TRUE,sep="\t")

# extract the sample names
csample <- unique(gsub(".{1}$",'',counts[,1]))

# get the max repeat count (median of all repeat counts)
# median is used to avoid some deviating counts of misnamed samples
repeats <- c()
for (sample in csample){
	rep <- length(grep(sample,counts[,1],fixed=TRUE))
	repeats <- c(repeats,rep)
}
repeats <- median(repeats)

# create an empty dataframe for the sample stat information
samplestat <- data.frame(matrix(NA,nrow=length(csample),
	ncol=(6*(repeats+2)+2)))

# fix the row and column names
# add the rownames based on the sample names
rownames(samplestat) <- csample

# create the column names based on the number of samples
tcolnames <- c()
for (cat in c('raw','prop_raw','prop_filt','prop_noniden',
	'prop_blacklist','prop_retained')){
	for (num in 1:repeats){
		tcolnames <- c(tcolnames, paste(cat,num,sep=''))
	}
	tcolnames <- c(tcolnames,paste(cat,'_avg',sep=''))
	tcolnames <- c(tcolnames,paste(cat,'_sd',sep=''))
}
tcolnames <- c(tcolnames,'avg_rep')
tcolnames <- c(tcolnames,'avg_filt_rep')

# add the column names
colnames(samplestat) <- tcolnames



#####################################################
# Add the raw counts data to the sample stat table, #
# both as raw counts and proportion of the sample.  #
#####################################################

# Parse through the sample names
for (sample in csample){

	# get the rows with data
	pos <- grep(sample,counts$Sample.name,fixed=TRUE)

	# get the row for the sample in the samplestat table
	samplepos <- grep(sample,rownames(samplestat),fixed=TRUE)

	# get the total raw read count for the sample
	totraw <- sum(counts$Read.count[pos])

	# create an empty vector that stores the repeat
	# proportions in order to calc the standard deviation
	sampleprop <- c()	

	# get the number of repeats for the sample, cap
	# at the expected repeat count in case of a sample
	# name mixup.
	if (length(pos) > repeats){
		srep <- repeats
	} else {
		srep <- length(pos)
	}
	
	# calculate the proportion of raw reads per repeat
	# and add it to the sample stat table, as well as
	# the raw number itself.
	for (rep in 1:srep){
		
		# calc the read proportion for the repeat
		prep <-	counts$Read.count[pos[rep]]/totraw

		# add the proportion to the sampleprop vector
		sampleprop <- c(sampleprop, prep)

		# add the raw repeat counts and proportional repeat
		# counts to the sample stat table
		samplestat[samplepos,grep(paste("raw",rep,sep=""),
			colnames(samplestat),fixed=TRUE)] <-
			counts$Read.count[pos[rep]]
		samplestat[samplepos,grep(paste("prop_raw",rep,sep=""),
			colnames(samplestat),fixed=TRUE)] <- prep

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
	obiclean <- grep(sub("sample","obiclean_status",
		colnames(tcombi)[samples[s]]), colnames(tcombi),fixed=TRUE)

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
	arct_info[seq,"count"] <- sum(tcombi[seq,samples])
}



################################################
# Remove sequences below the minimum total     #
# count or below the minimum number of repeats #
################################################

# remove the low total read count sequences
subset = which(arct_info[,"count"]>=min_total_reads)
tcombi=tcombi[subset,]
arct_info=arct_info[subset,]
ncbi_info=ncbi_info[subset,]
sequences=sequences[subset,]

# count the total number of repeats
arct_info$total_rep <- NA
for (seq in 1:nrow(tcombi)){
	arct_info[seq,"total_rep"] <- sum(tcombi[seq,samples]>0)
}

# remove the low total repeat sequences
subset = which(arct_info[,"total_rep"]>=min_total_rep)
tcombi=tcombi[subset,]
arct_info=arct_info[subset,]
ncbi_info=ncbi_info[subset,]
sequences=sequences[subset,]



##################################
# remove sequences with both     #
# identities below the threshold #
##################################

# copy the technical filtered combined sequence table
icombi <- tcombi

# find the rows that match the minimum identity
subset = which(arct_info[,"arct_iden"]>=min_iden|
	ncbi_info[,"ncbi_iden"]>=min_iden)

# apply the filtering to the data
icombi=icombi[subset,]
arct_info=arct_info[subset,]
ncbi_info=ncbi_info[subset,]
sequences=sequences[subset,,drop=FALSE]



##################################
# Apply  the blacklist filtering #
##################################

# copy the identity filtered combined sequence table
bcombi <- icombi

# find the rows that match sequences in the blacklist
# and select the non matching rows
subset = which(!sequences[,"sequence"] %in% blacklist[,1])

# apply the filtering to the data
bcombi=bcombi[subset,]
arct_info=arct_info[subset,]
ncbi_info=ncbi_info[subset,]
sequences=sequences[subset,,drop=FALSE]



##########################
# write the cleaned data #
##########################

# get the samples from the final filtered combi table, i.e. remove the
# obiclean information since this is no longer needed
clean_sample <- bcombi[,which(grepl("sample.",colnames(bcombi),fixed=TRUE))]

# get the final table
finaltable=cbind(arct_info,ncbi_info,sequences,clean_sample)

# resort the final table based on the total read count
finaltable <- finaltable[order(-finaltable[,"count"]),]

write.table(finaltable, file=paste(output_name,"_filtered.tsv",sep=""),
	quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)



############################################
# get the sample repeat count, read sums   #
# and the (weighted) proportion of repeats #
############################################

# copy the final filtered table for summarization
scombi <- bcombi

# get the unique sample names (by removing the last character
# from the sample name, which is presumed to the repeat number.)
usamples <- unique(gsub(".{1}$",'',colnames(scombi[samples])))

# loop through the samples
for (us in usamples){

	# get the relevant columns and column count
	pos <- grep(us,colnames(scombi),fixed=TRUE)
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
finaltable=cbind(arct_info,ncbi_info,sequences,summary)

# re-sort the final table based on the total read count
finaltable <- finaltable[order(-finaltable[,"count"]),]

# write the table
write.table(finaltable, file=paste(output_name,"_summary.tsv",sep=""),
	quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)



#############################################
# Calculate the proportion of reads removed #
# and retained for the filtering steps      #
#############################################

# loop through the samples
for (us in usamples){

	# clean the sample name so it matches the samplestat names
	cus <- sub("sample.","",us)

	# get the row in the samplestat table
	cusr <- grep(cus,rownames(samplestat),fixed=TRUE)
	
	# get the relevant columns and column count
	pos <- grep(us,colnames(rcombi),fixed=TRUE)
	lpos <- length(pos)	

	# create empty vectors for storing the proportional data
	# in order to calculate the mean and standard deviation
	tfilt <- c()
	ifilt <- c()
	bfilt <- c()
	retain <- c()

	# get the number of repeats for the sample, cap
	# at the expected repeat count in case of a sample
	# name mixup.
	if (length(pos) > repeats){
		srep <- repeats
	} else {
		srep <- length(pos)
	}

	# loop throug the repeats for a sample
	for (rep in 1:srep){

		# get the raw number of reads for the repeat
		raw <- samplestat[cusr,rep]

		# get the proportion of technical filtered reads
		pt <- (raw-(sum(tcombi[,pos[rep]])))/raw
		
		# get the proportion of identity filtered reads
		pi <- ((raw-(sum(icombi[,pos[rep]])))/raw)-pt
	
		# get the proportion of blacklist filtered reads
		pb <- ((raw-(sum(bcombi[,pos[rep]])))/raw)-(pt+pi)

		# get the proportion of retained reads
		pr <- sum(bcombi[,pos[rep]])/raw

		# add the proprotion of technical filtered reads
		# to the samplestat table and vector
		samplestat[cusr,(((repeats+2)*2)+rep)] <- pt
		tfilt <- c(tfilt,pt)

		# add the proprotion of identity filtered reads
		# to the samplestat table and vector
		samplestat[cusr,(((repeats+2)*3)+rep)] <- pi
		ifilt <- c(ifilt,pi)

		# add the proprotion of blacklist filtered reads
		# to the samplestat table and vector
		samplestat[cusr,(((repeats+2)*4)+rep)] <- pb
		bfilt <- c(bfilt,pb)

		# add the proprotion of retained reads to the
		# samplestat table and vector
		samplestat[cusr,(((repeats+2)*5)+rep)] <- pr
		retain <- c(retain,pr)

	}

	# calculate the mean and standard deviation for the
	# proportion of technical filtered reads and add them
	# to the samplestat table
	samplestat[cusr,(((repeats+2)*2)+(repeats+1))] <-
		mean(tfilt)
	samplestat[cusr,(((repeats+2)*2)+(repeats+2))] <-
		sd(tfilt)

	# the same for the identity filtered reads
	samplestat[cusr,(((repeats+2)*3)+(repeats+1))] <-
		mean(ifilt)
	samplestat[cusr,(((repeats+2)*3)+(repeats+2))] <-
		sd(ifilt)

	# and the blacklist filtered reads
	samplestat[cusr,(((repeats+2)*4)+(repeats+1))] <-
		mean(bfilt)
	samplestat[cusr,(((repeats+2)*4)+(repeats+2))] <-
		sd(bfilt)

	# and the proportion of retained reads
	samplestat[cusr,(((repeats+2)*5)+(repeats+1))] <-
		mean(retain)
	samplestat[cusr,(((repeats+2)*5)+(repeats+2))] <-
		sd(retain)

}



#######################################
# get the average repeats for the top #
# 10 filtered and unfiltered taxa     #
#######################################

# loop through the samples
for (us in usamples){

	# clean the sample name so it matches the samplestat names
	cus <- sub("sample.","",us)
	
	# get the relevant columns and column count
	pos <- grep(us,colnames(rcombi),fixed=TRUE)
	lpos <- length(pos)

	# create two subsets for the raw and filtered data
	rsubset <- rcombi[,pos]
	fsubset <- bcombi[,pos]
	
	# sort the subset columns based on the readcount
	# across the rows (note: just order by count
	# if there is only one column).
	if (lpos == 1){

		rsubset <- sort(rsubset, decreasing=TRUE)
		fsubset <- sort(fsubset, decreasing=TRUE)

		# calculate the average repeats for the
		# top 10 most abundant sequences
		ravgreps <- mean(rsubset[1:10]>0)
		favgreps <- mean(fsubset[1:10]>0)

	} else {

		rsubset <- rsubset[order(rowSums(-rsubset)),]
		fsubset <- fsubset[order(rowSums(-fsubset)),]

		# calculate the average repeats for the
		# top 10 most abundant sequences
		ravgreps <- mean(rowSums(rsubset[1:10,]>0))
		favgreps <- mean(rowSums(fsubset[1:10,]>0))

	}

	# add the average repeats to the samplestat table
	samplestat$avg_rep[grep(cus,rownames(samplestat),
		fixed=TRUE)] <- ravgreps
	samplestat$avg_filt_rep[grep(cus,rownames(samplestat),
		fixed=TRUE)] <- favgreps

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

