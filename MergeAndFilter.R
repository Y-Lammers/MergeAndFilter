#!/usr/bin/env Rscript

# R script for merging two sets of OBITools output. Applying some basic
# filters and producing a file that summarizes the PCR repeats.

# USAGE: MergeAndFilter.R [arctborbryo file] [ncbi file] [base output name] 
# [min iden] [min reads] [min total reads] [min total repeats]

# Contact: youri.lammers@gmail.com
# Version: 1.2.2

# set arguments
arct_name=commandArgs(trailingOnly = TRUE)[1]
ncbi_name=commandArgs(trailingOnly = TRUE)[2]
output_name=commandArgs(trailingOnly = TRUE)[3]

# minimum identity
if (is.na(commandArgs(trailingOnly = TRUE)[4])) {
	min_iden=1
} else
	min_iden=as.numeric(commandArgs(trailingOnly = TRUE)[4])

# minimum reads
if (is.na(commandArgs(trailingOnly = TRUE)[5])) {
	min_reads=3
} else
	min_reads=as.numeric(commandArgs(trailingOnly = TRUE)[5])

# minimum total reads
if (is.na(commandArgs(trailingOnly = TRUE)[6])) {
	min_total_reads=10
} else
	min_total_reads=as.numeric(commandArgs(trailingOnly = TRUE)[6])

# minimum total repeats
if (is.na(commandArgs(trailingOnly = TRUE)[7])) {
	min_total_rep=3
} else
	min_total_rep=as.numeric(commandArgs(trailingOnly = TRUE)[7])


#################################
# Values for testing or running #
# the R script manually         #
#################################

arct_name="AOHL3_8.ali.frm.uniq.tagswap.c2.cl.arctborbryo-iden.ann.sort.tsv"
ncbi_name="AOHL3_8.ali.frm.uniq.tagswap.c2.cl.NCBI-iden.ann.sort.tsv"
output_name="test"
min_iden=1
min_reads=3
min_total_reads=10
min_total_rep=3


###################
# Read the tables #
###################

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
sample <- arct[,which(grepl("sample.",colnames(arct)))]
clean <- arct[,which(grepl("obiclean_status.",colnames(arct)))]

combi=cbind(sample,clean)


##################################
# remove sequences with both     #
# identities below the threshold #
##################################

subset = which(arct_info[,"arct_iden"]>=min_iden|
	ncbi_info[,"ncbi_iden"]>=min_iden)
combi=combi[subset,]
arct_info=arct_info[subset,]
ncbi_info=ncbi_info[subset,]
sequences=sequences[subset,,drop=FALSE]


#######################
# Detect homopolymers #
#######################

# create two new rows for the homopolymer info
sequences$homopolymer <- NA
sequences$homopolymer_type <- NA

# go through the sequences
for (i in 1:nrow(sequences)){

	# search for homopolymers that are 5bp or longer
	rep <- regexpr("(.)\\1{4,}",sequences[i,"sequence"])
	if(rep!=-1){
		# if present: add the info to the columns
		sequences[i,"homopolymer"] <- attr(rep,"match.length")
		sequences[i,"homopolymer_type"] <- substr(
			sequences[i,"sequence"],rep,rep)
	}
}


########################################
# remove both internal and sequence    #
# below the minimum sequence threshold #
########################################

# detect the sample columns
samples <- grep("sample.",colnames(combi))

# for sample in the sample list
for (i in 1:length(samples)){

	# find the matching obiclean column
	obiclean <- grep(sub("sample","obiclean_status",
		colnames(combi)[samples[i]]), colnames(combi))

	# for sequence in the sequence list
	for (j in 1:nrow(combi)){

		# if the read number is below the
		# threshold, set it to zero
		if(combi[j,samples[i]]<min_reads){
			combi[j,samples[i]] <- 0
		} else
			# if the read number is above the
			# threshold, check if it is an internal
			# sequence, if so, set it to zero
			if(combi[j,obiclean]=="i"){
				combi[j,samples[i]] <- 0
			}
	}
}	
	
# recalculate the total count
for (i in 1:nrow(combi)){
	arct_info[i,"count"] <- sum(combi[i,samples])
}


################################################
# Remove sequences below the minimum total     #
# count or below the minimum number of repeats #
################################################

# remove the low total read count sequences
subset = which(arct_info[,"count"]>=min_total_reads)
combi=combi[subset,]
arct_info=arct_info[subset,]
ncbi_info=ncbi_info[subset,]
sequences=sequences[subset,]

# count the total number of repeats
arct_info$total_rep <- NA
for (i in 1:nrow(combi)){
	arct_info[i,"total_rep"] <- sum(combi[i,samples]>0)
}

# remove the low total repeat sequences
subset = which(arct_info[,"total_rep"]>=min_total_rep)
combi=combi[subset,]
arct_info=arct_info[subset,]
ncbi_info=ncbi_info[subset,]
sequences=sequences[subset,]


##########################
# write the cleaned data #
##########################

# get the samples from the combi table, i.e. remove the
# obiclean information since this is no longer needed
clean_sample <- combi[,which(grepl("sample.",colnames(combi)))]

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

# get the unique sample names (by removing the last character
# from the sample name, which is presumed to the repeat number.)
usample <- unique(gsub(".{1}$",'',colnames(combi[samples])))

# loop through the samples
for (u in usample){

	# get the relevant columns and column count
	pos <- grep(u,colnames(combi))
	lpos <- length(pos)

	# calculate the total sum of reads for the
	# sample across all sequences. If there is only one
	# repeat, just take the sum of that column.
	if (lpos > 1){
		pprop <- colSums(combi[,pos])/sum(combi[,pos])	
	} else {
		pprop <- 1
	}
	
	# create four new columns for total read sum, repeat count, 
	# proportion of repeats and weighted proportion of repeats
	combi[[paste("totread_",u,sep="")]] <- NA
	combi[[paste("avgread_",u,sep="")]] <- NA
	combi[[paste("sdread_",u,sep="")]] <- NA
	combi[[paste("totrep_",u,sep="")]] <- NA
	combi[[paste("proprep_",u,sep="")]] <- NA
	combi[[paste("weightrep_",u,sep="")]] <- NA


	# loop through the sequences in order to calculate
	# the total sum, count and proportional repeats per sequence
	for (i in 1:nrow(combi)){

		# get the total sum, count and prop count
		totsum <- sum(combi[i,pos])
		totrep <- sum(combi[i,pos]>0)
		proprep <- totrep/lpos

		# calculate the average read count and
		# standard deviation for samples with more than
		# one repeat
		if (length(combi[i,pos]) == 1){
			avgread <- combi[i,pos]
			sdread <- 0
		} else {
			avgread <- rowMeans(combi[i,pos])
			sdread <- sd(combi[i,pos])
		}

		# add the values to the new columns
		combi[i,paste("totread_",u,sep="")] <- totsum
		combi[i,paste("avgread_",u,sep="")] <- avgread
		combi[i,paste("sdread_",u,sep="")] <- sdread
		combi[i,paste("totrep_",u,sep="")] <- totrep
		combi[i,paste("proprep_",u,sep="")] <- proprep

	}

	# calculate the mean number of repeats for the sample
	# (ignoring 0 repeats) as well as the proportion of 
	# mean repeats
	if (sum(combi[,paste("proprep_",u,sep="")]) > 0){
		mean_rep <- mean(combi[combi[,paste("totrep_",u,sep="")]>0
			,paste("totrep_",u,sep="")], na.rm=TRUE)
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
		for (i in 1:nrow(combi)){

			# recalc the proportion based on the number of
			# reads in each repeat across the PCR replicates		
			wproprep <- sum((combi[i,pos]>0)*pprop)

			# add the weighted proportion to the dataframe
			combi[i,paste("weightrep_",u,sep="")] <- wproprep						
		}

	} else {

		# reuse the regular proprep if the ratio is too low
		combi[,paste("weightrep_",u,sep="")] <- 
			combi[,paste("proprep_",u,sep="")]

	}
}


################################
# Write the summarized results #
################################


# get the positions for the summarized columns
pos <- grep("totread_|avgread_|sdread_|totrep_|proprep_|weightrep_",colnames(combi))

# get the summary subset
summary <- combi[,pos]

# get the final table
finaltable=cbind(arct_info,ncbi_info,sequences,summary)

# resort the final table based on the total read count
finaltable <- finaltable[order(-finaltable[,"count"]),]

# write the table
write.table(finaltable, file=paste(output_name,"_summary.tsv",sep=""),
	quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
