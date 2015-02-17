#list all fastqc files in the directory and read in sample names
filelist = list.files(path = ".", pattern="fastqc_data")
samples = read.table("sample_names.txt", header=F, sep="\t",stringsAsFactors=F)

#create table for the values of interest
sample  			  = rep("N",2*length(samples[[1]]))
encoding 			  = rep("N",2*length(samples[[1]]))
GC_content 			  = rep("N",2*length(samples[[1]]))
overrepresented_seq   = rep("",2*length(samples[[1]]))
overrepresented_count = rep("",2*length(samples[[1]]))
overrepresented_prct  = rep("",2*length(samples[[1]]))
overrepresented_src   = rep("",2*length(samples[[1]]))
polyT_read			  = rep("N",2*length(samples[[1]]))
polyT_count			  = rep("",2*length(samples[[1]]))
polyT_prct			  = rep("",2*length(samples[[1]]))

#loop through the samples to pull the values of interest
for (g in seq(1,length(samples[[1]])))
{
	samp = samples[[1]][g]

	ind = grep("^N$",sample)[1]

	file_1 			= filelist[grep(samp,filelist)]
	qc_1   			= read.table(file_1,header=F, sep="\t", stringsAsFactors=F, blank.lines.skip = T, fill =T)
	sample[ind] 	= samp
	encoding[ind] 	= qc_1[grep("Encoding",qc_1[,1]), 2]
	GC_content[ind] = qc_1[grep("%GC",qc_1[,1]), 2]

	#iterate through all the Overrepresented Sequences
	for (i in seq((grep("Overrepresented",qc_1[,1])+1),(grep("Kmer",qc_1[,1])-3),2))
	{
		overrepresented_seq[ind] 	= paste(overrepresented_seq[ind],qc_1[i,1],",",sep="") #sequence
		overrepresented_count[ind]	= paste(overrepresented_count[ind], qc_1[i,2],",",sep="") #counts
		overrepresented_prct[ind]	= paste(overrepresented_prct[ind],qc_1[i,3],",",sep="") #percentage
		overrepresented_src[ind]	= paste(overrepresented_src[ind], qc_1[i+1,1],",",sep="") #source
	}

	#look for the presence of any specific adapter sequences that can contaminate results
	if (length(grep("AAAAAAATTTTTTCCCCGGGGGGG",qc_1)) > 0)
	{
		polyT_read[ind]  = "Y"
		polyT_count[ind] = qc_1[grep("AAAAAAATTTTTTCCCCGGGGGGG",qc_1[,1]),2]
		polyT_prct[ind]	 = qc_1[grep("AAAAAAATTTTTTCCCCGGGGGGG",qc_1[,1]),3]
	}

}

#combine variables; write out to text file
out <- cbind(sample,encoding,GC_content,overrepresented_seq,overrepresented_count,overrepresented_prct,overrepresented_src,polyT_read,polyT_count,polyT_prct)
write.table(out,"fastqc_summary.txt", row.names=F, sep="\t")
