read_fasta_aln<-function(f){
	con<-file(f)
	r<-readLines(con)
	close(con)

	name_line_numbers<-grep(">",r)
	starts<-name_line_numbers+1
	ends<-c(name_line_numbers-1,length(r))[-1]
	names<-r[name_line_numbers]	

	m<-mapply(function(x,y){paste(r[x:y],collapse="")},starts,ends)
	print(m)
	print(length(m))
	return(cbind(names,m))
}

#f<-"/home/koshka/work/DAVID/homologs_contact_map/SH2_pdb_70.fa"

f<-commandArgs(TRUE)[1]

aln<-read_fasta_aln(f)

aln_length<-nchar(aln[1,2])

gap_percent<-sapply(c(1:aln_length),function(x){
	sum(sapply(aln[,2],function(y,i=x){
		if (substr(y,i,i)=='-') return(1)
		else return(0)
	}))})

print(gap_percent)

print(aln_length)
print(length(gap_percent))
