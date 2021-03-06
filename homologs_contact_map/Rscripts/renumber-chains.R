library('hash')
letter_code<-hash(keys=c("ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"),values=c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"))
space_numbers<-c(6,5,1,4,1,3,1,1,4,1,3,8,8,8,6,6,6,4,2,2)
space_classes<-c("character","numeric",rep('character',6),"numeric",rep('character',2),rep('numeric',5),rep('character',4))
sum<-0;
space_numbers_starts<-array(dim=length(space_numbers)-1)
space_numbers_ends<-array(dim=length(space_numbers)-1)
k<-1
for(i in space_numbers){
	space_numbers_ends[k]<-sum+i
	if (k==1) space_numbers_starts[k]<-1
	else space_numbers_starts[k]<-space_numbers_ends[k-1]+1
	sum<-sum+i
	k<-k+1
}

write_pdb<-function(x,f){
	unlink(f)
	for(j in 1:nrow(x)){
		ss<-''
		l<-0
		for(col in 1:ncol(x)){
			s<-as.character(x[j,col])
			g<-space_numbers[col]-nchar(s)
			s2<-paste(rep(' ',g),collapse="")
			ss<-paste(ss,s2,s,sep='')
		}
		write(ss,f,append=T)				
	}
}

read_atom_pdb<-function(r){
	
	rr<-array(dim=length(r))
	for(k in 1:length(r)){
		x<-substring(r[k],space_numbers_starts,space_numbers_ends)
		y<-paste(x,collapse=";")
		rr[k]<-y
	}
	
	con<-textConnection(rr)
	t<-read.csv(con,sep=";",header=F,stringsAsFactors=F,colClasses=space_classes)
	return(t)
}

restrict_domain<-function(f_name,dir,dir2,ch,start,end){

	f<-paste(dir,f_name,".pdb",sep="")
	
	if(!file.exists(f)) return('')
	
	print(paste(f_name,ch,start,end,sep=" "))

	con <- file(f)
	on.exit(close(con))

	r<-readLines(con)
	e<-length(r)
	re<-grep("^ENDMDL",r)
	if (length(re)>0) e<-re[1]
	rr<-r[1:e]
	rr<-rr[grepl("^ATOM",rr)]

	t<-read_atom_pdb(rr)

	resi<-t[t[,8]==ch,]
	resi_sh2<-resi[(resi[,9]<=end)&(resi[,9]>=start),]
	if (nrow(resi_sh2)==0) return()
	resi_old<-sort(unique(resi_sh2[,9]))
	k<-1
	seq_a<-sapply(resi_old,function(x){
		resi_sh2[resi_sh2[,9]==x,9]<<-k
		k<<-k+1
		y<-resi_sh2[resi_sh2[,9]==(k-1),6][1]
		return(letter_code[[y]])})

	seq=paste(seq_a,collapse="")
	if(end-start+1!=nchar(seq)) {
		warning(paste("!not equal length", nchar(seq)-(end-start+1), pdb, chain,start,end))
		return()
	}
	if (nrow(resi_sh2[resi_sh2[,2]<=0,])>0) warning(paste("atom<=0",f_name))
	resi_sh2[,2]<-c(1:nrow(resi_sh2))
	resi_sh2<-resi_sh2[,1:16]
	f_i<-paste(dir2,f_name,ch,"_",start,"-",end,"-renum.pdb",sep="")
	print(f_i)
	write_pdb(resi_sh2,f_i)
	return()
}


mapping<-commandArgs(TRUE)[1]
map<-read.csv(mapping,sep="\t",stringsAsFactors=F)

pfam<-'PF00017'
map_pf<-map[grepl(pfam,map[,5]),]

#dir='/home/koshka/work/DAVID/homologs_contact_map/SH2_pdb_70'
#dir2='/home/koshka/work/DAVID/homologs_contact_map/SH2_pdb_70_renum'

dir=commandArgs(TRUE)[2]
dir2=commandArgs(TRUE)[3]

if(file.exists(dir2)) unlink(dir2, recursive=T)
dir.create(dir2)

for(i in 1:nrow(map_pf)){
	pdb<-map_pf[i,1]
	chain<-map_pf[i,2]
	start<-as.numeric(map_pf[i,3])
	end<-as.numeric(map_pf[i,4])
	seq<-restrict_domain(pdb,dir,dir2,chain,start,end)
}

