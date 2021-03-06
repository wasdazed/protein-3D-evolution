dir='/home/koshka/work/DAVID/homologs_contact_map/SH2_pdb_70/'
dir2='/home/koshka/work/DAVID/homologs_contact_map/SH2_pdb_70_renum/'
dir3='/home/koshka/work/DAVID/homologs_contact_map/SH2_pdb_70_beads/'
dir4='/home/koshka/work/DAVID/homologs_contact_map/SH2_pdb_70_contact_maps/'

mapping='/home/koshka/work/DAVID/homologs_contact_map/pdb_pfam_mapping.txt'
fasta='/home/koshka/work/DAVID/homologs_contact_map/SH2_pdb_70.fa'
pdf_file='/home/koshka/work/DAVID/homologs_contact_map/SH2_pdb_70_lendth.pdf'

mistral_path='/home/koshka/work/DAVID/homologs_contact_map/mistral3.6_exec/Linux_v3.6_32bit/./mistral.x'

rm -r $dir2
Rscript /home/koshka/work/DAVID/homologs_contact_map/Rscripts/renumber-chains.R $mapping $dir $dir2

rm -r $dir3
mkdir $dir3
rm -r $dir4
mkdir $dir4

for f in $(find $dir2 -name '*.pdb')
do
	echo $f
	n=$(echo $f | sed 's/.*\///')
	n=$(echo $n | sed 's/.pdb//')
	echo $n
	f2=$dir3$n.bead
	/home/koshka/work/DAVID/cass-1.1.1/./pdb2bead $f $f2
	FILESIZE=$(stat -c%s $f2)
	if [[ $FILESIZE = 0 ]]
	then 
	rm $f2
	rm $f
	continue
	else
	f3=$dir4$n.cm
	/home/koshka/work/DAVID/cass-1.1.1/./nadia_calc-contact-maps $f2 $f3
	fi
done


