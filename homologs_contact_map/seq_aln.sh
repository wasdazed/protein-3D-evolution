dir='/home/koshka/work/DAVID/homologs_contact_map/SH2_pdb_70/'
dir2='/home/koshka/work/DAVID/homologs_contact_map/SH2_pdb_70_renum/'
dir3='/home/koshka/work/DAVID/homologs_contact_map/SH2_pdb_70_beads/'
dir4='/home/koshka/work/DAVID/homologs_contact_map/SH2_pdb_70_contact_maps/'

mapping='/home/koshka/work/DAVID/homologs_contact_map/pdb_pfam_mapping.txt'
fasta='/home/koshka/work/DAVID/homologs_contact_map/SH2_pdb_70.fa'
fasta_aln='/home/koshka/work/DAVID/homologs_contact_map/SH2_pdb_70.fasta_aln'
pdf_file='/home/koshka/work/DAVID/homologs_contact_map/SH2_pdb_70_lendth.pdf'

mistral_path='/home/koshka/work/DAVID/homologs_contact_map/mistral3.6_exec/Linux_v3.6_32bit/./mistral.x'

#Rscript /home/koshka/work/DAVID/homologs_contact_map/Rscripts/make_fasta.R $dir3 $fasta $pdf_file

#t_coffee -in $fasta -output=aln,fa

Rscript /home/koshka/work/DAVID/homologs_contact_map/Rscripts/filter_aln.R $fasta_aln