dir='/home/koshka/work/DAVID/homologs_contact_map/SH2_pdb_70/'
dir2='/home/koshka/work/DAVID/homologs_contact_map/SH2_pdb_70_renum/'
dir3='/home/koshka/work/DAVID/homologs_contact_map/SH2_pdb_70_beads/'
dir4='/home/koshka/work/DAVID/homologs_contact_map/SH2_pdb_70_contact_maps/'

mistral_path='/home/koshka/work/DAVID/homologs_contact_map/mistral3.6_exec/Linux_v3.6_32bit/./mistral.x'


#structural aln

pdb_array=$(find $dir2 -name '*.pdb')

$mistral_path $padb_array


