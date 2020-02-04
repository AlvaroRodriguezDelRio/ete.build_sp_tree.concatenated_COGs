# Alvaro Rodriguez del Rio
# script for automatically calculating a species tree based on the 41 COGs concatenateed
# input: path to directory with one proteome (.faa extension) per species. The species name must be contained in the file name
# The file names will be the names in the branches of the final tree

import re
import sys
from os import listdir
from os.path import isfile, join
import subprocess
import os
from collections import Counter
import glob
from Bio import SeqIO
import operator
import argparse
import numpy as np

def run(cmd):
	try:
    		s = subprocess.call(cmd, shell=True, executable='/bin/bash')
	except:
        	print >>sys.stderr, '#### cmd failed', rescue[0], len(alg), cmd
        	print >>sys.stderr, '#### Error failed result:', rescue[0], len(alg)
                sys.stderr.flush()

def run_hmmer_search(proteome_fasta_file, species_name, hmmsearch_ouput_dir, COG_hmm, COG_name):
	print ("Running hmmer " + COG_name +" vs "+proteome_fasta_file)
	run('source "%s" && hmmsearch --tblout %s %s %s > %s' %\
       		 ("/scratch/alvaro/catalog/MGS_annotation/annotator/bin/config.sh", hmmsearch_ouput_dir+species_name+"."+COG_name+".hmm", COG_hmm, proteome_fasta_file, hmmsearch_ouput_dir+"tmp"))


def parse_hmmer_search_out(hmmsearch_ouput_dir,species_name,COG_name,min_eval):
	significant_hit_names_eval = dict()
        with open(hmmsearch_ouput_dir+species_name+"."+COG_name+".hmm") as table:
                for line in table:
			if not re.match("#",line):
				gene_name = re.split("\s+", line)[0]
				eval = float(re.split("\s+",line)[4])
				if eval < min_eval:
					significant_hit_names_eval[gene_name] = float(eval)

	if len(significant_hit_names_eval)>1:
		print (str(len(significant_hit_names_eval))+" hits found for "+COG_name+" in species "+species_name+", selecting best for building tree")
	elif len(significant_hit_names_eval)==1:
		print ("1 hit found for "+COG_name+" in species "+species_name)
	else:
		print ("No hit found for "+COG_name+" in species "+species_name)
	return(significant_hit_names_eval)


def get_best_hit(hits):
	min_eval_found = min(np.array(hits.values()))
	best_hit = [key for key in hits if hits[key] == min_eval_found][0]
	return(best_hit)

####
# main program
####

# parse options
parser = argparse.ArgumentParser(description='Calculates species tree after looking for the 41 bacterial COGs in the proteomes')
parser.add_argument("dir",help="directory where the proteomes are stored (.faa extension)",type=str)
parser.add_argument("-n", help="directory where the fasta file for each of the proteimes is (.fna extension)",type=str)
parser.add_argument("-e", help="minimum e-value for considering a hit as significant", type=float)

args = parser.parse_args()
faa_file_dir = args.dir
fna_file_dir = args.n
min_eval = args.e
if fna_file_dir:
	print("Reading nucleotide fastas in "+fna_file_dir)
	fna_files = glob.glob(fna_file_dir+'*.fna')
else:
	print("No nucleotide fasta files provided")

# variable definition
COG_hhm_dir = "/scratch/alvaro/ete.COG_species_tree/data/fetchMG/lib/"
hmmsearch_ouput_dir = "out_hmm/"
sep = "&"

# preparing results dir & getting input files
run('mkdir -p %s' %(hmmsearch_ouput_dir))
faa_files = glob.glob(faa_file_dir+'*.faa') 
COG_hmm_files = glob.glob(COG_hhm_dir+"*.hmm")

# prepare re and variables for getting COG and sp names
pattern_faa_name = re.compile("([A-Za-z0-9\_]+)\.faa")
species_genes_dict = dict()
pattern_COG_name = re.compile(".*\/(COG.*).hmm")
COG_names = set()

# run hmmer search for each species and COG
for file in faa_files:
	sp_name = pattern_faa_name.search(file).group(1)
	print (sp_name)
	species_genes_dict[sp_name] = []
	for COG_hmm in COG_hmm_files:
		COG_name = pattern_COG_name.search(COG_hmm).group(1)
		COG_names.add(COG_name)
		run_hmmer_search(file, sp_name, hmmsearch_ouput_dir, COG_hmm, COG_name)


# take hits for each species and COG, write to fasta and to cog file
gene_COG_sp = dict()
fasta_input_ete_build = open("proteome_seqs.faa","w")
for COG in COG_names:
	gene_COG_sp[COG] = dict()
	for species in species_genes_dict:
		
		# get hits by parsing hmm search output
		hits = parse_hmmer_search_out(hmmsearch_ouput_dir,species,COG,min_eval)
		if len(hits) == 0:
			continue

		# get best hit
		best_hit = get_best_hit(hits)
		gene_COG_sp[COG][species] = best_hit
		species_genes_dict[species].append(best_hit)

# add sequences of COGs in proteome of the species to fasta file
print ("adding sequences to faa file")
fasta_input_ete_build = open("proteome_seqs.faa","w")
for species in species_genes_dict:
	print (species)
	for seq in SeqIO.parse(faa_file_dir+species+".faa", "fasta"):
		if seq.id in species_genes_dict[species]:
			seq.id = species + sep + seq.id
			seq.description = ""
			SeqIO.write(seq, fasta_input_ete_build, "fasta")

if fna_file_dir:
	print ("add sequences to fna file")
	fna_input_ete_build = open("proteome_seqs.fna","w")
	for species in species_genes_dict:
		for seq in SeqIO.parse(fna_file_dir+species+".fna", "fasta"):
			if seq.id in species_genes_dict[species]:
				seq.id = species + sep + seq.id
				seq.description = ""
				SeqIO.write(seq, fna_input_ete_build, "fasta")


# create orthology file
cogs_file_input_ete_build = open("cog_file.txt","w")
for COG in gene_COG_sp:
	out_cogs_file = []
	for species in gene_COG_sp[COG]:
		out_cogs_file.append(species + sep + gene_COG_sp[COG][species])
	if len(out_cogs_file) > 0:
		cogs_file_input_ete_build.write('\t'.join(out_cogs_file)+"\n")

# run ete
print ("Running ete build...")
if fna_file_dir:
	run('export PATH=/scratch/alvaro/miniconda2/bin/:$PATH && ete3 build -n proteome_seqs.fna -a proteome_seqs.faa --cogs cog_file.txt --spname-delimiter "%s" -o tree -m sptree_fasttree_100 -w standard_fasttree --clearall --rename-dup-seqnames --noimg' %\
	(sep))	
else:
	run('export PATH=/scratch/alvaro/miniconda2/bin/:$PATH && ete3 build -a proteome_seqs.faa --cogs cog_file.txt --spname-delimiter "%s" -o tree -m sptree_fasttree_100 --nt-switch-threshold 0 -w standard_fasttree --clearall --rename-dup-seqnames --noimg' %\
	(sep))

print ("Done")
print ("Results are in " + "tree/cog_100-alg_concat_default-fasttree_full/")
