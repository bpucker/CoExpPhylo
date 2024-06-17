### Nele Fiene, Boas Pucker ###
### b.pucker@tu-braunschweig.de ###

### some functions taken from MYB_annotator.py and KIPEs ###
### some code blocks were generated using GPT 3.5 and GPT 4o ###

__version__ = "v1.17"

__usage__ = """
				python3 coexp_phylo.py
				--config <CONFIG_FILE>
				--out <OUTPUT_FOLDER>
				--mode <RUNNING_MODE>(batch|all)[batch]		# mode = all not implemented yet
				
				optional:
				----
				ANNOTATION
				----
				--anno <ANNOTATION_FILE>
				--araport <ARAPORT_PEP_FILE>
				--seqs_cluster_anno <PERCENTAGE_OF_SEQUENCES_PER_CLUSTER_USED_FOR_ANNOTATION>[50.0]

				----
				PER-SPECIES COEXPRESION ANALYSIS
				----
				--r <R_CUTOFF_FOR_CORRELATION_ANALYSIS>[0.7]
				--p <PVALUE_CUTOFF_FOR_CORRELATION_ANALYSIS>[0.05]
				--numcut <NUMBER_OF_COEXPRESSED_GENES_CONSIDERED>[100]
				--min_exp_cutoff <MIN_EXP_FOR_GENE_TO_CONSIDER>[30]
				--coexp_script_path <PATH_TO_HELPER_SCRIPT>[coexp_helper.py]

				----
				DIAMOND
				----
				--cpub <NUMBER_OF_CPUs_TO_USE_FOR_BLAST>[4]
				--batch <NUMBER_OF_BLASTP-JOBS_TO_RUN_IN_PARALLEL>[7]
				--scorecut <MIN_BLAST_SCORE_CUTOFF>[100]
				--simcut <MIN_BLAST_SIMILARITY_CUTOFF>[80]
				--lencut <MIN_BLAST_LENGTH_CUTOFF>[100]
				--evalue <E-VALUE_CUTOFF>

				----
				CLUSTERING
				----
				--mindetect <MIN_NUMBER_OF_COEXPED_BAITS>[1]
				--minseqcutoff <MIN_SEQ_CLUSTER_SIZE>[10]
				--mincoexpseqcutoff <MIN_SPECIES_WITH_COEXPED_GENES_FOR_CLUSTER>[3]
				--minintersec <INTERSEC_BETWEEN_CLUSTERS_FOR_MERGE>[1]
					
				----
				ALIGNMENT
				----
				--alnmethod <ALIGNMENT_ALGORITH>(mafft|muscle)[mafft]
				--mafft <MAFFT_PATH>[mafft]
				--muscle <MUSCLE_PATH>[muscle]

				----
				TREE CONSTRUCTION
				----
				--treemethod <TREE_ALGORITHM>(fasttree|raxml|iqtree)[fasttree]
				--raxml <RAXML_PATH>[raxml]
				--fasttree <FASTTREE_PATH>[FastTree]
				--iqtree <IQ-TREE_PATH>[iqtree]
				--cpur <CPUs_FOR_TREE_CONSTRUCTION(only for raxml and iqtree)[cpub]
					
				"""

import os, sys, subprocess, math, glob
from operator import itemgetter
import csv
import hashlib
import numpy as np
from scipy import stats
import datetime
import networkx as nx
import random


# --- end of imports --- #

def check_tools(name, docu_file):
	# check whether all tools and packages are available
	try:
		result = subprocess.run([name, '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
		version_info = result.stdout.strip() or result.stderr.strip()
		with open( docu_file, 'a') as f:
			f.write( str( name + ': ' + version_info + '\n--\n' ) )
		return True
	except subprocess.CalledProcessError:
		try:
			result = subprocess.run([name, '-version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
			version_info = result.stdout.strip() or result.stderr.strip()
			with open( docu_file, 'a') as f:
				f.write( str( name + ': ' + version_info + '\n--\n' ) )
			return True
		except subprocess.CalledProcessError:
			try:
				subprocess.run([name], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
				with open( docu_file, 'a') as f:
					f.write( str( name + ': no version number available \n--\n' ) )
				return True
			except subprocess.CalledProcessError:
				subprocess.run([name, 'echo', 'check'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
				with open( docu_file, 'a') as f:
					f.write( str( name + ': no version number available \n--\n' ) )
				return True
	except FileNotFoundError:
		return False

def calculate_md5sum( file ):
	hash_md5 = hashlib.md5()
	try:
		with open( file, "rb" ) as f:
			for chunk in iter( lambda: f.read( 4096 ), b"" ):
				hash_md5.update( chunk )
		return hash_md5.hexdigest()
	except FileNotFoundError:
		return str( file + "File not found" )

def generate_docu_file ( docu_file, config_file, mode, rcutoff, pvaluecutoff, scorecut, simcut, lencut, evalue, mindetection, minseqcutoff, mincoexpseqcutoff, min_exp_cutoff, min_intersection_cutoff, alnmethod, treemethod, anno_file, araport_seq_file, iqtree, seqs_cluster_anno ):
	with open( docu_file, 'w') as f:
		f.write('--------------\nDOCUMENTATION\n--------------\n\nSCRIPT VERSION\n--------------\n' )
		f.write( str ( __version__ ) + '\n\nCONFIG FILE\n--------------\n')
		f.write( str ( config_file + '\t\t' + str( calculate_md5sum( config_file ) ) + '\n\n' ) )
		f.write( 'PARAMETER\n--------------\n' )
		f.write( str ( 'mode = ' + mode + '\n' ) )	
		f.write( str ( 'r-cutoff = ' + str( rcutoff ) + '\n' ) )
		f.write( str ( 'p-value cutoff = ' + str( pvaluecutoff ) + '\n' ) )
		f.write( str ( 'scorecut = ' + str( scorecut ) + '\n' ) )	
		f.write( str ( 'simcut = ' + str( simcut ) + '\n' ) )
		f.write( str ( 'lencut = ' + str( lencut ) + '\n' ) )
		f.write( str ( 'evalue = ' + str( evalue ) + '\n' ) )
		f.write( str ( 'mindetect = ' + str( mindetection ) + '\n' ) )
		f.write( str ( 'minseqcutoff = ' + str( minseqcutoff ) + '\n' ) )
		f.write( str ( 'mincoexpseqcutoff = ' + str( mincoexpseqcutoff ) + '\n' ) )
		f.write( str ( 'min_exp_cutoff = ' + str( min_exp_cutoff ) + '\n' ) )
		f.write( str ( 'min_intersection_cutoff = ' + str( min_intersection_cutoff ) + '\n' ) )
		f.write( str ( 'alignment method = ' + alnmethod + '\n' ) )
		f.write( str ( 'tree method = ' + treemethod + '\n' ) )
		f.write( str ( 'anno = ' + anno_file + '\n' ) )
		f.write( str ( 'araport file = ' + araport_seq_file + '\n' ) )
		f.write( str ( 'seqs per cluster used for anno [%] = ' + str( seqs_cluster_anno ) + '\n' ) )
		f.write('\nTOOLS\n--------------\n')
	if treemethod != "iqtree":
		required_tools = [ 'nohup', 'parallel', 'diamond', alnmethod, treemethod ]
	else:
		required_tools = [ 'nohup', 'parallel', 'diamond', alnmethod, iqtree ]

    # Check each tool and handle the case where a tool is missing
	missing_tools = []
	for tool in required_tools:
		if not check_tools(tool, docu_file):
			missing_tools.append( tool )

	if missing_tools != []:
		print( 'Error: Required tools ', missing_tools , ' are not installed.\n Please install the required tools before starting again.' )
		sys.exit(1)

	with open( docu_file, 'a') as f:
		f.write( '\nDATASETS\t\t\t\tmd5sum\n--------------\n' )
	with open( config_file, "r" ) as config:
		config_content = list( csv.reader( config ) )

	with open( docu_file, 'a') as f:
		for i, row in enumerate( config_content ):
			f.write( row[0] + "\n" )		#writes species names

			for file_entry in row[1:]:
				md5sum = calculate_md5sum( file_entry )
				f.write( file_entry + "\t\t" + md5sum + "\n" )
			f.write( "---\n")
		f.write( "EOF\n" )

def translate( seqs ):
	"""! @brief translates the given nucleotide sequence into a peptide sequence
	@return dictionary with peptide sequences as values
	"""
	
	genetic_code = {	'CTT': 'L', 'ATG': 'M', 'AAG': 'K', 'AAA': 'K', 'ATC': 'I',
								'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T',
								'AGC': 'S', 'ACA': 'T', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N',
								'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H',
								'ACG': 'T', 'CCG': 'P', 'AGT': 'S', 'CAG': 'Q', 'CAA': 'Q',
								'CCC': 'P', 'TAG': '*', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C',
								'CGA': 'R', 'CCA': 'P', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R',
								'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': '*', 'GGA': 'G',
								'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'GAG': 'E', 'TCG': 'S',
								'TTA': 'L', 'GAC': 'D', 'TCC': 'S', 'GAA': 'E', 'TCA': 'S',
								'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A',
								'GTG': 'V', 'TTC': 'F', 'GTT': 'V', 'GCT': 'A', 'ACC': 'T',
								'TTG': 'L', 'CGT': 'R', 'TAA': '*', 'CGC': 'R'
							}
	
	final_peptide_seqs = {}
	for key in list( seqs.keys() ):	#iterate over all given nucleotide sequences
		seq = seqs[ key ].upper()
		peptide = []
		for i in range( int( len( seq ) / 3.0 ) ):	#loop to iterate over codons
			codon = seq[i*3:i*3+3]	#extracts a codon from the given nucleotide sequence
			try:
				peptide.append( genetic_code[ codon ] )
			except:
				peptide.append( "*" )
		final_peptide_seqs.update( { key: "".join( peptide ) } )	#store peptide sequence in dictionary
	return final_peptide_seqs


def load_sequences( fasta_file ):
	"""! @brief load sequences of given FASTA file into dictionary with sequence IDs as keys and sequences as values """
	
	sequences = {}
	with open( fasta_file ) as f:
		header = f.readline()[1:].strip()
		if " " in header:	#take only the space-free part of the header (if space present)
			header = header.split(' ')[0]
		if "\t" in header:	#take only the tab-free part of the header (if tab present)
			header = header.split('\t')[0]
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: "".join( seq ) } )
					header = line.strip()[1:]
					if " " in header:
						header = header.split(' ')[0]
					if "\t" in header:
						header = header.split('\t')[0]
					seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		sequences.update( { header: "".join( seq ) } )	
	return sequences


def load_baits( baitfile, spec_id_prefix ):
	"""! @brief load all bait IDs with the corresponding species prefix
	@warning this function requires adjustment if KIPEs input is modified
	"""
	
	IDs = []
	with open( baitfile, "r" ) as f:
		line = f.readline()
		#WARNING: adjustment required if modified in KIPEs
		if line == "ID	Gene	Similarity	ConservedResidues	ConservedRegions\n":	#remove first line if header produced by KIPEs
			line = f.readline()
		while line:
			if "\t" in line:
				parts = line.strip().split('\t')
				IDs.append( spec_id_prefix + "@" + parts[0] )
			else:
				IDs.append( line.strip() )
			line = f.readline()
	return IDs


def load_config_file_content( config_file, output_folder ):
	"""! @brief load content of config file, generate clean CDS file, and generate clean PEP file """
	
	infos = {}
	with open( config_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split(',')
				pep_file = output_folder + parts[0] + ".pep.fasta"	#name of new peptide sequence file
				cds_file = output_folder + parts[0] + ".cds.fasta"	#name of new CDS file
				cds_seqs = load_sequences( parts[2] )
				find_cds_seqs = {}
				with open( cds_file, "w" ) as out:	#generating clean CDS file
					for key in list( cds_seqs.keys() ):
						out.write( '>' + parts[0] + "@" + key + '\n' + cds_seqs[ key ] + "\n" )	#include species name in sequence name
						find_cds_seqs.update( { parts[0] + "@" + key: cds_seqs[ key ] } )
				pep_seqs = translate( find_cds_seqs )	#translation of all CDS
				with open( pep_file, "w" ) as out:	#generate clean PEP file
					for key in list( pep_seqs.keys() ):
						out.write( '>' + key + '\n' + pep_seqs[ key ] + "\n" )
				baits = load_baits( parts[3], parts[0] )
				#storing information in a dictionary under the species name
				infos.update( { parts[0]: { 'id': parts[0], 'tpm': parts[1], 'cds_file': cds_file, 'baits_file': parts[3], 'pep_file': pep_file, 'cds': find_cds_seqs, 'pep': pep_seqs, 'baits': baits } } )
			line = f.readline()
	return infos


def load_expression_values( filename, spec_prefix ):
	"""! @brief load all expression data per species
		@return dictionary with expression as value and species name + gene name as key
	"""
	
	expression_data = {}	#dictionary has species name + gene name as key; value is dictionary with accession as key and TPM as value
	with open( filename, "r" ) as f:
		header = f.readline().split('\t')
		if header[0] == '':		#check if first entry is empty and replace with 'gene'
			header[0] = 'gene'
		header[-1] = header[-1].strip() 	#remove \n at last entry of the header
		tissues = header[1:]	#list of accessions (first field can be empty string i.e. invisible in file)
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			expression = {}
			for idx, each in enumerate( parts[1:] ):
				expression.update( { tissues[  idx ] : float( parts[ idx+1 ] ) } )	#convert value into float when adding to dictionary
			line = f.readline()
			expression_data.update( { spec_prefix + "@" + parts[0]: expression } )
	return expression_data


def compare_candidates_against_all( candidate, spec_prefix, gene_expression, rcutoff, pvaluecutoff, coexpnumcutoff, min_exp_cutoff ):
	"""! @brief compare candidate gene expression against all genes to find co-expressed genes
		@return list of co-expressed genes that passed the filters
	"""
	
	tissues = sorted( list( gene_expression[ list( gene_expression.keys() )[0] ].keys() ) )	#accessions / sample names
	coexpressed_genes = []
	candidate = spec_prefix + "@" + candidate
	for i, gene2 in enumerate( list( gene_expression.keys() ) ):	#iterate over all genes
		if candidate != gene2:
			values = []
			total_expression = 0
			for tissue in tissues:
				try:
					x = gene_expression[ candidate ][ tissue ]
					y = gene_expression[ gene2 ][ tissue ]
					total_expression += y
					if not math.isnan( x ) and not math.isnan( y ) :	#exclude undefined values
						values.append( [ x, y ] )
				except KeyError:
					pass
			if total_expression > min_exp_cutoff:	#only consider genes with a minimal expression level
				if len( values ) > 0:
					r, p = stats.spearmanr( values )
					if not math.isnan( r ):
						if r > rcutoff and p < pvaluecutoff:
							coexpressed_genes.append( { 'id': gene2, 'correlation': r, 'p_value': p } )
				else:
					sys.stdout.write ( "WARNING: no expression detected - " + candidate + "\t" + gene2 + "\t" + ";".join( list( gene_expression.keys() )[:5] ) + "\n" )
					sys.stdout.flush()
	sys.stdout.flush()
	coexp_gene_IDs = []
	for gene in sorted( coexpressed_genes, key=itemgetter('correlation') )[::-1][:coexpnumcutoff]:
		coexp_gene_IDs.append( gene['id'] )
	return coexp_gene_IDs


def load_blast_hits( blast_result_file, scorecut, simcut, lencut, evalue):
	""""! @brief load BLAST hits above a certain score per bait
		@return list of sequence IDs that are considered good BLAST hits
	"""
	
	hits = []

	try:
		with open( blast_result_file, "r" ) as f:
			line = f.readline()
			while line:
				parts = line.strip().split('\t')
				if float( parts[-1] ) > scorecut:
					if float( parts[2] ) > simcut:
						if int( parts[3] ) > lencut:
							if float( parts[-2] ) <= evalue:
								hits.append( parts[1] )
				line = f.readline()

	except FileNotFoundError:
		pass

	return list( set( hits ) )	#reduce to list of unique elements


def load_hits_per_bait( blast_result_file, scorecut, simcut, lencut, evalue ):
	"""! @brief load BLAST hits per bait
	@return dictionary with baits as keys and lists of hits as values
	"""
	
	hits = {}
	with open( blast_result_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if float( parts[-1] ) > scorecut:	#BLAST hits are filtered with user defined criteria
				if float( parts[2] ) > simcut:
					if int( parts[3] ) > lencut:
						if float( parts[-2] ) <= evalue:
							try:
								hits[ parts[0] ].append( parts[1] )	#add to existing dictionary entry
							except KeyError:
								hits.update( { parts[0]: [ parts[1] ] } )	#generate new entry in dictionary
			line = f.readline()
	return hits


def load_alignment( aln_file ):
	"""! @brief load FASTA alignment into dictionary
	@return dictionary with sequence IDs as keys and corresponding alignment sequences as values
	"""
	
	sequences = {}
	with open( aln_file ) as f:
		header = f.readline()[1:].strip()
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: "".join( seq ) } )
					header = line.strip()[1:]
					seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		sequences.update( { header: "".join( seq ) } )	
	return sequences


def alignment_trimming( aln_file, cln_aln_file, occupancy ):
	"""! @brief remove all alignment columns with insufficient occupancy; produce clean file"""
	
	alignment = load_alignment( aln_file )	#alignment in FASTA format is loaded from file
	# --- if there is an alignment (expected case) 
	if len( list(alignment.keys()) ) > 0:
		# --- identify valid residues in aligned sequences (columns with sufficient occupancy) --- #
		valid_index = []
		for idx, aa in enumerate( list(alignment.values())[0] ):	#iterate over all positions in the aligned sequences
			counter = 0
			for key in list(alignment.keys()):
				if alignment[ key ][ idx ] != "-":	#check for match in other sequences (all sequences have same length)
					counter += 1
			if counter / float( len( list(alignment.keys()) ) ) > occupancy:	#collect positions of sufficient occupancy in list
				valid_index.append( idx )
		
		# --- generate new sequences --- #
		with open( cln_aln_file, "w" ) as out:
			for key in list(alignment.keys()):	#iterate over all sequences in the alignment
				seq = alignment[ key ]
				new_seq = []
				for idx in valid_index:	#only add residues of positions with sufficient occupancy
					new_seq.append( seq[ idx ] )
				out.write( ">" + key + '\n' + "".join( new_seq ) + '\n' )
	# --- just in case the alignment file is empyt (is this possible?) ---#
	else:
		with open( cln_aln_file, "w" ) as out:
			out.write( "" )


def tree_constructor( X_aln_input_file, treemethod, X_output_folder, Xname, alnmethod, tree_cmd_file, mafft, muscle, raxml, fasttree, iqtree, cpur ):
	"""! @brief handles the construction of alignments and phylogenetic tree
			@note second FASTA file can be an empty string to run this function just based on one FASTA file
	"""
	
	X_aln_file = X_aln_input_file + ".aln"
	X_cln_aln_file = X_aln_file + ".cln"
	
	# --- generate alignment --- #
	if not os.path.isfile( X_aln_file ):
		if alnmethod == "muscle":
			p = subprocess.Popen( args= muscle + " -align " + X_aln_input_file + " -output " + X_aln_file, shell=True )
			p.communicate()
		else:
			p = subprocess.Popen( args= mafft + " --quiet " + X_aln_input_file + " > " + X_aln_file, shell=True )
			p.communicate()
	
	# --- trim alignment i.e. remove positions with many gaps = low occupancy --- #
	if not os.path.isfile( X_cln_aln_file ):
		alignment_trimming( X_aln_file, X_cln_aln_file, occupancy=0.1 )
	
	if treemethod == "raxml":	#construct tree with RAxML
		prefix = X_output_folder + Xname + "RAxML_tree"
		tree_file = prefix + ".raxml.bestTree.tre"
		current_tree_command = raxml + " --all --threads " + str( cpur ) + " --model LG+G8+F --msa" + X_cln_aln_file + "--prefix" + prefix
	
	elif treemethod == "iqtree":	#construct tree with IQ-TREE2
		prefix = X_output_folder + Xname + "IQtree"
		tree_file = prefix + ".treefile"	#".treefile" is appended to provided alignment file name
		current_tree_command = iqtree + " -ntop " + str( cpur ) + " -alrt 1000 -bb 100 -s " + X_cln_aln_file + " --prefix " + prefix
	
	else:	#construct tree with FastTree2
		tree_file = X_output_folder  + Xname  + "FastTree_tree.tre"
		current_tree_command = fasttree + " -wag  -nopr -nosupport < " + X_cln_aln_file + " > " + tree_file
	
	with open( tree_cmd_file, 'a' ) as tree_command_out:
		tree_command_out.write( current_tree_command + "\n" )


def count_coexp_specs( seq_IDs_to_check ):
	"""! @brief count species with co-expressed sequences """
	
	spec_list = set()
	for each in seq_IDs_to_check:	
		if "_coexp" in each:
			species = each.split('@')[0]
			spec_list.add( species )
	return len( spec_list )


def load_annotation( anno_file ):
	"""! @brief load annotation from given file
	@return dictionary with sequence ID as key and corresponding annotation as value
	"""
	
	mapping_table = {}
	with open( anno_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')	#file needs to be TAB-separated
			mapping_table.update( { parts[0]: ". ".join( parts[1:] ) } )	#all columns following the first one are annotation
			line = f.readline()
	return mapping_table


def annotate_group( blast_result_file, anno_mapping_table ):
	#return string comprizing functional annotation if available
 
	best_hit = ""
	best_hit_score = 0
	with open( blast_result_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if float( parts[-1] ) > best_hit_score:
				best_hit_score = float( parts[-1] )
				best_hit = parts[1]
			line = f.readline()
	try:
		return best_hit + " - " + anno_mapping_table[ best_hit ]
	except KeyError:
		return best_hit + " - " + "n/a"


def get_groups(annotation_file, araport_seq_file, blast_results, anno_mapping_table, huge_pep_collection, aln_folder, mincoexpseqcutoff, minseqcutoff, tmp_folder, min_intersection_cutoff, cpub, batch, perc_cluster_anno):
	"""! @brief identify groups in huge cluster """

	# --- identify groups --- #
	if min_intersection_cutoff == 1:  # computationally easier
		G = nx.Graph()
		for key, values in blast_results.items():
			for value in values:
				G.add_edge(key, value)
			G.add_node(key)			# to ensure isolated nodes are also present
		final_groups = [list(component) for component in nx.connected_components(G)]	 #Identification of connected components

	else:   # needs to be impemented
		final_groups = []
		sys.stdout.write( str( datetime.datetime.now() ) + " CLUSTERING WITH MIN_INTERSECTION_CUTOFF > 1 NOT IMPLEMENTED YET!! \n" )
		sys.stdout.flush()

	sys.stdout.write( str( datetime.datetime.now() ) + " ... groups identified ... \n" )
	sys.stdout.flush()
	
	# --- prepare database for annotation --- #
	if len( araport_seq_file ) > 0:
		tmp_blast_folder = tmp_folder + "anno_blast/"
		if not os.path.exists( tmp_blast_folder ):
			os.makedirs( tmp_blast_folder )
		blast_db = tmp_blast_folder + "blastdb"
		p = subprocess.Popen( args= "diamond makedb --in " + araport_seq_file + " -d " + blast_db, shell=True )
		p.communicate()
	
	# --- write sequences into output files --- #
	if len( araport_seq_file ) > 0:
		anno_blast_cmd_file = tmp_blast_folder + "ANNO_BLAST_COMMANDS.sh"
		with open( anno_blast_cmd_file, "w" ) as anno_command_out:
			anno_command_out.write( str ( "#!/bin/bash\nnohup parallel -j " + str( batch ) + " <<EOF \n" ) )
	
	with open( annotation_file, "w" ) as anno_out:
		counter = 0
		final_cluster = [] 		#identify groups with sufficient coexpressed candidates 

		for fg in final_groups:
			if len( fg ) > minseqcutoff: 
				coexpseqcounter = count_coexp_specs( fg )
				if coexpseqcounter >= mincoexpseqcutoff:
					counter += 1
					final_cluster.append( fg )
					group_sequence_file = aln_folder + str( counter ).zfill( 3 ) + ".pep.fasta"	#save all sequenes of a group in a FASTA file
					group_anno_file = tmp_blast_folder + str( counter ).zfill( 3 ) + "_anno.pep.fasta"
					seq_collection_anno = []
					for gene in fg:
						current_seq = '>' + gene + "\n" + huge_pep_collection[ gene ] + "\n"
						with open( group_sequence_file, "a" ) as out:
							out.write( current_seq )
						if len( araport_seq_file ) > 0:
							seq_collection_anno.append( current_seq )
					
					seq_num = int( len ( seq_collection_anno ) * ( perc_cluster_anno / 100 ) )  #determine number of sequences for annotation
					final_seqs = random.sample( seq_collection_anno, seq_num )
					with open( group_anno_file, 'a') as out:
						for seq in final_seqs:
							out.write( seq )

					if len( araport_seq_file ) > 0:
						with open(anno_blast_cmd_file, "a" ) as anno_command_out:
							blast_result_file = tmp_blast_folder + group_sequence_file.split('/')[-1] + "_blast_results.txt"
							current_anno_command = "diamond blastp -q " + group_anno_file + " -d " + blast_db + " -o " + blast_result_file + " --threads " + str( cpub ) + "\n"
							anno_command_out.write( current_anno_command )

		if len( araport_seq_file ) > 0:
			with open( anno_blast_cmd_file, "a" ) as anno_command_out:
				anno_command_out.write( "EOF" )
			os.chmod( anno_blast_cmd_file, 0o755 )
			subprocess.run( ['nohup', '/bin/bash', anno_blast_cmd_file], check = True )
			sys.stdout.write( str( datetime.datetime.now() ) + " - DIAMOND blast for annotation completed! " + "\n" )
			sys.stdout.flush()

		counter = 0
		for fg in final_cluster:  # return best annotation
			counter += 1
			if len( araport_seq_file ) > 0: 
				group_sequence_file = aln_folder + str( counter ).zfill( 3 ) + ".pep.fasta"
				blast_result_file = tmp_blast_folder + group_sequence_file.split('/')[-1] + "_blast_results.txt"
				group_annotation = annotate_group( blast_result_file, anno_mapping_table )
			else:
				group_annotation = "n/a"
			anno_out.write( str( counter ).zfill(3) + "\t" + str( len( fg ) ) + "\t" + group_annotation + "\n" )


def main( arguments ):
	"""! @brief run everything """

	config_file = arguments[ arguments.index('--config')+1 ]
	#ID,tpm,cds,baits
	#ID = species ID/name; tpm= tpm file; cds = CDS file; baits = IDs of interest in file
	output_folder = arguments[ arguments.index('--out')+1 ]
	mode = arguments[ arguments.index('--mode')+1 ]	
	#batch: n jobs run in parallel
	#all: all jobs are send to the cluster at once
	
	if '--anno' in arguments:
		anno_file = arguments[ arguments.index('--anno')+1 ]
	else:
		anno_file = ""
	
	if '--r' in arguments:
		rcutoff = float( arguments[ arguments.index('--r')+1 ] )
	else:
		rcutoff = 0.7	#correlation coefficient cutoff for gene expression analysis
	
	if '--p' in arguments:
		pvaluecutoff = float( arguments[ arguments.index('--p')+1 ] )
	else:
		pvaluecutoff = 0.05	#p-value cutoff for co-expression analysis
	
	if '--numcut' in arguments:
		coexpnumcutoff = int( arguments[ arguments.index('--numcut')+1 ] )
	else:
		coexpnumcutoff = 100	#best X co-expressed genes to consider
	
	if '--cpub' in arguments:
		cpub = int( arguments[ arguments.index('--cpub')+1 ] )
	else:
		cpub = 4
	
	if '--batch'in arguments:
		batch = int( arguments[ arguments.index('--batch')+1 ] )
	else:
		batch = 7
	
	if '--scorecut' in arguments:	#score cutoff for BLAST hits
		scorecut = float( arguments[ arguments.index('--scorecut')+1 ] )
	else:
		scorecut = 100.0
	
	if '--simcut' in arguments:	#similarity cutoff for BLAST hits
		simcut = float( arguments[ arguments.index('--simcut')+1 ] )
	else:
		simcut = 80.0
	
	if '--lencut' in arguments:	#length cutoff for BLAST hits
		lencut = float( arguments[ arguments.index('--lencut')+1 ] )
	else:
		lencut = 100

	if '--evalue' in arguments:
		evalue = float( arguments[ arguments.index('--evalue')+1 ]  )
	else:
		evalue = 1e-5
	
	if '--mindetect' in arguments:
		mindetection = int( arguments[ arguments.index('--mindetect')+1 ] )
	else:
		mindetection = 1	#number of bait genes that a given sequence need to be co-expressed with to be considered (strict would be equal to number of baits)
	
	if '--minseqcutoff' in arguments:
		minseqcutoff = int( arguments[ arguments.index('--minseqcutoff')+1 ] )
	else:
		minseqcutoff = 10	#minimal number of sequences to compose a group as tree construction input
	
	if '--mincoexpseqcutoff' in arguments:
		mincoexpseqcutoff = int( arguments[ arguments.index('--minseqcutoff')+1 ] )
	else:
		mincoexpseqcutoff = 3	#minimal number of species with co-expressed sequences to compose a group as tree construction input
	
	if '--min_exp_cutoff' in arguments:
		min_exp_cutoff = int( arguments[ arguments.index('--min_exp_cutoff')+1 ] )
	else:
		min_exp_cutoff = 30	#minimal cumulative expression per gene to be considered in the co-expresssion analysis	
	
	if '--alnmethod' in arguments:	#alignment method
		alnmethod = arguments[ arguments.index('--alnmethod')+1 ]
	else:
		alnmethod = "mafft"
	
	if '--mafft' in arguments:	#path to MAFFT (alignment)
		mafft = arguments[ arguments.index('--mafft')+1 ]
	else:
		mafft = "mafft"
	
	if '--muscle' in arguments:	#path to MUSCLE (alignment)
		muscle = arguments[ arguments.index('--muscle')+1 ]
	else:
		muscle = "muscle"	
	
	if '--raxml' in arguments:	#path to RAxML (tree construction)
		raxml = arguments[ arguments.index('--raxml')+1 ]
	else:
		raxml = "raxml"
	
	if '--fasttree' in arguments:	#path to FastTree (tree construction)
		fasttree = arguments[ arguments.index('--fasttree')+1 ]
	else:
		fasttree = "FastTree"
	
	if '--iqtree' in arguments:	#path to FastTree (tree construction)
		iqtree = arguments[ arguments.index('--iqtree')+1 ]
	else:
		iqtree = "iqtree2"
	
	if '--cpur' in arguments:	#CPUs to use for tree
		cpur = int( arguments[ arguments.index('--cpur')+1 ] )
	else:
		cpur = min( [ 4, cpub ] )
	
	if '--treemethod' in arguments:	#tree construction methods
		treemethod = arguments[ arguments.index('--treemethod')+1 ]
		if treemethod not in [ "fasttree", "raxml", "iqtree" ]:
			treemethod = "fasttree"
	else:
		treemethod = "fasttree"
	
	if '--araport' in arguments:
		araport_seq_file = arguments[ arguments.index('--araport')+1 ]
	else:
		araport_seq_file = ""

	if '--seqs_cluster_anno' in arguments:
		seqs_cluster_anno = int( arguments[ arguments.index('--seqs_cluster_anno')+1 ] )
		if seqs_cluster_anno > 100.0:
			sys.exit( "ERROR: seqs_cluster_anno value must be <= 100" )
		if seqs_cluster_anno <= 0.0:
			sys.exit( "ERROR: seqs_cluster_anno value must be > 0" )
	else:
		seqs_cluster_anno = 50.0
	
	if '--minintersec' in arguments:
		min_intersection_cutoff = int( arguments[ arguments.index('--minintersec')+1 ] )
	else:
		min_intersection_cutoff = 1	#minimal intersection between two sequence clusters to trigger merging
	
	if '--coexp_script_path' in arguments:
		coexp_script_path = arguments[ arguments.index('--coexp_script_path')+1 ]
	else:
		coexp_script_path = "coexp_helper.py"	#needs coexp helper script path
	
	
	if output_folder[-1] != "/":
		output_folder += "/"
	
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	tmp_folder = output_folder + "tmp/"
	if not os.path.exists( tmp_folder ):
		os.makedirs( tmp_folder )

	### --- Documentation --- ###
	sys.stdout.write( str( datetime.datetime.now() ) + " - preparing documentation file ... \n" )
	docu_file = output_folder + datetime.date.today().strftime('%Y%m%d') + '_docu.txt'
	generate_docu_file ( docu_file, config_file, mode, rcutoff, pvaluecutoff, scorecut, simcut, lencut, evalue, mindetection, minseqcutoff, mincoexpseqcutoff, min_exp_cutoff, min_intersection_cutoff, alnmethod, treemethod, anno_file, araport_seq_file, iqtree, seqs_cluster_anno )
	sys.stdout.write( str( datetime.datetime.now() ) + " ... documentation file completed: " + docu_file + "\n" )
	sys.stdout.write( str( datetime.datetime.now() ) + " loading data ... \n" )
	config_data = load_config_file_content( config_file, output_folder )	#a species tag is included while loading the sequences to make them unique
	sys.stdout.write( str( datetime.datetime.now() ) + " ... loading data completed.\n" )
	sys.stdout.flush()
	
	if len( anno_file ) > 1:
		anno_mapping_table = load_annotation( anno_file )
	else:
		anno_mapping_table = {}
	
	### ---- BIG LOOP OVER ALL SPECIES --- ###
	huge_seq_collection_file = output_folder + "huge_seq_collection.cds.fasta"
	huge_seq_collection_file_pep = output_folder + "huge_seq_collection.pep.fasta"
	huge_seq_collection = {}

	# --- Preparation shell scripts --- #
	num_cores = os.cpu_count()
	if mode == "batch":
		coexp_cmd_file = output_folder + "COEXP_COMMANDS.sh"
		blast_dbs_file = output_folder + "BLAST_DBs.sh"
		blast_cmd_file = output_folder + "BLAST_COMMANDS.sh"
		tree_cmd_file = output_folder + "TREE_COMMANDS.sh"
		
		shell_header = "#!/bin/bash \n"
		batch_header = "nohup parallel -j " + str( batch ) + " <<EOF \n"
		batch_coexp = "nohup parallel -j " + str( int( batch/2 ) + 1 ) + " <<EOF \n"		#less jobs parallel due to high memory demand of coexpression analysis
		batch_tree = "nohup parallel -j " + str( num_cores - 1 ) + " <<EOF \n"
		with open( coexp_cmd_file, "w" ) as coexp_command_out:
			coexp_command_out.write( shell_header + batch_coexp )
		with open( blast_dbs_file, "w" ) as blast_dbs_out:
			blast_dbs_out.write( shell_header )			
		with open( blast_cmd_file, "w" ) as blast_command_out:
			blast_command_out.write( shell_header + batch_header )
		with open( tree_cmd_file, "w" ) as tree_command_out:
			tree_command_out.write( shell_header + batch_tree )
			
		# --- run co-expression per species --- #
		for spec in list( config_data.keys() ):
			info = config_data[ spec ]
			coexp_result_file = tmp_folder + info['id'] + "_coexp_results.txt"
			
			# --- write command for co-expression analysis and to build blast db --- #
			# modify gene names if there is a special character in gene name
			baits_modified = []
			for bait in ( info['baits'] ):
				baits_modified.append( bait.replace( "|", "\|" ).replace( "$", "\$" ).replace( "*", "\*" ).replace( ";", "\;" ).replace( "~", "\~" ).replace( "/", "\/" ) )

			with open( coexp_cmd_file, "a" ) as coexp_command_out:
				current_coexp_cmd = "".join( [ 	"python3 ", coexp_script_path, 
																	" --exp ", info['tpm'],
																	" --specid ", info['id'],
																	" --genes ", "@@@".join( baits_modified ),
																	" --rcutoff ", str( rcutoff ),
																	" --pvaluecutoff ", str( pvaluecutoff ),
																	" --coexpnumcutoff ", str( coexpnumcutoff ),
																	" --min_exp_cutoff ", str( min_exp_cutoff ),
																	" --mindetection ", str( mindetection ),
																	" --output ", coexp_result_file
																] )
				coexp_command_out.write( current_coexp_cmd + "\n" )
			
			blastdb = tmp_folder + spec + "_blastdb"
			with open( blast_dbs_file,  "a" ) as blast_dbs_out:
				current_blast_db_command = "diamond makedb --in " + config_data[ spec ]['pep_file'] + " -d " + blastdb
				blast_dbs_out.write( current_blast_db_command + "\n" + "sleep 2" + "\n" )				

		# --- execute coexpression and blast db scripts --- #			
		with open( coexp_cmd_file, "a" ) as coexp_command_out:
			coexp_command_out.write( "EOF" )
		os.chmod( coexp_cmd_file, 0o755 )
		subprocess.run( ['nohup', '/bin/bash', coexp_cmd_file], check = True )
		sys.stdout.write( str( datetime.datetime.now() ) + " - Coexpression analyses completed! " + "\n" )
		sys.stdout.flush()
		os.chmod( blast_dbs_file, 0o755 )
		subprocess.run( ['nohup', '/bin/bash', blast_dbs_file], check = True )
		sys.stdout.write( str( datetime.datetime.now() ) + " - All blast databases built! " + "\n" )
		sys.stdout.flush()		

	# --- loading of coexpression results --- #
	for spec in list( config_data.keys() ):
		info = config_data[ spec ]
		coexp_result_file = tmp_folder + info['id'] + "_coexp_results.txt"	
		with open( coexp_result_file, "r" ) as f:
			valid_coexp_genes = []
			for each in f.read().strip().split('\n'):
				if len( each ) > 1:
					valid_coexp_genes.append( each )
			config_data[ spec ].update( { 'coexp': valid_coexp_genes } )
			
		if os.path.isfile( coexp_result_file ):	#next step is only possible if co-expression results are available
			# --- run BLAST search of top X co-expressed genes against all other species --- #
			baits_file = tmp_folder + spec + ".baits.fasta"
			with open( baits_file, "w" ) as out:
				info = config_data[ spec ]
				for gene in info['coexp']:
					out.write( '>' + gene + "\n" + info['pep'][ gene ] + "\n" )
			
			for gene in valid_coexp_genes:	#add all sequences with co-expression to huge collection with a "_coexp" tag
				try:
					huge_seq_collection[ gene + "_coexp" ]
				except KeyError:
					huge_seq_collection.update( { gene + "_coexp": config_data[ spec ]['cds'][ gene ] } )	#only add new sequence if not present yet
			
			for spec2 in list( config_data.keys() ):
				if mode == "batch":
					if spec != spec2:
						blast_result_file = tmp_folder + spec + "_vs_" + spec2 + "blasthits.txt"	#run BLASTp search against each other species
						blastdb = tmp_folder + spec2 + "_blastdb"
						with open( blast_cmd_file, "a" ) as blast_command_out:
								current_search_cmd = "diamond blastp -q " + baits_file + " -d " + blastdb + " -o " + blast_result_file + " --threads " + str( cpub )
								blast_command_out.write( current_search_cmd + "\n" )
	
	# --- execute all blastp commands --- #
	if mode == "batch":
		with open( blast_cmd_file, "a" ) as blast_command_out:
			blast_command_out.write( "EOF" )
		sys.stdout.write( str( datetime.datetime.now() ) + " Starting DIAMOND blastp ... " + "\n" )
		os.chmod( blast_cmd_file, 0o755 )
		subprocess.run( ['/bin/bash', blast_cmd_file] )
		sys.stdout.write( str( datetime.datetime.now() ) + " ... DIAMOND blastp completed! " + "\n" )
		sys.stdout.flush()

	# --- add BLAST-based sequences to huge_seq_collection --- #
	for spec in list( config_data.keys() ):
		for spec2 in list( config_data.keys() ):
			if spec != spec2:
				blast_result_file = tmp_folder + spec + "_vs_" + spec2 + "blasthits.txt"			
				hits = load_blast_hits( blast_result_file, scorecut, simcut, lencut, evalue )	# lists of sequence IDs							
				for hit in hits:	#add all sequences from other species to huge collection
					try:
						huge_seq_collection[ hit + "_coexp"]	#check that this sequence was not coexpressed in another species
					except:
						try:
							huge_seq_collection[ hit ]
						except KeyError:
							huge_seq_collection.update( { hit: config_data[ spec2 ]['cds'][ hit ] } )
	
	# --- check if at least one non-coexpressed sequence is included in the collection --- #
	non_coexp_included = False
	for key in list( huge_seq_collection.keys() ):
		if not "_coexp" in key:
			non_coexp_included = True
			break
	if non_coexp_included:	#only continue if at BLAST-based sequence collection completed (at least one non-coexpressed sequence included)
		with open( huge_seq_collection_file, "w" ) as out:
			for key in huge_seq_collection.keys():
				out.write( '>' + key + '\n' + huge_seq_collection[ key ] + "\n" )
	
	### --- FINAL PART --- ###
	sys.stdout.write( str( datetime.datetime.now() ) + " - per species analyses completed.\n" )
	sys.stdout.flush()
	if os.path.isfile( huge_seq_collection_file ):
		# --- run BLAST all vs. all for huge sequence collection --- #
		huge_pep_collection = translate( huge_seq_collection )
		with open( huge_seq_collection_file_pep, "w" ) as out:
			for key in list( huge_pep_collection.keys() ):
				out.write( '>' + key + '\n' + huge_pep_collection[ key ] + "\n" )
		blastdb = tmp_folder + "huge_blastdb"
		blast_result_file = tmp_folder + "huge_blast_result_file.txt"
		if not os.path.isfile( blast_result_file ):
			sys.stdout.write( str( datetime.datetime.now() ) + " - starting DIAMOND BLAST for sequence clustering...\n" )
			sys.stdout.flush()
			p = subprocess.Popen( args= "diamond makedb --in " + huge_seq_collection_file_pep + " -d " + blastdb, shell=True )
			p.communicate()
			p = subprocess.Popen( args= "diamond blastp -q " + huge_seq_collection_file_pep + " -d " + blastdb + " -o " + blast_result_file + " --threads " + str( 27 ), shell=True )
			p.communicate()
			sys.stdout.write( str( datetime.datetime.now() ) + " ... DIAMOND blastp for sequence clustering completed.\n" )
			sys.stdout.flush()
		blast_results = load_hits_per_bait( blast_result_file, scorecut, simcut, lencut, evalue )
		
		tree_folder = output_folder + "trees/"
		if not os.path.exists( tree_folder ):
			os.makedirs( tree_folder )
		aln_folder = output_folder + "aln/"
		if not os.path.exists( aln_folder ):
			os.makedirs( aln_folder )
		
		
		# --- identify groups within the huge sequence collection --- #
		annotation_file = output_folder + "functional_annotation_of_clusters.txt"
		if not os.path.isfile( annotation_file ):
			sys.stdout.write( str( datetime.datetime.now() ) + " - constructing sequence clusters...\n" )
			sys.stdout.flush()
			get_groups( annotation_file, araport_seq_file, blast_results, anno_mapping_table, huge_pep_collection, aln_folder, mincoexpseqcutoff, minseqcutoff, tmp_folder, min_intersection_cutoff, cpub, batch, seqs_cluster_anno )
			sys.stdout.write( str( datetime.datetime.now() ) + " ... sequence clustering completed.\n" )
			sys.stdout.flush()
			
		# --- construct phylogenetic trees for each of them and highlight co-expressed genes --- #
		sys.stdout.write( str( datetime.datetime.now() ) + " - starting construction of phylogenetic trees...\n" )
		sys.stdout.flush()
		fasta_input_files = glob.glob( aln_folder + "*.pep.fasta" )
		for filename in fasta_input_files:
			name = filename.split('/')[-1].split('.')[0]
			tree_constructor( filename, treemethod, tree_folder, name, alnmethod, tree_cmd_file, mafft, muscle, raxml, fasttree, iqtree, cpur )
		
		with open( tree_cmd_file , "a" ) as tree_command_out:
			tree_command_out.write( "EOF" )
		os.chmod( tree_cmd_file, 0o755 )
		try:
			subprocess.run(['nohup', '/bin/bash', tree_cmd_file], check = True)
		except subprocess.CalledProcessError:
			sys.stdout.write( str( datetime.datetime.now() ) + " ERROR during Tree generation, please check tree-files ...\n" )
			sys.stdout.flush()		

		sys.stdout.write( str( datetime.datetime.now() ) + " ... construction of phylogenetic trees completed.\n" )
		sys.stdout.flush()
		
		# --- identify and define orthogroups --- #
		#find "coexp" sequences that are clustered = calculate distances between them and check for edges < number of species

if '--config' in sys.argv and '--out' in sys.argv and '--mode' in sys.argv:
	main( sys.argv )
elif '--version' in sys.argv:
	sys.exit( __version__ )
else:
	sys.exit( __usage__ )
