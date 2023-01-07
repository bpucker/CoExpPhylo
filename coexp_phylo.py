### Boas Pucker ###
### b.pucker@tu-braunschweig.de ###

### some functions taken from MYB_annotator.py and KIPEs ###

__version__ = "v0.11"

__usage__ = """
					python3 coexp_phylo.py
					--config <CONFIG_FILE>
					--out <OUTPUT_FOLDER>
					
					optional:
					--anno <ANNOTATION_FILE>
					--araport <ARAPORT_PEP_FILE>
					--r <R_CUTOFF_FOR_CORRELATION_ANALYSIS>
					--p <PVALUE_CUTOFF_FOR_CORRELATION_ANALYSIS>
					--numcut <NUMBER_OF_COEXPRESSED_GENES_CONSIDERED>
					--cpu <NUMBER_OF_CPUs_TO_USE>
					--scorecut <MIN_BLAST_SCORE_CUTOFF>
					--simcut <MIN_BLAST_SIMILARITY_CUTOFF>
					--lencut <MIN_BLAST_LENGTH_CUTOFF>
					--method <TREE_ALGORITHM>(fasttree|raxml|iqtree)[fasttree]
					--mafft <MAFFT_PATH>[mafft]
					--raxml <RAXML_PATH>[raxml]
					--fasttree <FASTTREE_PATH>[FastTree]
					--iqtree <IQ-TREE_PATH>[iqtree]
					--cpur <CPUs_FOR_TREE_CONSTRUCTION>[cpu]
					"""

import os, re, sys, subprocess, math, glob, time
from operator import itemgetter
import numpy as np
from scipy import stats

# --- end of imports --- #

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
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: "".join( seq ) } )
					header = line.strip()[1:]
					if " " in header:
						header = header.split(' ')[0]
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
		tissues = f.readline().strip().split('\t')[1:]	#list of accessions (first field can be empty string i.e. invisible in file)
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			expression = {}
			for idx, each in enumerate( parts[1:] ):
				expression.update( { tissues[  idx ] : float( parts[ idx+1 ] ) } )	#convert value into float when adding to dictionary
			line = f.readline()
			expression_data.update( { spec_prefix + "@" + parts[0]: expression } )
	return expression_data


def compare_candidates_against_all( candidate, gene_expression, rcutoff, pvaluecutoff, coexpnumcutoff ):
	"""! @brief compare candidate gene expression against all genes to find co-expressed genes
		@return list of co-expressed genes that passed the filters
	"""
	
	tissues = sorted( list( gene_expression[ list( gene_expression.keys() )[0] ].keys() ) )	#accessions / sample names
	coexpressed_genes = []
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
			if total_expression > 30:	#only consider genes with a minimal expression level
				if len( values ) > 0:
					r, p = stats.spearmanr( values )
					if not math.isnan( r ):
						if r > rcutoff and p < pvaluecutoff:
							coexpressed_genes.append( { 'id': gene2, 'correlation': r, 'p_value': p } )
				else:
					print ( "WARNING: no expression detected - " + candidate + "\t" + gene2 + "\t" + ";".join( list( gene_expression.keys() )[:5] ) )
	coexp_gene_IDs = []
	for gene in sorted( coexpressed_genes, key=itemgetter('correlation') )[::-1][:coexpnumcutoff]:
		coexp_gene_IDs.append( gene['id'] )
	return coexp_gene_IDs


def load_blast_hits( blast_result_file, scorecut, simcut, lencut):
	""""! @brief load BLAST hits above a certain score per bait
		@return list of sequence IDs that are considered good BLAST hits
	"""
	
	hits = []
	with open( blast_result_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if float( parts[-1] ) > scorecut:
				if float( parts[2] ) > simcut:
					if int( parts[3] ) > lencut:
						hits.append( parts[1] )
			line = f.readline()
	return list( set( hits ) )	#reduce to list of unique elements


def load_hits_per_bait( blast_result_file, scorecut, simcut, lencut ):
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


def tree_constructor( X_aln_input_file, method, X_output_folder, Xname, mafft, raxml, fasttree, iqtree, cpur ):
	"""! @brief handles the construction of alignments and phylogenetic tree
			@note second FASTA file can be an empty string to run this function just based on one FASTA file
	"""
	
	X_aln_file = X_aln_input_file + ".aln"
	X_cln_aln_file = X_aln_file + ".cln"
	
	# --- generate MAFFT alignment --- #
	if not os.path.isfile( X_aln_file ):
		p = subprocess.Popen( args= mafft + " --quiet " + X_aln_input_file + " > " + X_aln_file, shell=True )
		p.communicate()
	
	# --- trim alignment i.e. remove positions with many gaps = low occupancy --- #
	if not os.path.isfile( X_cln_aln_file ):
		alignment_trimming( X_aln_file, X_cln_aln_file, occupancy=0.1 )
	
	if method == "raxml":	#construct tree with RAxML
		prefix = X_output_folder + Xname + "RAxML_tree"
		tree_file = prefix + ".raxml.bestTree.tre"
		if not os.path.isfile( tree_file ):
			p = subprocess.Popen( args= " ".join( [ raxml, "--all --threads " + str( cpur ) + " --model LG+G8+F --msa", X_cln_aln_file, "--prefix", prefix ] ), shell=True )
			p.communicate()
	
	elif method == "iqtree":	#construct tree with IQ-TREE
		tree_file = X_cln_aln_file + ".treefile"	#".treefile" is appended to provided alignment file name
		if not os.path.isfile( tree_file ):
			p = subprocess.Popen( args= " ".join( [ iqtree, + " -nt " + str( cpur ) + " -alrt 1000 -bb 100 -s", X_cln_aln_file, "--prefix", prefix ] ), shell=True )
			p.communicate()
	
	else:	#construct tree with FastTree2
		tree_file = X_output_folder  + Xname  + "FastTree_tree.tre"
		if not os.path.isfile( tree_file ):
			p = subprocess.Popen( args= " ".join( [ fasttree, "-wag  -nopr -nosupport <", X_cln_aln_file, ">", tree_file ] ), shell=True )
			p.communicate()
	return tree_file


def count_exp_seqs( seq_IDs_to_check ):
	"""! @brief count co-expressed sequences """
	
	counter = 0
	for each in seq_IDs_to_check:
		if "_coexp" in each:
			counter += 1
	return counter


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


def annotate_group( group_sequence_file, tmp_blast_folder, blast_db, anno_mapping_table ):
	"""! @brief annotate given group of sequences based on best BLAST hit against Araport11
	@return string comprizing functional annotation if available
	"""
	
	blast_result_file = tmp_blast_folder + group_sequence_file.split('/')[-1] + "_blast_results.txt"
	p = subprocess.Popen( args= "blastp -query " + group_sequence_file+ " -db " + blast_db + " -out " + blast_result_file + " -outfmt 6 -evalue 0.001 ", shell=True )
	p.communicate()
	
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


def get_groups( annotation_file, araport_seq_file, blast_results, anno_mapping_table, huge_pep_collection, aln_folder, mincoexpseqcutoff, minseqcutoff, tmp_folder, min_intersection_cutoff ):
	"""! @brief identify groups in huge cluster """
	
	# --- identify groups --- #
	final_groups = []
	for key in list( blast_results.keys() ):
		new_group = blast_results[ key ] + [ key ]
		best_index = False
		best_intersection = min_intersection_cutoff + 0
		for group in final_groups:	#iterate over all 'final groups'
			intersection = len( list( set.intersection( set( new_group ), set( group ) ) ) )	#evaluate intersection of current group with all existing groups
			if intersection > best_intersection:
				best_intersection = intersection + 0	#any overlap with another group leads to merging
				best_index = final_groups.index( group )
		if best_index:
			final_groups[ best_index ] = list( set( new_group + final_groups[ best_index ] ) )	#merge with an existing group if overlap exists
		else:
			final_groups.append( new_group )
	
	# --- prepare database for annotation --- #
	if len( araport_seq_file ) > 0:
		tmp_blast_folder = tmp_folder + "anno_blast/"
		if not os.path.exists( tmp_blast_folder ):
			os.makedirs( tmp_blast_folder )
		blast_db = tmp_blast_folder + "blastdb"
		p = subprocess.Popen( args= "makeblastdb -in " + araport_seq_file + " -out " + blast_db + " -dbtype prot", shell=True )
		p.communicate()
	
	# --- write sequences into output files --- #
	with open( annotation_file, "w" ) as anno_out:
		counter = 0
		for fg in final_groups:
			print ( len( fg ) )
			if len( fg ) > minseqcutoff:
				coexpseqcounter = count_exp_seqs( fg )
				print ( "coexp: " + str( coexpseqcounter ) )
				if coexpseqcounter > mincoexpseqcutoff:
					counter += 1
					group_sequence_file = aln_folder + str( counter ).zfill( 3 ) + ".pep.fasta"	#save all sequenes of a group in a FASTA file
					with open( group_sequence_file, "w" ) as out:
						for gene in fg:
							out.write( '>' + gene + "\n" + huge_pep_collection[ gene ] + "\n" )
					if len( araport_seq_file ) > 0:
						group_annotation = annotate_group( group_sequence_file, tmp_blast_folder, blast_db, anno_mapping_table )
					else:
						group_annotation = "n/a"
					anno_out.write( str( counter ).zfill(3) + "\t" + str( len( fg ) ) + "\t" + group_annotation + "\n" )


def main( arguments ):
	"""! @brief run everything """
	
	config_file = arguments[ arguments.index('--config')+1 ]
	#ID,tpm,cds,baits
	#ID = species ID; tpm= tpm file; cds = CDS file; baits = IDs of interest in file
	output_folder = arguments[ arguments.index('--out')+1 ]
	
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
	if '--cpu' in arguments:
		cpu = int( arguments[ arguments.index('--cpu')+1 ] )
	else:
		cpu = 4
	if '--scorecut' in arguments:	#score cutoff for BLAST hits
		scorecut = float( arguments[ arguments.index('--scorecut')+1 ] )
	else:
		scorecut = 100.0
	if '--simcut' in arguments:	#similarity cutoff for BLAST hits
		simcut = float( arguments[ arguments.index('--simcut')+1 ] )
	else:
		simcut = 60.0
	if '--lencut' in arguments:	#length cutoff for BLAST hits
		lencut = float( arguments[ arguments.index('--lencut')+1 ] )
	else:
		lencut = 100
	
	mindetection = 1	#number of bait genes that a given sequence need to be co-expressed with to be considered (strict would be equal to number of baits)
	minseqcutoff = 10	#minimal number of sequences to compose a group as tree construction input
	mincoexpseqcutoff = 3	#minimal number of co-expressed sequences to compose a group as tree construction input
	
	if '--mafft' in arguments:	#path to MAFFT (alignment)
		mafft = arguments[ arguments.index('--mafft')+1 ]
	else:
		mafft = "mafft"
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
		iqtree = "iqtree"
	
	if '--cpur' in arguments:	#CPUs to use for tree
		cpur = int( arguments[ arguments.index('--cpur')+1 ] )
	else:
		cpur = min( [ 4, cpu ] )
	
	if '--method' in arguments:	#tree construction methods
		method = arguments[ arguments.index('--method')+1 ]
		if method not in [ "fasttree", "raxml", "iqtree" ]:
			method = "fasttree"
	else:
		method = "fasttree"
	
		
	if '--araport' in arguments:
		araport_seq_file = arguments[ arguments.index('--araport')+1 ]
	else:
		araport_seq_file = ""
	
	min_intersection_cutoff = 0	#minimal intersection between two sequence clusters that needs to be exceeded to trigger merging
	
	
	if output_folder[-1] != "/":
		output_folder += "/"
	
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	tmp_folder = output_folder + "tmp/"
	if not os.path.exists( tmp_folder ):
		os.makedirs( tmp_folder )
	
	config_data = load_config_file_content( config_file, output_folder )	#a species tag is included while loading the sequences to make them unique
	
	if len( anno_file ) > 1:
		anno_mapping_table = load_annotation( anno_file )
	else:
		anno_mapping_table = {}
	
	huge_seq_collection_file = output_folder + "huge_seq_collection.cds.fasta"
	huge_seq_collection_file_pep = output_folder + "huge_seq_collection.pep.fasta"
	if not os.path.isfile( huge_seq_collection_file ):
		huge_seq_collection = {}
		# --- run co-expression per species --- #
		for spec in list( config_data.keys() ):
			info = config_data[ spec ]
			exp = load_expression_values( info['tpm'], info['id'] )
			combined_coexp_genes = []
			for gene in info['baits']:
				combined_coexp_genes += compare_candidates_against_all( gene, exp, rcutoff, pvaluecutoff, coexpnumcutoff )
			unique_combined_coexp_genes = list( set( combined_coexp_genes ) )
			valid_coexp_genes = []
			for gene in unique_combined_coexp_genes:
				if combined_coexp_genes.count( gene ) >= mindetection:	#this slects genes co-expressed with entire pathway
					valid_coexp_genes.append( gene )
			config_data[ spec ].update( { 'coexp': valid_coexp_genes } )
		
			# --- run BLAST search of top X co-expressed genes against all other species --- #
			baits_file = tmp_folder + spec + ".baits.fasta"
			with open( baits_file, "w" ) as out:
				info = config_data[ spec ]
				for gene in info['coexp']:
					out.write( '>' + gene + "\n" + info['pep'][ gene ] + "\n"  )
			
			for gene in valid_coexp_genes:	#add all sequences with co-expression to huge collection
				try:
					del huge_seq_collection[ gene ]	#delete sequence entry that does not show coexpression (in other species)
				except:
					pass
				try:
					huge_seq_collection[ gene + "_coexp" ]
				except KeyError:
					huge_seq_collection.update( { gene + "_coexp": config_data[ spec ]['cds'][ gene ] } )	#only add new sequence if not present yet
			
			for spec2 in list( config_data.keys() ):
				if spec != spec2:
					blast_result_file = tmp_folder + spec + "_vs_" + spec2 + "blasthits.txt"	#run BLASTp search against each other species
					blastdb = tmp_folder + spec2 + "_blastdb"
					if not os.path.isfile( blast_result_file ):
						p = subprocess.Popen( args= "makeblastdb -in " + config_data[ spec2 ]['pep_file'] + " -out " + blastdb + " -dbtype prot", shell=True )
						p.communicate()
						p = subprocess.Popen( args= "blastp -query " + baits_file+ " -db " + blastdb + " -out " + blast_result_file + " -outfmt 6 -evalue 0.001 -num_threads " + str( cpu ), shell=True )
						p.communicate()
					hits = load_blast_hits( blast_result_file, scorecut, simcut, lencut )	# lists of sequence IDs
					
					for hit in hits:	#add all sequences from other species to huge collection
						try:
							huge_seq_collection[ hit + "_coexp"]	#check that this sequence was not coexpressed in another species
						except:
							try:
								huge_seq_collection[ hit ]
							except KeyError:
								huge_seq_collection.update( { hit: config_data[ spec2 ]['cds'][ hit ] } )
		with open( huge_seq_collection_file, "w" ) as out:
			for key in huge_seq_collection.keys():
				out.write( '>' + key + '\n' + huge_seq_collection[ key ] + "\n" )
	else:
		huge_seq_collection = load_sequences( huge_seq_collection_file )
	
	# --- run BLAST all vs. all for huge sequence collection --- #
	huge_pep_collection = translate( huge_seq_collection )
	with open( huge_seq_collection_file_pep, "w" ) as out:
		for key in list( huge_pep_collection.keys() ):
			out.write( '>' + key + '\n' + huge_pep_collection[ key ] + "\n" )
	blastdb = tmp_folder + "huge_blastdb"
	blast_result_file = tmp_folder + "huge_blast_result_file.txt"
	if not os.path.isfile( blast_result_file ):
		p = subprocess.Popen( args= "makeblastdb -in " + huge_seq_collection_file_pep + " -out " + blastdb + " -dbtype prot", shell=True )
		p.communicate()
		p = subprocess.Popen( args= "blastp -query " + huge_seq_collection_file_pep+ " -db " + blastdb + " -out " + blast_result_file + " -outfmt 6 -evalue 0.001 -num_threads " + str( cpu ), shell=True )
		p.communicate()
	blast_results = load_hits_per_bait( blast_result_file, scorecut, simcut, lencut )
	
	tree_folder = output_folder + "trees/"
	if not os.path.exists( tree_folder ):
		os.makedirs( tree_folder )
	aln_folder = output_folder + "aln/"
	if not os.path.exists( aln_folder ):
		os.makedirs( aln_folder )
	
	# --- identify groups within the huge sequence collection --- #
	annotation_file = output_folder + "functional_annotation_of_clusters.txt"
	if not os.path.isfile( annotation_file ):
		get_groups( annotation_file, araport_seq_file, blast_results, anno_mapping_table, huge_pep_collection, aln_folder, mincoexpseqcutoff, minseqcutoff, tmp_folder, min_intersection_cutoff )

	# --- construct phylogenetic trees for each of them and highlight co-expressed genes --- #
	fasta_input_files = glob.glob( aln_folder + "*.pep.fasta" )
	for filename in fasta_input_files:
		name = filename.split('/')[-1].split('.')[0]
		tree_constructor( filename, method, tree_folder, name, mafft, raxml, fasttree, iqtree, cpur )
	
	
	# --- identify and define orthogroups --- #
	#find "coexp" sequences that are clustered = calculate distances between them and check for edges < number of species



if '--config' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
