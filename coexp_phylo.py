### Boas Pucker ###
### b.pucker@tu-braunschweig.de ###

### some functions taken from MYB_annotator.py and KIPEs ###

__version__ = "v0.031"

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
					--mode <TREE_ALGORITHM>(fasttree|raxml)[fasttree]
					--mafft <MAFFT_PATH>
					--raxml <RAXML_PATH>
					--fasttree <FASTTREE_PATH>
					--cpur <CPUs_FOR_TREE_CONSTRUCTION>
					"""


import os, re, sys, subprocess, math, glob
from operator import itemgetter
import numpy as np
from scipy import stats

# --- end of imports --- #

def translate( seqs ):
	"""! @brief translates the given nucleotide sequence into peptide and splits at each star (stop codon) """
	
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
	for key in seqs.keys():
		seq = seqs[ key ].upper()
		peptide = []
		for i in range( int( len( seq ) / 3.0 ) ):
			codon = seq[i*3:i*3+3]
			try:
				peptide.append( genetic_code[ codon ] )
			except:
				peptide.append( "*" )
		final_peptide_seqs.update( { key: "".join( peptide ) } )
	return final_peptide_seqs


def load_sequences( fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	with open( fasta_file ) as f:
		header = f.readline()[1:].strip()
		if " " in header:
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
	"""! @brief load all bait IDs """
	
	IDs = []
	with open( baitfile, "r" ) as f:
		line = f.readline()
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
	"""! @brief load content of config file """
	
	infos = {}
	with open( config_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split(',')
				pep_file = output_folder + parts[0] + ".pep.fasta"
				cds_file = output_folder + parts[0] + ".cds.fasta"
				cds_seqs = load_sequences( parts[2] )
				find_cds_seqs = {}
				with open( cds_file, "w" ) as out:
					for key in list( cds_seqs.keys() ):
						out.write( '>' + parts[0] + "@" + key + '\n' + cds_seqs[ key ] + "\n" )
						find_cds_seqs.update( { parts[0] + "@" + key: cds_seqs[ key ] } )
				pep_seqs = translate( find_cds_seqs )
				with open( pep_file, "w" ) as out:
					for key in list( pep_seqs.keys() ):
						out.write( '>' + key + '\n' + pep_seqs[ key ] + "\n" )
				baits = load_baits( parts[3], parts[0] )
				infos.update( { parts[0]: { 'id': parts[0], 'tpm': parts[1], 'cds_file': cds_file, 'baits_file': parts[3], 'pep_file': pep_file, 'cds': find_cds_seqs, 'pep': pep_seqs, 'baits': baits } } )
			line = f.readline()
	return infos


def load_expression_values( filename, spec_prefix ):
	"""! @brief load all expression values """
	
	expression_data = {}
	with open( filename, "r" ) as f:
		tissues = f.readline().strip().split('\t')[1:]
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			expression = {}
			for idx, each in enumerate( parts[1:] ):
				expression.update( { tissues[  idx ] : float( parts[ idx+1 ] ) } )
			line = f.readline()
			expression_data.update( { spec_prefix + "@" + parts[0]: expression } )
	return expression_data


def compare_candidates_against_all( candidate, gene_expression, rcutoff, pvaluecutoff, coexpnumcutoff ):
	"""! @brief compare candidate gene expression against all genes to find co-expressed genes """
	
	tissues = sorted( list( gene_expression[ list( gene_expression.keys() )[0] ].keys() ) )
	coexpressed_genes = []
	for i, gene2 in enumerate( list( gene_expression.keys() ) ):
		if candidate != gene2:
			values = []
			total_expression = 0
			for tissue in tissues:
				try:
					x = gene_expression[ candidate ][ tissue ]
					y = gene_expression[ gene2 ][ tissue ]
					total_expression += y
					if not math.isnan( x ) and not math.isnan( y ) :
						values.append( [ x, y ] )
				except KeyError:
					pass
			if total_expression > 30:
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
	""""! @brief load BLAST hits above a certain score per bait """
	
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
	return list( set( hits ) )


def load_hits_per_bait( blast_result_file, scorecut, simcut, lencut ):
	"""! @brief load BLAST hits per bait """
	
	hits = {}
	with open( blast_result_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if float( parts[-1] ) > scorecut:
				if float( parts[2] ) > simcut:
					if int( parts[3] ) > lencut:
						try:
							hits[ parts[0] ].append( parts[1] )
						except KeyError:
							hits.update( { parts[0]: [ parts[1] ] } )
			line = f.readline()
	return hits


def load_alignment( aln_file, tmp_mapping ):
	"""! @brief load alignment and replace query IDs by real sequence names """
	
	sequences = {}
	
	with open( aln_file ) as f:
		header = f.readline()[1:].strip()
		try:
			header = tmp_mapping[ header ]
		except KeyError:
			pass
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: "".join( seq ) } )
					header = line.strip()[1:]
					try:
						header = tmp_mapping[ header ]
					except KeyError:
						pass
					seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		sequences.update( { header: "".join( seq ) } )	
	return sequences


def alignment_trimming( aln_file, cln_aln_file, occupancy ):
	"""! @brief remove all alignment columns with insufficient occupancy """
	
	alignment = load_alignment( aln_file, {} )
	# --- if there is an alignment (expected case) 
	if len( list(alignment.keys()) ) > 0:
		# --- identify valid residues in aligned sequences (columns with sufficient occupancy) --- #
		valid_index = []
		for idx, aa in enumerate( list(alignment.values())[0] ):
			counter = 0
			for key in list(alignment.keys()):
				if alignment[ key ][ idx ] != "-":
					counter += 1
			if counter / float( len( list(alignment.keys()) ) ) > occupancy:
				valid_index.append( idx )
		
		# --- generate new sequences --- #
		with open( cln_aln_file, "w" ) as out:
			for key in list(alignment.keys()):
				seq = alignment[ key ]
				new_seq = []
				for idx in valid_index:
					new_seq.append( seq[ idx ] )
				out.write( ">" + key + '\n' + "".join( new_seq ) + '\n' )
	# --- just in case the alignment file is empyt (is this possible?) ---#
	else:
		with open( cln_aln_file, "w" ) as out:
			out.write( "" )


def tree_constructor( X_aln_input_file, mode, X_output_folder, Xname, mafft, raxml, fasttree, cpur ):
	"""! @brief handles the construction of alignments and phylogenetic tree
			@note second FASTA file can be an empty string to run this function just based on one FASTA file
	"""
	
	X_aln_file = X_aln_input_file + ".aln"
	X_cln_aln_file = X_aln_file + ".cln"
	
	if not os.path.isfile( X_aln_file ):
		p = subprocess.Popen( args= mafft + " --quiet " + X_aln_input_file + " > " + X_aln_file, shell=True )
		p.communicate()
	
	if not os.path.isfile( X_cln_aln_file ):
		alignment_trimming( X_aln_file, X_cln_aln_file, occupancy=0.1 )
	
	if mode == "raxml":	#RAxML
		prefix = X_output_folder + Xname + "RAxML_tree"
		tree_file = prefix + ".raxml.bestTree.tre"
		if not os.path.isfile( tree_file ):
			p = subprocess.Popen( args= " ".join( [ raxml, "--all --threads " + str( cpur ) + " --model LG+G8+F --msa", X_cln_aln_file, "--prefix", prefix ] ), shell=True )
			p.communicate()
	else:	#FastTree2
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
	"""! @brief load annotation from given file """
	
	mapping_table = {}
	with open( anno_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			mapping_table.update( { parts[0]: ". ".join( parts[1:] ) } )
			line = f.readline()
	return mapping_table


#THIS FUNCTION IS DEPRECATED
def get_groups1( annotation_file, blast_results, anno_mapping_table, huge_pep_collection, aln_folder, mincoexpseqcutoff, minseqcutoff ):
	"""! @brief identify groups  """
	
	with open( annotation_file, "w" ) as anno_out:
		counter = 0
		while counter < len( list( blast_results.keys() ) ):
			first_key = list( blast_results.keys() )[0]
			if len( blast_results[ first_key ] ) > minseqcutoff:
				coexpseqcounter = count_exp_seqs( blast_results[ first_key ] )
				AGIs = re.findall( "AT\dG\d+", "".join( blast_results[ first_key ] ) )
				if len( AGIs ) > 0:
					anno = []
					for AGI in AGIs:
						try:
							anno.append( anno_mapping_table[ AGI ] )
						except KeyError:
							pass
					anno_out.write( str( counter ).zfill(3) + "\t" + ";".join( AGIs ) + "\t" + "; ".join( anno ) + "\n" )
				else:
					anno_out.wirte( str( counter ).zfill(3) + "\t" + "n/a\tn/a\n" )
				if coexpseqcounter >= mincoexpseqcutoff:
					output_file = aln_folder + str( counter ).zfill( 3 ) + ".pep.fasta"	#save all sequenes of a group in a FASTA file
					with open( output_file, "w" ) as out:
						for gene in blast_results[ first_key ]:
							out.write( '>' + gene + "\n" + huge_pep_collection[ gene ] + "\n" )
							try:
								del blast_results[ gene ]
							except KeyError:
								pass
			try:
				del blast_results[ first_key ]
			except KeyError:
				pass
			if len( list( blast_results.keys() ) ) == 0:	#exit loop once all sequences are assigned to a group
				break
			counter += 1
		print ( "number of groups: " + str( counter ) )


def annotated_group( group_sequence_file, tmp_blast_folder, blast_db, anno_mapping_table ):
	"""! @brief annotate group based on best BLAST hit against Araport11 """
	
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


def get_groups2( annotation_file, araport_seq_file, blast_results, anno_mapping_table, huge_pep_collection, aln_folder, mincoexpseqcutoff, minseqcutoff, tmp_folder ):
	"""! @brief identify groups in huge cluster """
	
	# --- identify groups --- #
	final_groups = []
	for key in list( blast_results.keys() ):
		new_group = blast_results[ key ] + [ key ]
		best_index = False
		best_intersection = 0
		for group in final_groups:
			intersection = len( list( set.intersection( set( new_group ), set( group ) ) ) )
			if intersection > best_intersection:
				best_intersection = intersection + 0
				best_index = final_groups.index( group )
		if best_index:
			final_groups[ best_index ] = list( set( new_group + final_groups[ best_index ] ) )
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
						group_annotation = annotated_group( group_sequence_file, tmp_blast_folder, blast_db, anno_mapping_table )
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
	if '--scorecut' in arguments:
		scorecut = float( arguments[ arguments.index('--scorecut')+1 ] )
	else:
		scorecut = 100.0
	if '--simcut' in arguments:
		simcut = float( arguments[ arguments.index('--simcut')+1 ] )
	else:
		simcut = 60.0
	if '--lencut' in arguments:
		lencut = float( arguments[ arguments.index('--lencut')+1 ] )
	else:
		lencut = 100
	
	mindetection = 1	#number of bait genes that a given sequence need to be co-expressed with to be considered (strict would be equal to number of baits)
	minseqcutoff = 10	#minimal number of sequences to compose a group as tree construction input
	mincoexpseqcutoff = 3	#minimal number of co-expressed sequences to compose a group as tree construction input
	
	if '--mafft' in arguments:
		mafft = arguments[ arguments.index('--mafft')+1 ]
	else:
		mafft = "mafft"
	if '--raxml' in arguments:
		raxml = arguments[ arguments.index('--raxml')+1 ]
	else:
		raxml = "raxml"
	if '--fasttree' in arguments:
		fasttree = arguments[ arguments.index('--fasttree')+1 ]
	else:
		fasttree = "FastTree"
	if '--cpur' in arguments:
		cpur = int( arguments[ arguments.index('--cpur')+1 ] )
	else:
		cpur = min( [ 4, cpu ] )
	if '--mode' in arguments:
		mode = arguments[ arguments.index('--mode')+1 ]
		if mode not in [ "fasttree", "raxml" ]:
			mode = "fasttree"
	else:
		mode = "fasttree"
	
	if '--araport' in arguments:
		araport_seq_file = arguments[ arguments.index('--araport')+1 ]
	else:
		araport_seq_file = ""
	
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
		#get_groups1( annotation_file, blast_results, anno_mapping_table, huge_pep_collection, aln_folder, mincoexpseqcutoff, minseqcutoff )
		get_groups2( annotation_file, araport_seq_file, blast_results, anno_mapping_table, huge_pep_collection, aln_folder, mincoexpseqcutoff, minseqcutoff, tmp_folder )

	# --- construct phylogenetic trees for each of them and highlight co-expressed genes --- #
	fasta_input_files = glob.glob( aln_folder + "*.pep.fasta" )
	for filename in fasta_input_files:
		name = filename.split('/')[-1].split('.')[0]
		tree_constructor( filename, mode, tree_folder, name, mafft, raxml, fasttree, cpur )
	
	
	# --- identify and define orthogroups --- #
	#find "coexp" sequences that are clustered = calculate distances between them and check for edges < number of species



if '--config' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
