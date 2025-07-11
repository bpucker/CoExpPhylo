### Nele Fiene, Boas Pucker ###
### b.pucker@tu-braunschweig.de ###

### some functions taken from MYB_annotator.py and KIPEs ###
### some code blocks were generated using GPT 3.5 and GPT 4o ###

__version__ = "v1.30"

__usage__ = """
				python3 coexp_phylo.py
				--config <CONFIG_FILE>
				--out <OUTPUT_FOLDER>
				
				optional:
				----
				ANNOTATION
				----
				--anno <ANNOTATION_FILE>
				--reference <REFERENCE_PEP_FILE>
				--seqs_cluster_anno <PERCENTAGE_OF_SEQUENCES_PER_CLUSTER_USED_FOR_ANNOTATION>[50.0]

				----
				PER-SPECIES COEXPRESION ANALYSIS
				----
				--r <R_CUTOFF_FOR_CORRELATION_ANALYSIS>[0.7]
				--p <PVALUE_CUTOFF_FOR_CORRELATION_ANALYSIS>[0.05]
				--numcut <NUMBER_OF_COEXPRESSED_GENES_CONSIDERED>[100]
				--min_exp_cutoff <MIN_EXP_FOR_GENE_TO_CONSIDER>[30
				--mindetect <MIN_NUMBER_OF_COEXPED_BAITS>[1]

				----
				DIAMOND
				----
				--cpub <NUMBER_OF_CPUs_TO_USE_FOR_BLAST>[4]
				--batch <NUMBER_OF_BLASTP-JOBS_TO_RUN_IN_PARALLEL>[7]
				--scorecut <MIN_BLAST_SCORE_CUTOFF>[100]
				--simcut <MIN_BLAST_SIMILARITY_CUTOFF>[80]
				--lencut <MIN_BLAST_LENGTH_CUTOFF>[100]
				--evalue <E-VALUE_CUTOFF>[1e-5]

				----
				CLUSTERING
				----
				--minseqcutoff <MIN_SEQ_CLUSTER_SIZE>[10]
				--mincoexpseqcutoff <MIN_SPECIES_WITH_COEXPED_GENES_FOR_CLUSTER>[3]
					
				----
				ALIGNMENT
				----
				--alnmethod <ALIGNMENT_ALGORITH>(mafft|muscle)[mafft]
				--mafft <MAFFT_PATH>[mafft]
				--muscle <MUSCLE_PATH>[muscle]
				--occupancy <OCCUPANCY_VALUE_FOR_ALIGNMENT>[0.1]

				----
				TREE CONSTRUCTION
				----
				--treemethod <TREE_ALGORITHM>(fasttree|raxml|iqtree)[fasttree]
				--raxml <RAXML_PATH>[raxml]
				--fasttree <FASTTREE_PATH>[FastTree]
				--iqtree <IQ-TREE_PATH>[iqtree]
				--clean <to clean gene ids for correct tree visualization>

				----
				BATCH UPLOAD ITOL
				----
				if batch upload is wanted, the following arguments are mandatory
				--API <API-KEY>
				--proj_name <PROJECT_NAME>(project must already exist in your iTOL account)
				optional:
				--upload_script <PATH_TO_iTOL_UPLOADER.PL>
					
				"""

import os, re, sys, subprocess, math, glob
from operator import itemgetter
import signal
import csv
import hashlib
from scipy import stats
import datetime
import networkx as nx
import random
import plotly.graph_objects as go
import numpy as np
from collections import Counter

# --- end of imports --- #

def check_tools( name, docu_file ):
	# check whether all tools and packages are available
	try:
		result = subprocess.run( [name, '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True )
		version_info = result.stdout.strip() or result.stderr.strip()
		with open( docu_file, 'a') as f:
			f.write( str( name + ': ' + version_info.splitlines()[0] + '\n--\n' ) )
		return True
	except subprocess.CalledProcessError:
		try:
			result = subprocess.run( [name, '-version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True )
			version_info = result.stdout.strip() or result.stderr.strip()
			with open( docu_file, 'a') as f:
				f.write( str( name + ': ' + version_info.splitlines()[0] + '\n--\n' ) )
			return True
		except subprocess.CalledProcessError:
			try:
				result = subprocess.run( [name, '-help'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True )
				version_info = result.stdout.strip() or result.stderr.strip()
				with open( docu_file, 'a') as f:
					f.write( str( name + ': ' + version_info.splitlines()[0] + '\n--\n' ) )
				return True
			except subprocess.CalledProcessError:
				try:
					subprocess.run( [name], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True )
					with open( docu_file, 'a') as f:
						f.write( str( name + ': no version number available \n--\n' ) )
					return True
				except subprocess.CalledProcessError:
					subprocess.run( [name, 'echo', 'check'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True )
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

def generate_docu_file( tmp_folder, docu_file, config_file, rcutoff, pvaluecutoff, scorecut, simcut, lencut, evalue, mindetection, minseqcutoff, mincoexpseqcutoff, min_exp_cutoff, alnmethod, treemethod, anno_file, ref_seq_file, muscle, mafft, raxml, iqtree, fasttree, seqs_cluster_anno, clean, api, proj_name, perl_path ):
	
	# --- Upload test tree to iTOL --- #
	if api != '':
		tree_folder = tmp_folder + "test_trees/"
		if not os.path.exists( tree_folder ):
			os.makedirs( tree_folder )
		test_tree = tree_folder + "upload_test.tree"
		with open( test_tree, 'w' ) as tree_file:
			tree_file.write( "((7h15:0.1,(Just ignore this tree:0.1),(15:0.09,(4:0.08,(7r33:0.07,(70:0.06,(7357:0.05,(7h3:0.04,(upl04d:0.03,(170l:0.02,(1n70:0.01):0.01):0.01):0.01):0.01):0.01):0.01):0.01):0.01):0.01):0.01);\n" )
		
		zip_folder = str(tree_folder + 'upload_test' + '.zip')
		subprocess.run(['zip', '-j', zip_folder, test_tree], stdout=subprocess.DEVNULL, check = True)
		# --- generate config file --- #
		upload_config = tree_folder + "upload_test" + "_iTOL.cfg"
		with open( upload_config, 'w' ) as cfg:
			cfg.write( "zipFile=" + zip_folder + "\nprojectName=" + proj_name + "\nAPIkey=" + api + "\n"+ "treeDescription=Test Upload :)\n")

		results = subprocess.run([perl_path, upload_config], capture_output = True, text = True)
		if 'error' in results.stdout.lower():
			sys.exit("ERROR with upload to iTOL: " + results.stdout + "\n" )
		if results.stderr:
			sys.exit("ERROR with upload to iTOL: " + results.stderr + "\n" )

	with open( docu_file, 'w') as f:
		f.write('--------------\nDOCUMENTATION\n--------------\n\nSCRIPT VERSION\n--------------\n' )
		f.write( str ( __version__ ) + '\n\nCONFIG FILE\n--------------\n')
		f.write( str ( config_file + '\t\t' + str( calculate_md5sum( config_file ) ) + '\n\n' ) )
		f.write( 'PARAMETER\n--------------\n' )
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
		f.write( str ( 'alignment method = ' + alnmethod + '\n' ) )
		f.write( str ( 'tree method = ' + treemethod + '\n' ) )
		f.write( str ( 'anno = ' + anno_file + '\n' ) )
		f.write( str ( 'reference peptide file = ' + ref_seq_file + '\n' ) )
		f.write( str ( 'seqs per cluster used for anno [%] = ' + str( seqs_cluster_anno ) + '\n' ) )
		if api != '':
			f.write( str ( 'API key for iTOL = ' + api + '\n' ) )
			f.write( str ( 'project name for iTOL upload = ' + proj_name + '\n' ) )		
		if clean:
			f.write( 'Warning: Due to activated clean option, Gene IDs might be changed (\':\' --> \'_\' to avoid incorrect representation in trees)\n')
		f.write('\nTOOLS\n--------------\n')

	missing_modules = []
	if api != '':
		required_modules = ['strict', 'LWP::UserAgent', 'HTTP::Request::Common', 'Config::File::Simple', 'Mozilla::CA']
		for module in required_modules:
			installed = subprocess.run( ['perl', '-M' + module, '-e', ''],  capture_output=True, text=True )
			if installed.returncode != 0:
				missing_modules.append( module )

	required_tools = ['nohup', 'parallel', 'diamond']
	if alnmethod == "mafft":
		required_tools.append( mafft )
	else: 
		required_tools.append( muscle )
	
	if treemethod == "raxml":
		required_tools.append( raxml )
	elif treemethod == "iqtree":
		required_tools.append( iqtree )
	else:
		required_tools.append( fasttree )

	missing_tools = []
	for tool in required_tools:
		if not check_tools( tool, docu_file ):
			missing_tools.append( tool )

	if missing_tools != [] or missing_modules != []:
		if missing_tools != []:
			missing_tools_string = ", ".join( missing_tools )
			sys.stdout.write( str( 'ERROR: Required tools ' + missing_tools_string + ' are not installed.\n Please install the required tools before starting again.\n' ) )
		if missing_modules != []:
			missing_modules_string = ", ".join( missing_modules )
			sys.stdout.write( str( 'ERROR: Required perl modules ' + missing_modules_string + ' are not installed but necessary for iTOL upload.\n Please install the required modules before starting again.\n' ) )
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

def handle_sighup(signum, frame):
	sys.stdout.write( str( datetime.datetime.now() ) + ' - ERROR: SIGHUP received. Please check if there is enough RAM available to conduct the analysis. Unfortunately, this error might occur. Please start the analysis again. \n' )
	sys.stdout.flush()
	raise SystemExit("Subprocess terminated due to SIGHUP")

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


def load_config_file_content( config_file, data_folder ):
	"""! @brief load content of config file, generate clean CDS file, and generate clean PEP file 
	"""
	infos = {}
	sample_counts = {}
	with open( config_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split(',')
				if len( parts ) > 5:
					sys.exit( 'Invalid structure of confgig file, requires:\nspecies,tpm,cds,baits(,pep) without header' )
				pep_file = data_folder + parts[0] + ".pep.fasta"	#name of new peptide sequence file
				cds_file = data_folder + parts[0] + ".cds.fasta"	#name of new CDS file
				cds_seqs = load_sequences( parts[2] )
				find_cds_seqs = {}
				with open( cds_file, "w" ) as out:	#generating clean CDS file
					for key in list( cds_seqs.keys() ):
						out.write( '>' + parts[0] + "@" + key + '\n' + cds_seqs[ key ] + "\n" )	#include species name in sequence name
						find_cds_seqs.update( { parts[0] + "@" + key: cds_seqs[ key ] } )

				find_pep_seqs = {}
				if len( parts ) == 5:	#if pep file is provided, use it directly
					pep_seqs = load_sequences( parts[4] )
				else:
					pep_seqs = translate( find_cds_seqs )	#translation of all CDS
				
				with open( pep_file, "w" ) as out:	#generate clean PEP file
					for key in list( pep_seqs.keys() ):
						if len( parts ) == 5:
							out.write( '>' + parts[0] + "@" + key + '\n' + pep_seqs[ key ] + "\n" )
							find_pep_seqs.update( { parts[0] + "@" + key: pep_seqs[ key ] } )
						else:
							out.write( '>' + key + '\n' + pep_seqs[ key ] + "\n" )
							find_pep_seqs = pep_seqs

				baits = load_baits( parts[3], parts[0] )

				with open( parts[1], "r" ) as tpm:
					samples = tpm.readline()
					sample_count = len( samples.strip().split( '\t' ) ) - 1		#subtract row-header (gene-ID)
					sample_counts[parts[0]] = ( sample_count )

				#storing information in a dictionary under the species name
				infos.update( { parts[0]: { 'id': parts[0], 'tpm': parts[1], 'cds_file': cds_file, 'baits_file': parts[3], 'pep_file': pep_file, 'cds': find_cds_seqs, 'pep': find_pep_seqs, 'baits': baits } } )
			
			line = f.readline()
	
	return infos, sample_counts

def generate_plot_exp( data, output_folder ):
	values = list(data.values())
	species = list(data.keys())

	# --- count height of each bin to determine shift of y-values --- #
	hist_counts, bin_edges = np.histogram(values, bins=10)	
	max_count = max(hist_counts)
	step = 1 / max_count
	bin_edges[-1] += 1e-6
	bin_indices = np.digitize(values, bin_edges) - 1
	y_positions = [0] * len(values)
	current_positions = [0] * (len(bin_edges) - 1)
	for i, val in enumerate(values):
		bin_idx = bin_indices[i]
		y_positions[i] = current_positions[bin_idx] * step
		current_positions[bin_idx] += 1

	# --- Plot generation --- #
	fig = go.Figure()
	fig.add_trace(go.Histogram(
		x=values,
		nbinsx=19,  # Specify the number of bins(somehow always one more than specified)
		name='',
		marker=dict( color='#40A578' ),
		opacity=0.75
	))

	fig.add_trace(go.Scatter(	#add individual data points
		x=values,
		y=y_positions,
		mode='markers',
		marker=dict(size=10, color='#006769'),
		text=[f'{species[i]}: {values[i]}' for i in range(len(values))],  # Hover text
		hoverinfo='text', 
		name='Single data points'
	))

	fig.update_layout(
		title='Sample counts',
		xaxis_title='Number of samples per species',
		yaxis_title='Count',
		showlegend = False,
		barmode='overlay'  # Overlay the histograms
	)

	histogram_file = output_folder + "species_count_histogram.html"
	fig.write_html( histogram_file )

	# --- warning if 20% of the species have less then 50 samples --- #
	total_values = len(values)
	low_sample_counts = sum(1 for value in values if value < 50)
	critical_percentage = total_values / 5

	perc_critical_species = low_sample_counts/total_values * 100
	if perc_critical_species >= critical_percentage:
		sys.stdout.write( "WARNING: " + str( round( perc_critical_species, 2 ) ) + "% of the analysed species have a low sample count number. Coexpression results might not be that powerful.\n" )


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


def compare_candidates_against_all( candidate, spec_prefix, gene_expression, rcutoff, pvaluecutoff, coexpnumcutoff, min_exp_cutoff, correlation_dict ):
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
					# Store the maximum correlation for each gene2 independent of whether its coexpressed or not
					if gene2 not in correlation_dict or r > correlation_dict[gene2]:
						correlation_dict[gene2] = r
				else:
					sys.stdout.write ( "WARNING: no expression detected - " + candidate + "\t" + gene2 + "\t" + ";".join( list( gene_expression.keys() )[:5] ) + "\n" )
					sys.stdout.flush()
	coexp_gene_IDs = []
	for gene in sorted( coexpressed_genes, key=itemgetter('correlation') )[::-1][:coexpnumcutoff]:
		coexp_gene_IDs.append( gene['id'] )
	return coexp_gene_IDs

def coexp_helper( arguments, __usage__ ):
	"""! @brief run coexpression analysis """

	mode = arguments[ arguments.index('--mode')+1 ]
	if mode != "coexp_helper":
		sys.exit( __usage__ )
	
	expression_file = arguments[ arguments.index('--exp')+1 ]
	species_id = arguments[ arguments.index('--specid')+1 ]
	genes = arguments[ arguments.index('--genes')+1 ].split('@@@')
	rcutoff = float( arguments[ arguments.index('--rcutoff')+1 ] )
	pvaluecutoff = float( arguments[ arguments.index('--pvaluecutoff')+1 ] )
	coexpnumcutoff = int( arguments[ arguments.index('--coexpnumcutoff')+1 ] )
	min_exp_cutoff = float( arguments[ arguments.index('--min_exp_cutoff')+1 ] )
	mindetection = int( arguments[ arguments.index('--mindetection')+1 ] )
	output_dir = arguments[ arguments.index('--output')+1 ]

	output_file = output_dir + species_id + '_coexp_results.txt'
	correlation_file = output_dir + 'correlation_coefficients.csv'			#one huge csv file for all sequences of all species
	
	
	exp = load_expression_values( expression_file, species_id )
	combined_coexp_genes = []
	correlation_dict = {}  # Store correlation coefficients

	for gene in genes:
		combined_coexp_genes += compare_candidates_against_all( gene, species_id, exp, rcutoff, pvaluecutoff, coexpnumcutoff, min_exp_cutoff, correlation_dict )
	unique_combined_coexp_genes = list( set( combined_coexp_genes ) )
	valid_coexp_genes = []
	for gene in unique_combined_coexp_genes:
		if combined_coexp_genes.count( gene ) >= mindetection:	#this selects genes co-expressed with entire pathway
			valid_coexp_genes.append( gene )
	
	# Write coexpressed candidates into file
	with open( output_file, "w" ) as out:
		out.write( "\n".join( valid_coexp_genes ) )

	# Write all correlation coefficients in csv file
	with open(correlation_file, "a", newline="") as cor_out:
		csv_writer = csv.writer(cor_out)
		for gene, correlation in correlation_dict.items():
			csv_writer.writerow([gene, correlation])

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

def aln_trimming( arguments, __usage__ ):
	"""! @brief trim alignment
	@can only be reached via terminal input
	"""
	mode = arguments[ arguments.index('--mode')+1 ]
	if mode != "aln_trimming":
		sys.exit( __usage__ )
	
	occupancy = arguments[ arguments.index('--occu')+1 ]

	aln_file = arguments[ arguments.index('--aln_file')+1 ]
	cln_aln_file = arguments[ arguments.index('--cln_aln_file')+1 ]

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
			if counter / float( len( list(alignment.keys()) ) ) > float(occupancy):	#collect positions of sufficient occupancy in list
				valid_index.append( idx )
		
		# --- generate new sequences --- #
		with open( cln_aln_file, "w" ) as out:
			for key in list(alignment.keys()):	#iterate over all sequences in the alignment
				seq = alignment[ key ]
				new_seq = []
				for idx in valid_index:	#only add residues of positions with sufficient occupancy
					new_seq.append( seq[ idx ] )
				out.write( ">" + key + '\n' + "".join( new_seq ) + '\n' )
	# --- just in case the alignment file is empty (is this possible?) ---#
	else:
		with open( cln_aln_file, "w" ) as out:
			out.write( "" )

def tree_constructor( X_aln_input_file, treemethod, X_output_folder, Xname, alnmethod, aln_tree_cmd_file, mafft, muscle, occupancy, raxml, fasttree, iqtree ):
	"""! @brief handles the construction of alignments and phylogenetic tree
		@note second FASTA file can be an empty string to run this function just based on one FASTA file
	"""
	
	X_aln_file = X_aln_input_file + ".aln"
	X_cln_aln_file = X_aln_file + ".cln"

	# --- generate alignment command --- #
	if alnmethod == "muscle":
		current_aln_command = muscle + " -in " + X_aln_input_file + " -out " + X_aln_file + " -quiet"
	elif alnmethod == "mafft":
		current_aln_command = mafft + " --quiet " + X_aln_input_file + " > " + X_aln_file

	# --- generate trim alignment-command i.e. remove positions with many gaps = low occupancy --- #
	current_aln_trimming_cmd = "python3 " + os.path.abspath(__file__) + " --mode aln_trimming --occu " + str(occupancy) + " --aln_file " + X_aln_file + " --cln_aln_file " + X_cln_aln_file
	
	# --- generate tree command
	if treemethod == "raxml":	#construct tree with RAxML
		prefix = X_output_folder + Xname + "RAxML_tree"
		current_tree_command = raxml  + " --all --threads 1 --model LG+G8+F --msa " + X_cln_aln_file + " --prefix " + prefix + " && cp " + prefix + ".treefile " + prefix + ".tree"
	
	elif treemethod == "iqtree":	#construct tree with IQ-TREE2
		prefix = X_output_folder + Xname + "IQtree"
		current_tree_command = iqtree + " -alrt 1000 -bb 1000 -s " + X_cln_aln_file + " --prefix " + tree_file + " && cp " + prefix + ".treefile " + prefix + ".tree"
	
	else:	#construct tree with FastTree2
		tree_file = X_output_folder  + Xname  + "FastTree_tree.tree"
		current_tree_command = fasttree + " -wag -nopr -nosupport < " + X_cln_aln_file + " > " + tree_file
	
	with open( aln_tree_cmd_file, 'a' ) as command_out:
		command_out.write( current_aln_command + " && " + current_aln_trimming_cmd + " && "+ current_tree_command + "\n" )


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


def annotate_group( blast_result_file, anno_mapping_table, seq_num ):
	"""! @brief annotate whole cluster: take annotation that occures the most
	@return string comprizing functional annotation if available and the confidence score
	"""

	blast_file_direct = blast_result_file.split("/")[-1]
	cluster_no = blast_file_direct.split(".")[0]

	seqs = {}				# dictionary{seq: bitscore}
	blast_hits = {}			# dictionary{seq: annotation}
	with open( blast_result_file, "r" ) as f:
		
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			new_seq = parts[0]
			if new_seq in seqs :
				if float( parts[-1] ) > seqs[new_seq]:		# only update if bitscore is higher --> only take best annotation per sequence
					seqs[new_seq] = float( parts[-1] )
					blast_hits[new_seq] = parts[1]
			else:
				seqs[new_seq] = float( parts[-1] )
				blast_hits[new_seq] = parts[1]

			line = f.readline()
	
	if blast_hits != {}:				# skip if no annotation was made
		if len( seqs ) <= 10:			# warning for low annotation number
			sys.stdout.write( str( 'WARNING: Low number of annotation matches for cluster  ' + cluster_no + '\n' ) )	
			sys.stdout.flush()

		counts = Counter( blast_hits.values() )			# determine occurence of each annotation
		max_count = max( counts.values() )		# get the number of the most common annotation
		best_hits = [key for key, value in counts.items() if value == max_count]		# collect each annotation that occurs the most

		if len( best_hits ) > 1:						#if multiple annos equally common: take best bit score
			best_hit = max(best_hits, key=lambda x: [seqs[seq_id] for seq_id, anno in blast_hits.items() if anno == x])
		else:
			best_hit = best_hits[0]

		# --- Determine confidence_score --- #
		confidence_score = ( max_count/seq_num ) * 100		# number of annotation found/all seqeunces used for annotation

	else:
		best_hit = "n/a"
		confidence_score = 0.0
		sys.stdout.write( str( 'WARNING: No annotation possible for cluster  ' + cluster_no + '\n' ) )

	try:
		return best_hit + " - " + anno_mapping_table[ best_hit ], confidence_score
	except KeyError:
		return best_hit + " - " + "n/a", confidence_score
	

def sort_groups( group ):
	"""! @brief sorts groups according to the presence of coexp-sequences """

	coexp_count = sum('coexp' in node for node in group)	# determine number of coexp seqs
	at_count = sum('@' in node for node in group)			# determine number of seqs in general
	if at_count == 0:
		return float('inf')
	return -(coexp_count / at_count)		# negative because of descanding sorting

def get_groups( annotation_file, ref_seq_file, blast_results, anno_mapping_table, huge_pep_collection, aln_folder, mincoexpseqcutoff, minseqcutoff, tmp_folder, cpub, batch, perc_cluster_anno, clean ):
	"""! @brief identify groups in huge cluster """

	# --- identify groups --- #
	G = nx.Graph()
	for key, values in blast_results.items():
		for value in values:
			G.add_edge(key, value)
		G.add_node(key)			# to ensure isolated nodes are also present

	num_nodes = G.number_of_nodes()		#determine number of sequences in general
	final_groups = [list(component) for component in nx.connected_components(G)]	 #Identification of connected components
	final_groups = sorted( final_groups, key = sort_groups )		# Sort groups according to ratio between coexp-seqs and number of seqs per cluster

	sys.stdout.write( str( datetime.datetime.now() ) + " ... groups identified ... \n" )
	sys.stdout.flush()
	
	# --- prepare database for annotation --- #
	if len( ref_seq_file ) > 0:
		tmp_blast_folder = tmp_folder + "anno_blast/"
		if not os.path.exists( tmp_blast_folder ):
			os.makedirs( tmp_blast_folder )
		blast_db = tmp_blast_folder + "blastdb"
		p = subprocess.Popen( args= "diamond makedb --in " + ref_seq_file + " -d " + blast_db, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL )
		p.communicate()
	
		# --- write sequences into output files --- #
		anno_blast_cmd_file = tmp_blast_folder + "ANNO_BLAST_COMMANDS.sh"
		with open( anno_blast_cmd_file, "w" ) as anno_command_out:
			anno_command_out.write( str ( "#!/bin/bash\nnohup parallel -j " + str( batch ) + " <<EOF \n" ) )
	
	with open( annotation_file, "w" ) as anno_out:
		anno_out.write( "Tree" + "\t" + "Nodes" + "\t" + "Confidence of annotation" + "\t" + "Annotation" + "\n" )
		counter = 0
		final_cluster = [] 		#identify groups with sufficient coexpressed candidates 
		seq_nums = {}			#store number of sequences per cluster

		for fg in final_groups:
			if len( fg ) > minseqcutoff: 
				coexpseqcounter = count_coexp_specs( fg )
				if coexpseqcounter >= mincoexpseqcutoff:
					counter += 1
					final_cluster.append( fg )
					group_sequence_file = aln_folder + str( counter ).zfill( 4 ) + ".pep.fasta"	#save all sequenes of a group in a FASTA file
					if len( ref_seq_file ) > 0:
						group_anno_file = tmp_blast_folder + str( counter ).zfill( 4 ) + "_anno.pep.fasta"
					seq_collection_anno = []
					for gene in fg:
						if clean:
							clean_gene = gene.replace( ":", "_" )
							current_seq = '>' + clean_gene + "\n" + huge_pep_collection[ gene ] + "\n"
						else:
							current_seq = '>' + gene + "\n" + huge_pep_collection[ gene ] + "\n"
						
						with open( group_sequence_file, "a" ) as out:
							out.write( current_seq )
						if len( ref_seq_file ) > 0:
							seq_collection_anno.append( current_seq )
					
					if len( seq_collection_anno ) > 20:
						seq_num = int( len ( seq_collection_anno ) * ( perc_cluster_anno / 100 ) )  #determine number of sequences for annotation
						if seq_num < 10:					#take at least 10 seqences for annotation
							seq_num = 10
						seq_nums[ counter ] = seq_num		#add dictionary entry to determine confidence score
						final_seqs = random.sample( seq_collection_anno, seq_num )
					else:
						seq_nums[ counter ] = len( seq_collection_anno )
						final_seqs = seq_collection_anno		#to avoid that less than 10 seqs are annotated
					
					if len( ref_seq_file ) > 0:
						with open( group_anno_file, 'w') as out:
							for seq in final_seqs:
								out.write( seq )
						with open(anno_blast_cmd_file, "a" ) as anno_command_out:
							blast_result_file = tmp_blast_folder + group_sequence_file.split('/')[-1] + "_blast_results.txt"
							current_anno_command = "diamond blastp -q " + group_anno_file + " -d " + blast_db + " -o " + blast_result_file + " --threads " + str( cpub ) + "\n"
							anno_command_out.write( current_anno_command )

		if len( ref_seq_file ) > 0:
			with open( anno_blast_cmd_file, "a" ) as anno_command_out:
				anno_command_out.write( "EOF" )
			os.chmod( anno_blast_cmd_file, 0o755 )
			subprocess.run( ['nohup', '/bin/bash', anno_blast_cmd_file], check = True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL )
			sys.stdout.write( str( datetime.datetime.now() ) + " - DIAMOND blast for annotation completed! " + "\n" )
			sys.stdout.flush()

		counter = 0
		for fg in final_cluster:  # return best annotation
			counter += 1
			if len( ref_seq_file ) > 0: 
				group_sequence_file = aln_folder + str( counter ).zfill( 4 ) + ".pep.fasta"
				blast_result_file = tmp_blast_folder + group_sequence_file.split('/')[-1] + "_blast_results.txt"
				seq_num = seq_nums[ counter ]
				group_annotation, score = annotate_group( blast_result_file, anno_mapping_table, seq_num )
			else:
				group_annotation = "n/a"
				score = 0
			anno_out.write( str( counter ).zfill( 4 ) + "\t" + str( len( fg ) ) + "\t" + str( round( score, 2 ) )+ "\t" + group_annotation + "\n" )

def interpolate_color( x ):
	"""! @brief determine corresponding color for max and min value for each iTOL tree to keep the color consitent over all trees
	indepentent from their highest and lowest values
	"""

	# colors as RGB (0-255)
	min_color = (255, 89, 78)  # #FF594E
	mid_color = (255, 255, 255)  # #FFFFFF
	max_color = (13, 116, 255)  # #0D74FF
	
	if x < 0:
		# Interpolate between red and white
		ratio = (x + 1)  # Map von [-1, 0] auf [0, 1]
		r = int(min_color[0] + (mid_color[0] - min_color[0]) * ratio)
		g = int(min_color[1] + (mid_color[1] - min_color[1]) * ratio)
		b = int(min_color[2] + (mid_color[2] - min_color[2]) * ratio)
	else:
		# Interpolate between white and blue
		ratio = x  # Map von [0, 1] auf [0, 1]
		r = int(mid_color[0] + (max_color[0] - mid_color[0]) * ratio)
		g = int(mid_color[1] + (max_color[1] - mid_color[1]) * ratio)
		b = int(mid_color[2] + (max_color[2] - mid_color[2]) * ratio)
	
	# Convert to hex
	return f"#{r:02X}{g:02X}{b:02X}"

def upload_itol( output_folder, tree_folder, api, proj_name, anno, perl_path, corr_file ):
	"""! @brief generate annotation file for each tree
	@zip tree-file + config_file
	@generate config file for each tree
	@conduct upload
	"""

	itol_folder = output_folder + "itol/"
	if not os.path.exists( itol_folder ):
		os.makedirs( itol_folder )

	iTOL_IDs = {}
	tree_files = sorted( glob.glob( tree_folder + "*.tree" ) )
	anno_dict = {}
	with open( anno, 'r' ) as anno_table:
		content = csv.reader( anno_table, delimiter='\t' )
		for row in content:
			anno_dict[str( row[0] ) ] = row[3]

	# Load correlation values
	correlation_dict = {}
	with open(corr_file, 'r') as corr_file:
		corr_reader = csv.reader(corr_file)
		for row in corr_reader:
			correlation_dict[row[0]] = float(row[1])

	for tree in tree_files:
		coexp_seqs = []
		all_seqs = []
		tree_no = os.path.basename( tree )
		key = tree_no[:4]		#identify number of tree
		anno_tree = anno_dict[ key ]		#collect annotation for tree description
		anno_coexp_file = itol_folder + key + "_iTOL_coexp_annotation.txt"
		gradient_file = itol_folder + key + "_iTOL_gradient_labels.txt"

		# --- identify coexp genes --- #
		with open( tree, 'r' ) as tree_file:
			tree_content = tree_file.read()
			for word in re.split( r'[(),:]', tree_content ):
				if 'coexp' in word:
					word = word.strip()
					coexp_seqs.append( word )
				all_seqs.append( word )

		with open( anno_coexp_file, 'w' ) as anno_itol:
			anno_itol.write( "DATASET_STYLE\nSEPARATOR COMMA\nDATASET_LABEL,coexp_genes\nCOLOR,#00ff00\nDATA")
			for coexp in coexp_seqs:
				anno_itol.write( "\n" + coexp + ",label,node,#000000,1,bold" )
				#anno_itol.write( "\n" + coexp + ",label,node,#000000,1,bold,#40A578" )
			
		# Collect all values
		all_values = {}
		for seq in all_seqs:
			if seq in coexp_seqs:
				base_name = seq.replace("_coexp", "")
			else:
				base_name = seq
			if base_name in correlation_dict:
				value = round(correlation_dict[base_name], 2)
				all_values[seq] = value

		# Determine max and min color
		min_value = min(all_values.values())
		min_col = interpolate_color( min_value )

		max_value = max(all_values.values())
		max_col = interpolate_color ( max_value )

		mid_value = max_value - 0.5 *( max_value - min_value )
		mid_color = interpolate_color( mid_value )


		with open(gradient_file, 'w') as gradient:
			gradient.write("DATASET_GRADIENT\nSEPARATOR COMMA\nDATASET_LABEL,Correlation coefficients,1\nCOLOR,#80B0F9\nSTRIP_WIDTH,100\nAUTO_LEGEND,1")
			gradient.write("\nCOLOR_MIN," + min_col )
			gradient.write("\nCOLOR_MAX," + max_col )
			gradient.write("\nUSE_MID_COLOR,1\nCOLOR_MID," + mid_color + "\nLABEL_ALIGN_TO_TREE,0\nDATA\n")
			for seq_id, value in all_values.items():
				gradient.write(seq_id + "," + str(value) + "\n")			
		
		zip_folder = str(itol_folder + key + '.zip')
		subprocess.run(['zip', '-j', zip_folder, tree, anno_coexp_file, gradient_file], stdout=subprocess.DEVNULL, check = True)
		# --- generate config file --- #
		upload_config = itol_folder + key + "_iTOL.cfg"
		with open( upload_config, 'w' ) as cfg:
			cfg.write( "zipFile=" + zip_folder + "\nprojectName=" + proj_name + "\nAPIkey=" + api + "\n")
			if anno_table != 'n/a':
				cfg.write( "treeDescription=" + anno_tree + "\n" )
	
		results = subprocess.run([perl_path, upload_config], capture_output = True, text = True)
		if 'error' in results.stdout.lower():
			sys.stdout.write( str( datetime.datetime.now() ) + "ERROR with uploading" + tree + " : " + results.stdout + "\n" )
			sys.stdout.flush()
			iTOL_IDs[tree_no] = "Upload failed"
		if results.stderr:
			sys.stdout.write( str( datetime.datetime.now() ) + "ERROR with uploading" + tree + " : " + results.stderr + "\n" )
			sys.stdout.flush()
			iTOL_IDs[tree_no] = "Upload failed"
		else:
			id = results.stdout.strip().split('\n')[-1]
			iTOL_IDs[tree_no] = id
	
	id_file = tree_folder + "docu_iTOL_IDs.txt"
	with open( id_file, 'w' ) as out:
		out.write( "Your trees are accessible using the following iTOL tree IDs:\n" )
		for tree_id in iTOL_IDs:
			out.write( tree_id + "\t" + iTOL_IDs[tree_id] + "\n" )

	sys.stdout.write( str( datetime.datetime.now() ) + " Upload finished. Documentation at: " + id_file + "\n" )
	sys.stdout.flush()



def main( arguments ):
	"""! @brief run everything """

	config_file = arguments[ arguments.index('--config')+1 ]
	#ID,tpm,cds,baits
	#ID = species ID/name; tpm= tpm file; cds = CDS file; baits = IDs of interest in file
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
	if alnmethod != "mafft" and alnmethod != "muscle":
		sys.exit( "ERROR: alignment tool must be either mafft or muscle!" )
	
	
	if '--mafft' in arguments:	#path to MAFFT (alignment)
		mafft = arguments[ arguments.index('--mafft')+1 ]
	else:
		mafft = "mafft"
	
	if '--muscle' in arguments:	#path to MUSCLE (alignment)
		muscle = arguments[ arguments.index('--muscle')+1 ]
	else:
		muscle = "muscle"	
	
	if '--occupancy' in arguments:
		occupancy = float( arguments[ arguments.index('--occupancy')+1 ]  )
	else:
		occupancy = 0.1
	
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
	
	if '--treemethod' in arguments:	#tree construction methods
		treemethod = arguments[ arguments.index('--treemethod')+1 ]
		if treemethod not in [ "fasttree", "raxml", "iqtree" ]:
			treemethod = "fasttree"
	else:
		treemethod = "fasttree"
	
	if '--reference' in arguments:
		ref_seq_file = arguments[ arguments.index('--reference')+1 ]
	else:
		ref_seq_file = ""

	if '--seqs_cluster_anno' in arguments:
		seqs_cluster_anno = int( arguments[ arguments.index('--seqs_cluster_anno')+1 ] )
		if seqs_cluster_anno > 100.0:
			sys.exit( "ERROR: seqs_cluster_anno value must be <= 100" )
		elif seqs_cluster_anno <= 0.0:
			sys.exit( "ERROR: seqs_cluster_anno value must be > 0" )
		elif seqs_cluster_anno < 10.0:
			seqs_cluster_anno = 10.0
			sys.stdout.write( "WARNING: 10% of the sequences per cluster should be annotated. seqs_cluster_anno was automatically set to 10. \n" )
			sys.stdout.flush()
	else:
		seqs_cluster_anno = 50.0

	if '--clean' in arguments:	#replace colons in gene IDs before tree generation
		clean = True
	else:
		clean = False

	if '--API' in arguments:
		api =  arguments[ arguments.index('--API')+1 ] 
	else:
		api = ''
	
	if '--proj_name' in arguments:
		proj_name = arguments[ arguments.index('--proj_name')+1 ] 
	else:
		proj_name = ''

	if ( api != '' and proj_name == '' ) or ( api == '' and proj_name != '' ):
		sys.exit( "ERROR: Batch upload to iTOL requires API key (--API) AND project name (--proj_name). Please provide both parameters or none if no iTOL batch upload is wanted.\n")

	if api != '':
		if '--upload_script' in arguments:
			perl_path = arguments[ arguments.index('--upload_script')+1 ] 
			if not os.path.isfile( perl_path ):
				sys.exit( "ERROR: iTOL_uploader.pl not found at " + perl_path + "\n" )
		else:
			perl_path = str ( os.path.dirname( __file__ )+ "/iTOL_uploader.pl" )
			if not os.path.isfile( perl_path ):
				perl_path2 = str ( os.path.dirname( __file__ )+ "/batch_scripts/iTOL_uploader.pl" )
				if not os.path.isfile( perl_path2 ):
					sys.exit( "ERROR: iTOL_uploader.pl not found at " + perl_path + "\nPlease provide a path to the upload script.\n" )
		os.chmod( perl_path, 0o755 )
	else:
		perl_path = ''
	
	if output_folder[-1] != "/":
		output_folder += "/"
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	tmp_folder = output_folder + "tmp/"
	if not os.path.exists( tmp_folder ):
		os.makedirs( tmp_folder )

	data_folder = output_folder + "data/"
	if not os.path.exists( data_folder ):
		os.makedirs( data_folder )

	shell_folder = output_folder + "shell/"
	if not os.path.exists( shell_folder ):
		os.makedirs( shell_folder )

	### --- Documentation --- ###
	sys.stdout.write( str( datetime.datetime.now() ) + " - preparing documentation file ... \n" )
	sys.stdout.flush()
	docu_file = output_folder + datetime.date.today().strftime('%Y%m%d') + '_docu.txt'
	generate_docu_file ( tmp_folder, docu_file, config_file, rcutoff, pvaluecutoff, scorecut, simcut, lencut, evalue, mindetection, minseqcutoff, mincoexpseqcutoff, min_exp_cutoff, alnmethod, treemethod, anno_file, ref_seq_file, muscle, mafft, raxml, iqtree, fasttree, seqs_cluster_anno, clean, api, proj_name, perl_path )
	sys.stdout.write( str( datetime.datetime.now() ) + " ... documentation file completed: " + docu_file + "\n" )
	sys.stdout.write( str( datetime.datetime.now() ) + " loading data ... \n" )
	sys.stdout.flush()
	config_data, sample_counts = load_config_file_content( config_file, data_folder )	#a species tag is included while loading the sequences to make them unique
	sys.stdout.write( str( datetime.datetime.now() ) + " ... loading data completed.\n" )
	sys.stdout.write( str( datetime.datetime.now() ) + " Generating plot to display sample counts per species ...\n" )
	sys.stdout.flush()
	generate_plot_exp( sample_counts, output_folder )
	sys.stdout.write( str( datetime.datetime.now() ) + " ... Expression data plot stored at " + output_folder + "species_count_histogram.html\n" )
	sys.stdout.flush()
	
	if len( anno_file ) > 1:
		anno_mapping_table = load_annotation( anno_file )
	else:
		anno_mapping_table = {}
	
	# --- Determine maximal number of jobs with one core to run in parallel --- #
	max_cores = cpub * batch
	if max_cores > 1:
		max_cores -= 1

	### ---- BIG LOOP OVER ALL SPECIES --- ###
	huge_seq_collection_file = output_folder + "huge_seq_collection.cds.fasta"
	huge_seq_collection_file_pep = output_folder + "huge_seq_collection.pep.fasta"
	huge_seq_collection = {}

	# --- Preparation shell scripts --- #
	coexp_cmd_file = shell_folder + "COEXP_COMMANDS.sh"
	blast_dbs_file = shell_folder + "BLAST_DBs.sh"
	blast_cmd_file = shell_folder + "BLAST_COMMANDS.sh"
	aln_tree_cmd_file = shell_folder + "ALIGNMENT_TREE_COMMANDS.sh"
	
	shell_header = "#!/bin/bash \n"
	batch_header = "nohup parallel -j " + str( batch ) + " <<EOF \n"
	batch_coexp = "nohup parallel -j " + str( int( batch/2 ) + 1 ) + " <<EOF \n"		#less jobs parallel due to high memory demand of coexpression analysis
	batch_tree = "nohup parallel -j " + str( max_cores - 1 ) + " <<EOF \n"
	with open( coexp_cmd_file, "w" ) as coexp_command_out:
		coexp_command_out.write( shell_header + batch_coexp )
	with open( blast_dbs_file, "w" ) as blast_dbs_out:
		blast_dbs_out.write( shell_header )			
	with open( blast_cmd_file, "w" ) as blast_command_out:
		blast_command_out.write( shell_header + batch_header )
	with open( aln_tree_cmd_file, "w" ) as aln_tree_command_out:
		aln_tree_command_out.write( shell_header + batch_tree )

	# --- SIGHUP iinterruption handling --- #
	signal.signal(signal.SIGHUP, handle_sighup)
		
	# --- run co-expression per species --- #
	sys.stdout.write( str( datetime.datetime.now() ) + " - Starting coexpression analyses ... " + "\n" )
	sys.stdout.flush()
	for spec in list( config_data.keys() ):
		info = config_data[ spec ]
		coexp_result_file = tmp_folder + info['id'] + "_coexp_results.txt"
		if not os.path.isfile( coexp_result_file ):
		
			# --- write command for co-expression analysis and to build blast db --- #
			# modify gene names if there is a special character in gene name
			baits_modified = []
			for bait in ( info['baits'] ):
				baits_modified.append( bait.replace( "|", "\|" ).replace( "$", "\$" ).replace( "*", "\*" ).replace( ";", "\;" ).replace( "~", "\~" ).replace( "/", "\/" ) )

			with open( coexp_cmd_file, "a" ) as coexp_command_out:
				current_coexp_cmd = "".join( [ "python3 " + os.path.abspath(__file__),
																	" --mode coexp_helper",
																	" --exp ", info['tpm'],
																	" --specid ", info['id'],
																	" --genes ", "@@@".join( baits_modified ),
																	" --rcutoff ", str( rcutoff ),
																	" --pvaluecutoff ", str( pvaluecutoff ),
																	" --coexpnumcutoff ", str( coexpnumcutoff ),
																	" --min_exp_cutoff ", str( min_exp_cutoff ),
																	" --mindetection ", str( mindetection ),
																	" --output ", tmp_folder
																] )
				coexp_command_out.write( current_coexp_cmd + "\n" )
			
		blastdb = tmp_folder + spec + "_blastdb"
		blastdb_file = blastdb + ".dmnd"		# to check whether file already exists
		if not os.path.isfile( blastdb_file ):
			with open( blast_dbs_file,  "a" ) as blast_dbs_out:
				current_blast_db_command = "diamond makedb --in " + config_data[ spec ]['pep_file'] + " -d " + blastdb
				blast_dbs_out.write( current_blast_db_command + "\n" + "sleep 2" + "\n" )				

	# --- execute coexpression and blast db scripts --- #			
	with open( coexp_cmd_file, "a" ) as coexp_command_out:
		coexp_command_out.write( "EOF" )
	os.chmod( coexp_cmd_file, 0o755 )
	try:
		subprocess.run( ['nohup', '/bin/bash', coexp_cmd_file], check = True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE )
		sys.stdout.write( str( datetime.datetime.now() ) + " ... coexpression analyses completed! " + "\n" )
		sys.stdout.write( str( datetime.datetime.now() ) + " - Building BLASTDBs " + "\n" )
		sys.stdout.flush()
	except subprocess.CalledProcessError as e:
		error_file = output_folder + "error_coexp.log"
		with open ( error_file, 'wb' ) as error_log:
			error_log.write( e.stderr )
		sys.exit( "ERROR during coexpression analysis, please check " + error_file + "\n" )

	os.chmod( blast_dbs_file, 0o755 )
	try:
		subprocess.run( ['nohup', '/bin/bash', blast_dbs_file], check = True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE )
		sys.stdout.write( str( datetime.datetime.now() ) + " - All blast databases built! " + "\n" )
		sys.stdout.flush()
	except subprocess.CalledProcessError as e:
		error_file = output_folder + "error_blastdbs.log"
		with open ( error_file, 'wb' ) as error_log:
			error_log.write( e.stderr )
		sys.exit( "ERROR during blast DB generation analysis, please check " + error_file + "\n" )

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
				if spec != spec2:
					blast_result_file = tmp_folder + spec + "_vs_" + spec2 + "blasthits.txt"	#run BLASTp search against each other species
					if not os.path.isfile( blast_result_file ):
						blastdb = tmp_folder + spec2 + "_blastdb"
						with open( blast_cmd_file, "a" ) as blast_command_out:
								current_search_cmd = "diamond blastp -q " + baits_file + " -d " + blastdb + " -o " + blast_result_file + " --threads " + str( cpub )
								blast_command_out.write( current_search_cmd + "\n" )
	
	# --- execute all blastp commands --- #
	with open( blast_cmd_file, "a" ) as blast_command_out:
		blast_command_out.write( "EOF" )
	sys.stdout.write( str( datetime.datetime.now() ) + " - Starting DIAMOND blastp ... " + "\n" )
	sys.stdout.flush()
	os.chmod( blast_cmd_file, 0o755 )
	subprocess.run( ['/bin/bash', blast_cmd_file], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL )
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

	# --- run BLAST all vs. all for huge sequence collection --- #
	huge_pep_collection = translate( huge_seq_collection )
	with open( huge_seq_collection_file_pep, "w" ) as out:
		for key in list( huge_pep_collection.keys() ):
			out.write( '>' + key + '\n' + huge_pep_collection[ key ] + "\n" )
	blastdb = tmp_folder + "huge_blastdb"
	blast_result_file = tmp_folder + "huge_blast_result_file.txt"
	if not os.path.isfile( blast_result_file ):
		sys.stdout.write( str( datetime.datetime.now() ) + " - Starting DIAMOND BLAST for sequence clustering...\n" )
		sys.stdout.flush()
		p = subprocess.Popen( args= "diamond makedb --in " + huge_seq_collection_file_pep + " -d " + blastdb, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL )
		p.communicate()
		p = subprocess.Popen( args= "diamond blastp -q " + huge_seq_collection_file_pep + " -d " + blastdb + " -o " + blast_result_file + " --threads " + str( max_cores ), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL )
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
	sys.stdout.write( str( datetime.datetime.now() ) + " - constructing sequence clusters...\n" )
	sys.stdout.flush()
	get_groups( annotation_file, ref_seq_file, blast_results, anno_mapping_table, huge_pep_collection, aln_folder, mincoexpseqcutoff, minseqcutoff, tmp_folder, cpub, batch, seqs_cluster_anno, clean )
	sys.stdout.write( str( datetime.datetime.now() ) + " ... sequence clustering completed.\n" )
	sys.stdout.flush()
		
	# --- construct phylogenetic trees for each of them and highlight co-expressed genes --- #
	sys.stdout.write( str( datetime.datetime.now() ) + " - Starting construction of phylogenetic trees...\n" )
	sys.stdout.flush()
	fasta_input_files = glob.glob( aln_folder + "*.pep.fasta" )
	for filename in fasta_input_files:
		name = filename.split('/')[-1].split('.')[0]
		tree_constructor( filename, treemethod, tree_folder, name, alnmethod, aln_tree_cmd_file, mafft, muscle, occupancy, raxml, fasttree, iqtree )
	
	with open( aln_tree_cmd_file , "a" ) as aln_tree_command_out:
		aln_tree_command_out.write( "EOF" )
	os.chmod( aln_tree_cmd_file, 0o755 )
	try:
		subprocess.run(['nohup', '/bin/bash', aln_tree_cmd_file], check = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE )
	except subprocess.CalledProcessError as e:
		error_file = output_folder + "error_alignment_tree.log"
		with open ( error_file, 'wb' ) as error_log:
			if e.stderr:
				error_log.write( e.stderr )
			elif e.stdout:			# only wirte stdout if no stderr was captured
				error_log.write( e.stdout )

		sys.exit( "ERROR during alignment or tree generation, please check " + error_file )

	sys.stdout.write( str( datetime.datetime.now() ) + " ... construction of phylogenetic trees completed.\n" )
	sys.stdout.flush()

	if api != '':
		sys.stdout.write( str( datetime.datetime.now() ) + " Preparing upload to iTOL ...\n" )
		sys.stdout.flush()
		corr_file = tmp_folder + 'correlation_coefficients.csv'
		upload_itol( output_folder, tree_folder, api, proj_name, annotation_file, perl_path, corr_file )

		sys.stdout.write( str( datetime.datetime.now() ) + " Trees uploaded to project \"" + proj_name + "\" in iTOL.\n" )
	
	sys.stdout.write( str( datetime.datetime.now() ) + " FINISHED!\n" )
	sys.stdout.flush()
		
	# --- identify and define orthogroups --- #
	#find "coexp" sequences that are clustered = calculate distances between them and check for edges < number of species

if '--config' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
elif '--mode' in sys.argv and '--aln_file' in sys.argv:
	aln_trimming( sys.argv, __usage__ )
elif '--mode' in sys.argv and '--exp' in sys.argv:
	coexp_helper( sys.argv, __usage__ )
elif '--version' in sys.argv or '-v' in sys.argv:
	sys.exit( __version__ )
else:
	sys.exit( __usage__ )
