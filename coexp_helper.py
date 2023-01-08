### Boas Pucker ###
### b.pucker@tu-braunschweig.de ###

### WARNING: THIS SCRIPT IS PARTH OF CoExpPhylo - DO NOT RUN IT ON ITS OWN ###
### https://github.com/bpucker/CoExpPhylo ###

import os, re, sys, subprocess, math, glob, time
from operator import itemgetter
import numpy as np
from scipy import stats

# --- end of imports --- #

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


def compare_candidates_against_all( candidate, gene_expression, rcutoff, pvaluecutoff, coexpnumcutoff, min_exp_cutoff ):
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
			if total_expression > min_exp_cutoff:	#only consider genes with a minimal expression level
				if len( values ) > 0:
					r, p = stats.spearmanr( values )
					if not math.isnan( r ):
						if r > rcutoff and p < pvaluecutoff:
							coexpressed_genes.append( { 'id': gene2, 'correlation': r, 'p_value': p } )
				else:
					sys.stdout.write ( "WARNING: no expression detected - " + candidate + "\t" + gene2 + "\t" + ";".join( list( gene_expression.keys() )[:5] ) + "\n" )
					sys.stdout.flush()
	if len( coexpressed_genes ) > 0:
		coexp_gene_IDs = []
		for gene in sorted( coexpressed_genes, key=itemgetter('correlation') )[::-1][:coexpnumcutoff]:
			coexp_gene_IDs.append( gene['id'] )
		return coexp_gene_IDs
	else:
		return []


def main( arguments ):
	"""! @brief run everything """
	
	expression_file = arguments[ arguments.index('--exp')+1 ]
	species_id = arguments[ arguments.index('--specid')+1 ]
	genes = arguments[ arguments.index('--genes')+1 ].split('@@@')
	rcutoff = float( arguments[ arguments.index('--rcutoff')+1 ] )
	pvaluecutoff = float( arguments[ arguments.index('--pvaluecutoff')+1 ] )
	coexpnumcutoff = int( arguments[ arguments.index('--coexpnumcutoff')+1 ] )
	min_exp_cutoff = float( arguments[ arguments.index('--min_exp_cutoff')+1 ] )
	mindetection = int( arguments[ arguments.index('--mindetection')+1 ] )
	output_file = arguments[ arguments.index('--output')+1 ]
	
	exp = load_expression_values( expression_file, species_id )
	combined_coexp_genes = []
	for gene in genes:
		combined_coexp_genes += compare_candidates_against_all( gene, exp, rcutoff, pvaluecutoff, coexpnumcutoff, min_exp_cutoff )
	unique_combined_coexp_genes = list( set( combined_coexp_genes ) )
	valid_coexp_genes = []
	for gene in unique_combined_coexp_genes:
		if combined_coexp_genes.count( gene ) >= mindetection:	#this slects genes co-expressed with entire pathway
			valid_coexp_genes.append( gene )
	
	with open( output_file, "w" ) as out:
		out.write( "\n".join( valid_coexp_genes ) )

main( sys.argv )
