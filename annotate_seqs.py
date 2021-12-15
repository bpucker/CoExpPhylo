### Boas Pucker ###
### b.pucker@tu-braunschweig.de ###

### some functions taken from MYB_annotator.py and KIPEs ###

__version__ = "v0.01"

__usage__ = """
					python3 annotate_seqs.py
					--in <INPUT_FILE> | --in <INPUT_FOLDER>
					--ref <REF_FASTA_FILE>
					--out <OUTPUT_FOLDER>
					
					optional:
					--anno <ANNOTATION_FILE>
					"""

import os, re, sys, subprocess, glob
from operator import itemgetter

# --- end of imports --- #

def load_annotation( anno_file ):
	"""! @brief load annotation from given file """
	
	anno = {}
	with open( anno_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			anno.update( { parts[0]: ";".join( parts[1:] ) } )
			line = f.readline()
	return anno


def load_best_hit_per_query( blast_result_file ):
	"""! @brief load best hit per query """
	
	best_hits = {}
	with open( blast_result_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				data = best_hits[ parts[0] ]
				if float( parts[-1] ) > data['score']:
					del best_hits[ parts[0] ]
					best_hits.update( { parts[0]: { 'score': float( parts[-1] ), 'subject': parts[1] } } )
			except:
				best_hits.update( { parts[0]: { 'score': float( parts[-1] ), 'subject': parts[1] } } )
			line = f.readline()
	final_best_hits = {}
	for each in list( best_hits.keys() ):
		final_best_hits.update( { each: best_hits[ each ]['subject'] } )
	return final_best_hits


def find_most_abundant_hit( best_hit_per_query ):
	"""! @brief find most abundant hit """
	
	hits = list( best_hit_per_query.values() )
	abundances = []
	for hit in list( set( hits ) ):
		abundances.append( { 'id': hit, 'value': hits.count( hit ) } )
	sorted_hits = sorted( abundances, key=itemgetter('value') )[::-1]
	return sorted_hits[0]['id']


def main( arguments ):
	"""! @brief run everything """
	
	user_input = arguments[ arguments.index('--in')+1 ]
	if os.path.isfile( user_input ):
		input_files = [ user_input ]
	else:
		input_files = sorted( glob.glob( user_input + "*.fasta" ) )
	ref_file = arguments[ arguments.index('--ref')+1 ]
	output_folder = arguments[ arguments.index('--out')+1 ]
	
	if '--anno' in arguments:
		anno_file = arguments[ arguments.index('--anno')+1 ]
		anno = load_annotation( anno_file )
	else:
		anno = {}
	
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	ref_db = output_folder + "ref_blastdb"
	p = subprocess.Popen( args= "makeblastdb -in " + ref_file + " -out " + ref_db + " -dbtype prot", shell=True )
	p.communicate()
	
	summary_file = output_folder + "summary.txt"
	with open( summary_file, "w" ) as out:
		out.write( "FileID\tBestHit\tFunctionalAnnotation\n" )
		for idx, input_file in enumerate( input_files ):
			print( str( idx+1 ) + "/" + str( len( input_files ) ) )
			ID = input_file.split('/')[-1].split('.')[0]
			blast_result_file = output_folder + ID + ".txt"
			if not os.path.isfile( blast_result_file ):
				p = subprocess.Popen( args= "blastp -query " + input_file + " -db " + ref_db + " -out " + blast_result_file + " -outfmt 6 -evalue 0.001 -num_threads 4", shell=True )
				p.communicate()
			best_hit_per_query = load_best_hit_per_query( blast_result_file )
			most_abundant_best_hit = find_most_abundant_hit( best_hit_per_query )
			try:
				out.write( ID + "\t" + most_abundant_best_hit + "\t" + anno[ most_abundant_best_hit ] + "\n" )
			except KeyError:
				out.write( ID + "\t" + most_abundant_best_hit + "\tn/a\n" )


if '--in' in sys.argv and '--ref' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
