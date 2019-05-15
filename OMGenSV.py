import sys, argparse, os
from os import path
from argparse import RawTextHelpFormatter
import subprocess

def main():
	parser = argparse.ArgumentParser(description="Script for using single-molecule optical maps to genotype structural variation\n", formatter_class=RawTextHelpFormatter)
	parser.add_argument('sample_list', help='list of samples to be processed, one sample name per line')
	parser.add_argument('molecules_dir', help='directory containing CMAP files of single molecules for each sample')
	parser.add_argument(metavar='critical_regions_bedfile', dest='hapinfo', help='bedfile containing the following info in tab-delimited format, one line per haplotype:\n  haplotype_cmap_ID critical_region_start critical_region_end haplotype_name')
	parser.add_argument('out_dir', help='output directory')
	parser.add_argument('-a', dest='run_alignment', help='Tells script to run alignments with OMTools. If not set, script will look for alignments in the standard location (output_directory/samplename.oma)', action='store_true')
	parser.add_argument('--suffix', dest='molecules_file_suffix', help='Suffix for molecules files, to be combined with molecules_dir and sample name for the full molecule filename, e.g.:\n  molecules_dir=/home/user/local_molecules\n  Sample list contains sample "NA12878"\n  molecules_file_suffix=_local_molecules.cmap\n  Final path to molecules file: /home/user/local_molecules/NA12878_local_molecules.cmap\nDefault: "_local_molecules.cmap"', default="_local_molecules.cmap")
	parser.add_argument('--jar', dest='OMTools_jar', help='Location of OMTools.jar file (only required if -a is set)')
	parser.add_argument('--refs', dest='haplotypes_cmap', help='CMAP file containing the haplotypes to which molecules will be aligned')
	parser.add_argument('-t', dest='threads', help='Number of threads to use with OMTools aligner. Default is 1.', default=1, type=int)
	parser.add_argument('-j', help='Optional: location of java to use with OMTools.jar (default is just "java" assumed to be in the user\'s path)', dest="java", default="java")
	parser.add_argument('--min_perc', dest='min_perc', help='If a sample has an allele with molecule support below %%min_perc of its best-supported allele, molecules supporting the former will be flagged for manual evaluation (range: 0-1, default: 0.25)', type=float, default=0.25)
	parser.add_argument('--min_mols', dest='min_mols', help='An allele with molecule support below min_mols will have its molecules flagged for manual evaluation (default: 3)', type=int, default=3)
	parser.add_argument('--min_score', dest='min_score', help='Minimum cumulative alignment score for reporting hits. Default: 0', type=float, default=0.0)
	args = parser.parse_args()

	check_parameters(args)

	# Read initial input files
	samples = read_samples(args.sample_list)
	hapinfo = read_hapinfo(args.hapinfo)

	# Get hits for each sample
	log = open("%s/alignments.log" % args.out_dir,'a')
	hits = {} #sample-->haplotype ID-->list of molecules supporting that haplotype
	for sample in samples:
		sys.stderr.write("Processing sample %s\n" % sample)
		hits = get_hits_wrapper(sample, hits, args, log, hapinfo)
	log.close()

	# Print output files
	print_hits(hits, '%s/hits' % args.out_dir)
	haps = [x[3] for x in hapinfo]
	print_results(hits, haps, '%s/results' % args.out_dir)

	# Find samples that should be manually evaluated and write files to visualize them
	min_perc = 0.25 # If minor allele has fewer supporting molecules than min_perc * major allele, evaluate the former manually
	min_mols = 3 # If allele has fewer supporting molecules than min_mols, evaluate them manually
	to_evaluate = select_molecules_to_evaluate(hits, min_perc, min_mols) #dict of hap-->list of samples
	create_evaluation_files(to_evaluate, args, hits)

	sys.stderr.write("Finished!\n")

def check_parameters(args):
	if args.run_alignment and not args.OMTools_jar:
		sys.stderr.write("Error: -a is set to run alignments, but the location of OMTools.jar isn't specified. Please add '--jar path_to_OMTools.jar' to run alignments\n")
		sys.exit(1)  

	if args.run_alignment and not args.haplotypes_cmap:
		sys.stderr.write("Error: -a is set to run alignments, but the location of the haplotypes CMAP isn't specified. Please add '--refs path_to_haplotypes.cmap' to run alignments\n")
		sys.exit(1)

	if not os.path.exists(args.molecules_dir):
		sys.stderr.write("Error: molecules directory does not exist: %s\n" % args.molecules_dir)
		sys.exit(1)

	if not os.path.exists(args.out_dir):
		sys.stderr.write("Error: output directory does not exist: %s\n" % args.out_dir)
		sys.exit(1)

	if args.min_perc < 0 or args.min_perc > 1:
		sys.stderr.write("Error: min_perc must be between 0 and 1.\n")
		sys.exit(1)

	# Clean up directory names
	args.out_dir = args.out_dir.rstrip('/')
	args.molecules_dir = args.molecules_dir.rstrip('/')

def create_evaluation_files(to_evaluate, args, hits):
	#renumbering the molecules so that there are no repeats when combining different samples
	crossref = {} #sample-->mol ID-->renumbered ID
	new_ID_counter = 1
	for hap in to_evaluate:
		oma_outfile = open('%s/%s_hits_to_evaluate.oma' % (args.out_dir, hap), 'w')
		qcmap_outfile = open('%s/%s_hits_to_evaluate_q.cmap' % (args.out_dir, hap), 'w')
		print_headers = True #on first pass, need to print headers
		samples = to_evaluate[hap]

		for sample in samples:
			# create cross-referenced IDs
			if sample not in crossref:
				crossref[sample] = {}
			for mol in hits[sample][hap]:
				crossref[sample][mol] = '%d' % new_ID_counter
				new_ID_counter += 1

			write_subset_of_queries_in_cmap(hits[sample][hap], '%s/%s%s' % (args.molecules_dir, sample, args.molecules_file_suffix), qcmap_outfile, crossref[sample], print_headers)
			write_subset_of_queries_in_oma(hits[sample][hap], '%s/%s_best.oma' % (args.out_dir, sample), oma_outfile, crossref[sample], print_headers)
			print_headers = False
		oma_outfile.close()
		qcmap_outfile.close()
		
	# print crossref info
	crossref_filename='%s/to_evaluate_mol_IDs_crossreference' % args.out_dir
	crossref_outfile = open(crossref_filename,'w')
	crossref_outfile.write("Sample\toriginal_mol_ID\tnew_mol_ID\n")
	for sample in crossref:
		for mol in crossref[sample]:
			crossref_outfile.write("%s\t%s\t%s\n" % (sample, mol, crossref[sample][mol]))
	crossref_outfile.close()
	subprocess.call(['sort', '-k3,3n', crossref_filename, '-o', crossref_filename], stdout=None, stderr=None, stdin=None)

def write_subset_of_queries_in_oma(queries, oma_filename, out_file, crossref_mols, print_headers):
	oma_file = open(oma_filename, 'r')
	for line in oma_file:
		if line.startswith('#'):
			if print_headers:
				out_file.write('%s' % line)
			continue
		line = line.strip().split('\t')
		if line[0] in queries:
			line[0] = crossref_mols[line[0]]
			out_file.write('%s\n' % '\t'.join(line))
	oma_file.close()

def write_subset_of_queries_in_cmap(queries, cmap_filename, out_file, crossref_mols, print_headers):
	cmap_file = open(cmap_filename,'r')
	for line in cmap_file:
		if line.startswith('#'):
			if print_headers:
				out_file.write('%s' % line)
			continue
		line = line.strip().split('\t')
		if line[0] in queries:
			line[0] = crossref_mols[line[0]]
			out_file.write('%s\n' % '\t'.join(line))
	cmap_file.close()

def select_molecules_to_evaluate(hits, min_perc, min_mols):
	to_evaluate = {} # hap-->samples
	for sample in hits:
		current_hits = hits[sample] #dict of hap-->list of molecules
		major = max([len(x) for x in current_hits.values()]) # get the number of molecules for the haplotype with the most hits for this sample
		if not major:
			continue
		for hap in current_hits:
			num_mols = len(current_hits[hap])
			if num_mols==0:
				continue
			if num_mols < min_mols or float(num_mols)/major < min_perc:
				if hap not in to_evaluate:
					to_evaluate[hap] = [sample]
				else:
					to_evaluate[hap].append(sample)
	return to_evaluate

def run_alignment(log, args, sample):
	log.write("########## %s ##########\n" % sample)
	log.flush()
	try:
		subprocess.check_call([args.java, '-jar', args.OMTools_jar, 'OMBlastMapper', '--refmapinformat', '6', '--optmapinformat', '6', '--refmapin', '%s' % args.haplotypes_cmap, '--optmapin', '%s/%s%s' % (args.molecules_dir, sample, args.molecules_file_suffix), '--optresout', '%s/%s.oma' % (args.out_dir, sample), '--thread', '%d' % args.threads, '--optresoutformat', '0', '--writeunmap', 'false', '--alignmentjoinmode', '1', '--filtermode', '1', '--trimmode', '1', '--minconf', '0', '--minjoinscore', '0', '--maxalignitem', '2', '--meas', '600', '--ear', '0.5'], stdout=log, stderr=None, stdin=None)
	except subprocess.CalledProcessError as exc:
		sys.stderr.write("OMTools returned an error code:\n%s\n" % exc.returncode)
		sys.exit(1)

def print_hits(hits, outfilename):
	outfile = open(outfilename, 'w') # overwrites previous file
	for sample in hits:
		for hap in hits[sample]:
			outfile.write("%s\t%s\t%s\n" % (sample, hap, ','.join(hits[sample][hap])))
	outfile.close()

def print_results(hits, haps, outfilename):
	outfile = open(outfilename, 'w') # overwrites previous file
	outfile.write("Sample\t%s\n" % '\t'.join(haps)) #header
	for sample in hits:
		outfile.write("%s" % sample)
		for hap in haps:
			outfile.write("\t%d" % len(hits[sample][hap]))
		outfile.write('\n')
	outfile.close()

def get_hits(oma_filename, qcoords, chrom, start, end):
	mols = {}
	indels = {} #by mol, list indels where the size change between mol and ref is larger than max_size_change
	max_size_change=50000

	omafile = open(oma_filename, 'r')
	for line in omafile:
		if line.startswith("#"):
			continue
		line = line.strip().split()
		#QueryID    QuerySeg    QuerySegInfo    RefID   Strand  Score   Confidence  RefSegStart RefSegStop  QuerySegStart   QuerySegStop    RefStartCoord   RefStopCoord    Cigar
		if line[3] != chrom:
			continue
		if line[0] not in mols:
			qe_seg = max(int(line[9]), int(line[10]))
			qe = qcoords[line[0]][qe_seg-1]
			mols[line[0]] = {'start':float(line[11]), 'end':float(line[12]), 'conf':float(line[6]), 'qe':qe}
		else:
			### indels present ###
			rs, re = float(line[11]), float(line[12]) # these are always lower to higher regardless of orientations
			qs_seg = min(int(line[9]), int(line[10]))
			qs = qcoords[line[0]][qs_seg-1]
			strand=line[4]
			if line[0] not in indels:
				indels[line[0]] = [{'qlen':qs-mols[line[0]]['qe']}]
			else:
				indels[line[0]].append({'qlen':qs-mols[line[0]]['qe']})
			if strand=='+':
				indels[line[0]][-1]['rs'] = mols[line[0]]['end']
				indels[line[0]][-1]['re'] = rs
			else:
				indels[line[0]][-1]['rs'] = re
				indels[line[0]][-1]['re'] = mols[line[0]]['start']

			if rs < mols[line[0]]['start']:
				mols[line[0]]['start'] = rs
			if re > mols[line[0]]['end']:
				mols[line[0]]['end'] = re
			qe_seg = max(int(line[9]), int(line[10]))
			mols[line[0]]['qe'] = qcoords[line[0]][qe_seg-1]

	omafile.close()

	final_mols = []
	flagged = 0
	for mol in mols:
		skip=False
		if mols[mol]['start'] <= start+1 and mols[mol]['end'] >= end-1: #adding 1bp buffer to avoid rounding error
			if mol in indels:
				for indel in indels[mol]:
					if abs(indel['qlen']-(indel['re']-indel['rs'])) > max_size_change:
						if indel['rs'] < end and indel['re'] > start:
							skip=True
							break
			if not skip and mols[mol]['conf']>0:
				final_mols.append(mol)
			elif not skip and mols[mol]['conf']==0:
				flagged += 1
	if flagged:
		sys.stderr.write("WARNING: file %s had %d molecules that matched the following critical region but were filtered out for having >1 equally good alignments. Suggestion: check that this critical region isn't also present elsewhere in the reference cmap: %s\t%d\t%d\n" % (oma_filename, flagged, chrom, start, end))
	return final_mols

def get_hits_wrapper(sample, hits, args, log, hapinfo):
	if sample in hits:
		sys.stderr.write("Warning: samples list contains %s twice (skipping second one)\n" % sample)
		return
	hits[sample] = {}

	if args.run_alignment:
		run_alignment(log, args, sample)

	# Get best hit for each molecule, filtering out ones with multiple best hits
	process_oma_get_best('%s/%s.oma' % (args.out_dir, sample), args.min_score, '%s/%s_best.oma' % (args.out_dir, sample))

	# Get hits for each haplotype that span the critical region(s)
	qcoords = process_query_file('%s/%s%s' % (args.molecules_dir, sample, args.molecules_file_suffix))
	for cmap_ID, CR_start, CR_end, hapname in hapinfo:
		current_hits = get_hits('%s/%s_best.oma' % (args.out_dir, sample), qcoords, cmap_ID, CR_start, CR_end)
		if hapname not in hits[sample]:
			hits[sample][hapname] = current_hits
		else:
			hits[sample][hapname] += current_hits
			hits[sample][hapname] = list(set(hits[sample][hapname]))
	return hits

def process_query_file(query_filename):
	qcmapfile = open(query_filename,'r')
	qcoords = {}
	for line in qcmapfile:
		if line.startswith('#'):
			continue
		line=line.strip().split('\t')
		if line[0] not in qcoords:
			qcoords[line[0]] = [float(line[5])]
			continue
		if line[4] == '0':
			continue
		qcoords[line[0]].append(float(line[5]))
	qcmapfile.close()
	return qcoords

def read_samples(samples_filename):
	samples = []
	try:
		samplefile=open(samples_filename,'r')
	except:
		sys.stderr.write("Error: unable to open samples file: %s\n" % samples_filename)
		sys.exit(1)
	for line in samplefile:
		samples.append(line.strip())
	samplefile.close()
	return samples

def read_hapinfo(hapinfo_filename):
	hapinfo = [] #list of lists
	try:
		hapfile=open(hapinfo_filename,'r')
	except:
		sys.stderr.write("Error: unable to open haplotype info file: %s\n" % hapinfo_filename)
		sys.exit(1)
	for line in hapfile:
		line=line.strip().split('\t')
		hapinfo.append([line[0], float(line[1]), float(line[2]), line[3]])
	hapfile.close()
	hapinfo.sort(key=lambda x: x[3]) #sort of hap name
	return hapinfo

def get_best_hit(hits, confs, minscore, outfile):
	# Find the best hit for each molecule, if there is one
	flag_dup = False
	if len(list(set(confs)))==1:
		if len(hits)>1:
			flag_dup = True
		else:
			return

	indices_to_print = [counter for counter, value in enumerate(confs) if value==max(confs)]

	# calculate score of best hit
	high_score=0
	for index in indices_to_print:
		high_score += float(hits[index][5])
	high_score -= 10*(len(indices_to_print)-1)

	# print best hit if its score is above the minimum
	if high_score >= minscore:
		for index in indices_to_print:
			if flag_dup:
				hits[index][6] = '0' # to mark an alignment with 2 equally good alignments, set its confidence to 0
			outfile.write("%s\n" % '\t'.join(hits[index]))

def process_oma_get_best(omapath, minscore, outfilename):
	if not os.path.isfile(omapath):
		sys.stderr.write("Error: alignment file %s not found. Use -a option to run alignment(s).\n" % omapath)
		sys.exit(1)
	omafile = open(omapath,'r')
	outfile = open(outfilename, 'w')
	
	hits = []
	confs = []
	query = None

	for line in omafile:
		if line.startswith("#"):
			outfile.write("%s" % line)
			continue
		line = line.strip().split()
		#columns are:
		#QueryID    QuerySeg    QuerySegInfo    RefID   Strand  Score   Confidence  RefSegStart RefSegStop  QuerySegStart   QuerySegStop    RefStartCoord   RefStopCoord    Cigar
		if not query:
			query = line[0]
		elif query != line[0]:
			get_best_hit(hits, confs, minscore, outfile)
			query = line[0]
			hits, confs = [], []
		hits.append(line)
		confs.append(float(line[6]))

	get_best_hit(hits, confs, minscore, outfile) #process the last entry
	omafile.close()
	outfile.close()

if __name__=="__main__": main()
