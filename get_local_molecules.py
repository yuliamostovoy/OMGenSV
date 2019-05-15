import sys, argparse, os
from os import path
import subprocess
import pybedtools
import datetime

def main():
	parser = argparse.ArgumentParser(description="Given an assembly output directory from Bionano Saphyr, produce a cmap file containing molecules that are local to one or more regions")
	parser.add_argument('assembly_dir', help="Path to Bionano Saphyr assembly directory (one level above the 'contigs' directory)")
	parser.add_argument('coords_bedfile', help="Bedfile containing region(s) from which to extract local molecules. If multiple regions are defined, the results will all be grouped together. Note that chromosome IDs must match the IDs used in the reference cmap (i.e. they must be numerical)")
	parser.add_argument('output_dir', help="Output directory")
	parser.add_argument('OMTools_jar', help='Location of OMTools.jar file')
	parser.add_argument('-j', help='Location of java to use with OMTools.jar (default is just "java" assumed to be in the user\'s path)', dest="java", default="java")
	parser.add_argument('-m', help='Directory to which molecules.tar.gz was extracted (will speed up script)', dest="mol_aln_dir")
	parser.add_argument('-p', dest='prefix', help="Prefix for output cmap file", default='')
	args = parser.parse_args()

	check_parameters(args)

	coords = pybedtools.BedTool(args.coords_bedfile)

	molecules = set()
	contig_regions = [] #list of lists with chr, start, end
	# get local molecules for each region
	sys.stderr.write("Getting molecules aligned to ref %s\n" % datetime.datetime.now())
	for chrom, start, end in coords:
		start, end = int(start), int(end)
		# define the contig regions that align to the specific region on the reference
		contig_regions += crossref_xmap_ref_to_query("%s/contigs/exp_refineFinal1/alignref_final" % args.assembly_dir, chrom, start, end) #list of lists with chr, start, end
		contig_regions += crossref_xmap_ref_to_query("%s/contigs/exp_refineFinal1_sv" % args.assembly_dir, chrom, start, end, full=True) #list of lists with chr, start, end

		# get molecules aligned directly to the specified region on the reference
		molecules = molecules.union(get_overlapping_queries("%s/contigs/alignmolvref/merge/alignmolvref_contig%s.xmap" % (args.assembly_dir, chrom), chrom, start, end))

	sys.stderr.write("Getting molecules aligned to contigs %s\n" % datetime.datetime.now())

	if not args.mol_aln_dir:
		# extract the molecule-to-contig alignment files for the relevant contigs
		contigs = [x[0] for x in contig_regions]
		mol_aln_dir = "%s/molecule_contig_alignments" % args.output_dir
		extract_molecule_alignments_to_contigs(args.assembly_dir, contigs, mol_aln_dir)
		sys.stderr.write("Done extracting %s\n" % datetime.datetime.now())
	else:
		mol_aln_dir = args.mol_aln_dir

	# get molecules aligned to contig regions that aligned to the specified regions on the reference
	contig_regions = merge_bed_regions(contig_regions)
	for contig, contig_start, contig_end in contig_regions:
		molecules = molecules.union(get_overlapping_queries("%s/exp_refineFinal1_contig%s.xmap" % (mol_aln_dir, contig), contig, int(contig_start), int(contig_end)))

	# write output bnx
	sys.stderr.write("Writing BNX file... %s\n" % datetime.datetime.now())
	extract_molecules_from_bnx("%s/contigs/auto_noise/autoNoise1_rescaled.bnx" % args.assembly_dir, molecules, "%s/%slocal_molecules.bnx" % (args.output_dir, args.prefix))
	# convert to cmap format
	sys.stderr.write("Converting to CMAP file...\n")
	subprocess.call([args.java, '-jar', args.OMTools_jar, 'DataTools', '--optmapin', "%s/%slocal_molecules.bnx" % (args.output_dir, args.prefix), '--optmapout', "%s/%slocal_molecules.cmap" % (args.output_dir, args.prefix)], stdout=None, stderr=None, stdin=None)
	sys.stderr.write("Finished!\n")

def extract_molecules_from_bnx(bnx_filename, molecules, output_filename):
	bnx_file = open(bnx_filename, 'r')
	output_file = open(output_filename, 'w')

	to_print = False
	for line in bnx_file:
		if line.startswith('#'):
			output_file.write("%s" % line)
			continue
		if line.startswith('0'):
			line = line.strip().split('\t')
			ID = line[1]
			if ID in molecules:
				to_print = True
				output_file.write("%s\n" % '\t'.join(line))
			else:
				to_print = False
		else:
			if to_print:
				output_file.write("%s" % line)
	bnx_file.close()
	output_file.close()

def get_overlapping_queries(xmap_filename, rchr, rstart, rend):
	xmap_file = open(xmap_filename, 'r')
	queries = set()

	for line in xmap_file:
		if line.startswith('#'):
			continue
		line = line.strip().split('\t')
		if line[2]==rchr and float(line[5])<rend and float(line[6])>rstart:
			queries.add(line[1])

	return queries

def merge_bed_regions(contig_regions):
	# convert list of lists to a format readable by pybedtools
	input_string = ""
	for chrom, start, end in contig_regions:
		input_string += "%s %.0f %.0f\n" % (chrom, start, end)

	# create bed file and merge it
	bed = pybedtools.BedTool(input_string, from_string=True)
	bed = bed.sort()
	bed = bed.merge()

	return bed	

def crossref_xmap_ref_to_query(directory, rchr, rstart, rend, full=False):
	# given reference positions and an alignment, return the aligned query positions
	if full:
		xmap_file = open("%s/EXP_REFINEFINAL1_full.xmap" % directory)
	else:
		xmap_file = open("%s/EXP_REFINEFINAL1.xmap" % directory)
	q_cmap_file = open("%s/EXP_REFINEFINAL1_q.cmap" % directory)
	r_cmap_file = open("%s/EXP_REFINEFINAL1_r.cmap" % directory)

	# convert ref position to label
	rstart_label, rend_label = convert_pos_to_labels(rchr, rstart, rend, r_cmap_file)
	# get the corresponding query labels
	query_aligned_labels = crossref_labels_from_xmap(rchr, rstart, rend, rstart_label, rend_label, xmap_file)
	# get the corresponding query positions
	query_positions = convert_labels_to_pos_batch(query_aligned_labels, q_cmap_file)

	xmap_file.close()
	q_cmap_file.close()
	r_cmap_file.close()

	return query_positions

def convert_labels_to_pos_batch(aligned_labels, cmap):
	# aligned_labels is list of tuples containing (contig, label_start, label_end)

	#parse cmap
	label_to_pos = {} #contig_ID-->label-->pos (dict of dicts)
	for line in cmap:
		if line.startswith("#"):
			continue
		line = line.strip().split('\t')
		if line[4]=="0": #label channel
			continue
		contig_ID, label, pos = line[0], int(line[3]), int(round(float(line[5])))
		if contig_ID not in label_to_pos:
			label_to_pos[contig_ID] = {label:pos}
		else:
			label_to_pos[contig_ID][label] = pos

	aligned_coords = []
	for contig, label_start, label_end in aligned_labels:
		if contig not in label_to_pos:
			sys.stderr.write("Warning: aligned contig ID %d not found in provided cmap\n" % contig)
			continue
		if label_start not in label_to_pos[contig]:
			sys.stderr.write("Error: label number %d not found for contig %d in provided cmap\n" % (label_start, contig))
			sys.exit(1)
		if label_end not in label_to_pos[contig]:
			sys.stderr.write("Error: label number %d not found for contig %d in provided cmap\n" % (label_end, contig))
			sys.exit(1)
		aligned_coords.append([contig, label_to_pos[contig][label_start], label_to_pos[contig][label_end]])
	return aligned_coords

def crossref_labels_from_xmap(chrom, rstart, rend, rstart_label, rend_label, xmap, strict=False):
	q_data = [] #list of tuples (query, start_label, end_label)
	q_data_tmp = {} #dict of start, end labels where keys are query IDs
	for line in xmap:
		if line.startswith("#"):
			continue
		line = line.strip().split('\t')
		line_rchr, line_rstart, line_rend, ori = line[2], float(line[5]), float(line[6]), line[7]
		if line_rchr != chrom:
			continue
		if not (line_rend > rstart and line_rstart < rend):
			continue
		labels_ref_to_query = process_label_string(line[13])
		if strict:
			if rstart_label in labels_ref_to_query:
				qstart_label = labels_ref_to_query[rstart_label]
			else:
				continue
		else:
			if min(labels_ref_to_query.keys()) > rstart_label:
				qstart_label = labels_ref_to_query[min(labels_ref_to_query.keys())]
			elif rstart_label in labels_ref_to_query:
				qstart_label = labels_ref_to_query[rstart_label]
			elif rstart_label-1 in labels_ref_to_query:
				qstart_label = labels_ref_to_query[rstart_label-1]
			elif rstart_label-2 in labels_ref_to_query:
				qstart_label = labels_ref_to_query[rstart_label-2]
			elif rstart_label+1 in labels_ref_to_query:
				qstart_label = labels_ref_to_query[rstart_label+1]
			else:
				continue

		if strict:
			if rend_label in labels_ref_to_query:
				qend_label = labels_ref_to_query[rend_label]
			else:
				continue
		else:
			if max(labels_ref_to_query.keys()) < rend_label:
				qend_label = labels_ref_to_query[max(labels_ref_to_query.keys())]
			elif rend_label in labels_ref_to_query:
				qend_label = labels_ref_to_query[rend_label]
			elif rend_label+1 in labels_ref_to_query:
				qend_label = labels_ref_to_query[rend_label+1]
			elif rend_label+2 in labels_ref_to_query:
				qend_label = labels_ref_to_query[rend_label+2]
			elif rend_label-1 in labels_ref_to_query:
				qend_label = labels_ref_to_query[rend_label-1]
			else:
				continue

		if ori=='-':
			qstart_label, qend_label = qend_label, qstart_label

		q_data.append((line[1], min(qstart_label), max(qend_label)))

	return q_data

def process_label_string(label_str):
	label_str=label_str.rstrip(')').replace('(','').split(')')#list of 'x,y' strings, format is ref,query, where ref is always ascending but query is descending if - orientation
	labels_ref_to_query={}
	for entry in label_str:
		entry=entry.split(',')
		r,q = int(entry[0]), int(entry[1])
		if r in labels_ref_to_query:
			labels_ref_to_query[r].append(q)
		else:
			labels_ref_to_query[r] = [q]
	return labels_ref_to_query

def convert_pos_to_labels(chrom, start, end, cmap):
	if start > end:
		start, end = end, start
	start_label, end_label = None, None
	i=0

	for line in cmap:
		if line.startswith("#"):
			continue
		line = line.strip().split('\t')
		cmap_chr, pos = line[0], float(line[5])
		if cmap_chr != chrom:
			continue
		if start_label == None and pos >= start:
			start_label = i
		if end_label == None:
			if pos == end:
				end_label = i
			elif pos > end:
				end_label = i-1
		i+=1

	return start_label+1, end_label+1 #counting 0-based, returning results as 1-based (same as in cmap)

def extract_molecule_alignments_to_contigs(assembly_dir, contigs, mol_aln_dir):
	# extract the molecule-to-contig alignment files for a given list of contigs
	if not os.path.exists("%s" % mol_aln_dir):
		os.mkdir("%s" % mol_aln_dir)
	for contig in contigs:
		try:
			subprocess.check_call(['tar', 'zxf', '%s/molecules.tar.gz' % assembly_dir, '-C', '%s' % mol_aln_dir, 'exp_refineFinal1_contig%s*map' % contig], stdout=None, stderr=None, stdin=None)
		except:
			sys.exit(1)

def check_parameters(args):
	# Clean up directory names
	args.assembly_dir = args.assembly_dir.rstrip('/')
	args.output_dir = args.output_dir.rstrip('/')
	args.mol_aln_dir = args.mol_aln_dir.rstrip('/')

	# check paths
	if not os.path.exists(args.assembly_dir):
		sys.stderr.write("Error: assembly directory does not exist: %s\n" % args.assembly_dir)
		sys.exit(1)
	if not os.path.exists('%s/contigs' % args.assembly_dir):
		sys.stderr.write("Error: assembly directory does not contain the 'contigs' directory, please check the path: %s\n" % args.assembly_dir)
		sys.exit(1)
	if not os.path.exists(args.output_dir):
		sys.stderr.write("Error: output directory does not exist: %s\n" % args.output_dir)
		sys.exit(1)

	if args.prefix:
		args.prefix += '_'


if __name__=="__main__": main()
