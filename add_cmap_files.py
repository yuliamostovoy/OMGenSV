import sys
import argparse

def main():
	parser = argparse.ArgumentParser(description="Add two cmap entries together (must be one per file) with specified gap size and orientations\n")
	parser.add_argument('cmapfile1',help="name of first cmap file")
	parser.add_argument('cmapfile2',help="name of second cmap file")
	parser.add_argument('-g', '--gap_length', dest='gapsize', help="Length of gap between cmap1 and cmap2 (default 0)", type=float)
	parser.add_argument('-r','--reverse', nargs='+', dest="reverse", metavar="cmap_to_reverse", help="Indicate which cmaps to reverse (1 and/or 2) (default is none)", type=int, default=[])
	parser.add_argument('-o', '--new_id', dest='new_id', type=int, help="New ID to assign to output cmap (default is 1)")

	args = parser.parse_args()

	cmapfile1=open(args.cmapfile1,'r')
	cmapfile2=open(args.cmapfile2,'r')
	if args.gapsize:
		gapsize=float(args.gapsize)
	else:
		gapsize=0

	if args.new_id:
		new_id=args.new_id
	else:
		new_id=1

	if args.reverse:
		if len(args.reverse)>2:
			sys.stderr.write("Cmaps to reverse can be at most two (1, 2)\n")
			sys.exit(1)
	for i in args.reverse:
		if i not in [1,2]:
			sys.stderr.write("Error: -r parameter must be '1' and/or '2'\n")
			sys.exit(1)

	for line in cmapfile1:
		if line.startswith("#"):
			print line.strip()
		else:
			break

	line1_c1=line.strip().split('\t')
	length_c1=float(line1_c1[1])
	num_labels_c1=int(line1_c1[2])

	for line in cmapfile2:
		if line.startswith("#"):
			continue
		else:
			break

	line1_c2=line.strip().split('\t')
	length_c2=float(line1_c2[1])
	num_labels_c2=int(line1_c2[2])

	#CMapId	ContigLength(bp)	NumSites	SiteID	LabelChannel	Position	StdDev	Coverage	Occurrence
	firstcols="%s\t%.1f\t%d" % (new_id, length_c1+gapsize+length_c2, num_labels_c1+num_labels_c2)
	addtosecond=gapsize+length_c1
	cmap1=[float(line1_c1[5])]
	for line in cmapfile1:
		line=line.strip().split('\t')
		cmap1.append(float(line[5]))	

	if 1 in args.reverse:
		label = print_reverse(cmap1, extra=0, firstcols=firstcols, label=1)
	else:
		label = print_simple(cmap1, extra=0, firstcols=firstcols, label=1)

	cmap2=[float(line1_c2[5])]
	for line in cmapfile2:
		line=line.strip().split('\t')
		cmap2.append(float(line[5]))

	if 2 in args.reverse:
		print_reverse(cmap2, extra=addtosecond, firstcols=firstcols, label=label, second=True)
	else:
		print_simple(cmap2, extra=addtosecond, firstcols=firstcols, label=label, second=True)

def print_simple(positions, extra, firstcols, label, second=False):
	if extra and not second:
		sys.stderr.write("Doesn't make sense: padding on first cmap\n")
		sys.exit(1)
	for pos in positions[:-1]: #skip last one for now
		print "%s\t%d\t1\t%.1f\t1\t1\t1" % (firstcols, label, pos + extra)
		label += 1
	if second:
		print "%s\t%d\t0\t%.1f\t0\t0\t0" % (firstcols, label, positions[-1] + extra)
	return label

def print_reverse(positions, extra, firstcols, label, second=False):
	if extra and not second:
		sys.stderr.write("Doesn't make sense: padding on first cmap\n")
		sys.exit(1)
	positions.reverse()
	prev_pos=positions[0]
	prev_total=0
	for pos in positions[1:]:
		print "%s\t%d\t1\t%.1f\t1\t1\t1" % (firstcols, label, prev_pos-pos+prev_total+extra)
		label += 1
		prev_total=prev_pos-pos+prev_total
		prev_pos=pos
	if second:
		print "%s\t%d\t0\t%.1f\t0\t0\t0" % (firstcols, label, prev_total + positions[-1] + extra)
	return label

if __name__=="__main__": main()
