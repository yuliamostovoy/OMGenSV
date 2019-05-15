# OMGenSV: a toolkit for using optical maps to genotype structural variation

## Prerequisites:

Python 2.7.x

OMTools >=v1.3 (https://github.com/TF-Chan-Lab/OMTools)

Python packages: 
- argparse
- pybedtools (https://daler.github.io/pybedtools/main.html)

## Usage
### Step 1: Identify local configurations and create corresponding cmaps
First, you'll need to identify structural variation at your locus of interest by examining assembled contigs that align to the locus.

Next, create accurate and complete CMAP representations of each variant, using the reference CMAP and assembled contig CMAPs as templates to be manipulated with the following tools:

get_cmap_subsets.py - Create a new CMAP from a specified region of a larger CMAP

add_cmap_files.py   - Add two CMAP files together into a larger CMAP, with options to reverse one or both of the input CMAPs and to add any amount of label-less space between the two CMAPs.

Usage example: You find an assembled contig that has an inversion at 1:520,000-550,000. The contig is short and spans from 1:510,000-560,000, so in order to include sufficient flanking labels for molecules to align to, you run the following commands to expand the contig using the reference genome and label the final contig with the CMAP ID '1001':

```
echo $'1\t310000\t509999' > upstream_coords.bed
echo $'1\t560001\t760000' > downstream_coords.bed

python get_cmap_subsets.py hg38_chr1.cmap upstream_coords.bed > upstream.cmap
python get_cmap_subsets.py hg38_chr1.cmap downstream_coords.bed > downstream.cmap

python add_cmap_files.py upstream.cmap inversion_contig.cmap > inversion_with_upstream.cmap
python add_cmap_files.py inversion_with_upstream.cmap downstream.cmap -o 1001 > inversion_with_flanking.cmap
```

### Step 2: Create a single reference CMAP file
Decide which CMAPs will be analyzed together and put them into a single CMAP file to serve as the reference for molecule alignment. Some guidelines:
  - Each entry in the file should have a different CMAP ID
  - Focus on a single region at a time. If the locus being studied is long (i.e. several hundred kb, or more than the length of a typical single molecule in your data), it may be best to analyze different parts of it separately. For example, a strategy that we've used successfully for regions with long segmental duplications has been to identify variants anchored in either the upstream or downstream unique flanking areas, and analyze those two sets separately.
  - Each analysis must include at least two structurally different versions of the same region. The simplest case would be the reference configuration vs. an alternate configuration. 
  - If the region has very high sequence identity to another region elsewhere, consider including the duplicate region as a decoy to prevent its molecules from aligning to your region of interest and giving false positive results.

### Step 3: Identify the 'critical regions' for each variant that you want to genotype
The critical regions should be recorded in BED format (with columns CMAP ID, start, end). This is the minimal region that molecules must be aligned to in order to support a given variant. These regions can be modified as needed later. Some guidelines:

  - The label pattern represented by each critical region should be present only once within your reference CMAP file. 
  - Any decoy entries don't need to be represented. 
  - The longer the critical region, the fewer hits you are likely to get, especially when it exceeds the length of the average molecule in your dataset. For a dataset with an average molecule length of ~250kb, we've found that critical regions of up to 300kb are ideal, while exceeding 400kb significantly reduces the number of hits. Keep this issue in mind when comparing the number of hits for critical regions of different lengths.

### Step 4: Isolate the local molecules from your region of interest for each sample
Pull out the molecules that are local to the region of interest from the Bionano _de novo_ assembly output directory:

a) Define the local region in a 3-column BED file. We recommend padding your locus with a few hundred kb of flanking regions to be comprehensive.

b) To speed up this step, you may choose to extract the molecules.tar.gz file from the assembly output directory. Otherwise, the script will extract only the required files, which will save disk space but may take longer.

c) Run the extraction command:
```
python get_local_molecules.py assembly_dir coords_bedfile output_dir path_to/OMTools.jar [-j java_location -m molecules_extraction_folder -p output_file_prefix]
```
To run on many samples, you can do something like this to write the results to the same output folder, with the assembly directory names as prefixes:

assembly_dirs = list of all the assembly directory names
```
for i in `cat assembly_dirs`; do python get_local_molecules.py $i coords_bedfile output_dir path_to/OMTools.jar -p $i; done
```

### Step 5: Run the genotyping script
```
python OMGenSV.py sample_list local_molecules_dir critical_regions_bedfile output_dir -a --jar path_to/OMTools.jar -t 16 --refs configurations.cmap
```

This command will create a [sample].oma file and a [sample]_best.oma file for each sample, as well as a joint "results" file showing the number of supporting molecules for each critical region in the critical_regions_bedfile for each sample. The IDs of the supporting molecules are written to "hits".

The script assumes a specific convention for the local molecules CMAP files, namely that they are all present in the same directory (local_molecules_dir), with the same name format: [sample]_local_molecules.cmap, which is the output format of get_local_molecules.py. If your files have a different suffix, that can be specified with --suffix.

It's recommended that manual evaluation be done for poorly-supported configurations in each sample. OMGenSV defines these as configurations that are supported in a sample by fewer than 3 molecules, or that are supported by fewer than 25% of the number of molecules that support the best-supported configuration in the sample. You can change these defaults with the --min_mols and --min_perc flags. For each configuration, the script prepares files for manual evaluation by pooling the supporting molecules from each sample where the configuration was poorly-supported. These are written to [configuration_ID]_hits_to_evaluate.oma and [configuration_ID]_hits_to_evaluate_q.cmap and can be visualized as follows:

```
java -jar path_to/OMTools.jar --viewrefin configurations.cmap --viewmapin [configuration_ID]_hits_to_evaluate_q.cmap --viewresin [configuration_ID]_hits_to_evaluate.oma --viewannoin critical_regions_bedfile --viewregion [configuration_ID]:[start]-[end]
```

Make a note of which molecules have strong alignments to the configuration and which look questionable or poor. Because these evaluation files potentially compile molecules from multiple samples, which might have repeated molecule IDs, the script renumbers the molecules and creates a cross-reference file at to_evaluate_mol_IDs_crossreference. You can use this file to match molecules to sample IDs and then modify the results file accordingly.

To rerun OMGenSV.py without rerunning the alignments (e.g. you may want to experiment with different versions of critical_regions_bedfile or change parameters like --min_perc, --min_mols, or --min_score), just omit the -a option. The script will look for completed [sample].oma alignment files in the designated output directory. The following flags are only used for alignments and will be ignored if -a is not set: -j, -t, --refs, --jar

## Author
Yulia Mostovoy, in Pui-Yan Kwok's lab at the University of California, San Francisco
