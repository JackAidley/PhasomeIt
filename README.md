# PhasomeIt
PhasomeIt is a tool for locating, classfying, and grouping phase variable SSRs in small genomic collections. See the paper cited below for a full description.

At the time I wrote it for my PhD, there were less than a hundred Campylobacter genomes publicly available and the software was designed for collections on this sort of scale. The largest collection I used it with contained somewhat over 200 genomes. If you have access to a computing cluster it will handle larger numbers than this but the running time increases at a rate larger than the square of the number of genomes, which means it will quickly exhaust available processing and memory limits for large collections. Not only that but the output of the software is designed to be intelligable for the now relatively small collections available at the time so you may find that it is difficult to interpret if used with much larger genome collections.

Installation
------------

SYSTEM REQUIREMENTS

PhasomeIt is written in Python and thus platform-independent, however it is reliant on Bossref (included) which will only run on 64-bit Intel-compatible processors. The memory requirement will depend on the dataset you are processing but >4Gb is recommended for moderate size datasets (20-70 genomes) and increases from there. PhasomeIt has been tested under Windows and Linux (64-bit versions required). It is also expected to work under Mac OS X but this has not been tested.

INSTALL PRE-REQUISITES

1. Required: Python 3, version 3.9 or later (it may function with earlier versions)
	a. Go to http://python.org and download the latest version of Python 3 for your system
	b. Install according to the instructions on that site

2. Required: The Python libraries BioPython and NatSort
	a. To install BioPython, use Pip with the command line instruction 'pip install biopython'
	b. To install NatSort, use Pip with the instruction 'pip install NatSort'
	(If you don't know how to use pip, see https://packaging.python.org/tutorials/installing-packages/

3. Require: Install the BLAST+ suite
	a. Go to the NCBI site, download and install according to their instructions - https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download
	b. Ensure that the BLAST+ suite is added to the system path so that it can be invoked from anywhere

INSTALL PHASOMEIT

PhasomeIt does not require any special installation steps, merely unzip to your desired location and it's ready to go. You can also acquire the files from github, at //http://github.com/JackAidley/PhasomeIt and obtain them using git. If you don't know how to do this already I suggest you use the .zip release version

For Linux users: you may need to restore the executable bit on the 'bossref' executable, navigate to the PhasomeIt directory from the shell and type 'chmod +x bossref'

Usage
-----

PREPARE GENOMES

Genomes should be prepared as GenBank flat format files, with all contigs for a single genome in one file (if not complete). These files should have the extension '.gb' ('.gbk' and '.gbf' are also accepted). PhasomeIt will create temporary files within this folder and output files within a subfolder called 'summary_tracts'. The total size of these files is likely to be several times the size of the original genome so ensure that there is enough disk space available. Processing time will be faster if these files are stored locally on a fast disk.

INVOKE PHASOMEIT

PhasomeIt is invoked by calling python on the script from the command line, navigate to the folder containing phasomeit and invoke the command below. Depending on how you have set up your system you may need to use 'python3' instead of 'python':

python phasomeit -t *target-folder* -c *cutoffs* [-f *filter*]

Where: *target-folder* is the path containing your genomes
       *cutoffs* are the minimum repeat lengths to accept at each repeat length (see below)
       *filter* is used to remove some results and is optional (see below)

HOW TO CHOOSE CUTOFFS

Each cutoff is set for a different tract length. So, for example if you set it to '5 4 2' then it will accept GGGGG or TTTTTTTTT but not CCCC or AAA; ATATATAT but not GAGAGA; and ATCATC but not TTA. Repeats are included under the lower length that they are accepted on, so GAGAGAGA would be accepted as 4 GA repeats not 2 GAGA repeats in the preceeding example. However if cutoffs of 5 5 2 were set then it would be accepted as a 2 GAGA repeat. Any value of 0 means that length is ignored, so '5 0 2' would omit dinucleotide repeats from the preceeding example.

There are no exact guidelines for what cutoffs to choose, it will depend on your organism and you will need the domain knowledge to choose them correctly but if you set them too low then the program may not be able to cope with the number of false positives produced. For single tracts in particular values less than 7-10 may be problematic depending on genome size and G/C composition. Generally you want lower cutoffs for longer repeat units sizes. For our Campylobacter analysis values of '7 6 0 5 5' were used. This excluded trinucleotide repeats as these cannot cause phase variation by frame-shift mutation, however it is known that variation in these tracts can cause variation of function in some species.

FILTERING RESULTS

With species with a low (or high) G/C content, chance is likely to skew the results, especially as A/T tracts are less variable and need to be longer to produce a phase variable effect. Filtering allows extraneous results to be removed. Currently only a single filter can be set for any run. Filters are specified as {repeat-unit}{maximum-length-to-ignore} (no space between the two) and repeat-unit types support the normal four DNA letters plus W and S, so - for example - adding '-f W9' would remove all A/T tracts of length 9 or less from the analysis (i.e. treat the minimum A/T length as 10). A filter cutoff of 'W9' was used in our analysis.

OTHER FUNCTIONS

See 'python phasomeit --help' for a full list of commands.

# Support and citation

PhasomeIt is described in the 2018 publication _PhasomeIt: an ‘omics’ approach to cataloguing the potential breadth of phase variation in the genus Campylobacter_  https://doi.org/10.1099/mgen.0.000228

Please cite this publication if you use this program in your research.
