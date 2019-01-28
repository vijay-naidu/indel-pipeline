"# indel-pipeline" 

Start with fasta files as downloaded from ENA website. These are assumed to be paired-end reads. The examples included here can be found using study accession number PRJEB4239. 

A description of the process inclusing software dependencies is contained in the file PilonPipeline.txt. 

First trim the input files:
scripts/indeltrim.sh
Switch to trim's output sub-directory '/trim'.

Then run the main script:
scripts/indelscript.sh
This will output vcf files for both snps and indels.

There is a python 3 script, scripts/process_indels_only.py, which runs through the vcf files and produsces summary information including a matrix of indel difference counts and locations of all the indels. In the script there is a depth threshold which by default is set to 10.
Another script, scripts/process_snps.py, outputs a SNP difference matrix.

Directory paths will need to be altered as appropriate, depending on where you run this from.
