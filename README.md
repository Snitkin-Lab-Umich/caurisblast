# caurisblast - python scripts to make custom BLAST easier

## Really Quick Start

This script currently runs best on the Great Lakes HPC cluster. Log into the cluster and run the following lines to set everything up:

```

git clone https://github.com/Snitkin-Lab-Umich/caurisblast
module load Bioinformatics blast-plus/2.16.0

```

Next, run your desired BLAST:

```

python blast.py --query [your_query] --subject [your_subject] --type [search_type] --name [run_name]

```

[your_query] should be your query sequences. These are what you are searching for. This should be a .fasta file or a directory of containing .fasta files.
[your_subject] should be your subject sequences. These are what you are searching through for hits. This should also be a .fasta file or a directory of containing .fasta files.
[search_type] should indicate if you are doing a nucleotide or protein BLAST.
[run_name] should be a unique name for your run.

Running this should place the outcome of your BLAST into the results/ directory.

## Detailed Explanation

This python script is a simple wrapper for BLAST+, the command line utility for running standard BLASTs with custom databases. It aims to slightly simplify using BLAST+ by automatically combining .fasta files, generating custom databases as needed, performing repeat masking, and comparing annotations to BLAST results. An explanation of each input is below.

### --query

These are the query sequences you want to search for in your BLAST database. Currently, these must be provided in the .fasta format. You can provide a single input file or a directory containing multiple .fasta files. If multiple files are provided, they will be concatenated into single input file for the BLAST search. This will be located at [run_name]_combined_query.fasta in the main directory.

### --subject

This is either your existing BLAST database or a set of .fasta files you want to make into a BLAST database. If you are using an existing BLAST database, provide the path to the .fasta file here and add the --database flag. If you want to make a new BLAST database, provide a single input file or a directory containing multiple files. As with --query, these currently must be in the .fasta format. If multiple files are provided, they will be concatenated into db/[run_name]_blastdb.fasta. 

### --database

Add this flag if you are using an existing BLAST database. If you have previously run this script, the database will be at db/[run_name]_blastdb.fasta.

### --type

This is the type of BLAST you are doing, as well as the type of BLAST database you are making. The options are 'nucleotide' or 'protein'. The type of database must match the search you are performing.

### --name

This is the run name, which will prefix most output files.

### --format

(Optional) This is the output format. By default, this is format 10, which is a .csv file. See https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.options_common_to_all_blast/ for all options.



Additional options are currently in progress!


