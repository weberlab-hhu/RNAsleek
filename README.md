# RNAsleek
Semi-automated RNAseq pipeline for processing public RNAseq samples on a PBSpro cluster.

SleekRNAseq aims to get you both mapped reads and some nice QC for a whole bunch of
RNAseq samples with some safety checks built in to make sure everything is 
working as expected.

**fair warning, this is more documentation of what we've done than anything
designed to be used by others**

## Prep work
To get started you need a few things

### All the tools
#### python + packages
It's only been tested under python3.6, hopefully works forward at least

You will also need the packages listed in requirements.txt, e.g.
```
pip install -r requirements.txt
```

#### Bioinformatics tools
Each step has certain dependencies. The versions 
listed in parentheses
indicate what we tested it with/used.

##### Wget or Fetch
In practice Wget was much less likely to throw an error for
the download, but both of the above still need SRA Toolkit
for `fastq-dump` (2.8.2).

##### Trimmomatic
Current non-dynamic code expects the jar
and adapters to be found within `$HOME/extra_programs/Trimmomatic-0.36/`.
Configurability is also on the todo list (0.36).

##### Fastqc
FastQC (v0.11.5)

##### Hisat
samtools, assumed to be available with `module load` (1.6)
hisat2 (2.1.0)

##### CollectRNAseqMetrics
Picard Tools (52.0)
assumed to be available at `$HOME/extra_programs/picard.jar`

### RunInfo file
That is a SraRunInfo.csv (or file with the exact same columns) 
for all the samples you wish to process.
You generally get one of these by searching for whatever you're interested in 
on SRA (https://www.ncbi.nlm.nih.gov/sra), then in the upper right you 
- click on "Send to:"
- select "File" under Choose Destination
- change Format to "RunInfo"
- click "Create File"

### config file
This specifies species info, what jobs you wish to run and any customization.
In particular, the scientific name, taxid and for mapping the species (sp)
must be specified.

See example.ini

### Prepped genome information
For now this is awkward, inflexilble, and manual... 
cleaning it up is on the todo list...

You only need this for steps that require a reference genome 
(currently Hisat and CollectRNAseqMetrics), and
this needs to be setup independent of this code base.

In a folder in the same directory as you will be running these
analyses, you will need a folder named 'genomes' which should
contain a folder with the species name specified in the config
under 'sp' (e.g. `sp = example_species`). Continuing with this
example you will want the genomic fasta file to be found
under `genomes/example_species/example_species.fa`, the gff3 annotation
file to be found under `genomes/example_species/example_species.gff3`
and the hisat2 indexes (for Hisat only) to be found under 
`genomes/example_species/example_species.\*.ht2`

## Setup
To setup all the scripts and qsub files you will just need to run
```
python <path_to>/RNAsleek/rnasleek.py <project_directory> <RunInfo_file> -c <config.ini>
```

## Step wise
Once you have the steps you can cd into your chosen 'project_directory'
and qsub the files once their dependencies are met. 
The qsub script will setup a job array with all the
samples.

The Job dependencies are as follows
- Wget or Fetch: None
- Trimming: Wget or Fetch
- Fastqc: Trimming
- Hisat: Trimming
- CollectRNAseqMetrics: Hisat

After each job finishes you should run 
```
python <path_to>/RNAsleek/rnasleek.py <project_directory> <RunInfo_file> -c <config.ini> --check_output
```
This will produce the file `project_directory/output_report.txt`. This files shows
any and all errors found in the stderr output files as well as any deviations from 
expected output for each step. You'll probably want to use `grep` to look at just
the step you just ran. Errors have to be fixed manually. Often it just requires increasing
the requested memory in the qsub file; except when it doesn't, which is why it's hard to code.

Once you are happy with the output of a step, qsub the next one until you're done

## That's a wrap
If you've ran all the steps from example.ini, you can also get a 
nice output summary via multiQC and some plotting here.

```
python <path_to>/RNAsleek/rnasleek.py <project_directory> <RunInfo_file> -c <config.ini> --prep_multiqc
cd <project_directory>/multiqc/
multiqc .
cd ../multiqc_untrimmed
multiqc .
cd ../..
python <path_to>/RNAsleek/viz/summarizer.py <project_directory> <RunInfo_file> -o <output.pdf>
```
