# RNAsleek
Semi-automated RNAseq pipeline for processing public RNAseq samples on a cluster with PBSpro. 
Actually... it's really not configurable enough to be used for anything but the HHU HPC cluster, just right now.

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
pip install .  # for contributers/if code might change: pip install -e .
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

##### BWA
bwa 2.1 (`module load bwa-mem2/2.1`)

##### CollectRNAseqMetrics
Picard Tools (52.0), which is
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

## Directory Organization at start
I would recommend using:
(SraRunInfo.csv and config.ini don't have to be there, it just makes it easier
to keep track of)

```
<my directory>/SraRunInfo.csv
<my directory>/config.ini
genomes  # this will be autocreated
```

### Prepped genome information

If you haven't setup a genome with indexes for this species yet,
then you should run:

```
rnasleek -d <project_directory> -s <RunInfo_file> -c <config.ini> genome \
  -f <genome.fa> -g <genome.gff3> -s <SpeciesName>
```

This will copy the files `<genome.fa>` and `<genome.gff3>` into the 
directory `genomes/<SpeciesName>` where `genomes` is located sister 
to the project directory.

`<SpeciesName>` must exactly match that listed under `sp = ` in mapping and
other genome-requiring jobs in the config file.

To run the indexing, from the project directory run `qsub qsubs/genome_prep.qsub`

## Setup

To setup all the scripts and qsub files you will just need to run
```
rnasleek -d <project_directory> -s <RunInfo_file> -c <config.ini> setup
```


## Step wise
Once you have the steps you can `cd` into your chosen _project_directory_
and qsub the files once their dependencies are met (e.g. `qsub qsubs/wget.qsub`)
The qsub script will setup a job array with all the
samples.

The Job dependencies are as follows
- Wget or Fetch: None
- Trimming: Wget or Fetch
- Fastqc: Trimming
- Hisat: Trimming
- BWA: Trimming
- CollectRNAseqMetrics: Hisat

After each job finishes you should run 
```
rnasleek -d <project_directory> -s <RunInfo_file> -c <config.ini> check
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
rnasleek -d <project_directory> -s <RunInfo_file> -c <config.ini> multiqc
cd <project_directory>/multiqc/
multiqc .
cd ../multiqc_untrimmed
multiqc .
cd ../..
python <path_to>/RNAsleek/viz/summarizer.py <project_directory> <RunInfo_file> -o <output.pdf>
```

## Work arounds for the lack of internet access (e.g. on our HPC)
- run the first setup on a machine _with_ internet access
- run the `wget` on the same machine.
```
# needs slight customization because each script calls the variable `$PBS_O_WORKDIR`
# so make this variable (obviously, this assumes your wd is the project directory,
# and if not, `pwd` can be replaced with the full path to the project directory)
export PBS_O_WORKDIR=`pwd`
# run the download scripts
ls scripts/wgetSRS*|xargs -n1 -P2 -I% bash %
# copy the result to the hpc (or modify for your internet-free machine)
cd ..
rsync -ravu <project_dir> <user_name>@storage.hpc.rz.uni-duesseldorf.de:/gpfs/project/<user_name>/
```

## Thanks
Thank you to @danidey for some code and a whole lot of inspiration and organization
