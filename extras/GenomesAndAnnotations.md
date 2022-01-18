# Genome & Annotation

## pre-prep
First you want to setup a `genomes/` directory 
(sister to your project directories).

For each species you will then want to setup
a subdirectory in the following style

```
genomes/
└── Mesculenta
    ├── Mesculenta.fa
    └── Mesculenta.gff3
```

Where the species name (in this example is `Mesculenta`) of
course matches the actual species in question.

Symlinks to genome & annotation files as downloaded are 
both allowed and encouraged.

## processing
For both Hisat2 and CollectRnaSeqMetrics you will
need not the raw files, but an index, and the annotation
specifically processed (to refflat3), so the final
structure would look like 

```
genomes/
└── Mesculenta
    ├── Mesculenta.1.ht2
    ├── Mesculenta.2.ht2
    ├── Mesculenta.3.ht2
    ├── Mesculenta.4.ht2
    ├── Mesculenta.5.ht2
    ├── Mesculenta.6.ht2
    ├── Mesculenta.7.ht2
    ├── Mesculenta.8.ht2
    ├── Mesculenta.fa
    ├── Mesculenta.gff3
    ├── Mesculenta.gtf
    ├── Mesculenta.refflat2
    └── Mesculenta.refflat3.gz
```

At some point this will be automated and added 
to RNAsleek. Until then, please see the example
qsub files that would take care of this `template_anno_build.qsub`
and `template_hisat_build.qsub`. You just need 
to replace (by the metold of your choice)
"SPECIES" with "<your_species_name>" and submit, e.g. 
for the example above you might run 

```
# create species-specific 
cat template_anno_build.qsub | sed 's/SPECIES/Mesculenta/g' > anno_build_Mesculenta.qsub
# and from the parent directory of `genomes` and any <project directories>
# submit with
qsub anno_build_Mesculenta.qsub
# and check output
tree genomes/Mesculenta
```
