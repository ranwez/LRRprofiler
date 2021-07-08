# LRRprofiler

LRRprofiler is an automated program to detect and annotate LRR (Leucine-Rich Repeat) containing plant receptors.

Please cite : 
 A New Comprehensive Annotation of Leucine-Rich Repeat-Containing Receptors in Rice. Gottin C., Diévart A., Summo M., Droc G., Périn C., Ranwez V. and Chantret N.
 preprint on bioRxiv. doi: https://doi.org/10.1101/2021.01.29.428842

The program was developed on Linux os.

The program uses several external tools and databases:
- iTAK (v1.7)
- HMMER (v3.1b2)
- MAFFT (v7.271)
- TMHMM (v2.0c) [academic use only]
- SMART database

## Using Singularity (.sif) container

A singularity container for the LRRprofiler program can be download directly in your workspace with:
```
singularity pull LRRprofiler.sif library://cgottin/default/lrr_profiler:0.1
```
or from the [Sylabs cloud](https://cloud.sylabs.io/library/_container/600ea381517f0358917abf0a) .

The file LRRprofiler_v0.1_sing_3.3.def provide the corresponding singularity recipe.

### Usage
The program can be run with the command line :
```
singularity run lrr_profiler0.1.sif --in_proteome <fastafile> --name <jobname> [--dev]
```

--in_proteome: Path of the proteome fasta file. The program will work if all sequence headers are parsed without description (i.e. ">OS01g10200" and not ">OS01g10200 expressed protein")

--name: Character string use for output directory and file names

--dev: If provided, the pipeline will retain the working directory containing temporary files

Using the example file :
```
singularity run lrr_profiler0.1.sif --in_proteome Arabidopsis_Thaliana_reviewed_proteom_SwissProt_05-2020.fasta --name ARATH
```


