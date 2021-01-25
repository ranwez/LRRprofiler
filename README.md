# LRRprofiler

LRRprofiler is an automated program to detect and annotate LRR (Leucine-Rich Repeat) containing plant receptor.

Please cite : ...

The program was developt for Linux os.

The program uses several external tools and database:
- iTAK 
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

### Running LRRprofiler
The program can be run with the command line :
singularity run lrr_profiler0.1.sif --in_proteome <fastafile> --name <jobname>

The program will work if all sequence headers are parsed without description (i.e. ">OS01g10200" and not ">OS01g10200 expressed protein")

Using the example file :
singularity run lrr_profiler0.1.sif --in_proteome Arabidopsis_Thaliana_reviewed_proteom_SwissProt_05-2020.fasta --name ARATH