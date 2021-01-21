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

A singularity container for the LRRprofiler program can be download from the Sylabs cloud (link). 
The file LRRprofiler_v0.1_sing_3.3.def provide the corresponding singularity recipe.

### Running LRRprofiler
The program can be run with the command line :
singularity run LRRprofiler _v0.1_sing_3.3.sif --in_proteome <fastafile> --name <jobname>

Using the example file :
singularity run LRRprofiler _v0.1_sing_3.3.sif --in_proteome Arabidopsis_Thaliana_reviewed_proteom_SwissProt_05-2020.fasta --name ARATH