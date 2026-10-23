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

## Using Singularity/Apptainer (.sif) containers

LRRprofiler Singularity/Apptainer containers are available on the [GitHub Container Registry (GHCR)](https://github.com/JohGi?tab=packages). They can be downloaded directly using Apptainer or Singularity with the `oras://` protocol.

**Latest version (0.24):**

```bash
apptainer pull LRRprofiler_0.24.sif oras://ghcr.io/johgi/lrrprofiler:0.24
```

### Usage

LRRprofiler can be run using Apptainer (or Singularity):

```bash
apptainer run LRRprofiler_0.24.sif \
    --in_proteome <fastafile> \
    --name <jobname> \
    [--dev] [--nobuild]
```

For example, using the provided Arabidopsis thaliana proteome:

```bash
apptainer run LRRprofiler_0.24.sif \
    --in_proteome /opt/LRRprofiler/example/Arabidopsis_Thaliana_reviewed_proteom_SwissProt_05-2020.fasta \
    --name ARATH
```

#### Pipeline options
--in_proteome: (mandatory) Path of the proteome fasta file. Currently, the program will work if all sequence headers are parsed without description (i.e. ">OS01g10200" and not ">OS01g10200 expressed protein")

--name: (mandatory) Character string used for output directory and file names.

--dev: (optional) If provided, the pipeline will retain the working directory containing temporary files.

--nobuild: (optional) If provided, skip the profile refinement process and use exclusively profiles from the HMM_lib folder. This option can be useful to obtain homogeneous and reproductible annotation when comparing several proteomes.
