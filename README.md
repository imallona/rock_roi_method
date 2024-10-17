# Aim

This repository belongs to the rock and roi project from the University of Zurich. You might want to browse https://github.com/imallona/rock_roi_paper too.

This repository includes a Snakemake workflow to automate data processing from raw reads to count tables (and R `singleCellExperiment` objects) listing both on-target TSO and off-target WTA readouts. The software stack needed to run the method is containerized using Docker.

If you are planning to analyze BD Rhapsody data without RoCK nor RoCK+ROI, that is, just WTA plus sample tags, you might want to check https://github.com/imallona/rhapsodist (under development) instead.

# Preprint

- [RoCK and ROI: Single-cell transcriptomics with multiplexed enrichment of selected transcripts and region-specific sequencing](https://www.biorxiv.org/content/10.1101/2024.05.18.594120v2) (2024)
  
# Components

To analyze their data, users need to provide their sequencing files in fastq format (one for the cell barcode plus UMI; and another for the cDNA) and a configuration file specifying the experiment characteristics and extra information, including:

1. A genome (fasta) to align the genome to (i.e. hg38, mm10 etc). The genome needs to contain all (ontarget) captured sequences, so if these do not belong to the standard genome (i.e. GFP, tdTomato) it needs to be updated to append the extra sequences.
2. A gene annotation (GTF) whose features are quantified separately for WTA and TSO. It is expected to contain a whole transcriptome gene annotation (i.e. Gencode, Refseq etc) as well as an explicit definition of the rock and/or roi targets captured by the TSO. Instructions to build this GTF are included within the software’s documentation.
3. A set of cell barcode whitelists (standard BDRhapsody cell barcodes are included within this repository)
4. Parameters to fine tune CPU and memory usage.
5. Parameters to fine tune the expected number of cells and other EmptyDrops parameters.

# Workflow layout

The workflow follows these steps:

1. Index the reference genome with STAR. Please finetune the `sjdbOverhang` (`config.yaml`) accordingly (i.e. cDNA length - 1). 
2. Subset reads matching the WTA cell barcodes and map those to the transcriptome (genome plus GTF) using STARsolo. Detected cell barcodes (cells) are filtered in at two levels: first, by matching to the user-provided cell barcode whitelist; and second, by applying the EmptyDrops algorithm to discard empty droplets. We report two outputs from this step: the filtered-in cells according to the aforementioned filters; and the unbiased, whole-transcriptome WTA count table
3. Subset reads matching both the TSO CB structure and the filtered in cell barcodes and map those to the transcriptome. Our reasoning is that the expected TSO transcriptional complexity is undefined and not usable to tell apart cells from empty droplets, so we borrow the filtered-in cells from the EmptyDrops results from the WTA analysis.
4. (optional) Count on-target features in a more lenient way, filtering in multioverlapping and multimapping reads. This run mode is recommended when the captured regions target non unique loci (i.e. repetitive sequences).

Hence, our workflow always reports a WTA count table with as many genes as on-target and off-target gene features in the GTF, and per filtered-in cell barcode. As for the TSO, we offer these run modes:

- `tso off- and ontarget unique`: generates a count table for TSO reads from filtered-in cells; this count table has the same dimensions as the WTA.
- `tso ontarget multi`: creates a count table for TSO reads from filtered-in cells for only on-target features while allowing for multioverlapping and multimapping alignments.
- `all`: produces both `tso off- and ontarget unique` and `tso ontarget multi` outputs.

Finally, we generate an R SingleCellExperiment object with the aforementioned count tables and the following structure:
- `wta` assay: raw counts from the WTA analysis (unique reads).
- (optional) `tso_off_and_ontarget_unique` assay: raw counts from the `tso off- and ontarget unique` or `all` run modes.
- (optional) `tso_ontarget_multi` altExp alternative experiment: raw counts from the `tso ontarget multi` run mode.  A complementary altExp built on WTA data, named `wta_ontarget_multi`, quantifies multioverlapping and multimapping reads to the on-target regions in WTA data. 

We also provide a simulation runmode to showcase the method, where raw reads (fastqs), genome and GTF are generated for three on-target features and one off-target features across hundreds of cells before running the method.

# Repository structure

- `data`: Contains BD Rhapsody whitelists and a bash script to generate fake data for run the method (in bash).
- `Snakefile`, main workflow file.
- `Dockerfile`, docker file.
- `src`: scripts called by `Snakefile`.
- `env`: conda environment.

# Snakefile layout

## Without simulation

- Rulegraph: [pdf](./docs/rulegraph.pdf), [png](./docs/rulegraph.png).
- DAG for a run with four samples: [txt gz](./docs/dag.txt.gz).

As generated by:

```
snakemake --configfile ~/src/rock_roi_paper/00_mixing_experiment/mixing_conf.yaml --dag | \
   gzip -c > docs/dag.txt.gz
snakemake --configfile ~/src/rock_roi_paper/00_mixing_experiment/mixing_conf.yaml --rulegraph | \
   grep -v "^\[" | \
   dot -Tpdf > docs/rulegraph.pdf
snakemake --configfile ~/src/rock_roi_paper/00_mixing_experiment/mixing_conf.yaml --rulegraph | \
   grep -v "^\[" | \
   dot -Tpng > docs/rulegraph.png
```


## With simulation

- Rulegraph: [pdf](./docs/simul.pdf), [png](./docs/simul.png)
- DAG: [txt gz](.docs/dag_simulation.gz) (gzip compressed- hundreds of cells)

As generated by:

```
snakemake --configfile config.yaml --dag | gzip -c > docs/dag_simulation.gz

snakemake --configfile config.yaml --rulegraph | \
   grep -v "^\[" | \
   dot -Tpdf > docs/simul.pdf
snakemake --configfile config.yaml --rulegraph | \
   grep -v "^\[" | \
   dot -Tpng > docs/simul.png
```

# Installs and/or running the method

We asume a GNU/Linux system. We haven't tested it on Mac and it won't run on Windows.

We provide a Dockerfile to define a suitable GNU/Linux environment to run our method.

## Using docker

Assuming a single-threaded execution (--cores 1):

```
docker build . -t rock && docker run -it --entrypoint /bin/bash rock

## once inside the `rock` container

snakemake -p --cores 1 --configfile config.yaml

```

Data access to the host (i.e. mounting) is needed.

## Using conda

In a GNU/Linux install conda, i.e. miniconda, with something along these lines:

```
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3

```

Activate conda:

```
source ~/miniconda3/bin/activate
```

Create an env containing snakemake and mamba (to speed up installs):

```
conda install -c conda-forge mamba
mamba create -c conda-forge -c bioconda -n rock snakemake

```

Activate the env with snakemake and run the method (in this case, a simulation).

```
conda activate rock
snakemake -s Snakefile --configfile config.yaml --use-conda --cores 2
```

## Compiling and installing dependencies directly

### STAR (STARsolo)

```
mkdir -p ~/soft/star
cd $_

wget https://github.com/alexdobin/STAR/archive/2.7.10b.tar.gz
tar -xzf 2.7.10b.tar.gz
cd STAR-2.7.10b

cd source
make
make install

export PATH="$HOME"/soft/star/STAR-2.7.10b/source:$PATH

```

### Subread

```
mkdir -p ~/soft/subread
cd $_
wget 'https://sourceforge.net/projects/subread/files/subread-2.0.6/subread-2.0.6-source.tar.gz/download' 
tar xzvf download 
cd subread-2.0.6-source/src 
make -f Makefile.Linux 

## assuming ~/.local/bin/ is part of the $PATH
ln -s  ~/soft/subread/subread-2.0.6-source/bin/featureCounts ~/.local/bin/featureCounts
```

### Python deps

Including snakemake, pandas and deeptools.

Caution this will use the system's pip and current pythonpath; using a virtualenv is advised.

```
pip install snakemake pandas

```

### BEDtools

```
mkdir -p ~/soft/bedtools
cd $_
wget https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz
tar -zxvf bedtools-2.29.1.tar.gz
cd bedtools2
make

## assuming ~/.local/bin/ is part of the $PATH
cp bin/* ~/.local/bin/
```

### Kent utils

This recipe is for Linux x86 users. For other OSes and/or architectures, select the relevant binaries from https://hgdownload.soe.ucsc.edu/admin/exe/ .

```
cd ~/.local/bin # or somewhere in your $PATH where you can write

wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v385/faSize
chmod +x faSize

wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v385/bedGraphToBigWig
chmod +x bedGraphToBigWig
```

# Running the method (without conda)

```
snakemake --cores 50 --configfile config.yaml -p
```

## Usage tips

### Writing a config file for snakemake

Please check [our example config yaml](config.yaml).

### Simulation

We provide a data simulation procedure to be able to run the workflow on them. In brief, we aim to quantify in 2 cells a set of features: one is `offtarget` and captured by the WTA assay; and others are (multioverlapping, multimapper) `ontarget` readouts captured by the TSO assay.

### Picking up the right whitelist

- 96x3: Enhanced beads, 2022 and 2023
- 384x3: Latest enhanced beads batches, 2023

The 96x3 list is a subset of the 384x3 list. You can browse these allowedlists at [https://github.com/imallona/rock_roi_method/tree/main/data](https://github.com/imallona/rock_roi_method/tree/main/data).

### TSO counting only target features/regions

To subset the GTF features to count during the custom TSO counting process:

1. Write the starting GTF (`config['gtf']`) so it contains an indentifier in column2, i.e. you could flag the features that are being 'rock'-ed or 'roi'-ed with a searchable `source` [GTF column 2](https://www.ensembl.org/info/website/upload/gff.html). 

For instance, with `captured` in column2.

```
offtarget       ERCC    exon    1       1061    .       +       .       gene_id "offtarget"; transcript_id "offtarget_1";
ontarget        captured        exon    1       1023    .       +       .       gene_id "ontarget_1"; transcript_id "ontarget_1";
ontarget        captured        exon    100     800     .       +       .       gene_id "ontarget_1b"; transcript_id "ontarget_1b";
ontarget        captured        exon    1030    1090    .       +       .       gene_id "ontarget_2"; transcript_id "ontarget_2";
```

Then, edit the `config.yaml` file so the `capture_gtf_column_2_pattern` reads `captured`.

# License

GPLv3

# Extra reading

- https://github.com/imallona/rock_roi_paper Actual analysis using this method
- (deprecated) https://gitlab.uzh.ch/izaskun.mallona/ebrunner_spectral

## Contributors

- Izaskun Mallona 
- We reuse and adapt tools from STAR, subread (featurecounts), samtools and others; we are extremely grateful to their contributors for their unvaluable resources free and openly provided to the community

## Contact

izaskun.mallona at mls.uzh.ch, Mark D. Robinson lab
https://www.mls.uzh.ch/en/research/robinson.html
