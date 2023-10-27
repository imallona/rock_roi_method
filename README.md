# Aim

This repository belongs to the rock and roi project from the University of Zurich. You might want to browse https://github.com/imallona/rock_roi_paper too.

This repository includes a Snakemake workflow to automate data processing from raw reads to count tables (and R `singleCellExperiment` objects) listing both on-target TSO and off-target WTA readouts. The software stack needed to run the method is containerized using Docker.

To analyze their data, users need to provide their sequencing files in fastq format (one for the cell barcode plus UMI; and another for the cDNA) and a configuration file specifying the experiment characteristics and extra information, including:

1. A genome (fasta) to align the genome to (i.e. hg38, mm10 etc). The genome needs to contain all (ontarget) captured sequences, so if these do not belong to the standard genome (i.e. GFP, tdTomato) it needs to be updated to append the extra sequences.
2. A gene annotation (GTF) whose features are quantified separately for WTA and TSO. It is expected to contain a whole transcriptome gene annotation (i.e. Gencode, Refseq etc) as well as an explicit definition of the rock and/or roi targets captured by the TSO. Instructions to build this GTF are included within the software’s documentation.
3. A set of cell barcode whitelists following BDRhapsody’s standards (standard BDRhapsody cell barcodes are included within the software)
4. Parameters to fine tune CPU and memory usage.
5. Parameters to fine tune the expected number of cells and other EmptyDrops parameters.

The workflow follows these steps:

1. Index the reference genome with STAR.
2. Subset reads matching the WTA cell barcodes and map those to the transcriptome (genome plus GTF) using STARsolo. Detected cell barcodes (cells) are filtered in at two levels: first, by matching to the user-provided cell barcode whitelist; and second, by applying the EmptyDrops algorithm to discard empty droplets. We report two outputs from this step: the filtered-in cells according to the aforementioned filters; and the unbiased, whole-transcriptome WTA count table
3. Subset reads matching both the TSO CB structure and the filtered in cell barcodes and map those to the transcriptome. Our reasoning is that the expected TSO transcriptional complexity is undefined and not usable to tell apart cells from empty droplets, so we borrow the filtered-in cells from the EmptyDrops results from the WTA analysis.
4. (optional) Count on-target features in a more lenient way, filtering in multioverlapping and multimapping reads. This run mode is recommended when the captured regions target non unique loci (i.e. repetitive sequences).

Hence, our workflow always reports a WTA count table with as many genes as on-target and off-target gene features in the GTF, and per filtered-in cell barcode. As for the TSO, we offer these run modes:

- `tso off- and ontarget unique`: generates a count table for TSO reads from filtered-in cells; this count table has the same dimensions as the WTA.
- `tso ontarget multi`: creates a count table for TSO reads from filtered-in cells for only on-target features while allowing for multioverlapping and multimapping alignments.
- `all`: produces both `tso off- and ontarget unique` and `tso ontarget multi` outputs.

Finally, we generate an R SingleCellExperiment object with the aforementioned count tables and the following structure:
- `wta` assay: raw counts from the WTA analysis.
- (optional) `tso_off_and_ontarget_unique` assay: raw counts from the `tso off- and ontarget unique` or `all` run modes.
- (optional) `tso_ontarget_multi` altExp alternative experiment: raw counts from the `tso ontarget multi` run mode.

We also provide a simulation runmode to showcase the method, where raw reads (fastqs), genome and GTF are generated for three on-target features and one off-target features across hundreds of cells before running the method.

## Repository structure

- `data`: Contains BD Rhapsody whitelists and a bash script to generate fake data for run the method (in bash).
- `main`: Contains the Snakemafile and python module to align and count both targeted and untargeted gene expression.

## Snakefile layout

### Without simulation

```
## i.e. 
snakemake -s main/Snakefile --configfile ~/src/rock_roi_paper/00_mixing_experiment/mixing_conf.yaml  --dag
```

- Rulegraph: [png](./docs/rulegraph.png), [pdf](./docs/rulegraph.pdf), [specs](./docs/rulegraph).
- DAG for a run with four samples: [png](./docs/dag.png), [pdf](./docs/rulegraph.pdf), [specs](./docs/rulegraph).

### With simulation


- Rulegraph: [png](./docs/simul.png), [specs](./docs/rulegraph_simulation)
- DAG: [specs](.docs/dag_simulation.gz) (gzip compressed- hundreds of cells)

## Installs and/or running the method

We asume a GNU/Linux system. We haven't tested it on Mac and it won't run on Windows.

We provide a Dockerfile to define a suitable GNU/Linux environment to run our method.

### Using docker

Assuming a single-threaded execution:

```
docker build . -t rock && docker run -it --entrypoint /bin/bash rock

## once inside the `rock` container

cd main
snakemake -p --cores 1 --configfile config.yaml

```

### Installing and user the system's dependencies

#### STAR (STARsolo)

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

#### Subread

```
mkdir -p ~/soft/soft/subread
cd $_
wget 'https://sourceforge.net/projects/subread/files/subread-2.0.6/subread-2.0.6-source.tar.gz/download' 
tar xzvf download 
cd subread-2.0.6-source/src 
make -f Makefile.Linux 

## assuming ~/.local/bin/ is part of the $PATH
ln -s  ~/soft/soft/subread/subread-2.0.6-source/bin/featureCounts ~/.local/bin/featureCounts
```

#### Python dependencies

Including snakemake, pandas and deeptools.

Caution this will use the system's pip and current pythonpath; using a virtualenv is advised.

```
cd /home/rock/main/module && \
    pip install -r rock_n_roi_requirements.txt && \
    pip install snakemake pandas deeptools

```

#### Running the method

```
cd main
snakemake --cores 50 --configfile config.yaml

```

## Usage tips

### Writing a config file for snakemake

Please check [our example config yaml](./main/config.yaml). **Important** This file must be named `config.yaml`.

### Simulation

We provide a data simulation procedure to be able to run the workflow on them. In brief, we aim to quantify in 2 cells a set of features: one is `offtarget` and captured by the WTA assay; and others are (multioverlapping, multimapper) `ontarget` readouts captured by the TSO assay.

The rationale behind is: lorem ipsum.

### Picking up the right whitelist

- 96x3: Enhanced beads, 2022 and 2023
- 384x3: Latest enhanced beads batches, 2023

### TSO counting only target features/regions

To subset the GTF features to count during the custom TSO counting process:

1. Write the starting GTF (`config['gtf']`) so it contains an indentifier in column2, i.e. you could flag the features that are being 'rock'-ed or 'roi'-ed with a searchable `source` (https://www.ensembl.org/info/website/upload/gff.html)[GTF column 2]. 

For instance, with `captured` in column2.

```
offtarget       ERCC    exon    1       1061    .       +       .       gene_id "offtarget"; transcript_id "offtarget_1";
ontarget        captured        exon    1       1023    .       +       .       gene_id "ontarget_1"; transcript_id "ontarget_1";
ontarget        captured        exon    100     800     .       +       .       gene_id "ontarget_1b"; transcript_id "ontarget_1b";
ontarget        captured        exon    1030    1090    .       +       .       gene_id "ontarget_2"; transcript_id "ontarget_2";
```

Then, edit the `config.yaml` file so the `capture_gtf_column_2_pattern` reads `captured`.

## License

GPLv3

## Contributors

- Izaskun Mallona 
- Nidhi Agrawal
- We reuse and adapt tools from STAR, subread (featurecounts), samtools and others; we are extremely grateful to their contributors for their unvaluable resources free and openly provided to the community

## Contact

izaskun.mallona at gmail dot com
