## Repository structure

- `data`: Contains the scripts to generate fake data for testing the pipeline
- `main`: Contains the scripts and python module to align and count both targeted and untargeted gene expression.

## Installs / requirements

We asume a GNU/Linux system. We haven't tested it on Mac and it won't run on Windows.

We provide a Dockerfile to define a suitable GNU/Linux environment to run our method.

### Using docker

Within the root of this directory (so the `Dockerfile` is there):

```
docker build . -t rock && docker run -it --entrypoint /bin/bash rock
```

### Manually

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

## Usage example

We provide (within the `Dockerfile`) a data simulation plus analysis workflow. In brief, we aim to quantify in 2 cells a set of features: one is `offtarget` and captured by the WTA assay; and others are (multioverlapping, multimapper) `ontarget` readouts captured by the TSO assay.

The rationale behind is: lorem ipsum.

## License

GPLv3

## Contributors

- Izaskun Mallona 
- Nidhi Agrawal
- We reuse and adapt tools from STAR, subread (featurecounts), samtools and others; we are extremely grateful to their contributors for their unvaluable resources free and openly provided to the community
