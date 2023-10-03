## Repository structure

- `data`: Contains the scripts to generate fake data for testing the pipeline
- `main`: Contains the scripts and python module to align and count both targeted and untargeted gene expression.

## Installs / requirements

We asume a GNU/Linux system. We haven't tested it on Mac and it won't run on Windows.

We provide a Dockerfile to define a suitable GNU/Linux environment to run our method.

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
to-do
```

## License

GPLv3

## Contributors

- Izaskun Mallona 
- Nidhi Agrawal
- We reuse and adapt tools from STAR, subread (featurecounts), samtools and others; we are extremely grateful to their contributors for their unvaluable resources free and openly provided to the community
