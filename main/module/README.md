### About
This repository encompasses the source code designed for the purpose of appending RG (Read Group) tags to BAM (Binary Alignment Map) files, as well as the subsequent recovery of CB (Cell Barcode) tags. The specific operations conducted within this codebase include:

1. Inclusion of RG tags to individual reads: The RG tag is composed as the concatenation of the 'bam-name' and the 'cb-tag'.
2. Creation of a sequence logo representing the TSO tail: This entails a visual representation of sequence motifs.
3. Derivation of canonical sequence from the generated logo: Canonical sequences are identified as the primary representatives of specific motifs.
4. Matching of dervied canonical sequence with the tail sequences of reads lacking CB tags and those that deviate from the provided linker patterns.
5. Determination of CB tags using positional information, for alignments matched with the canonical sequence from step 4.
6. Validation of the obtained CB codes against a set of permissible CB codes, extracted from the input BAM file.
7. Merging of all input BAM files into a singular BAM file.
8. Appending RG tags to the BAM header, facilitating enhanced organization and categorization.


### Setup and execution on Terminal
1. Navigate to directory combined_rock_roi using command:
    $ cd combined_rock_roi

2. Virtual environment establishment and activation:
    a. Create conda environment named 'rock_n_roi' using rock_roi_conda_env.yml file:
        $ conda env create -f rock_n_roi_conda_env.yml
    b. Activate this environment using the following command:
        $ conda activate rock_n_roi

    Note: Name of the environment can be modified and tailored as per requirement in rock_roi_conda_env.yml
    
2. Navigate to directory combined_rock_roi/src using command:
    combined_rock_roi$ cd src

3. Configuration of input parameters:
    Customize input parameters such as path of files, chromosome names, logs folder names etc. within src/utilities/config_input_parameters.yml file.

4. Execution of the entire pipeline: 
    Invoke src/run.sh file using command:
        src$ ./run.sh

5. Output storage and location: 
    The generated output is stored either within the 'logs' folder or within the designated path for 'log_folder' as stipulated in the 'config_input_parameters.yml' file.


