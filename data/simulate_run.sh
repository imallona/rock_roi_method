#!/bin/bash
##
## Simulates some WTA and TSO reads showcasing the rock/roi workflow
##

mkdir -p ~/simulate
cd $_

NTHREADS=5

## these are whitelists (BDRhapsody)
wget https://teichlab.github.io/scg_lib_structs/data/BD_CLS1.txt
wget https://teichlab.github.io/scg_lib_structs/data/BD_CLS2.txt
wget https://teichlab.github.io/scg_lib_structs/data/BD_CLS3.txt

## fake, short genome
cat << EOF > genome.fa
>offtarget ERCC-00002
TCCAGATTACTTCCATTTCCGCCCAAGCTGCTCACAGTATACGGGCGTCG
GCATCCAGACCGTCGGCTGATCGTGGTTTTACTAGGCTAGACTAGCGTAC
GAGCACTATGGTCAGTAATTCCTGGAGGAATAGGTACCAAGAAAAAAACG
AACCTTTGGGTTCCAGAGCTGTACGGTCGCACTGAACTCGGATAGGTCTC
AGAAAAACGAAATATAGGCTTACGGTAGGTCCGAATGGCACAAAGCTTGT
TCCGTTAGCTGGCATAAGATTCCATGCCTAGATGTGATACACGTTTCTGG
AAACTGCCTCGTCATGCGACTGTTCCCCGGGGTCAGGGCCGCTGGTATTT
GCTGTAAAGAGGGGCGTTGAGTCCGTCCGACTTCACTGCCCCCTTTCAGC
CTTTTGGGTCCTGTATCCCAATTCTCAGAGGTCCCGCCGTACGCTGAGGA
CCACCTGAAACGGGCATCGTCGCTCTTCGTTGTTCGTCGACTTCTAGTGT
GGAGACGAATTGCCAGAATTATTAACTGCGCAGTTAGGGCAGCGTCTGAG
GAAGTTTGCTGCGGTTTCGCCTTGACCGCGGGAAGGAGACATAACGATAG
CGACTCTGTCTCAGGGGATCTGCATATGTTTGCAGCATACTTTAGGTGGG
CCTTGGCTTCCTTCCGCAGTCAAAACCGCGCAATTATCCCCGTCCTGATT
TACTGGACTCGCAACGTGGGTCCATCAGTTGTCCGTATACCAAGACGTCT
AAGGGCGGTGTACACCCTTTTGAGCAATGATTGCACAACCTGCGATCACC
TTATACAGAATTATCAATCAAGCTCCCCGAGGAGCGGACTTGTAAGGACC
GCCGCTTTCGCTCGGGTCTGCGGGTTATAGCTTTTCAGTCTCGACGGGCT
AGCACACATCTGGTTGACTAGGCGCATAGTCGCCATTCACAGATTTGCTC
GGCAATCAGTACTGGTAGGCGTTAGACCCCGTGACTCGTGGCTGAACGGC
CGTACAACTCGACAGCCGGTGCTTGCGTTTTACCCTTAAAAAAAAAAAAA
AAAAAAAAAAA
>ontarget ERCC-00003 dimer
CAGCAGCGATTAAGGCAGAGGCGTTTGTATCTGCCATTATAAAGAAGTTT
CCTCCAGCAACTCCTTTCTTAATTCCAAACTTAGCTTCAGTTATAAATTC
CCCTCCCATGATTGGGATTTTATAAACTTTTCTTCCATATAATTCATCTT
TCTTCTCATAACCGTCTCCGAAAAACTTCAACTTAAATCCAACCTTTAAC
TGCTCATCAGCCATGTCTCCCACAGCATCAAAAATAGCAGTTGTTGGACA
TGTTAAGACACACTGCCCCAATCTCTCTAACATTTGATGCTCTAACTCTG
ACTTTTTAGGGTGGCATATCTGTATTATAAATCCTGGTCTTCCATCTGGT
GTTTTTGATGGAGGGACATATTTCTCAATTCCTGCTTCTGCTGGACACAT
TATAACTGAACAACCAAAACCTGTTGCCTCTGTAGCTGCAATCTTAGCCC
ACTTCTTTGTAGCTGCTGTTATTAAAACTCTTGAAACCCATATTGGGAAT
GCTTCTGCAAATGTATCTTCAATATATACTCCATTTATTTCCATAGTTTC
CCTCCATTAAGATTTTAACAATTATAGTTTATCTTAGGGGCTATTAATAT
CTTATCATTTGGTTTTTAATATTCGATAAATCCATAAATAAAAATATATC
AACAATAATTTTAAATAATCTAAGTATAGGTAATATAACAATTAAAAAGA
TTTAGAGGGATAGAATTGAACGGCATTAGGAGAATTGTTTTAGATATATT
GAAGCCGCATGAGCCAAAAATAACAGATATGGCATTAAAATTAACATCAT
TATCAAACATTGATGGGGTTAATATTACAGTCTATGAAATAGATAAAGAG
ACTGAGAATGTTAAAGTTACAATTGAAGGGAATAATTTAGATTTTGATGA
GATTCAGGAAATTATTGAAAGTTTGGGAGGGACTATTCACAGTATAGATG
AGGTTGTTGCAGGTAAAAAGATTATTGAAGAGTTAGAACACCACAAGATA
AAAAAAAAAAAAAAAAAAAAAAACAGCAGCGATTAAGGCAGAGGCGTTTGTATCTGCCATTATAAAGAAGTTT
CCTCCAGCAACTCCTTTCTTAATTCCAAACTTAGCTTCAGTTATAAATTC
CCCTCCCATGATTGGGATTTTATAAACTTTTCTTCCATATAATTCATCTT
TCTTCTCATAACCGTCTCCGAAAAACTTCAACTTAAATCCAACCTTTAAC
TGCTCATCAGCCATGTCTCCCACAGCATCAAAAATAGCAGTTGTTGGACA
TGTTAAGACACACTGCCCCAATCTCTCTAACATTTGATGCTCTAACTCTG
ACTTTTTAGGGTGGCATATCTGTATTATAAATCCTGGTCTTCCATCTGGT
GTTTTTGATGGAGGGACATATTTCTCAATTCCTGCTTCTGCTGGACACAT
TATAACTGAACAACCAAAACCTGTTGCCTCTGTAGCTGCAATCTTAGCCC
ACTTCTTTGTAGCTGCTGTTATTAAAACTCTTGAAACCCATATTGGGAAT
GCTTCTGCAAATGTATCTTCAATATATACTCCATTTATTTCCATAGTTTC
CCTCCATTAAGATTTTAACAATTATAGTTTATCTTAGGGGCTATTAATAT
CTTATCATTTGGTTTTTAATATTCGATAAATCCATAAATAAAAATATATC
AACAATAATTTTAAATAATCTAAGTATAGGTAATATAACAATTAAAAAGA
TTTAGAGGGATAGAATTGAACGGCATTAGGAGAATTGTTTTAGATATATT
GAAGCCGCATGAGCCAAAAATAACAGATATGGCATTAAAATTAACATCAT
TATCAAACATTGATGGGGTTAATATTACAGTCTATGAAATAGATAAAGAG
ACTGAGAATGTTAAAGTTACAATTGAAGGGAATAATTTAGATTTTGATGA
GATTCAGGAAATTATTGAAAGTTTGGGAGGGACTATTCACAGTATAGATG
AGGTTGTTGCAGGTAAAAAGATTATTGAAGAGTTAGAACACCACAAGATA
AAAAAAAAAAAAAAAAAAAAAAA
EOF

## fake, short canonical transcripts
echo -e 'offtarget\tERCC\texon\t1\t1061\t.\t+\t.\tgene_id "offtarget"; transcript_id "offtarget_1";
ontarget\tERCC\texon\t1\t1023\t.\t+\t.\tgene_id "ontarget_1"; transcript_id "ontarget_1";
ontarget\tERCC\texon\t100\t800\t.\t+\t.\tgene_id "ontarget_1b"; transcript_id "ontarget_1b";
ontarget\tERCC\texon\t1030\t1090\t.\t+\t.\tgene_id "ontarget_2"; transcript_id "ontarget_2";' > genome.gtf

## index
STAR --runThreadN $NTHREADS \
     --runMode genomeGenerate \
     --sjdbGTFfile genome.gtf \
     --genomeDir genome \
     --genomeSAindexNbases 4 \
     --sjdbOverhang 100 \
     --genomeFastaFiles genome.fa
     
## generate 100 R1s (cDNA) and R2s (CB + UMI) WTA for offtargets
## same UMI 100 times
seq 1 100 | awk ' {print "@"$0"\nTACTGGACTCGCAACGTGGGTCCATCAGTTGTCCGTATACCAAGACGTCTAAGGGCGGTGTACACCCTTTTGAGCAATGATTGCACAACCTGCGATCACCTTATACAGAATTATCAATCAAGCTCCCCGAGGAGCGGACTTGTAAGGACCGCCGCTTTCGCTCGGGTCTGCGG\n+\nFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF";}' | gzip -c > r1.fq.gz

seq 1 100 | awk ' {print "@"$0"\nCTTGTACTAGTGATGTTCTCCAGACAGGCTACAGATTTGATGGTTTTTTTTTTTTTTTTT\n+\nFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF";}' | gzip -c > r2.fq.gz

## similarly for region 1000-1100 for ontargets
## same UMI 100 times
seq 101 200 | awk ' {print "@"$0"\nATTATTGAAGAGTTAGAACACCACAAGATAAAAAAAAAAAAAAAAAAAAAAAACAGCAGCGATTAAGGCAGAGGCGTTTGTATCTGCCATTATAAAGAAGTTTCCTCCAGCAACTCCTTTCTTAATTCCAAACTTAGC\n+\nFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF";}' | gzip -c >> r1.fq.gz

seq 101 200 | awk ' {print "@"$0"\nCTTGTACTAAATGTGTTCTCCACCACGGCTACAGATTTGATGGTTTTTTTTTTTTTTTTT\n+\nFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF";}' | gzip -c >> r2.fq.gz

## another ontarget region - TSO barcode
## same UMI 100 times
seq 201 300 | awk ' {print "@"$0"\nATTCCAAACTTAGCTTCAGTTATAAATTCCCCTCCCATGATTGGGATTTTATAAACTTTTCTTCCATATAATTCATCTTTCTTCTCATAACCGTCTCCGAAAAACTTCAACTTAAATCCAACCTTTAACTGCTCATCA\n+\nFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF";}' | gzip -c >> r1.fq.gz

seq 201 300 | awk ' {print "@"$0"\nCTTGTACTAAATGTGTTCTCCACCACGGCTACAGATTTGATGGTTTTTTTTTTTTTTTTT\n+\nFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF";}' | gzip -c >> r2.fq.gz


## run starsolo for WTA
STAR --runThreadN $NTHREADS \
     --genomeDir genome \
     --readFilesCommand zcat \
     --outFileNamePrefix mapping_wta/ \
     --readFilesIn  r1.fq.gz r2.fq.gz  \
     --soloType CB_UMI_Complex \
     --soloAdapterSequence GTGANNNNNNNNNGACA \
     --soloCBposition 2_-9_2_-1 2_4_2_12 2_17_2_25 \
     --soloUMIposition 3_10_3_17 \
     --soloCBwhitelist BD_CLS1.txt BD_CLS2.txt BD_CLS3.txt \
     --soloCBmatchWLtype 1MM \
     --soloCellFilter None \
     --soloCellReadStats Standard \
     --sjdbOverhang 100 \
     --outSAMattributes sS sM CB UB CR CY UR UY \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode GeneCounts \
     --soloUMIlen 8 \
     --sjdbGTFfile genome.gtf

## same for TSO
STAR --runThreadN $NTHREADS \
     --genomeDir genome \
     --readFilesCommand zcat \
     --outFileNamePrefix mapping_tso/ \
     --readFilesIn  r1.fq.gz r2.fq.gz  \
     --soloType CB_UMI_Complex \
     --soloAdapterSequence AATGANNNNNNNNCCAC \
     --soloCBposition 2_-9_2_-1 2_4_2_12 2_17_2_25 \
     --soloUMIposition 3_10_3_17 \
     --soloCBwhitelist BD_CLS1.txt BD_CLS2.txt BD_CLS3.txt \
     --soloCBmatchWLtype 1MM \
     --soloCellFilter None \
     --soloCellReadStats Standard \
     --sjdbOverhang 100 \
     --outSAMattributes sS sM CB UB CR CY UR UY \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode GeneCounts \
     --soloUMIlen 8 \
     --sjdbGTFfile genome.gtf


## merge them and run an extra multifeature, multioverlap count for TSO data
##

cd ~/module/src

cp ./utilities/config_input_parameters.yml{,.backup}
cat <<EOF > ./utilities/config_input_parameters.yml

# input file paths
input_bam_file_paths: 
  - /home/rock/simulate/mapping_wta/Aligned.sortedByCoord.out.bam
  - /home/rock/simulate/mapping_tso/Aligned.sortedByCoord.out.bam
unique_cb_bamfile_list:               # for valid cell bar codes extraction
  - /home/rock/simulate/mapping_wta/Aligned.sortedByCoord.out.bam
gtf_file: /home/rock/simulate/genome.gtf
featureCounts_path: featureCounts

# input parameters
linker_pattern_for_seq_logo: ^([A-Z]{9})(AATG[A-Z]{9}CCAC) # tso pattern to get tso tail seq logo
linker_pattern_list:                  # list of linkers to be matched with alignment. If no match then CB tag is recovered where possible.
  - ^([A-Z]{9})(AATG[A-Z]{9}CCAC)     # tso linkers
  - ^(A|GT|TCA)?([A-Z]{9})(GTGA[A-Z]{9}GACA)  # wta linkers
seq_logo_filename: tso_seq_logo.png
cb1_len: 9                            # length of first part of CB tag
cb2_len: 9                            # length of second part of CB tag
cb3_len: 9                            # length of third part of CB tag
linker_left_len: 4                    # length of linker on left
linker_right_len: 4                   # length of linker on right
umi_len: 8                            # length of UMI in sS tag
chunk_genomic_area: 1000              # number of base pairs to be processed by each process in multiprocessing

# NOTE: number of rg tag names must be equal to the number of input bam files and in the order of input bam files
rg_tag_names:
  - wta
  - tso 

chromosomes:                          # chromosomes to be processed
  - offtarget
  - ontarget

subset_gtf: 0                         # subset gtf file by pattern (1: yes, 0: no)
subset_gtf_output_file: alien.gtf     # output file name of subset of gtf file
subset_gtf_pattern: 'ERCC'            # filter gtf file by this pattern if subset_gtf = 1
write_final_bam_to_csv: 1             # write final merged bam file to a csv file(1: yes, 0: no)
write_final_bam_header_to_txt: 1      # write final merged bam header to text file (1: yes, 0: no)
nProcessors: 4                        # for multiprocessing
nthreads: 4                           # for featureCounts

# output file names and paths
output_folder: ../simulated_output/
final_merged_file: merged.bam
featureCounts_output_file: featurecounted
log_file: log.txt
error_log: error_log.txt
EOF

python -m  main.combined_rock_roi

ls -l ../simulated_output
