#!/usr/bin/env snakemake -s
##
## Snakefile to process rock/roi data (general method)
##
## Started 11th Oct 2023
##
## Izaskun Mallona
## GPLv3

import os.path as op

configfile: "config.yaml"
           
if config['simulate']:
    config['gtf'] = op.join(config['working_dir'], 'data', 'genome.gtf')
    config['genome'] = op.join(config['working_dir'], 'data', 'genome.fa')
    config['samples'] = [{'name': 'simulated',
                          'uses': {
                              'cdna_fq': op.join(config['working_dir'], 'data', 'simulated', 'r1.fq.gz'),
                              'cb_umi_fq': op.join(config['working_dir'], 'data', 'simulated', 'r2.fq.gz'),
                              'whitelist': '96x3',
                              'expected_cells': 500}},
                         {'name': 'simulated2', 'uses': {
                             'cdna_fq': op.join(config['working_dir'], 'data', 'simulated', 'r1.fq.gz'),
                             'cb_umi_fq': op.join(config['working_dir'], 'data', 'simulated', 'r2.fq.gz'),
                             'whitelist': '96x3', 'expected_cells': 500}}]
    config['capture_gtf_column_2_pattern'] = 'captured'
    config['run_mode'] = 'all'


## to ease whitelists symlinking
if not op.isabs(config['rock_method_path']):
    config['rock_method_path'] = op.join(workflow.basedir, config['rock_method_path'])

def get_sample_names():
    return([x['name'] for x in config['samples']])

def get_cbumi_by_name(name):
    for i in range(len(config['samples'])):
        if config['samples'][i]['name'] == name:
             return(config['samples'][i]['uses']['cb_umi_fq'])

def get_cdna_by_name(name):
    for i in range(len(config['samples'])):
        if config['samples'][i]['name'] == name:
             return(config['samples'][i]['uses']['cdna_fq'])

def get_expected_cells_by_name(name):
    for i in range(len(config['samples'])):
        if config['samples'][i]['name'] == name:
             return(config['samples'][i]['uses']['expected_cells'])

def get_barcode_whitelist_by_name(name):
    for i in range(len(config['samples'])):
        if config['samples'][i]['name'] == name:
             return(config['samples'][i]['uses']['whitelist'])

print(get_sample_names())
         
## canonical CB 9-mers or just A{9}
## (used for simulations)
cb1s = ['GTCGCTATA','CTTGTACTA','CTTCACATA','ACACGCCGG','CGGTCCAGG','AATCGAATG','CCTAGTATA']
cb2s = ['TACAGGATA','CACCAGGTA','TGTGAAGAA','GATTCATCA','CACCCAAAG','CACAAAGGC','GTGTGTCGA']
cb3s = ['AAGCCTTCT','ATCATTCTG','CACAAGTAT','ACACCTTAG','GAACGACAA','AGTCTGTAC','AAATTACAG', 'AAAAAAAAA']
umis = ['AACCTTGG', 'CCGGTTAA', 'TTGGCCAA', 'GACATAGG']

rule all:
    input:
        op.join(config['working_dir'], 'multimodal', 'descriptive_report.html'),
        expand(op.join(config['working_dir'], 'multimodal', '{sample}', '{sample}_sce.rds'),
               sample = get_sample_names())#,
        # expand(op.join(config['working_dir'], 'align_{modality}', '{sample}', '{sample}_{modality}_coverage.bw'),
        #        modality = ['tso', 'wta'],
        #        sample = get_sample_names())
    
rule index:
    conda:
        "envs/all_in_one.yaml"
    input:
        gtf = config['gtf'],
        fa = config['genome']
    output:
        index_path =  op.join(config['working_dir'] , 'data', 'index', 'SAindex')
    threads:
        config['nthreads']
    params:
        simulate = config['simulate'],
        processing_path = op.join(config['working_dir'], 'data'),
        nthreads = config['nthreads'],
        star = config['STAR'],
        sjdbOverhang = config['sjdbOverhang']
    log:
        op.join(config['working_dir'], 'data', 'indexing.log')
    shell:
      """
    mkdir -p {params.processing_path}
    # cd {params.processing_path}

    ({params.star} --runThreadN {params.nthreads} \
     --runMode genomeGenerate \
     --sjdbGTFfile {input.gtf} \
     --genomeDir {params.processing_path}/index \
     --genomeSAindexNbases 4 \
     --sjdbOverhang {params.sjdbOverhang} \
     --genomeFastaFiles {input.fa} ) 2> {log}
        """

rule simulate_genome:
    output:
        op.join(config['working_dir'], 'data', 'genome.fa')
    shell:
        """
cat << EOF > {output}
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
        """

## notice the `captured` flag to focus on the last three GTF records during the custom featurecounting
rule simulate_gtf:
    output:
        op.join(config['working_dir'], 'data', 'genome.gtf')
    shell:
        """
echo -e 'offtarget\tERCC\texon\t1\t1061\t.\t+\t.\tgene_id "offtarget"; transcript_id "offtarget_1";
ontarget\tcaptured\texon\t1\t1023\t.\t+\t.\tgene_id "ontarget_1"; transcript_id "ontarget_1";
ontarget\tcaptured\texon\t100\t800\t.\t+\t.\tgene_id "ontarget_1b"; transcript_id "ontarget_1b";
ontarget\tcaptured\texon\t1030\t1090\t.\t+\t.\tgene_id "ontarget_2"; transcript_id "ontarget_2";' > {output}
        """

rule simulate_fastqs:
    conda:
        "envs/all_in_one.yaml"
    input:
        part_r1 = expand(op.join(config['working_dir'], 'data', '{sample}',
                                 'part_{cb1}_{cb2}_{cb3}_{umi}_r1.fq.gz'),
                         sample = get_sample_names(),
                         cb1 = cb1s, cb2 = cb2s, cb3 = cb3s, umi = umis),
        part_r2 = expand(op.join(config['working_dir'], 'data', '{sample}',
                                 'part_{cb1}_{cb2}_{cb3}_{umi}_r2.fq.gz'),
                         sample = get_sample_names(),
                         cb1 = cb1s, cb2 = cb2s, cb3 = cb3s, umi = umis)
    params:
        path = op.join(config['working_dir'], 'data', '{sample}')
    output:
        r1 = op.join(config['working_dir'], 'data', '{sample}', 'r1.fq.gz'),
        r2 = op.join(config['working_dir'], 'data', '{sample}', 'r2.fq.gz'),
        r1_extra = temp(op.join(config['working_dir'], 'data', '{sample}', 'r1_extra.fq.gz')),
        r2_extra = temp(op.join(config['working_dir'], 'data', '{sample}', 'r2_extra.fq.gz')),

    shell:
        """
# cd {params.path}

seq 1 100 | awk ' {{print "@extra"$0"\\nTACTGGACTCGCAACGTGGGTCCATCAGTTGTCCGTATACCAAGACGTCTAAGGGCGGTGTACACCCTTTTGAGCAATGATTGCACAACCTGCGATCACCTTATACAGAATTATCAATCAAGCTCCCCGAGGAGCGGACTTGTAAGGACCGCCGCTTTCGCTCGGGTCTGCGG\\n+\\nFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF";}}' | gzip -c > {output.r1_extra}

seq 1 100 | awk ' {{print "@extra"$0"\\nCTTGTACTAGTGATGTTCTCCAGACAGGCTACAGATTTGATGGTTTTTTTTTTTTTTTTT\\n+\\nFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF";}}' | gzip -c > {output.r2_extra}

cat {input.part_r1} {output.r1_extra} > {output.r1}
cat {input.part_r2} {output.r2_extra} > {output.r2}
      """

## notice we have a ~10x duplication rate (UMIs are reused multiple times)      
rule simulate_fastqs_from_a_cell:
    output:
        part_r1 =  temp(op.join(config['working_dir'], 'data', '{sample}',
                                'part_{cb1}_{cb2}_{cb3}_{umi}_r1.fq.gz')),
        part_r2 =  temp(op.join(config['working_dir'], 'data', '{sample}',
                                'part_{cb1}_{cb2}_{cb3}_{umi}_r2.fq.gz'))
    params:
        path = op.join(config['working_dir'], 'data', '{sample}'),
        tso_fix1 = "AATG",
        tso_fix2 = "CCAC",
        wta_fix1 = "GTGA",
        wta_fix2 = "GACA",
    threads:
        config['nthreads']
    shell:
        """
mkdir -p {params.path}

seq 1 10 | awk  '{{print "@wta"$0"_{wildcards.cb1}_{wildcards.cb2}_{wildcards.cb3}__{wildcards.umi}\\nTACTGGACTCGCAACGTGGGTCCATCAGTTGTCCGTATACCAAGACGTCTAAGGGCGGTGTACACCCTTTTGAGCAATGATTGCACAACCTGCGATCACCTTATACAGAATTATCAATCAAGCTCCCCGAGGAGCGGACTTGTAAGGACCGCCGCTTTCGCTCGGGTCTGCGG\\n+\\nFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF";}}' | gzip -c > {output.part_r1}

seq 1 10 | awk ' {{print "@wta"$0"_{wildcards.cb1}_{wildcards.cb2}_{wildcards.cb3}__{wildcards.umi}\\nA{wildcards.cb1}{params.wta_fix1}{wildcards.cb2}{params.wta_fix2}{wildcards.cb3}{wildcards.umi}AGATTTGATGGTTTTT\\n+\\nFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF";}}' | gzip -c > {output.part_r2}

## similarly for region 1000-1100 for ontargets/tso

seq 11 20 | awk ' {{print "@tso"$0"_{wildcards.cb1}_{wildcards.cb2}_{wildcards.cb3}__{wildcards.umi}\\nATTATTGAAGAGTTAGAACACCACAAGATAAAAAAAAAAAAAAAAAAAAAAAACAGCAGCGATTAAGGCAGAGGCGTTTGTATCTGCCATTATAAAGAAGTTTCCTCCAGCAACTCCTTTCTTAATTCCAAACTTAGC\\n+\\nFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF";}}' | gzip -c >> {output.part_r1}

seq 11 20 | awk  ' {{print "@tso"$0"_{wildcards.cb1}_{wildcards.cb2}_{wildcards.cb3}__{wildcards.umi}\\nT{wildcards.cb1}{params.tso_fix1}{wildcards.cb2}{params.tso_fix2}{wildcards.cb3}{wildcards.umi}TTTTTTTTTTTTTTTT\\n+\\nFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF";}}' | gzip -c >> {output.part_r2}

## another ontarget region - TSO barcode
## same UMI 100 times
seq 21 30 | awk ' {{print "@tso"$0"_{wildcards.cb1}_{wildcards.cb2}_{wildcards.cb3}__{wildcards.umi}\\nATTCCAAACTTAGCTTCAGTTATAAATTCCCCTCCCATGATTGGGATTTTATAAACTTTTCTTCCATATAATTCATCTTTCTTCTCATAACCGTCTCCGAAAAACTTCAACTTAAATCCAACCTTTAACTGCTCATCA\\n+\\nFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF";}}' | gzip -c >> {output.part_r1}

seq 21 30 | awk ' {{print "@tso"$0"_{wildcards.cb1}_{wildcards.cb2}_{wildcards.cb3}__{wildcards.umi}\\nG{wildcards.cb1}{params.tso_fix1}{wildcards.cb2}{params.tso_fix2}{wildcards.cb3}{wildcards.umi}TTTTTTTTTTTTTTTT\\n+\\nFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF";}}' | gzip -c >> {output.part_r2}

## let's add one read for ontargets on WTA

seq 31 32 | awk ' {{print "@wta"$0"_{wildcards.cb1}_{wildcards.cb2}_{wildcards.cb3}__{wildcards.umi}\\nATTCCAAACTTAGCTTCAGTTATAAATTCCCCTCCCATGATTGGGATTTTATAAACTTTTCTTCCATATAATTCATCTTTCTTCTCATAACCGTCTCCGAAAAACTTCAACTTAAATCCAACCTTTAACTGCTCATCA\\n+\\nFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF";}}' | gzip -c >> {output.part_r1}

seq 31 32 | awk ' {{print "@wta"$0"_{wildcards.cb1}_{wildcards.cb2}_{wildcards.cb3}__{wildcards.umi}\\nG{wildcards.cb1}{params.wta_fix1}{wildcards.cb2}{params.wta_fix2}{wildcards.cb3}{wildcards.umi}TTTTTTTTTTTTTTTT\\n+\\nFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF";}}' | gzip -c >> {output.part_r2}
      """

## bd offers a couple of sets of whitelists, so we fetch the right one according to the config.yaml file
def symlink_whitelist(sample):
    os.makedirs(op.join(config['working_dir'], 'align_wta'), exist_ok = True)
                
    if get_barcode_whitelist_by_name(name = sample) == '96x3':
        for x in ['BD_CLS1.txt', 'BD_CLS2.txt', 'BD_CLS3.txt']:
            try:
                os.symlink(src = op.join(config['rock_method_path'], 'data', 'whitelist_96x3', x),
                           dst = op.join(config['working_dir'], 'align_wta', sample, 'whitelists', x))
            except FileExistsError:
                break
    elif get_barcode_whitelist_by_name(name = sample) == '384x3':
        for x in ['BD_CLS1.txt', 'BD_CLS2.txt', 'BD_CLS3.txt']:
            try:
                os.symlink(src = op.join(config['rock_method_path'], 'data', 'whitelist_384x3', x),
                           dst = op.join(config['working_dir'], 'align_wta', sample, 'whitelists', x))
            except FileExistsError:
                break

rule prepare_whitelists:
    # conda:
    #     "envs/all_in_one.yaml"
    input:
        cbumi = lambda wildcards: get_cbumi_by_name(wildcards.sample),
        cdna = lambda wildcards: get_cdna_by_name(wildcards.sample)
        # r1 = op.join(config['working_dir'], 'data', "{sample}", 'r1.fq.gz'),
        # r2 = op.join(config['working_dir'], 'data', "{sample}", 'r2.fq.gz')
    output:
        cb1 = op.join(config['working_dir'], 'align_wta', "{sample}", 'whitelists', 'BD_CLS1.txt'),
        cb2 = op.join(config['working_dir'], 'align_wta', "{sample}", 'whitelists', 'BD_CLS2.txt'),
        cb3 = op.join(config['working_dir'], 'align_wta', "{sample}", 'whitelists', 'BD_CLS3.txt')
    run:
        sample = wildcards.sample
        symlink_whitelist(sample)


rule align_wta:
    conda:
        "envs/all_in_one.yaml"
    input:
        # r1 = op.join(config['working_dir'], 'data', "{sample}", 'r1.fq.gz'),
        # r2 = op.join(config['working_dir'], 'data', "{sample}", 'r2.fq.gz'),
        cdna = lambda wildcards: get_cdna_by_name(wildcards.sample),
        cbumi = lambda wildcards: get_cbumi_by_name(wildcards.sample),
        index_flag = op.join(config['working_dir'] , 'data', 'index', 'SAindex'),
        gtf = config['gtf'],
        cb1 = op.join(config['working_dir'], 'align_wta', "{sample}",  'whitelists', 'BD_CLS1.txt'),
        cb2 = op.join(config['working_dir'], 'align_wta', "{sample}", 'whitelists', 'BD_CLS2.txt'),
        cb3 = op.join(config['working_dir'], 'align_wta', "{sample}", 'whitelists', 'BD_CLS3.txt')
    output:
        bam = op.join(config['working_dir'], 'align_wta', '{sample}', 'Aligned.sortedByCoord.out.bam'),
        filtered_barcodes = op.join(config['working_dir'], 'align_wta', '{sample}', 'Solo.out', 'Gene',
                                    'filtered', 'barcodes.tsv')
    threads: workflow.cores
    params:
        threads = min(10, workflow.cores),
        path = op.join(config['working_dir'], 'align_wta', "{sample}/"),
        index_path = op.join(config['working_dir'] , 'data', 'index'),
        STAR = config['STAR'],
        # num_cells = get_expected_cells_by_name("{sample}"),
        tmp = op.join(config['working_dir'], 'tmp_align_wta_{sample}'),
        maxmem = config['max_mem_mb'] * 1024 * 1024,
        sjdbOverhang = config['sjdbOverhang']
    shell:
        """
   rm -rf {params.tmp}
   mkdir -p {params.path} 

   {params.STAR} --runThreadN {params.threads} \
     --genomeDir {params.index_path} \
     --readFilesCommand zcat \
     --outFileNamePrefix {params.path} \
     --readFilesIn  {input.cdna} {input.cbumi}  \
     --soloType CB_UMI_Complex \
     --soloAdapterSequence GTGANNNNNNNNNGACA \
     --soloCBposition 2_-9_2_-1 2_4_2_12 2_17_2_25 \
     --soloUMIposition 3_10_3_17 \
     --soloCBwhitelist {input.cb1} {input.cb2} {input.cb3} \
     --soloCBmatchWLtype 1MM \
     --soloCellFilter EmptyDrops_CR \
     --outSAMattributes NH HI AS nM NM MD jM jI MC ch CB UB gx gn sS CR CY UR UY \
     --soloCellReadStats Standard \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode GeneCounts \
     --soloUMIlen 8 \
     --sjdbGTFfile {input.gtf} \
     --outTmpDir {params.tmp} \
     --sjdbOverhang {params.sjdbOverhang} \
     --limitBAMsortRAM {params.maxmem}

    rm -rf {params.tmp}
        """

        
# rule generate_tso_whitelist_from_wta_filtered:
#     conda:
#        "envs/all_in_one.yaml"
#    input:
#         filtered_barcodes = op.join(config['working_dir'], 'align_wta', '{sample}', 'Solo.out', 'Gene',
#                                     'filtered', 'barcodes.tsv')
#     output:
#         tso_whitelist = op.join(config['working_dir'], 'data', '{sample}', 'tso_whitelist.txt')
#     params:
#         path = op.join(config['working_dir'], 'data', '{sample}')
#     shell:
#         """
#         mkdir -p {params.path}
#         awk -F "_" '{{print $1"AATG"$2"CCAC"$3}}' {input.filtered_barcodes} > {output.tso_whitelist}
#         """
        
rule align_tso:
    conda:
        "envs/all_in_one.yaml"
    input:
        cdna = lambda wildcards: get_cdna_by_name(wildcards.sample),
        cbumi = lambda wildcards: get_cbumi_by_name(wildcards.sample),
        index_flag = op.join(config['working_dir'] , 'data', 'index', 'SAindex'),
        gtf = config['gtf'],
        # filtered_barcodes = op.join(config['working_dir'], 'align_wta', '{sample}', 'Solo.out', 'Gene',
        #                             'filtered', 'barcodes.tsv'),
        cb1 = op.join(config['working_dir'], 'align_wta', "{sample}",  'whitelists', 'BD_CLS1.txt'),
        cb2 = op.join(config['working_dir'], 'align_wta', "{sample}", 'whitelists', 'BD_CLS2.txt'),
        cb3 = op.join(config['working_dir'], 'align_wta', "{sample}", 'whitelists', 'BD_CLS3.txt')
    output:
        #bam = temp(op.join(config['working_dir'], 'align_tso', '{sample}', 'Aligned.sortedByCoord.out.bam'))
        bam = op.join(config['working_dir'], 'align_tso', '{sample}', 'Aligned.sortedByCoord.out.bam')
    threads:  workflow.cores
    params:
        threads = min(10, workflow.cores),
        path = op.join(config['working_dir'], 'align_tso', "{sample}/"),
        index_path = op.join(config['working_dir'] , 'data', 'index'),
        STAR = config['STAR'],
        tmp = op.join(config['working_dir'], 'tmp_align_tso_{sample}'),
        maxmem = config['max_mem_mb'] * 1024 * 1024, 
        sjdbOverhang = config['sjdbOverhang']
    shell:
        """
   rm -rf {params.tmp}
   mkdir -p {params.path} 

        ## with this whitelist the BAM looks good in CR terms but the count tables are broken (!)
        ## not reporting the issue 
        # #--soloCBwhitelist input.tso_whitelist
        #      # --soloCBmatchWLtype Exact 

        {params.STAR} --runThreadN {params.threads} \
        --genomeDir {params.index_path} \
        --readFilesCommand zcat \
        --outFileNamePrefix {params.path} \
        --readFilesIn  {input.cdna} {input.cbumi}  \
        --soloType CB_UMI_Complex \
         --soloAdapterSequence AATGNNNNNNNNNCCAC \
         --soloCBposition 2_-9_2_-1 2_4_2_12 2_17_2_25 \
         --soloUMIposition 3_10_3_17 \
         --soloUMIlen 8 \
         --soloCellReadStats Standard \
         --soloCBwhitelist {input.cb1} {input.cb2} {input.cb3} \
         --soloCBmatchWLtype 1MM \
         --soloCellFilter None \
         --outSAMattributes NH HI AS nM NM MD jM jI MC ch CB UB gx gn sS CR CY UR UY\
         --outSAMtype BAM SortedByCoordinate \
         --quantMode GeneCounts \
        --sjdbGTFfile {input.gtf} \
        --outTmpDir {params.tmp} \
        --sjdbOverhang {params.sjdbOverhang} \
        --limitBAMsortRAM {params.maxmem}

        rm -rf {params.tmp}
        """

        
## this is a dirty workaround to reduce mem usage by selecting which chromosomes have 'captured' features
## @todo filter the coordinates, not only the chromosomes
rule subset_chromosomes_for_custom_counting:
    conda:
        "envs/all_in_one.yaml"
    input:
        bam = op.join(config['working_dir'], 'align_{modality}', '{sample}', 'cb_ub_filt.bam')
    output:
        red_bam = op.join(config['working_dir'], 'align_{modality}', '{sample}', 'subset_{modality}.bam'),
        chrs = temp(op.join(config['working_dir'], 'align_{modality}', '{sample}', 'captured_chrs.txt'))        
    threads: 10        
    params:
        gtf = config['gtf'],
        run_mode = config['run_mode'],
        pattern = config['capture_gtf_column_2_pattern'],
        temp_red_bam = op.join(config['working_dir'], 'align_{modality}', '{sample}',
                               'tmp_subset_{modality}.bam')
    # run:
    #     if params.run_mode in ['all', 'tso ontarget multi']:
    #         shell ("""
    #         grep {params.pattern} {params.gtf} | cut -f1 | sort | uniq  > {output.chrs}
    #         captured=$(cat {output.chrs} | tr '\\n' ' ')

    #         samtools index -@ {threads} {input.bam}

    #         samtools view -h {input.bam}  $captured -@ {threads} | \
    #             samtools view -Sb > {output.red_bam}
    #         """)
    shell:
        """
        run_mode={params.run_mode}

        if [ "$run_mode" = "all" ] || [ "$run_mode" = "tso_ontarget_multi" ]
        then
            grep {params.pattern} {params.gtf} | cut -f1 | sort | uniq  > {output.chrs}
            captured=$(cat {output.chrs} | tr '\\n' ' ')

            samtools index -@ {threads} {input.bam}

            samtools view -h {input.bam}  $captured -@ {threads} | \
                samtools view -Sb > {output.red_bam}
        fi
        """

rule subset_gtf_for_custom_counting:
    conda:
        "envs/all_in_one.yaml"
    input:
        gtf = config['gtf']
    output:
        subset_gtf = temp(op.join(config['working_dir'], 'multimodal', 'subset.gtf'))
    threads: 10
    params:
        pattern = config['capture_gtf_column_2_pattern']
    shell:
        """
        grep {params.pattern} {input.gtf} > {output.subset_gtf}
        """

rule add_readgroups_to_bam:
    conda:
        "envs/all_in_one.yaml"
    input:
        without_rgs = op.join(config['working_dir'], 'align_{modality}', '{sample}', 'subset_{modality}.bam'),
        filtered_barcodes = op.join(config['working_dir'], 'align_wta', '{sample}', 'Solo.out', 'Gene',
                                    'filtered', 'barcodes.tsv')
    output:
        header = temp(op.join(config['working_dir'], 'align_{modality}',
                                       '{sample}', 'header')),
        header_barcodes = temp(op.join(config['working_dir'], 'align_{modality}',
                                       '{sample}', 'header_rgs')),
        with_rgs_no_header = temp(op.join(config['working_dir'], 'align_{modality}', '{sample}',
                                          'subset_{modality}_rg_no_header.bam')),
        with_rgs = temp(op.join(config['working_dir'], 'align_{modality}', '{sample}',
                                'subset_{modality}_rg.bam'))
    threads: config['nthreads']
    shell:
        """
        ## replace the CBs by RGs
        samtools view -h {input.without_rgs} -@ {threads} | \
             sed 's/CB:Z:/RG:Z:tso_/g' | samtools view -h -Sb >  {output.with_rgs_no_header}

        ## get the old bam header
        samtools view -H {input.without_rgs} > {output.header}

        ## generate a RG header
        TAB="$(printf '\t')"

        while read -r barcode
        do
           # echo "@RG${{TAB}}ID:wta_${{barcode}}"
           echo "@RG${{TAB}}ID:tso_${{barcode}}"

        done < {input.filtered_barcodes} > {output.header_barcodes}

        samtools view {output.with_rgs_no_header} | \
          cat {output.header} {output.header_barcodes} - | \
          samtools view -Sb -@ {threads} | \
          samtools sort -@ {threads}  > {output.with_rgs}

        """

rule count_custom_regions_tso_no_module:
    conda:
        "envs/all_in_one.yaml"
    input:
        bam = op.join(config['working_dir'], 'align_tso', '{sample}', 'subset_tso_rg.bam'),
        gtf = op.join(config['working_dir'], 'multimodal', 'subset.gtf')
    output:
        fc = op.join(config['working_dir'], 'multimodal', '{sample}', 'featurecounted')
    threads: config['nthreads']
    resources:
         mem_mb=config['max_mem_mb']
    log: op.join(config['working_dir'], 'multimodal', '{sample}', 'featurecounts.log')
    params:
        featureCounts = config['featureCounts'],
        run_mode = config['run_mode'],
        overallconfig = 'config.yaml',
        max_mem = lambda wildcards, resources: resources.mem_mb * 1024,
        t = config['featurecounts_t'],
        g = config['featurecounts_g']
    # run:
    #     if params.run_mode in ['all', 'tso ontarget multi']:
    #         shell ("""
    #         ulimit -v {params.max_mem}
            
    #         ## featurecounts, notice the -M and -T and --fraction
    #         {params.featureCounts} \
    #              -a {input.gtf} \
    #              -o {output.fc} \
    #              {input.bam} \
    #              -F GTF \
    #              -t {params.t} \
    #              -g {params.g} \
    #              -f \
    #              -O \
    #              -M  \
    #              -T {threads} \
    #              --fraction \
    #              --byReadGroup &> {log}
    #         """)
    shell:
        """
        run_mode={params.run_mode}

        if [ "$run_mode" = "all" ] || [ "$run_mode" = "tso_ontarget_multi" ]
        then
            ulimit -v {params.max_mem}
            
            ## featurecounts, notice the -M and -T and --fraction
            {params.featureCounts} \
                 -a {input.gtf} \
                 -o {output.fc} \
                 {input.bam} \
                 -F GTF \
                 -t {params.t} \
                 -g {params.g} \
                 -f \
                 -O \
                 -M  \
                 -T {threads} \
                 --fraction \
                 --byReadGroup &> {log}
        fi
        """
        
rule retrieve_genome_sizes:
    conda:
        "envs/all_in_one.yaml"
    input:
        fa = config['genome']
    params:
        faSize = config['faSize']
    output:
        op.join(config['working_dir'], 'data', 'chrom.sizes')
    shell:
        """
        {params.faSize} -detailed -tab {input.fa} > {output}
        """
        
rule create_deduped_coverage_tracks_all_filtered_in_cbs:
    conda:
        "envs/all_in_one.yaml"
    input:
        bam = op.join(config['working_dir'], 'align_{modality}', '{sample}', 'Aligned.sortedByCoord.out.bam'),
        valid_barcodes = op.join(config['working_dir'], 'align_wta', '{sample}', 'Solo.out', 'Gene',
                                 'filtered', 'barcodes.tsv'),
        chromsizes = op.join(config['working_dir'], 'data', 'chrom.sizes')
    threads:
        config['nthreads']
    params:
        # bamCoverage = config['bamCoverage'],
        binSize = 10,
        bedtools = config['bedtools'],
        bedGraphToBigWig = config['bedGraphToBigWig']        
    output:
        cb_bam = temp(op.join(config['working_dir'], 'align_{modality}', '{sample}', 'cb_filt.bam')),
        cb_ub_bam = temp(op.join(config['working_dir'], 'align_{modality}', '{sample}', 'cb_ub_filt.bam')),
        header  = temp(op.join(config['working_dir'], 'align_{modality}', '{sample}', 'header')),
        cb_ub_bg = temp(op.join(config['working_dir'], 'align_{modality}', '{sample}', 'cb_ub_filt.bw')),
        bw = op.join(config['working_dir'], 'align_{modality}', '{sample}', '{sample}_{modality}_coverage.bw')
    shell:
        """
        ## this is unrelated to the bamgeneration; fixes starsolo's default permissions
        chmod -R ug+rwX $(dirname {input.bam})

        ## first filter in 'valid' CBs

        ## CB are the error-corrected barcodes
        samtools view -h -@ {threads} {input.bam} -D CB:{input.valid_barcodes} \
           -o {output.cb_bam}

        ## second, deduplicate by UB (error corrected barcodes)
        ##  so it gets only one alignment per UB:locus:CB combination
        # 27 CB (error corrected CB)
        # 28 UB (error corrected UMI)
        # mind the file is coordinate sorted already - so we deduplicate based on columns 27 and 28
        samtools view -H {output.cb_bam} > {output.header}

        samtools view {output.cb_bam} -@ {threads} | uniq -f26 | cat {output.header} - | \
          samtools view -Sb -@ {threads} > {output.cb_ub_bam}

        samtools index -@ {threads} {output.cb_ub_bam}

        {params.bedtools} genomecov -ibam {output.cb_ub_bam} \
            -bg -split | LC_COLLATE=C sort -k1,1 -k2,2n > {output.cb_ub_bg}

        ## bedgraph to bigwig
        {params.bedGraphToBigWig} {output.cb_ub_bg} {input.chromsizes} {output.bw}
        
        """

# rule create_deduped_coverage_tracks_single_cb:
#     conda:
#        "envs/all_in_one.yaml"
#    input:
#         bam = op.join(config['working_dir'], 'align_{modality}', '{sample}', 'Aligned.sortedByCoord.out.bam'),
#         valid_barcodes = op.join(config['working_dir'], 'align_wta', '{sample}', 'Solo.out', 'Gene',
#                                  'filtered', 'barcodes.tsv')
#     threads:
#         config['nthreads']
#     params:
#         bamCoverage = config['bamCoverage'],
#         binSize = 10
#     output:
#         cb_bam = temp(op.join(config['working_dir'], 'align_{modality}', '{sample}', 'cb_{cell}_filt.bam')),
#         cb_ub_bam = temp(op.join(config['working_dir'], 'align_{modality}', '{sample}',
#                                  'cb_{cell}_ub_filt.bam')),
#         header  = temp(op.join(config['working_dir'], 'align_{modality}', '{sample}', '{cell}_header')),
#         bw = op.join(config['working_dir'], 'align_{modality}', '{sample}', '{sample}_cb_{cell}_coverage.bw')
#     shell:
#         """
#         ## first filter in the query CB
#         samtools view -h -@ {threads} {input.bam} -d CB:{wildcards.cell} \
#            -o {output.cb_bam}

#         ## second, deduplicate by UB (error corrected barcodes)
#         ##  so it gets only one alignment per UB:locus:CB combination
#         # 27 CB (error corrected CB; the only CB being queried)
#         # 28 UB (error corrected UMI)
#         # mind the file is coordinate sorted already - so we deduplicate based on columns 27 and 28
#         samtools view -H {output.cb_bam} > {output.header}

#         samtools view -h {output.cb_bam} -@ {threads} | uniq -f26 | cat {output.header} - | \
#           samtools view -Sb -@ {threads} > {output.cb_ub_bam}

#         samtools index -@ {threads} {output.cb_ub_bam}

#         # deeptools coverage it, no normalization
#         {params.bamCoverage} --numberOfProcessors {threads} \
#           -b {output.cb_ub_bam} \
#           --normalizeUsing None \
#           -o {output.bw} --binSize {params.binSize}
        
#         """

rule install_r_deps:
    conda:
        "envs/all_in_one.yaml"
    input:
        script = op.join(config['rock_method_path'], 'installs.R')
    output:
        log = op.join(config['working_dir'], 'log', 'installs.log')
    params:
        run_mode = config['run_mode'],
        working_dir = config['working_dir'],
        sample = "{wildcards.sample}",
        Rbin = config['Rbin'],
        log_path = op.join(config['working_dir'], 'log')
    shell:
        """
        mkdir -p {params.log_path}

        {params.Rbin} -q --no-save --no-restore --slave \
             -f {input.script} &> {output.log}
        """
        
rule generate_sce:
    conda:
        "envs/all_in_one.yaml"
    input:
        # tso_bam = op.join(config['working_dir'], 'align_tso', '{sample}', 'Aligned.sortedByCoord.out.bam'),
        # wta_bam = op.join(config['working_dir'], 'align_wta', '{sample}', 'Aligned.sortedByCoord.out.bam'),
        fc = op.join(config['working_dir'], 'multimodal', '{sample}', 'featurecounted'),
        # config = op.join(config['working_dir'], 'multimodal', '{sample}', 'config.yaml'),
        gtf = config['gtf'],
        script = op.join(config['rock_method_path'], 'generate_sce_object.R'),
        installs = op.join(config['working_dir'], 'log', 'installs.log')
    output:
        sce = op.join(config['working_dir'], 'multimodal', '{sample}', '{sample}_sce.rds')
    params:
        run_mode = config['run_mode'],
        working_dir = config['working_dir'],
        sample = "{wildcards.sample}",
        Rbin = config['Rbin']
    # run:
    #     if params.run_mode in ['all', 'tso ontarget multi']:
    #         shell ("""
    #         echo 'all or multi' > {output.sce}
    #         """)
    #     else:
    #         shell(" echo 'something else' > {output.sce}")
    shell:
        """
        {params.Rbin} -q --no-save --no-restore --slave \
             -f {input.script} --args --sample {wildcards.sample} \
             --run_mode {params.run_mode} \
             --working_dir {params.working_dir} \
             --output_fn {output.sce}
        """


rule render_descriptive_report:
    conda:
        "envs/all_in_one.yaml"
    input:
        gtf = config['gtf'],
        script = op.join(config['rock_method_path'], 'process_sce_objects.Rmd'),
        sces = expand(op.join(config['working_dir'], 'multimodal', '{sample}', '{sample}_sce.rds'),
               sample = get_sample_names()),
        installs = op.join(config['working_dir'], 'log', 'installs.log')
    output:
        html = op.join(config['working_dir'], 'multimodal', 'descriptive_report.html'),
        cache = temp(op.join(config['rock_method_path'], 'process_sce_objects_cache')),
        cached_files = temp(op.join(config['rock_method_path'], 'process_sce_objects_files'))
    log: op.join(config['working_dir'], 'multimodal', 'descriptive_report.log')
    params:
        multimodal_path = op.join(config['working_dir'], 'multimodal'),
        run_mode = config['run_mode'],
        working_dir = config['working_dir'],
        sample = "{wildcards.sample}",
        Rbin = config['Rbin']
    shell:
        """
        {params.Rbin} -e 'rmarkdown::render(\"{input.script}\", 
          output_file = \"{output.html}\", 
          params = list(multimodal_path = \"{params.multimodal_path}\", 
                        run_mode = \"{params.run_mode}\"))' &> {log}
        """