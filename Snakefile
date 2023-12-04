#!/usr/bin/env snakemake -s
##
## Snakefile to process rock/roi data (general method)
##
## Started 11th Oct 2023
##
## Izaskun Mallona
## GPLv3

import os.path as op

include: op.join('src', 'simulate.snmk')
    
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

print(get_sample_names())

rule all:
    input:
        op.join(config['working_dir'], 'multimodal', 'descriptive_report.html'),
        expand(op.join(config['working_dir'], 'multimodal', '{sample}', '{sample}_sce.rds'),
               sample = get_sample_names()),
        expand(op.join(config['working_dir'], 'align_{modality}', '{sample}', '{sample}_{modality}_coverage.bw'),
               modality = ['tso', 'wta'],
               sample = get_sample_names())
    
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

    samtools index -@ {threads} {output.bam} 

    rm -rf {params.tmp}
        """

        
rule generate_tso_whitelist_from_wta_filtered:
    conda:
        "envs/all_in_one.yaml"
    input:
        filtered_barcodes = op.join(config['working_dir'], 'align_wta', '{sample}', 'Solo.out', 'Gene',
                                    'filtered', 'barcodes.tsv')
    output:
        tso_w1 = op.join(config['working_dir'], 'data', '{sample}', 'tso_whitelist_1.txt'),
        tso_w2 = op.join(config['working_dir'], 'data', '{sample}', 'tso_whitelist_2.txt'),
        tso_w3 = op.join(config['working_dir'], 'data', '{sample}', 'tso_whitelist_3.txt')
    params:
        path = op.join(config['working_dir'], 'data', '{sample}')
    shell:
        """
        mkdir -p {params.path}
        awk -F "_" '{{print $1}}' {input.filtered_barcodes} > {output.tso_w1}
        awk -F "_" '{{print $2}}' {input.filtered_barcodes} > {output.tso_w2}
        awk -F "_" '{{print $3}}' {input.filtered_barcodes} > {output.tso_w3}
        """
        
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
        # cb1 = op.join(config['working_dir'], 'align_wta', "{sample}",  'whitelists', 'BD_CLS1.txt'),
        # cb2 = op.join(config['working_dir'], 'align_wta', "{sample}", 'whitelists', 'BD_CLS2.txt'),
        # cb3 = op.join(config['working_dir'], 'align_wta', "{sample}", 'whitelists', 'BD_CLS3.txt')
        cb1 = op.join(config['working_dir'], 'data', '{sample}', 'tso_whitelist_1.txt'),
        cb2 = op.join(config['working_dir'], 'data', '{sample}', 'tso_whitelist_2.txt'),
        cb3 = op.join(config['working_dir'], 'data', '{sample}', 'tso_whitelist_3.txt')
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

        samtools index -@ {threads} {output.bam} 

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

rule count_custom_regions_no_module:
    conda:
        "envs/all_in_one.yaml"
    input:
        bam = op.join(config['working_dir'], 'align_{modality}', '{sample}', 'subset_{modality}_rg.bam'),
        gtf = op.join(config['working_dir'], 'multimodal', 'subset.gtf')
    output:
        fc = op.join(config['working_dir'], 'multimodal', '{sample}', '{modality}_featurecounted')
    threads: min(64, config['nthreads']) # featurecounts has a max 64 for -T
    resources:
         mem_mb=config['max_mem_mb']
    log: op.join(config['working_dir'], 'multimodal', '{sample}', '{modality}_featurecounts.log')
    params:
        featureCounts = config['featureCounts'],
        run_mode = config['run_mode'],
        overallconfig = 'config.yaml',
        max_mem = lambda wildcards, resources: resources.mem_mb * 1024,
        t = config['featurecounts_t'],
        g = config['featurecounts_g']
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
        
checkpoint retrieve_genome_sizes:
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

rule index_bam:
    conda:
        "envs/all_in_one.yaml"
    input:
        bam = op.join(config['working_dir'], 'align_{modality}', '{sample}', 'Aligned.sortedByCoord.out.bam')
    output:
        bai = op.join(config['working_dir'], 'align_{modality}', '{sample}',
                      'Aligned.sortedByCoord.out.bam.bai')
    threads: workflow.cores
    shell:
        """
        samtools index -@ {threads} {input.bam}     
        """
        
## by chrom
rule split_by_chr:
    conda:
        "envs/all_in_one.yaml"
    input:
        bam = op.join(config['working_dir'], 'align_{modality}', '{sample}', 'Aligned.sortedByCoord.out.bam'),
        bai = op.join(config['working_dir'], 'align_{modality}', '{sample}',
                      'Aligned.sortedByCoord.out.bam.bai'),
        valid_barcodes = op.join(config['working_dir'], 'align_wta', '{sample}', 'Solo.out', 'Gene',
                                 'filtered', 'barcodes.tsv')
    output:
        mini = temp(op.join(config['working_dir'], 'align_{modality}', '{sample}',
                            '{chrom}_subset.bam'))
        # sorted = temp(op.join(config['working_dir'], 'align_{modality}', '{sample}',
        #                     '{chrom}_cb_umi_sorted.bam'))
    threads: min(10, workflow.cores * 0.1) # workaround to reduce IO bottleneck
    shell:
        """
        samtools view -h -b {input.bam} {wildcards.chrom} > {output.mini}      
        """

## by chrom
rule dedup_by_cb_umi_gx:
    conda:
        "envs/all_in_one.yaml"
    input:
        # bam = op.join(config['working_dir'], 'align_{modality}', '{sample}', 'Aligned.sortedByCoord.out.bam'),
        mini = op.join(config['working_dir'], 'align_{modality}', '{sample}',
                            '{chrom}_subset.bam'),
        chromsizes = op.join(config['working_dir'], 'data', 'chrom.sizes')
    output:
        cb_ub_bam = temp(op.join(config['working_dir'], 'align_{modality}', '{sample}',
                            '{chrom}_cb_umi_deduped.bam')),
        header = temp(op.join(config['working_dir'], 'align_{modality}', '{sample}',
                            '{chrom}_cb_umi_deduped_header.txt'))
    threads: min(10, workflow.cores * 0.1) # workaround to reduce IO bottleneck
    shell:
        """
        samtools view -H {input.mini} > {output.header}

        # 20 gx (gene id)
        # 27 CB (error corrected CB)
        # 28 UB (error corrected UMI)
        # the chromosome is implicit - from the per-chr run
        samtools view {input.mini} |  sort -k27 -k28 -k 20 -u | cat {output.header} - | \
           samtools view -Sbh > {output.cb_ub_bam}  

        """

rule merge_deduped_bams:
    conda:
        "envs/all_in_one.yaml"
    input:
        # chromsizes = op.join(config['working_dir'], 'data', 'chrom.sizes'),
        # bam = op.join(config['working_dir'], 'align_{modality}', '{sample}', 'Aligned.sortedByCoord.out.bam'),
        bams = lambda wildcards: [op.join(config['working_dir'], 'align_{modality}', '{sample}', x) for x in list_by_chr_dedup_bams(wildcards)]
    output:
        unsorted = temp(op.join(config['working_dir'], 'align_{modality}', '{sample}',
                                'cb_ub_filt_unsorted.bam')),
        merged = op.join(config['working_dir'], 'align_{modality}', '{sample}', 'cb_ub_filt.bam')
    threads:
        config['nthreads']
    shell:
        """
        samtools merge -@ {threads} -o {output.unsorted} {input.bams}
        samtools sort -@ {threads} {output.unsorted} -o {output.merged}

        """

rule create_deduped_coverage_tracks_all_filtered_in_cbs:
    conda:
        "envs/all_in_one.yaml"
    input:
        bam = op.join(config['working_dir'], 'align_{modality}', '{sample}', 'cb_ub_filt.bam'),
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
        cb_ub_bam = op.join(config['working_dir'], 'align_{modality}', '{sample}', 'cb_ub_filt_twice.bam'),
        bw = op.join(config['working_dir'], 'align_{modality}', '{sample}', '{sample}_{modality}_coverage.bw'),
        cb_ub_bg = temp(op.join(config['working_dir'], 'align_{modality}', '{sample}', 'cb_ub_filt.bw'))
    shell:
        """
        ## this is unrelated to the bamgeneration; fixes starsolo's default permissions
        chmod -R ug+rwX $(dirname {input.bam})

        ## first filter in 'valid' CBs

        ## CB are the error-corrected barcodes
        samtools view -h -@ {threads} {input.bam} -D CB:{input.valid_barcodes} \
           -o {output.cb_ub_bam}
        
        samtools index -@ {threads} {output.cb_ub_bam}

        {params.bedtools} genomecov -ibam {output.cb_ub_bam} \
            -bg -split | LC_COLLATE=C sort -k1,1 -k2,2n > {output.cb_ub_bg}

        ## bedgraph to bigwig
        {params.bedGraphToBigWig} {output.cb_ub_bg} {input.chromsizes} {output.bw}
        
        """
        
# rule create_deduped_coverage_tracks_all_filtered_in_cbs:
#     conda:
#         "envs/all_in_one.yaml"
#     input:
#         bam = op.join(config['working_dir'], 'align_{modality}', '{sample}', 'Aligned.sortedByCoord.out.bam'),
#         valid_barcodes = op.join(config['working_dir'], 'align_wta', '{sample}', 'Solo.out', 'Gene',
#                                  'filtered', 'barcodes.tsv'),
#         chromsizes = op.join(config['working_dir'], 'data', 'chrom.sizes')
#     threads:
#         config['nthreads']
#     params:
#         # bamCoverage = config['bamCoverage'],
#         binSize = 10,
#         bedtools = config['bedtools'],
#         bedGraphToBigWig = config['bedGraphToBigWig']        
#     output:
#         cb_bam = temp(op.join(config['working_dir'], 'align_{modality}', '{sample}', 'cb_filt.bam')),
#         cb_ub_bam = temp(op.join(config['working_dir'], 'align_{modality}', '{sample}', 'cb_ub_filt.bam')),
#         header  = temp(op.join(config['working_dir'], 'align_{modality}', '{sample}', 'header')),
#         cb_ub_bg = temp(op.join(config['working_dir'], 'align_{modality}', '{sample}', 'cb_ub_filt.bw')),
#         bw = op.join(config['working_dir'], 'align_{modality}', '{sample}', '{sample}_{modality}_coverage.bw')
#     shell:
#         """
#         ## this is unrelated to the bamgeneration; fixes starsolo's default permissions
#         chmod -R ug+rwX $(dirname {input.bam})

#         ## first filter in 'valid' CBs

#         ## CB are the error-corrected barcodes
#         samtools view -h -@ {threads} {input.bam} -D CB:{input.valid_barcodes} \
#            -o {output.cb_bam}

#         ## second, deduplicate by UB (error corrected barcodes)
#         ##  so it gets only one alignment per UB:locus:CB combination
#         # 27 CB (error corrected CB)
#         # 28 UB (error corrected UMI)
#         # mind the file is coordinate sorted already - so we deduplicate based on columns 27 and 28
#         samtools view -H {output.cb_bam} > {output.header}

#         samtools view {output.cb_bam} -@ {threads} | uniq -f26 | cat {output.header} - | \
#           samtools view -Sb -@ {threads} > {output.cb_ub_bam}

#         samtools index -@ {threads} {output.cb_ub_bam}

#         {params.bedtools} genomecov -ibam {output.cb_ub_bam} \
#             -bg -split | LC_COLLATE=C sort -k1,1 -k2,2n > {output.cb_ub_bg}

#         ## bedgraph to bigwig
#         {params.bedGraphToBigWig} {output.cb_ub_bg} {input.chromsizes} {output.bw}
        
#         """

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

## yes the log is considered an output - to pass as a flag
## R_LIBS are conda's if run in conda, but /home/rock/R_LIBs if run in docker, and user's if run directly
rule install_r_deps:
    conda:
        "envs/all_in_one.yaml"
    input:
        script = op.join(config['rock_method_path'], 'src', 'installs.R')
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
        tso_fc = op.join(config['working_dir'], 'multimodal', '{sample}', 'tso_featurecounted'),
        wta_fc = op.join(config['working_dir'], 'multimodal', '{sample}', 'wta_featurecounted'),
        gtf = config['gtf'],
        script = op.join(config['rock_method_path'], 'src', 'generate_sce_object.R'),
        installs = op.join(config['working_dir'], 'log', 'installs.log'),
        subset_gtf = op.join(config['working_dir'], 'multimodal', 'subset.gtf')
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
             --output_fn {output.sce} \
             --captured_gtf {input.subset_gtf}
        """


rule render_descriptive_report:
    conda:
        "envs/all_in_one.yaml"
    input:
        mapping_report = op.join(config['working_dir'], 'multimodal', 'mapping_summary.txt'),
        gtf = config['gtf'],
        script = op.join(config['rock_method_path'], 'src', 'generate_descriptive_singlecell_report.Rmd'),
        sces = expand(op.join(config['working_dir'], 'multimodal', '{sample}', '{sample}_sce.rds'),
               sample = get_sample_names()),
        installs = op.join(config['working_dir'], 'log', 'installs.log')
    output:
        html = op.join(config['working_dir'], 'multimodal', 'descriptive_report.html'),
        # cache = temp(op.join(config['rock_method_path'], 'process_sce_objects_cache')),
        # cached_files = temp(op.join(config['rock_method_path'], 'process_sce_objects_files'))
    log: op.join(config['working_dir'], 'multimodal', 'descriptive_report.log')
    params:
        multimodal_path = op.join(config['working_dir'], 'multimodal'),
        run_mode = config['run_mode'],
        working_dir = config['working_dir'],
        sample = "{wildcards.sample}",
        Rbin = config['Rbin'],
        simulate = config['simulate']
    shell:
        """
        simulate={params.simulate}

        if [ "$simulate" = "False" ]
        then

        {params.Rbin} --vanilla -e 'rmarkdown::render(\"{input.script}\", 
          output_file = \"{output.html}\", 
          params = list(multimodal_path = \"{params.multimodal_path}\", 
                        run_mode = \"{params.run_mode}\"))' &> {log}
        else
          echo "no report - that just just a simulation; but SCE objects are ready" > {output.html}
          # touch [output.cache]          
          # touch [output.cached_files]
        fi

        """

rule generate_mapping_report:
    conda:
        "envs/all_in_one.yaml"
    input:
        gtf = config['gtf'],
        script = op.join(config['rock_method_path'], 'src', 'generate_mapping_report.R'),
        sces = expand(op.join(config['working_dir'], 'multimodal', '{sample}', '{sample}_sce.rds'),
                      sample = get_sample_names()),
        installs = op.join(config['working_dir'], 'log', 'installs.log')
    output:
        summary = op.join(config['working_dir'], 'multimodal', 'mapping_summary.txt')
    log: op.join(config['working_dir'], 'multimodal', 'mapping_report.log')
    params:
        path = config['working_dir'],
        Rbin = config['Rbin'],
        simulate = config['simulate']
    shell:
        """
        simulate={params.simulate}

        if [ "$simulate" = "False" ]
        then

        {params.Rbin} --no-echo --no-restore --file={input.script} \
           --args --path {params.path} > {output.summary}
        else
            echo "no report - that just just a simulation; but SCE objects are ready" > {output.summary}
        fi

        """
