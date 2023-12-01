#!/usr/bin/env python
##
## Functions process rock/roi data (general method)
##
## Started 11th Oct 2023
##
## Izaskun Mallona
## GPLv3

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

def get_chromosomes(wildcards):
    # with open(op.join(config['working_dir'], 'data', 'chrom.sizes')) as fh:
    # with open(chromsizes_fn) as fh:
    fn = checkpoints.retrieve_genome_sizes.get(**wildcards).output[0]
    with open(fn) as fh:
        return(list(line.strip().split('\t')[0] for line in fh))

def list_by_chr_dedup_bams(wildcards):
    chroms = get_chromosomes(wildcards)
    return(chrom + '_cb_umi_deduped.bam' for chrom in chroms)

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
