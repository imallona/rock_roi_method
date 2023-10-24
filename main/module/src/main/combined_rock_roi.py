import pysam
import pandas as pd
import pathlib
import time
import os
import sys
from datetime import datetime
import re
import yaml
import warnings
warnings.filterwarnings('ignore')

from utilities.reference_length import get_reference_length
from utilities.iterable import get_iterable
from utilities.format_time import get_format_time
from utilities.bam_to_df import get_bam_to_df
from utilities.seq_logo import get_seq_logo
from utilities.alignments_multiprocessing import alignments_multiprocessing
from utilities.seq_logo_multiprocessing import seq_logo_multiprocessing

import argparse

parser = argparse.ArgumentParser(description='')
parser.add_argument('--config', metavar='config', type=str,
                    help='full path to the config file')
cmdargs = parser.parse_args()

def main():
    try:            
        start_time = time.time()

        ######################### Start: Setting up input and output files/parameters #########################
        
        print('Setting up input and output parameters...')

        config_path = cmdargs.config #'./utilities/config_input_parameters.yml'
        with open(config_path) as cf_file:
            args = yaml.safe_load(cf_file.read())

        input_bam_file_paths = args['input_bam_file_paths']
        unique_cb_bamfile_list = args['unique_cb_bamfile_list']
        gtf_file = args['gtf_file']
        featureCounts_path = args['featureCounts_path']
        linker_pattern_for_seq_logo = args['linker_pattern_for_seq_logo']
        linker_pattern_list = args['linker_pattern_list']
        seq_logo_filename = args['seq_logo_filename']
        cb1_len = args['cb1_len']
        cb2_len = args['cb2_len']
        cb3_len = args['cb3_len']
        linker_left_len = args['linker_left_len']
        linker_right_len = args['linker_right_len']
        umi_len = args['umi_len']
        chunk_genomic_area = args['chunk_genomic_area']
        rg_tag_names = args['rg_tag_names']
        chromosomes = args['chromosomes']
        subset_gtf = args['subset_gtf']
        output_folder = args['output_folder']
        subset_gtf_output_file = output_folder + args['subset_gtf_output_file']
        subset_gtf_pattern = args['subset_gtf_pattern']
        write_final_bam_to_csv = args['write_final_bam_to_csv']
        write_final_bam_header_to_txt = args['write_final_bam_header_to_txt']
        nProcessors = args['nProcessors']
        nthreads = args['nthreads']
        final_merged_file = output_folder + args['final_merged_file']
        featureCounts_output_file = output_folder + args['featureCounts_output_file']
        log_file = output_folder + args['log_file']
        error_log = output_folder + args['error_log']

        pathlib.Path(output_folder).mkdir(parents=True, exist_ok=True)
        print('Done!\n')

        # if rg_tag_names is not provided then get rg_tag_names from input bam file names and add '_' to its elements
        print('Setting up rg tag names...')
        if rg_tag_names is None:
            rg_tag_names = []
            for input_bam_file_path in input_bam_file_paths:
                rg_tag_names.append(os.path.basename(input_bam_file_path).split('.')[0] + '_')
        else:
            rg_tag_names = [rg_tag + '_' for rg_tag in rg_tag_names]
        print('Done!\n')

        # sort and index input bam files if not already sorted and indexed
        for file in input_bam_file_paths:
            bam = pysam.AlignmentFile(input_bam_file_paths[0], "rb")
            header = bam.header
            if header['HD']['SO'] != 'coordinate':
                print('Sorting and indexing bam file: {}...'.format(file))
                pysam.sort("-o", file, file)
                pysam.index(file)
            elif not os.path.exists(file + '.bai'):
                print('Indexing bam file: {}...'.format(file))
                pysam.index(file)
            else:
                print('Bam file {} is already sorted and indexed'.format(file))
            bam.close()
        print('Done!\n')
        
        # check if chromosomes is provided else get all chromosomes from first input bam file header
        # open first input bam file
        print('Setting up chromosome names to be processed...')
        input_bam = pysam.AlignmentFile(input_bam_file_paths[0], "rb")
        if chromosomes is None:
            header = input_bam.header
            chromosomes = [ref['SN'] for ref in header['SQ']]
        print('Done!\n')

        ######################### End: setting up input and output files/parameters #########################

        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        f = open(log_file, 'w')
        f.write('\n\nTimestamp: {}'.format(timestamp))
        f.write('\nInput bam files:\n{}'.format('\n'.join(input_bam_file_paths)))
        f.close()
        
        # get genomic length of chromosome from first input bam file header
        print('Getting genomic length of chromosome...')
        chr_length_dict = {}

        ## patch to retrievel chr_names when inputted as `None` within the `chromosomes` yaml field
        
        # print(chromosomes)
        # print(input_bam.header["SQ"])
        chromosomes = [x['SN'] for x in input_bam.header["SQ"]]
        ## we also subset for those within the 'chromosomes' yaml file
        if args['chromosomes'] is not None:
            chromosomes = [x for x in chromosomes if x in args['chromosomes']]
        print(chromosomes)
        ## patch end
        
        for chr_name in chromosomes:
            chr_length_dict[chr_name] = get_reference_length(input_bam.header["SQ"], chr_name)
        # print(f"\nchromosomes: chr_length\n{chr_length_dict}\n")
        print('Done!\n')

        header_dict = input_bam.header.to_dict()
        input_bam.close()
        
        # Get genomic area as iterable for multiprocessing for each chromosome
        print('Getting genomic area as iterable for multiprocessing for each chromosome...')
        iterable_dict = {}
        
        for chr_name in chromosomes:
            end_pos = chr_length_dict[chr_name]
            start_pos = 0
            iterable = get_iterable(chunk_genomic_area, end_pos)
            iterable_dict[chr_name] = iterable
            # print(f"\n{chr_name}, length: {end_pos}\nstart and end pos for each process:\n{iterable}")
        print('Done!\n')

        # get sequence logo of the tail
        print('Generating sequence logo of the tail of reads...')
        unique_cb_df = pd.DataFrame(columns=['CB_tag', 'chromosome'])
        try:
            tail_list, unique_cb_list, seq_logo_mult_time = seq_logo_multiprocessing(chromosomes, input_bam_file_paths, unique_cb_bamfile_list, linker_pattern_for_seq_logo, umi_len, iterable_dict, nProcessors)
            canonical_seq = get_seq_logo(tail_list, 'Tail Seq Logo', output_folder + seq_logo_filename)
            df = pd.DataFrame(unique_cb_list, columns=['CB_tag', 'chromosome'])
            unique_cb_df = df.drop_duplicates(subset=['CB_tag', 'chromosome'])
            unique_cb_df.to_csv(output_folder+'unique_cb_df.csv', index = False)
        except:
            print('Error in generating canonical_seq. Setting canonical_seq as empty string.')
            canonical_seq = ''
        f = open(log_file, 'a+')
        f.write(f'\ncanonical_seq: {canonical_seq}\n')
        f.close()
        print('canonical_seq: ', canonical_seq)
        print('Done!\n')

        cb_linker_umi_len_list = [cb1_len, cb2_len, cb3_len, linker_left_len, linker_right_len, umi_len]

        # process alignments using multiprocessing
        print('Processing alignments using multiprocessing...')
        prefinal_headers_list, final_seen_alignments_set, final_alignments_list, final_complete_new_cb_linker_umi_list, total_count, total_count_recovered, mult_time = alignments_multiprocessing(
                                                                                                                                                                                chromosomes, 
                                                                                                                                                                                input_bam_file_paths, 
                                                                                                                                                                                iterable_dict, 
                                                                                                                                                                                rg_tag_names, 
                                                                                                                                                                                cb_linker_umi_len_list, 
                                                                                                                                                                                unique_cb_df, 
                                                                                                                                                                                canonical_seq, 
                                                                                                                                                                                linker_pattern_list, 
                                                                                                                                                                                nProcessors, 
                                                                                                                                                                                log_file
                                                                                                                                                                                )
        print('Done!\n')

        # make dataframe for final_complete_new_cb_linker_umi_list
        df_new_cb_linker_umi = pd.DataFrame(final_complete_new_cb_linker_umi_list)
        df_new_cb_linker_umi.to_csv(output_folder+'new_cb_linker_umi.csv', index=False)

        # make the header unique
        def make_hashable(d): return frozenset(d.items())
        final_headers_list = list(map(dict, set(map(make_hashable, prefinal_headers_list))))

        # write final bam file
        print('Final Bam file...')
        merge_count = 0
        header_dict['RG'] = final_headers_list
        new_header = pysam.AlignmentHeader.from_dict(dict(header_dict))
        merged_bam = pysam.AlignmentFile(final_merged_file, "wb", header=new_header)
        for alignment in final_alignments_list:
            aln = pysam.AlignedSegment.fromstring(alignment, header=merged_bam.header)
            unique_tag = aln.get_tag('RG') + '_' + aln.get_tag('UB')
            if unique_tag in final_seen_alignments_set:
                merge_count += 1
                merged_bam.write(aln)
                final_seen_alignments_set.remove(unique_tag)
        merged_bam.close()
        pysam.sort("-o", final_merged_file, final_merged_file)
        pysam.index(final_merged_file)
        print('- Done!')

        f = open(log_file, 'a+')
        f.write(f'\n{merge_count} total count of unique RG+CB+UB tags')
        f.write(f'\n{round((merge_count/total_count)*100, 2)}% of total alignments are unique RG+CB+UB tags')
        f.write(f'\n{round((total_count_recovered/total_count)*100, 2)}% of total alignments are recovered alignments')
        f.close()

        # write final bam file as df and its header
        if write_final_bam_to_csv == 1:
            print('\nWrite final bam file as csv...')
            final_bam = pysam.AlignmentFile(final_merged_file, "rb")
            get_bam_to_df(final_bam, output_folder + 'final_bam_df.csv')
            final_bam.close()
            print('- Done!')

        if write_final_bam_header_to_txt == 1:
            print('\nWrite final bam file header...')
            final_bam = pysam.AlignmentFile(final_merged_file, "rb")
            header = str(final_bam.header)
            with open(output_folder + 'final_bam_header.txt', "w") as outfile: outfile.write(header)
            final_bam.close()
            print('- Done!')

        # gtf file: get subset of gtf corresponding to chromosomes of interest
        if subset_gtf == 1:
            if subset_gtf_pattern != None:
                print('\nGtf file: get only subset of chromosomes...')
                with open(gtf_file, 'r') as input_file:
                    with open(subset_gtf_output_file, 'w') as output_file:
                        matching_lines = filter(lambda line: re.search(subset_gtf_pattern, line), input_file)
                        output_file.writelines(matching_lines)
                print('- Done!')
            else:
                subset_gtf_output_file = gtf_file
                print('\nsubset_gtf_pattern not provided, utilizing complete gtf file for featureCounts.')    
        else: 
            subset_gtf_output_file = gtf_file
            print('\nUtilizing complete gtf file for featureCounts.')

        # featureCounts
        print('\nFeatureCounts...')
        featurecounts_cmd = f"{featureCounts_path} \
                                -a {subset_gtf_output_file} \
                                -o {featureCounts_output_file} \
                                {final_merged_file} \
                                -F GTF \
                                -t exon \
                                -g gene_id \
                                -f \
                                -O \
                                -M  \
                                -T {nthreads} \
                                --byReadGroup"
        os.system(featurecounts_cmd)
        print('- Done!\n')

        end_time = time.time()
        hh, mm, ss = get_format_time(start_time, end_time, 'Finished complete pipeline: ')

        hh_seq_mult, mm_seq_mult, ss_seq_mult = seq_logo_mult_time
        hh_mult, mm_mult, ss_mult = mult_time

        f = open(log_file, 'a+')
        f.write("\n\nTime taken for seq logo multiprocessing: {:0>2} hours {:0>2} min {:05.2f} sec".format(hh_seq_mult, mm_seq_mult, ss_seq_mult))
        f.write("\nTime taken for alignments multiprocessing: {:0>2} hours {:0>2} min {:05.2f} sec".format(hh_mult, mm_mult, ss_mult))
        f.write("\nTime taken for complete pipeline: {:0>2} hours {:0>2} min {:05.2f} sec\n\n".format(hh, mm, ss))
        f.close()

    except Exception as e:
        print('Check log file for error: {}'.format(e))
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        f = open(error_log, 'a+')
        f.write('\n\nTimestamp:{}\nFile:{}\nLine:{} Error:{}\n\n'.format(timestamp, fname, exc_tb.tb_lineno, e))
        # print(exc_type, fname, exc_tb.tb_lineno)
        f.close()
        sys.exit(1)



if __name__ == '__main__':
    main()


