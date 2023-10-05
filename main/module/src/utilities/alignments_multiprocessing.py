import time
from multiprocessing import Pool
import functools
from tqdm import tqdm


from utilities.process_alignments_chunksize import process_alignments_chunksize
from utilities.format_time import get_format_time

def alignments_multiprocessing(chromosomes, input_bam_file_paths, iterable_dict, rg_tag_names, cb_linker_umi_len_list, unique_cb_df, canonical_seq, linker_pattern_list, nProcessors, log_file):
    complete_headers_list = []
    complete_seen_alignments_list = []
    complete_alignment_list = []
    complete_aln_count_list = []
    complete_new_cb_linker_umi_list = []
    complete_count_recovered_cb_list = []
    start_time_mult = time.time()
    for chr_name in chromosomes:
        for i in range(len(input_bam_file_paths)):
            # print(f'\n\n################## Start: bamfile: {input_bam_file_paths[i]} - chr: {chr_name} ##################')
            print('\n\n', chr_name, input_bam_file_paths[i])
            try:
                f = open(log_file, 'a+')
                f.write(f'\n\nInput bam file: {input_bam_file_paths[i]}, Chr: {chr_name}')
                f.write(f'\nstart pos, end pos, total_alignments, count_no_cb_tag, count_unique_tag, count_new_cb_linker_umi, count_none_new_cb_linker_umi, count_recovered_cb')
                f.close()
                new_header_list = []
                seen_alignments_set = set()
                alignment_list = []
                aln_count_list = []
                new_cb_linker_umi_list = []
                count_recovered_cb_list = []
                pool = Pool(processes=nProcessors)
                seen_alignments_set, new_header_list, alignment_list, aln_count_list, new_cb_linker_umi_list, count_recovered_cb_list = zip(
                                                                                                                                            *pool.map(functools.partial(
                                                                                                                                            process_alignments_chunksize, 
                                                                                                                                            chr_name = chr_name,
                                                                                                                                            input_bam_path = input_bam_file_paths[i],
                                                                                                                                            rg_tag_name = rg_tag_names[i],
                                                                                                                                            seen_alignments_set = seen_alignments_set,
                                                                                                                                            new_header_list = new_header_list,
                                                                                                                                            alignment_list = alignment_list,
                                                                                                                                            new_cb_linker_umi_list = new_cb_linker_umi_list,
                                                                                                                                            cb_linker_umi_len_list = cb_linker_umi_len_list,
                                                                                                                                            unique_cb_df = unique_cb_df,
                                                                                                                                            canonical_seq = canonical_seq,
                                                                                                                                            linker_pattern_list = linker_pattern_list,
                                                                                                                                            log_file = log_file,
                                                                                                                                            aln_count_list = aln_count_list,
                                                                                                                                            count_recovered_cb_list = count_recovered_cb_list
                                                                                                                                        ),
                                                                                                                                            tqdm([iterable for iterable in iterable_dict[chr_name]]))
                                                                                                                                        )
                pool.close()
                pool.join()
            except KeyboardInterrupt:
                print("Caught KeyboardInterrupt, terminating workers")
                pool.terminate()
            finally:
                pool.terminate()

            # Merge list and sets for all processes of all chunks of one input file
            header_list_pool = [e for elem in new_header_list for e in elem]
            seen_alignments_pool = set([e for elem in seen_alignments_set for e in elem])
            alignments_pool = [e for elem in alignment_list for e in elem]
            count_list = [e for elem in aln_count_list for e in elem]
            new_cb_linker_umi_list_pool = [e for elem in new_cb_linker_umi_list for e in elem]
            count_recovered_list = [e for elem in count_recovered_cb_list for e in elem]
            
            count = sum(count_list)
            count_recovered = sum(count_recovered_list)
            
            f = open(log_file, 'a+')
            f.write(f'\nTotal count of alignments: {count}')
            f.write(f'\nTotal count of recovered alignments: {count_recovered}')
            f.close()

            complete_headers_list.append(header_list_pool)
            complete_seen_alignments_list.append(seen_alignments_pool)
            complete_alignment_list.append(alignments_pool)
            complete_aln_count_list.append(count)
            complete_new_cb_linker_umi_list.append(new_cb_linker_umi_list_pool)
            complete_count_recovered_cb_list.append(count_recovered)

            # print(f'################## End: bamfile: {input_bam_file_paths[i]} - chr: {chr_name} ##################\n\n')

    end_time_mult = time.time()
    mult_time = get_format_time(start_time_mult, end_time_mult, 'Multiprocessing')

    # Final merger of list and sets for all input files and chromosomes
    prefinal_headers_list = [e for elem in complete_headers_list for e in elem]
    final_seen_alignments_set = set([e for elem in complete_seen_alignments_list for e in elem])
    final_alignments_list = [e for elem in complete_alignment_list for e in elem]
    final_complete_new_cb_linker_umi_list = [e for elem in complete_new_cb_linker_umi_list for e in elem]

    total_count = sum(complete_aln_count_list)
    total_count_recovered = sum(complete_count_recovered_cb_list)
    f = open(log_file, 'a+')
    f.write(f'\n\nIn all input bam files for all chromosomes of interest:')
    f.write(f'\n{total_count} total count of alignments')
    f.write(f'\n{total_count_recovered} total count of recovered alignments')
    f.close()

    return prefinal_headers_list, final_seen_alignments_set, final_alignments_list, final_complete_new_cb_linker_umi_list, total_count, total_count_recovered, mult_time



