import time
from multiprocessing import Pool
import functools
from tqdm import tqdm

from utilities.format_time import get_format_time
from utilities.tail_n_unique_cb import get_tail_n_uniue_cb


def seq_logo_multiprocessing(chromosomes, input_bam_file_paths, unique_cb_bamfile_list, linker_pattern, umi_len, iterable_dict, nProcessors):
    complete_tail_list = []
    complete_unique_cb_list = []
    start_time_mult = time.time()
    for chr_name in chromosomes:
        for i in range(len(input_bam_file_paths)):
            print('\n\n', chr_name, input_bam_file_paths[i])
            try:
                tail_list = []
                unique_cb_list = []
                pool = Pool(processes=nProcessors)
                tail_list,  unique_cb_list = zip(
                                                    *pool.map(functools.partial(
                                                    get_tail_n_uniue_cb, 
                                                    chr_name = chr_name,
                                                    input_bam_path = input_bam_file_paths[i],
                                                    unique_cb_bamfile_list = unique_cb_bamfile_list,
                                                    linker_pattern = linker_pattern,
                                                    umi_len =umi_len,
                                                    tail_list = tail_list,
                                                    unique_cb_list = unique_cb_list
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
            tail_list_pool = [e for elem in tail_list for e in elem]
            unique_cb_list_pool = [e for elem in unique_cb_list for e in elem]

            complete_tail_list.append(tail_list_pool)
            complete_unique_cb_list.append(unique_cb_list_pool)

    end_time_mult = time.time()
    seq_logo_mult_time = get_format_time(start_time_mult, end_time_mult, 'Seq Logo Multiprocessing')

    # Final merger of list and sets for all input files and chromosomes
    final_tail_list = [e for elem in complete_tail_list for e in elem]
    final_complete_unique_cb_list = [e for elem in complete_unique_cb_list for e in elem]

    return final_tail_list, final_complete_unique_cb_list, seq_logo_mult_time

