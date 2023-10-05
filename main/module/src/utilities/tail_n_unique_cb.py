import pysam
import re
from tqdm import tqdm


def get_tail_n_uniue_cb(iterable, chr_name, input_bam_path, unique_cb_bamfile_list, linker_pattern, umi_len, tail_list, unique_cb_list):
    start_pos = iterable[0]
    end_pos = iterable[1]
    input_bam = pysam.AlignmentFile(input_bam_path, "rb")
    alignments = input_bam.fetch(chr_name, start = start_pos, end = end_pos, multiple_iterators=True)
    for alignment in tqdm(alignments):
        sS_string = alignment.get_tag('sS')
        cb_tag = alignment.get_tag('CB')
        if cb_tag != '-':
            match = re.search(linker_pattern, sS_string)
            if match: 
                # case1: linker pattern matched in sS_string and CB tag available
                end_pos_linker = match.end()
                tail = sS_string[end_pos_linker + 9 + umi_len:] # 9(random) + 8(UMI/UB tag)
                tail_list.append(tail)
            if input_bam_path in unique_cb_bamfile_list:
                chromosome = input_bam.get_reference_name(alignment.reference_id)
                unique_cb_list.append([cb_tag, chromosome])
    input_bam.close()
    return tail_list, unique_cb_list

