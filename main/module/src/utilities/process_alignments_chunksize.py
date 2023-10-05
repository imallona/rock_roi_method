import pysam
import re
from tqdm import tqdm

from utilities.canonical_seq_match import get_canonical_seq_match

def process_alignments_chunksize(iterable, chr_name, input_bam_path, rg_tag_name, seen_alignments_set, new_header_list, alignment_list, new_cb_linker_umi_list, cb_linker_umi_len_list, unique_cb_df, canonical_seq, linker_pattern_list, log_file, aln_count_list, count_recovered_cb_list):
    start_pos = iterable[0]
    end_pos = iterable[1]
    # print(f'Started a process from pool for start and end postions: {iterable} for chr {chr_name}')
    input_bam = pysam.AlignmentFile(input_bam_path, "rb")
    alignments = input_bam.fetch(chr_name, start = start_pos, end = end_pos, multiple_iterators=True)
    total_count = 0
    count_cb_tag = 0
    count_unique_tag = 0
    count_new_cb_linker_umi = 0
    count_none_new_cb_linker_umi = 0
    count_recovered_cb = 0

    for alignment in tqdm(alignments):
        total_count += 1
        cb_tag = alignment.get_tag('CB')
        if '-' not in cb_tag:
            count_cb_tag += 1
            ub_tag = alignment.get_tag('UB')
            new_rg_tag = rg_tag_name + cb_tag
            unique_tag = new_rg_tag + '_' + ub_tag
            if unique_tag not in seen_alignments_set:
                count_unique_tag += 1
                alignment.set_tag('RG', new_rg_tag)
                new_header_list.append({'ID':new_rg_tag})
                seen_alignments_set.add(unique_tag)
                alignment_list.append(alignment.to_string())
        else:
            sS_string = alignment.get_tag('sS')
            for linker_pattern in linker_pattern_list:
                match = re.search(linker_pattern, sS_string)
                if match: break
            if not match: # case: CB tag = '-' and linker patterns do not match
                new_cb_linker_umi_dict = get_canonical_seq_match(sS_string, canonical_seq, unique_cb_df, cb_linker_umi_len_list)
                if new_cb_linker_umi_dict != None:
                    new_cb_linker_umi_list.append(new_cb_linker_umi_dict)
                    count_new_cb_linker_umi += 1
                    if new_cb_linker_umi_dict['matched_cb'] == 1:
                        count_recovered_cb += 1
                        cb_tag = new_cb_linker_umi_dict['new_cb_tag']
                        ub_tag = new_cb_linker_umi_dict['new_umi_tag']
                        alignment.set_tag('CB', cb_tag)
                        alignment.set_tag('UB', ub_tag)
                        count_cb_tag += 1
                        new_rg_tag = rg_tag_name + cb_tag
                        unique_tag = new_rg_tag + '_' + ub_tag
                        if unique_tag not in seen_alignments_set:
                            count_unique_tag += 1
                            alignment.set_tag('RG', new_rg_tag)
                            new_header_list.append({'ID':new_rg_tag})
                            seen_alignments_set.add(unique_tag)
                            alignment_list.append(alignment.to_string())
            else: count_none_new_cb_linker_umi += 1

    count_no_cb_tag = total_count - count_cb_tag
    # print(f'Total alignments: {total_count}, for start and end pos: {iterable} for chr {chr_name}\n')
    f = open(log_file, 'a+')
    f.write(f'\n{iterable[0]}\t{iterable[1]}\t{total_count}\t{count_no_cb_tag}\t{count_unique_tag}\t{count_new_cb_linker_umi}\t{count_none_new_cb_linker_umi}\t{count_recovered_cb}')
    f.close()
    aln_count_list.append(total_count)
    count_recovered_cb_list.append(count_recovered_cb)
    input_bam.close()
    # print(f'\nCompleted process for start and end postions: {iterable} for chr {chr_name}')
    return seen_alignments_set, new_header_list, alignment_list, aln_count_list, new_cb_linker_umi_list, count_recovered_cb_list

