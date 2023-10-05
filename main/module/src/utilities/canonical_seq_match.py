import re
from utilities.linkers_cb_umi_from_pos import get_linkers_cb_umi_from_pos


def match_canonical_seq(sS_string, canonical_seq, cb_linker_umi_len_list):
    seq_len = sum(cb_linker_umi_len_list)
    match_pattern =  r"[A-Z]{" + str(seq_len) + "}" + re.escape(canonical_seq)
    match = re.search(match_pattern, sS_string)
    if match:
        start_pos=0
        new_cb_tag, left_linker, right_linker, new_umi_tag = get_linkers_cb_umi_from_pos(sS_string, start_pos, cb_linker_umi_len_list)
        return new_cb_tag, left_linker, right_linker, new_umi_tag
    else: 
        return None, None, None, None
    

def get_canonical_seq_match(sS, canonical_seq, unique_cb_chr_df, cb_linker_umi_len_list):
    new_cb_tag, left_linker, right_linker, new_umi_tag = match_canonical_seq(sS, canonical_seq, cb_linker_umi_len_list)
    if new_cb_tag == None:
        return None
    else:
        if new_cb_tag in unique_cb_chr_df['CB_tag'].values: 
            matched_cb = 1
            # chrom = unique_cb_chr_df[unique_cb_chr_df['CB_tag']==new_cb_tag]['chromosome'].values
            chrom = unique_cb_chr_df.loc[unique_cb_chr_df['CB_tag'] == new_cb_tag, 'chromosome'].values
        else: 
            matched_cb = 0
            chrom = 'no_match'
        new_cb_linker_umi_dict = {
                                    'new_cb_tag' : new_cb_tag,
                                    'left_linker' : left_linker,
                                    'right_linker' : right_linker,
                                    'new_umi_tag' : new_umi_tag,
                                    'matched_cb' : matched_cb,
                                    'matched_chr' : chrom
                                }
        return new_cb_linker_umi_dict

