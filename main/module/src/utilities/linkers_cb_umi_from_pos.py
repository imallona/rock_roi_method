def get_linkers_cb_umi_from_pos(sS_string, start_pos, cb_linker_umi_len_list):
    cb1_len, cb2_len, cb3_len, linker_left_len, linker_right_len, umi_len = cb_linker_umi_len_list
    start_pos_1 = start_pos
    end_pos_1 = start_pos_1+cb1_len
    start_pos_2 = end_pos_1+linker_left_len
    end_pos_2 = start_pos_2+cb2_len
    start_pos_3 = end_pos_2+linker_right_len
    end_pos_3 = start_pos_3+cb3_len
    new_cb_tag = sS_string[start_pos_1:end_pos_1] + '_' + sS_string[start_pos_2:end_pos_2] + '_' + sS_string[start_pos_3:end_pos_3]
    left_linker = sS_string[end_pos_1:start_pos_2]
    right_linker = sS_string[end_pos_2:start_pos_3]
    new_umi_tag = sS_string[end_pos_3:end_pos_3+umi_len]
    return new_cb_tag, left_linker, right_linker, new_umi_tag



