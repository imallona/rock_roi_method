def get_reference_length(SQ_header_list, chr_name):
    for header_line in SQ_header_list:
        if header_line["SN"] == chr_name:
            chr_length = header_line["LN"]
            break
    print(f"Chromosome {chr_name} has length {chr_length}")
    return chr_length
