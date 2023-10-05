import pandas as pd

def get_bam_to_df(bamfile, save_df_path, chr=None):
    rows = []
    alignments_list = []
    alignment = 0
    for alignment in bamfile.fetch(chr):
        row = {
            "tag": alignment.get_tags(),
            "query_name": alignment.query_name,
            "query_length": alignment.query_length,
            "query_sequence": alignment.query_sequence,
            "query_qualities": alignment.query_qualities,
            "reference_name": alignment.reference_name,
            "reference_id": alignment.reference_id,
            "reference_start": alignment.reference_start,
            "mapping_quality": alignment.mapping_quality,
            "pos": alignment.pos,
            "flag": alignment.flag,
            "rname": alignment.rname,
            "mapq": alignment.mapq,
            "rnext": alignment.rnext,
            "pnext": alignment.pnext,
            "tlen": alignment.tlen,
            'cigar': alignment.cigarstring,
        }
        rows.append(row)
        alignments_list.append(alignment)
    df = pd.DataFrame(rows)
    df.to_csv(save_df_path, index=False)
    return

