def get_iterable(chunk_genomic_area, end_pos, start_pos=0):
    next_pos = chunk_genomic_area
    iterable = []
    while start_pos!=end_pos:
            iterable.append((start_pos, next_pos))
            start_pos = next_pos
            next_pos = next_pos + chunk_genomic_area
            if next_pos > end_pos: next_pos = end_pos
    return iterable
