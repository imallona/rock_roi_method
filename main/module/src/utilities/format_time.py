def get_format_time(start, end, name):
    total_time = end - start
    hours, rem = divmod(total_time, 3600)
    minutes, seconds = divmod(rem, 60)
    print("\n{} total time taken:\n{:0>2} hours {:0>2} min {:05.2f} sec".format(name, int(hours),int(minutes),seconds))
    return int(hours),int(minutes),seconds
