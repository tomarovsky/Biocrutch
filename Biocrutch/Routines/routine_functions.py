import bz2
import gzip


def metaopen(filename, flags, buffering=None):
    '''function for open file, file.gz, file.bz2'''
    if not isinstance(filename, str):
        return filename
    elif filename[-3:] == ".gz":
        return gzip.open(filename, flags)
    elif filename[-4:] == ".bz2":
        return bz2.open(filename, flags)
    else:
        if buffering is not None:
            return open(filename, flags, buffering=buffering)
        else:
            return open(filename, flags)