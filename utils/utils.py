from __future__ import print_function
from sys import getsizeof, stderr
from itertools import chain
from collections import deque
from reprlib import repr
import numpy as np
import os
import shutil
from zipfile import ZipFile

# ########################################################################################
def total_size(o, handlers={}, verbose=False):
    """ Returns the approximate memory footprint an object and all of its contents.
        in bytes.

    Automatically finds the contents of the following builtin containers and
    their subclasses:  tuple, list, deque, dict, set and frozenset.
    To search other containers, add handlers to iterate over their contents:

        handlers = {SomeContainerClass: iter,
                    OtherContainerClass: OtherContainerClass.get_elements}
                    
    https://code.activestate.com/recipes/577504/

    """
    dict_handler = lambda d: chain.from_iterable(d.items())
    all_handlers = {tuple: iter,
                    list: iter,
                    deque: iter,
                    dict: dict_handler,
                    set: iter,
                    frozenset: iter,
                   }
    all_handlers.update(handlers)     # user handlers take precedence
    seen = set()                      # track which object id's have already been seen
    default_size = getsizeof(0)       # estimate sizeof object without __sizeof__

    def sizeof(obj):
        if id(obj) in seen:       # do not double count the same object
            return 0
        seen.add(id(obj))
        s = getsizeof(obj, default_size)

        if verbose:
            print(s, type(obj), repr(obj), file=stderr)

        for typ, handler in all_handlers.items():
            if isinstance(obj, typ):
                s += sum(map(sizeof, handler(obj)))
                break
        return s

    return sizeof(o)

# ########################################################################################
def test_total_size():
    """
    Test the total_size function with some data
     
    """

    print ("=========== Dictionary ==============")
    d = dict(a=1, b=2, c=3, d=[4,5,6,7], e='a string of chars')
    # print(total_size(d, verbose=False))
    # print(getsizeof(d))
    print("Size calculated with total size function      : {} bytes".format(total_size(d, verbose=False)))
    print("Size calculated with getsizeof native function: {} bytes".format(getsizeof(d)))

    print ("=========== List ==============")
    l = [10,100, 200, 300]
    print("Size calculated with total size function      : {} bytes".format(total_size(l, verbose=False)))
    print("Size calculated with getsizeof native function: {} bytes".format(getsizeof(l)))

    print ("=========== Set ==============")
    l = (10,100, 200, 300)
    print("Size calculated with total size function      : {} bytes".format(total_size(l, verbose=False)))
    print("Size calculated with getsizeof native function: {} bytes".format(getsizeof(l)))

    print ("=========== File ==============")
    f = open("utils.py", 'r')
    print("Size calculated with total size function      : {} bytes".format(total_size(f, verbose=False)))
    print("Size calculated with getsizeof native function: {} bytes".format(getsizeof(f)))

    print ("=========== Numpy ==============")
    arr = np.zeros((100,100), dtype=float)
    print("Size calculated with total size function      : {} bytes".format(total_size(arr, verbose=False)))
    print("Size calculated with getsizeof native function: {} bytes".format(getsizeof(arr)))

    print ("=========== String ==============")
    s = "El Quijote de la Mancha"
    print("Size calculated with total size function      : {} bytes".format(total_size(s, verbose=False)))
    print("Size calculated with getsizeof native function: {} bytes".format(getsizeof(s)))

# ########################################################################################
def compress_files(dirbase="./", zipname="gaussian_inputs.zip"):

    # Empty path file paths list
    file_paths = []

    # Walking through directory and subidrectories
    for root, directories, files in os.walk(dirbase):
        for filename in files:
            # join the two strings in order to form the full filepath.
            filepath = os.path.join(root, filename)
            file_paths.append(filepath)

    # writing files to a zipfile
    with ZipFile(zipname, 'w') as zip:
        #writing each file
        for file in file_paths:
            zip.write(file)

# ########################################################################################
def delete_folder(folder):

    shutil.rmtree(folder, ignore_errors=False)

# ########################################################################################
def get_name_file_ext(fullpath):

    """
    This function returns the name of the path, the name of the file and the extension
    of an arbitrary file.

    Warning: This function only works in Linux
    TODO: Generalize to work in Windows

    On input:

    | ``fullnametrj``: The full path to the trajectory

    Return:
    | ``filename`` : Name of the file
    | ``path``     : Full Path
    | ``typetrj``  : Type of trajectory
    | ``isfile``   : True if fullpath is a file
    """

    if os.path.isfile(fullpath):

        filename = os.path.basename(fullpath)
        path = os.path.dirname(os.path.realpath(fullpath))
        typetrj = filename.split(".")[-1]
        isfile = True

        return filename, path, typetrj, isfile

    else:

        isfile = False
        return None, None, None, isfile

# ########################################################################################
def padding_list(l, fillval=np.nan):

    """

    Efficiently convert uneven list of lists to minimal containing
    array padded with nan (or other value passed through fillval).

    Solution from:
    https://stackoverflow.com/questions/40569220/efficiently-convert-uneven-list-of-lists-to-minimal-containing-array-padded-with
    :param l: List
    :param fillval: Value to pad the empty places
    :return: a padded np.array
    """

    lens = np.array([len(item) for item in l])
    mask = lens[:,None] > np.arange(lens.max())
    out = np.full(mask.shape,fillval)
    out[mask] = np.concatenate(l)
    return out