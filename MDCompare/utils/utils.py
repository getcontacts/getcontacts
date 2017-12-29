from __future__ import division
import os
import errno


def append_to_filename(filename, some_string):
    before, after = os.path.splitext(filename)
    return "%s_%s%s" % (before, some_string, after)
    # return '/'.join(filename.split('/')[:-1] + ["generic_" + filename.split('/')[-1]])


def clean_path(path):
    if path[-1] != '/':
        path += '/'
    return path


def get_file_descriptor(full_filename):
    return os.path.splitext(os.path.basename(full_filename))[0]


def open_dir(dir_name):
    try:
        os.makedirs(dir_name)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


def listdir_files(path):
    all_files = [path + filename for filename in os.listdir(path) if os.path.isfile(path + filename)]
    for filex in all_files:
        if os.path.basename(filex)[0] == '.':
            all_files.remove(filex)
    return all_files


def listdir_dirs(path):
    return [clean_path(path + dirname) for dirname in os.listdir(path) if os.path.isdir(path + dirname)]


def correct_interaction_list(int_list):
    if "hb" in int_list:
        int_list.remove("hb")
        int_list.append("hbss")
        int_list.append("hbsb")
        int_list.append("hbbb")
        # int_list.append("vdw2")
        int_list.append("wb")
        int_list.append("wb2")
    # if "vdw" in int_list:

    return int_list
