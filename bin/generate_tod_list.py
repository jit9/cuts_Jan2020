#!/usr/bin/env python

"""This script takes a pickle file and generate three lists of
TODs corresponding to train, validate and test for developing
machine learning pipeline.
"""
import os
import cPickle
import random

#########################
# define run parameters #
#########################

pickle_file = "/mnt/act3/users/lmaurin/work/pickle_cuts/mr3_pa3_s16_results.pickle"
output_dir = "data"
train_fname = "mr3_pa3_s16_train.txt"
validate_fname = "mr3_pa3_s16_validate.txt"
test_fname = "mr3_pa3_s16_test.txt"

n_train = 60
n_validate = 20
n_test = 20

#########
# main  #
#########

if not os.path.exists(output_dir):
    print("Folder %s doesn't exist, creating..." % output_dir)
    os.makedirs(output_dir)

# load pickle file
print("Loading pickle file: %s" % pickle_file)
with open(pickle_file, "r") as f:
    data = cPickle.load(f)

# random select required number of tods
total_tods = n_train + n_validate + n_test
tod_list = random.sample(data['name'], total_tods)

# split into train, validate and test
train_list = tod_list[:n_train]
validate_list = tod_list[n_train:n_train+n_validate]
test_list = tod_list[n_train+n_validate:]

# write list to file
def write_to_file(outfile, lst):
    print("Writing to %s" % outfile)
    with open(outfile, "w") as f:
        for l in lst:
            f.write("%s\n" % l)

# output filenames
outfile_train = os.path.join(output_dir, train_fname)
outfile_validate = os.path.join(output_dir, validate_fname)
outfile_test = os.path.join(output_dir, test_fname)

# write to files
write_to_file(outfile_train, train_list)
write_to_file(outfile_validate, validate_list)
write_to_file(outfile_test, test_list)

print("Done!")

