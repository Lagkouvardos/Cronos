import os
import io

# get files in current directory
curent_directory = os.getcwd()
files_of_current_directory = [f for f in os.listdir(curent_directory) if os.path.isfile(os.path.join(curent_directory, f))]

sample_names = []

# create list of sample names
for filename in files_of_current_directory:

    splits = filename.split('@')

    if len(splits) != 2:
        continue

    sample_name = splits[0]

    if sample_name in sample_names:
        continue

    sample_names.append(sample_name)


create_file = True

new_sample_names = []
for sample_name in sample_names:
	new_sample_name = sample_name.split(".")[0].zfill(3) + sample_name.split(".")[1].zfill(2)
	new_sample_names.append(new_sample_name)

f = None
f_full_path_name = os.path.join(curent_directory, "sample-to-file-mappings.tab")

for new_sample_name in new_sample_names.sort():

    r_filename = new_sample_name + '@R.fastq'
    f_filename = new_sample_name + '@F.fastq'

    if r_filename in files_of_current_directory:
        if f_filename in files_of_current_directory:
            if create_file:
                f = open(f_full_path_name, "w")
                create_file = False
            f.write('{}\t{}\t{}\n'.format(new_sample_name, f_filename, r_filename))

if f is not None:
    f.close()








