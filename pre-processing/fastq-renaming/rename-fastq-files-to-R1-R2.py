import os
import io

curent_directory = os.getcwd()
files_of_current_directory = [f for f in os.listdir(curent_directory) if os.path.isfile(os.path.join(curent_directory, f))]

for file_of_current_directory in files_of_current_directory:
    if file_of_current_directory.endswith('.fastq'):
        
        full_path_name = os.path.join(curent_directory, file_of_current_directory)

        newFileName = None
        if '@F' in file_of_current_directory:
            newFileName = file_of_current_directory.replace('@F', '_R1')

        if '@R' in file_of_current_directory:
            newFileName = file_of_current_directory.replace('@R', '_R2')

        splits = newFileName.split('.')
        if len(splits) != 3:
            continue

        newFileName = splits[0] + '-' + splits[1] + '.' + splits[2]

        os.rename(full_path_name, full_path_name.replace(file_of_current_directory, newFileName))