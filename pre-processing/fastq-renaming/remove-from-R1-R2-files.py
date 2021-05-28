import os
import io

curent_directory = os.getcwd()
files_of_current_directory = [f for f in os.listdir(curent_directory) if os.path.isfile(os.path.join(curent_directory, f))]

for file_of_current_directory in files_of_current_directory:
    if file_of_current_directory.endswith('.fastq'):
        
        full_path_name = os.path.join(curent_directory, file_of_current_directory)

        splits = file_of_current_directory.split('.')
        if len(splits) != 2:
            continue

        name = splits[0]
        extension = splits[1]

        parts = name.split('_')
        if len(parts) != 2:
            continue
        
        sampleAndTime = parts[0] # 86-9 or 105-8
        sequenceType = parts[1] # R1 or R2

        sample = sampleAndTime.split('-')[0]
        time = sampleAndTime.split('-')[1]

        sample = sample.zfill(3)
        time = time.zfill(2)

        newFileName =  sample + time + '_' + sequenceType + '.' + extension

        os.rename(full_path_name, full_path_name.replace(file_of_current_directory, newFileName))