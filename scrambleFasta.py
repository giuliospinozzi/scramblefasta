# -*- coding: utf-8 -*-
"""
Created on Thu Aug  7 17:18:38 2014

@author: stefano
"""

import os
import sys
import random
from numpy import ceil

import scrambleModule



################ FILES, PATHS AND NAMES ####################################################################################################


# File filter per extension
fasta_files_extension = 'fa'
                        # 'fa' or whatever, it determines the files on which this program will run, by means of the extension

#INPUT
input_path = None
             # If a folder, the program works with all *.fasta_files_extension in the folder.
             # If None, the program gets cwd as folder.
             # If a file, the program works on that file only.

# OUTPUT             
output_path = None
             # If None, the program gets cwd as folder.
             # Else, it must be a valid folder path: if the folder doesn't exist, will be created.

# ID STRING for out files labelling             
ID_for_scrambled_files = 'scrambled'
                         # ID string to discern scrambled files from original ones
                         # If None or an empty string, check path variable... ORIGINAL FILES MAY BE OVERWRITTEN!
                         
                         
                         
################ SETTINGS #################################################################################################################


#-------- Configure deletions --------#

# del_type_extractor FUNCTION
# IN: available_del_types tuple; OUT: one del type
# Example: implementation of random del type --> del_type_extractor = lambda available_del_types: random.choice(available_del_types)
# Actual implementation: static choice = 'Nbp'
del_type_extractor = lambda available_del_types: 'Nbp'

# del_calls_per_sequence FUNCTION
# IN: len_input_string, ratio between total nucleotides and deletion to perform; OUT: n of deletions per sequence (n of doDeletion calls)
# Actual implementation: random n of deletions per sequence (between 1 and del-to-do-computed-through-one_each_n_neclotides-var)
del_calls_per_sequence = lambda len_input_string, one_each_n_neclotides: random.randint(1,int(ceil(len_input_string/float(one_each_n_neclotides))))

deletions = {'do': True,
             'available_del_types': ('from_last', 'from_first', 'Nbp'),
             'del_type_extractor': del_type_extractor,
             'del_calls_per_sequence': del_calls_per_sequence,
             'one_each_n_neclotides': 300, # support variable 'del_calls_per_sequence' function
             'doDeletion': scrambleModule.doDeletion,
             'del_args': scrambleModule.getArgs # Go there to see how arguments are passed to scrambleModule.doDeletion
             }


#-------- Configure insertions --------#

# ins_type_extractor FUNCTION
# IN: available_ins_type tuple; OUT: one ins type
# Example: implementation of random ins type --> ins_type_extractor = lambda available_ins_type: random.choice(available_ins_type)
# Actual implementation: static choice = 'Nbp'
ins_type_extractor = lambda available_ins_type: 'Nbp'

# ins_calls_per_sequence FUNCTION
# IN: len_input_string, ratio between total nucleotides and insertions to perform; OUT: n of insertions per sequence (n of doInsertion calls)
# Actual implementation: random n of insertions per sequence (between 1 and ins-to-do-computed-through-one_each_n_neclotides-var)
ins_calls_per_sequence = lambda len_input_string, one_each_n_neclotides: random.randint(1,int(ceil(len_input_string/float(one_each_n_neclotides))))
          
insertions = {'do': True,
              'available_ins_types': ('from_last', 'from_first', 'Nbp'),
              'ins_type_extractor': ins_type_extractor,
              'ins_calls_per_sequence': ins_calls_per_sequence,
              'one_each_n_neclotides': 300, # support variable for 'ins_calls_per_sequence'
              'doInsertion': scrambleModule.doInsertion,
              'ins_args': scrambleModule.getArgs # Go there to see how arguments are passed to scrambleModule.getArgs
              }


#-------- Configure mutations --------#

# mut_type_extractor FUNCTION
# IN: available_mut_type tuple; OUT: one mut type
# Example: implementation of random mut type --> mut_type_extractor = lambda available_mut_type: random.choice(available_mut_type)
# Actual implementation: static choice = 'Nbp'
mut_type_extractor = lambda available_mut_type: 'Nbp'

# mut_calls_per_sequence FUNCTION
# IN: len_input_string, ratio between total nucleotides and mutations to perform; OUT: n of mutations per sequence (n of doMutation calls)
# Actual implementation: random n of mutations per sequence (between 1 and mut-to-do-computed-through-one_each_n_neclotides-var)
mut_calls_per_sequence = lambda len_input_string, one_each_n_neclotides: random.randint(1,int(ceil(len_input_string/float(one_each_n_neclotides))))
             
mutations = {'do': True,
             'available_mut_types': ('Nbp',),
             'mut_type_extractor': mut_type_extractor,
             'mut_calls_per_sequence': mut_calls_per_sequence,
             'one_each_n_neclotides': 300, # support variable for 'mut_calls_per_sequence'
             'doMutation': scrambleModule.doMutation,
             'mut_args': scrambleModule.getArgs # Go there to see how arguments are passed to scrambleModule.getArgs
             }



################ CODE #####################################################################################################################


print "\n [SCRAMBLE SEQUENCES IN FASTA FILES]"


# Get fasta file list (complete paths) looking at input_path --> fasta_files list
fasta_files = []
if input_path is None: # default, search for fasta files in cwd and append to fasta_files
    fasta_files = [os.path.normpath(os.path.join(os.getcwd(),filename)) for filename in os.listdir(os.getcwd()) if (os.path.isfile(os.path.normpath(os.path.join(os.getcwd(),filename))) and filename.split('.')[-1] == fasta_files_extension)]    
else:
    try:
        input_path = os.path.normpath(input_path)
    except:
        sys.exit("\n [ERROR] '{input_path}' is not formatted as a valid input path!\n\n [QUIT]\n\n".format(input_path=input_path))
    if os.path.isfile(input_path): # check if input_path is a *.fa file and append to fasta_files
        if input_path.split('.')[-1] == fasta_files_extension:
            fasta_files.append(input_path)
        else:
            sys.exit("\n [ERROR] *.fa file is expected!\n\n [QUIT]\n\n")        
    else: # search for fasta files in input_path and append to fasta_files
        if os.path.exists(input_path):
            fasta_files = [os.path.normpath(os.path.join(input_path,filename)) for filename in os.listdir(input_path) if (os.path.isfile(os.path.normpath(os.path.join(input_path,filename))) and filename.split('.')[-1] == fasta_files_extension)]
        else:
            sys.exit("\n [ERROR] '{wd}' is not a valid path or you don't have sufficient privileges to access!\n\n [QUIT]\n\n".format(wd=input_path))
# Quit if fasta_files list is empty      
if len(fasta_files) < 1:
    if input_path is None:
        input_path = os.getcwd()
    sys.exit("\n [ERROR] Can't find any fasta file in '{wd}'!\n\n [QUIT]\n\n".format(wd=input_path))
    

# Check output_path and create it if doesn't exist
if (output_path is None or output_path == ''):
    output_path = os.getcwd()
else:
    try:
        output_path = os.path.normpath(output_path)
    except:
        sys.exit("\n [ERROR] '{output_path}' is not formatted as valid output path!\n\n [QUIT]\n\n".format(output_path=output_path))
    if os.path.isfile(output_path):
        sys.exit("\n [ERROR] Output path must be a folder path, not a file! Your input: '{output_path}'\n\n [QUIT]\n\n".format(output_path=output_path))
    if not os.path.exists(output_path):
        os.makedirs(output_path)
# Prepare output file names through ID_for_scrambled_files
if ID_for_scrambled_files is None:
    ID_for_scrambled_files = ''
else:
    ID_for_scrambled_files = "_"+ID_for_scrambled_files


# ---> Put here a launch config summary: cwd, fasta files found, action to perform...


# For each fasta file
for fasta_filepath in fasta_files:
    
    # Starting message
    sys.stdout.write("\n * Processing {} ... ".format(os.path.basename(fasta_filepath)))
    sys.stdout.flush()
    
    # Read file and scramble sequences, then appended to scrambled_file_lines
    scrambled_file_lines = []
    with open(fasta_filepath, 'r') as source_file:
        
        # Load file as a list --> source_file_lines
        source_file_lines = source_file.readlines()
        # Get file features / quick check
        n_of_seq = len(source_file_lines)
        if n_of_seq == 0:
            print " [WARNING] this file turned out to be empty!"
            print "   [SKIP]"
            continue
        elif n_of_seq%2 != 0:
            print " [WARNING] this file has an odd number of lines! Check it!"
            print "   [SKIP]"
            continue
        else:
            n_of_seq = int(n_of_seq/2.0)
        
        # Loop over source file lines --> fill scrambled_file_lines
        for source_file_line in source_file_lines:
            if '>' in source_file_line:
                # Here you can act on sequence header
                header = source_file_line
                scrambled_file_lines.append(header)
            else:
                # Here you can act on sequence string
                sequence = source_file_line.rstrip()
                len_seq = len(sequence)                
                # deletions
                if deletions['do'] is True:
                    for i in range(0, deletions['del_calls_per_sequence'](len_seq, deletions['one_each_n_neclotides'])):
                        del_type = deletions['del_type_extractor'](deletions['available_del_types'])
                        del_args = deletions['del_args']('doDeletion', sequence, del_type)
                        sequence = deletions['doDeletion'](*del_args)
                # insertions
                if insertions['do'] is True:
                    for i in range(0, insertions['ins_calls_per_sequence'](len_seq, insertions['one_each_n_neclotides'])):
                        ins_type = insertions['ins_type_extractor'](insertions['available_ins_types'])
                        ins_args = insertions['ins_args']('doInsertion', sequence, ins_type)
                        sequence = insertions['doInsertion'](*ins_args)
                # mutations
                if mutations['do'] is True:
                    for i in range(0, mutations['mut_calls_per_sequence'](len_seq, mutations['one_each_n_neclotides'])):
                        mut_type = mutations['mut_type_extractor'](mutations['available_mut_types'])
                        mut_args = mutations['mut_args']('doMutation', sequence, mut_type)
                        sequence = mutations['doMutation'](*mut_args)
                # Newline
                sequence = sequence + "\n"
                # Final append
                scrambled_file_lines.append(sequence)
    
    # set-up a name for the scrambled version of the file
    scrambled_filename = os.path.basename(fasta_filepath)
    scrambled_filename = ''.join(scrambled_filename.split('.')[:-1] + [ID_for_scrambled_files, "."+fasta_files_extension])
    scrambled_filepath = os.path.join(output_path, scrambled_filename)
        
    # Write the scrambled file
    with open(scrambled_filepath, 'w') as scrambled_file:
        scrambled_file.writelines(scrambled_file_lines)
        
    # Ending message
    sys.stdout.write("done!")
    sys.stdout.flush()
    if ID_for_scrambled_files != '':
        print "\n   outfile name: {scrambled_filename}".format(scrambled_filename=scrambled_filename)


print "\n [FINISHED] \n"