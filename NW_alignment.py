#!/usr/bin/python3

import numpy as np
import argparse
import json
import re

parser = argparse.ArgumentParser(description='NW alignment. You have to specify two files with sequences in FASTA format. Results are save in ./results directory. Example: python3 NW_alignment.py -a one.txt -b two.txt')
parser.add_argument("-a", "--seq_a", type=str, action='store', help="Please specify the sequence a file", required=True)
parser.add_argument("-b", "--seq_b", type=str, action='store', help="Please specify the sequence b file.", required=True)

args = parser.parse_args()

config_path = 'config.json'

def read_config():
    """
    Loads parameters from a json file from config_path
    """
    with open(config_path) as json_file:
        try:
            data = json.load(json_file)
            gap = data['config']['penalties']['gap']
            same = data['config']['penalties']['same']
            diff = data['config']['penalties']['diff']
            max_seq_length = data['config']['max_seq_length']
            max_number_paths = data['config']['max_number_paths']
            output_dir = data['config']['output_dir']
        except:
            raise Exception("The conflig file is corrupted.")
        return data, gap, same, diff, max_seq_length, max_number_paths, output_dir

    
data, gap, same, diff, max_seq_length, max_number_paths, output_dir = read_config()
    
def read_sequence(file_path):
    """
    return: string with a sequence from a FASTA text file
    """
    with open(file_path, "r") as f:
        f.readline()
        fasta_data = f.read()
        converted_data = re.sub("\n", "", fasta_data)
        return converted_data

def algorithm(seq_a, seq_b):
    def matrix_fill(k, l):
        """
        Fill matrix with scores
        k,l - upper left index
        """
        up_left_val = mt[k][l] + (same if seq_a[k] == seq_b[l] else diff)
        left_val = mt[k+1][l] + gap
        up_val = mt[k][l+1] + gap
        
        values = np.array([up_left_val, left_val, up_val])
        arg_max = np.argmax(values)
        
        mt[k+1][l+1] = values[arg_max]
        
        mt_paths[k+1][l+1] = list(np.transpose((values[arg_max] == values).nonzero()).flatten())
        return
    def generate_paths():
        """
        finds paths/alignments on the base of mt
        """
        def recurrent_path(fun_params):
            nonlocal number_paths
            nonlocal result_file
            seq_align_a = fun_params[0]
            seq_align_b = fun_params[1]
            k = fun_params[2]
            l = fun_params[3]
            def alignment_state(n_possibility = 0):
                #on the base of the arrows append the sequence alignment
#                 nonlocal seq_a
#                 nonlocal seq_b
                if mt_paths[k][l][n_possibility] == 0:
                    return seq_align_a + seq_a[k-1], seq_align_b + seq_b[l-1], k-1, l-1
                elif mt_paths[k][l][n_possibility] == 1:
                    return seq_align_a + "-", seq_align_b + seq_b[l-1], k, l-1
                elif mt_paths[k][l][n_possibility] == 2:
                    return seq_align_a + seq_a[k-1], seq_align_b + "-", k-1, l
            
            #reccurence stop condition
            if k == 0 and l == 0:
                #specific indexation stems from that the matrix has more rows and columns than the length of seq_a, seq_b
                result_file += seq_align_a[len(seq_align_a)-2::-1] + "\n" + seq_align_b[len(seq_align_b)-2::-1] + "\n\n"
                return
            #new paths
            elif len(mt_paths[k][l]) == 1:
                recurrent_path(alignment_state())
            elif number_paths == max_number_paths:
                recurrent_path(alignment_state())
            elif len(mt_paths[k][l]) == 2:
                number_paths += 1
                recurrent_path(alignment_state(0))
                recurrent_path(alignment_state(1))
            elif number_paths + 2 == max_number_paths:
                #mamy trzy i dodajemy trzy
                number_paths += 2
                recurrent_path(alignment_state(0))
                recurrent_path(alignment_state(1))
                recurrent_path(alignment_state(2))
            else:
                #mamy trzy ale mozemy dodac tylko jedna
                number_paths += 1
                recurrent_path(alignment_state(0))
                recurrent_path(alignment_state(1))
        recurrent_path(("", "", m, n))
    
    #sequence change to get coherent indexation in matrix and in sequences
    seq_a = " " + seq_a
    seq_b = " " + seq_b
    
    #variable initialization
    m = len(seq_a)
    n = len(seq_b)
    mt = np.ndarray((m + 1,n + 1), np.int16)
    mt_paths = [[[-1]]*(n+1) for m in range(m+1)]
    result_file = ""
    number_paths = 1
    
    #matrix obvious values
    mt[0][0] = 0
    for i in range(1, m + 1):
        mt[i][0] = gap * i
        mt_paths[i][0] = [2]
    for i in range(1, n + 1):
        mt[0][i] = gap * i
        mt_paths[0][i] = [1]
    
    #matrix non-obvious values
    for i in range(m):
        for j in range(n):
            matrix_fill(i, j)
    #score saved to file
    result_file += str(mt[m][n]) + "\n\n"
    
    generate_paths()
    
    return result_file


seq_a = read_sequence(args.seq_a)
seq_b = read_sequence(args.seq_b)

if(len(seq_a) >= max_seq_length or len(seq_b) >= max_seq_length):
    print("Too long sequence. Max length: ", max_seq_length)
    exit()
else:
    result_file = algorithm(seq_a, seq_b)
    name_a = re.sub("\.[A-z]+", "", args.seq_a)
    name_b = re.sub("\.[A-z]+", "", args.seq_b)
    
    with open(output_dir + "result_" + name_a + "_" + name_b + "_.txt", "w+") as f:
        f.write(result_file)
    print("FINISHED")