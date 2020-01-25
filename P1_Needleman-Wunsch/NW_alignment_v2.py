#!/usr/bin/python3

import numpy
import argparse
import json
import re

from NW_alignment_algorithm import Needleman_Wunsch

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

seq_a = read_sequence(args.seq_a)
seq_b = read_sequence(args.seq_b)

if(len(seq_a) >= max_seq_length or len(seq_b) >= max_seq_length):
    print("Too long sequence. Max length: ", max_seq_length)
    exit()
else:
    algo = Needleman_Wunsch(gap=gap, same=same, diff=diff, max_seq_length=max_seq_length, max_number_paths=max_number_paths)
    algo.set_input_sequences(seq_a, seq_b)
    result_file = algo.export_results()

    name_a = re.sub("\.[A-z]+", "", args.seq_a)
    name_b = re.sub("\.[A-z]+", "", args.seq_b)
    
    with open(output_dir + "result_" + name_a + "_" + name_b + "_.txt", "w+") as f:
        f.write(result_file)
    print("FINISHED")