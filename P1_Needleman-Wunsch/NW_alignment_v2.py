#!/usr/bin/python3

import numpy
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

class Needleman_Wunsch:
    """
    Neddleman_Wunsch algorithm
    """
    #----------------------------------------------------
    # algorithm properties
    gap: int
    same: int
    diff: int
    max_seq_length: int
    max_number_paths: int
    #----------------------------------------------------
    # algorithm state - matrix
    is_data_initilized: bool
    is_matrix_computed: bool
    seq_a: str
    seq_b: str
    mt: numpy.ndarray
    m: int
    n: int
    #----------------------------------------------------
    # algorithm state - matrix
    is_path_computed: bool
    mt_paths: numpy.ndarray
    number_paths: int
    paths: str
    #----------------------------------------------------
    
    def __init__(self, gap=-2, same=5, diff=-5, max_seq_length=100, max_number_paths=5):
        self.gap = gap
        self.same = same
        self.diff = diff
        self.max_seq_length = max_seq_length
        self.max_number_paths = max_number_paths
        
        self.is_data_initialized = False
        self.is_matrix_computed = False
        self.is_path_computed = False
        
    
    def set_input_sequences(self, seq_a, seq_b):
        self.seq_a = ' ' + seq_a
        self.seq_b = ' ' + seq_b
        self.m = len(self.seq_a)
        self.n = len(self.seq_b)

        
        self.is_data_initialized = True
        self.is_matrix_computed = False
        self.is_path_computed = False
        
    def compute_algorithm_matrix(self):
        if not self.is_data_initialized:
            print("There is no input sequences. Use set_input_sequences method first.")
            return None
        if self.is_matrix_computed:
            return self.mt
        
        #matrix initialization
        self.__init_input_matrix__()
        
        #matrix obvious values
        self.mt[0][0] = 0
        for i in range(1, self.m):
            self.mt[i][0] = self.gap * i
            self.mt_paths[i][0] = [2]
        for i in range(1, self.n):
            self.mt[0][i] = self.gap * i
            self.mt_paths[0][i] = [1]
    
        #matrix non-obvious values
        for i in range(self.m-1):
            for j in range(self.n-1):
                self.__matrix_element_fill__(i, j)
                
        self.is_matrix_computed = True
        
        return self.mt
                
    def compute_score(self):
        if not self.is_matrix_computed:
            self.compute_algorithm_matrix()
        return self.mt[self.m-1][self.n-1]
    
    def compute_paths(self):
        if not self.is_data_initialized:
            print("There is no input sequences. Use set_input_sequences method first.")
            return None
        if self.is_path_computed:
            return self.paths
        if not self.is_matrix_computed:
            self.compute_algorithm_matrix()
        
        self.__recurrent_path__(("", "", self.m-1, self.n-1))
        self.is_path_computed = True
        print(self.paths)
    
    def export_results(self):
        if self.is_data_initialized is True:
            self.compute_paths()
            return str(self.compute_score()) + "\n\n" + self.paths
        else:
            print("There is no input sequences. Use set_input_sequences method first.")
        
    def __init_input_matrix__(self):
        self.mt = numpy.ndarray((self.m, self.n), numpy.int16)
        self.mt_paths = [[[-1]]*(self.n) for i in range(self.m)]
        self.number_paths = 1
        self.paths = ""
        
    def __matrix_element_fill__(self, k, l):
        """
        Fill one element of the lagorithm matrix with score
        k,l - upper left index
        """
        up_left_val = self.mt[k][l] + (self.same if self.seq_a[k+1] == self.seq_b[l+1] else self.diff)
        left_val = self.mt[k+1][l] + self.gap
        up_val = self.mt[k][l+1] + self.gap

        values = numpy.array([up_left_val, left_val, up_val])
        arg_max = numpy.argmax(values)

        self.mt[k+1][l+1] = values[arg_max]

        self.mt_paths[k+1][l+1] = list(numpy.transpose((values[arg_max] == values).nonzero()).flatten())
        
    def __recurrent_path__(self, fun_params):
        seq_align_a = fun_params[0]
        seq_align_b = fun_params[1]
        k = fun_params[2]
        l = fun_params[3]
        
        def alignment_state(n_possibility = 0):
            #on the base of the arrows append the sequence alignment

            if self.mt_paths[k][l][n_possibility] == 0:
                return seq_align_a + self.seq_a[k], seq_align_b + self.seq_b[l], k-1, l-1
            elif self.mt_paths[k][l][n_possibility] == 1:
                return seq_align_a + "-", seq_align_b + self.seq_b[l], k, l-1
            elif self.mt_paths[k][l][n_possibility] == 2:
                return seq_align_a + self.seq_a[k], seq_align_b + "-", k-1, l

        #reccurence stop condition
        if k == 0 and l == 0:
            #specific indexation stems from that the matrix has more rows and columns than the length of seq_a, seq_b
            self.paths += seq_align_a[len(seq_align_a)-1::-1] + "\n" + seq_align_b[len(seq_align_b)-1::-1] + "\n\n"
            return
        #new paths
        elif len(self.mt_paths[k][l]) == 1:
            self.__recurrent_path__(alignment_state())
        elif self.number_paths == self.max_number_paths:
            self.__recurrent_path__(alignment_state())
        elif len(self.mt_paths[k][l]) == 2:
            self.number_paths += 1
            self.__recurrent_path__(alignment_state(0))
            self.__recurrent_path__(alignment_state(1))
        elif self.number_paths + 2 == self.max_number_paths:
            #mamy trzy i dodajemy trzy
            self.number_paths += 2
            self.__recurrent_path__(alignment_state(0))
            self.__recurrent_path__(alignment_state(1))
            self.__recurrent_path__(alignment_state(2))
        else:
            #mamy trzy ale mozemy dodac tylko jedna
            self.number_paths += 1
            self.__recurrent_path__(alignment_state(0))
            self.__recurrent_path__(alignment_state(1))

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