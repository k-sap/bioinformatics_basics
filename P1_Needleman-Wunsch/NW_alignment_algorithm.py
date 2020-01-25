import numpy

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
            """
            defines in which way should algorithm look for next path element
            0 - on diagonal
            1 - above
            2 - on left

            returns reccurency state base on this decision as a tuple
            """

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