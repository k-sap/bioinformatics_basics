from NW_alignment_algorithm import Needleman_Wunsch

import unittest

class Needleman_Wunsch_Tests(unittest.TestCase):
    def check_score(self, algo, correct, seq_a, seq_b):
        algo.set_input_sequences(seq_a, seq_b)
        self.assertEqual(algo.compute_score(), correct)
        
    def test_score(self):
        algo = Needleman_Wunsch()
        self.check_score(algo, 23, "Radom", "Radoom")
        self.check_score(algo, 9, "Alan", "Dolan")
        
        algo = Needleman_Wunsch(gap=-7, same=6)
        self.check_score(algo, 23, "Radom", "Radoom")
        self.check_score(algo, 6, "Alan", "Dolan")

    def check_matrix_initialization(self, correct_shape, first_row, first_column, seq_a, seq_b):
        algo = Needleman_Wunsch()
        algo.set_input_sequences(seq_a, seq_b)
        algo.compute_algorithm_matrix()
        self.assertEqual(algo.mt.shape, correct_shape)
        self.assertCountEqual(algo.mt[0,:], first_row)
        self.assertCountEqual(algo.mt[:,0], first_column)
        
    def test_matrix_initialization(self):
        self.check_matrix_initialization((6,7), [0, -2, -4, -6, -8, -10, -12], [0, -2, -4, -6, -8, -10], "Radom", "Radoom")
        self.check_matrix_initialization((5,6), [0, -2, -4, -6, -8, -10], [0, -2, -4, -6, -8], "Alan", "Dolan")
        
    def check_status_init(self, algo):
        self.assertEqual(algo.is_data_initialized, False)
        self.assertEqual(algo.is_matrix_computed, False)
        self.assertEqual(algo.is_path_computed, False)
        
    def check_status_flow(self, algo, seq_a, seq_b):        
        algo.set_input_sequences(seq_a, seq_b)
        self.assertEqual(algo.is_data_initialized, True)
        self.assertEqual(algo.is_matrix_computed, False)
        self.assertEqual(algo.is_path_computed, False)
        
        algo.compute_algorithm_matrix()
        self.assertEqual(algo.is_matrix_computed, True)
        self.assertEqual(algo.is_data_initialized, True)
        self.assertEqual(algo.is_path_computed, False)
        
        algo.compute_paths()
        self.assertEqual(algo.is_matrix_computed, True)
        self.assertEqual(algo.is_data_initialized, True)
        self.assertEqual(algo.is_path_computed, True)
        
    def test_status_flow(self):
        algo = Needleman_Wunsch()
        self.check_status_flow(algo, "Radom", "Radoom")
        self.check_status_flow(algo, "Alan", "Dolan")
        self.check_status_flow(algo, "Mienio", "Wiecio")

if __name__ == '__main__':
    unittest.main()