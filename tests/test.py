import unittest
import numpy as np
import sys
import os

current_path = os.path.dirname(os.path.abspath(__file__))
relative_path = os.path.normpath(os.path.join(current_path, '../'))
sys.path.append(relative_path)
from duplex_prob import calc_duplex_prob
from secondary_structure import filter_secondary_structure

class TestDuplexProb(unittest.TestCase):
    def setUp(self):
        self.valid_sam_filename = "tests/valid_test.sam"
        self.valid_bed_filename = "tests/valid_test.bed"
        self.temp = 37

    # if it executes successfully we should return an array
    def test_calc_duplex_prob_valid(self):
        result = calc_duplex_prob(self.valid_sam_filename, self.valid_bed_filename, self.temp)
        self.assertIsInstance(result, np.ndarray)

    # test with invalid temperature input
    def test_calc_duplex_prob_invalid_temp(self):
        invalid_temp = 100
        with self.assertRaises(ValueError):
            calc_duplex_prob(self.valid_sam_filename, self.valid_bed_filename, invalid_temp)

class TestSecondaryStructure(unittest.TestCase):
    def setUp(self):
        self.valid_bed_filename = "tests/valid_filtered.bed"
        self.filter_MFE = 0

    def test_filter_secondary_structure_valid(self):
        filter_secondary_structure(self.valid_bed_filename, self.filter_MFE)

    def test_filter_secondary_structure_invalid_MFE(self):
        invalid_MFE = "not a float"
        with self.assertRaises(TypeError):
            filter_secondary_structure(self.valid_bed_filename, invalid_MFE)

if __name__ == '__main__':
    unittest.main()
