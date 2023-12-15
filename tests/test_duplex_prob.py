import unittest
import numpy as np
import sys
import os

current_path = os.path.dirname(os.path.abspath(__file__))
relative_path = os.path.normpath(os.path.join(current_path, '../DNA-Probe-Designer'))
sys.path.append(relative_path)
from duplex_prob import calc_duplex_prob, filter_duplex_prob

# test the calc_duplex_prob function using test input files
# mostly checking the sam file here since its used only in this function
class TestCalcDuplexProb(unittest.TestCase):
    # valid files should create the correct output
    # invalid files should toss errors
    def setUp(self):
        self.valid_sam_filename = "tests/files/valid_test.sam"
        self.valid_bed_filename = "tests/files/valid_test.bed"

        self.invalid_sam_filename = "tests/files/invalid_test.sam" 

        self.valid_temp = 37
        self.invalid_temp = 100

    # if it executes successfully we should return an array with probs between 0 and 1
    def test_calc_duplex_prob_valid(self):
        dupe_probs = calc_duplex_prob(self.valid_sam_filename, self.valid_bed_filename, self.valid_temp)
        self.assertTrue(all(0 <= p <= 1 for p in dupe_probs))
        self.assertIsInstance(dupe_probs, np.ndarray)

    # test with invalid temperature input
    def test_calc_duplex_prob_invalid_temp(self):
        with self.assertRaises(ValueError):
            calc_duplex_prob(self.valid_sam_filename, self.valid_bed_filename, self.invalid_temp)

    # bad sam file
    def test_calc_duplex_prob_wrong_format_sam(self):
        with self.assertRaises(ValueError):
            calc_duplex_prob(self.invalid_sam_filename, self.valid_bed_filename, self.valid_temp)

        
class TestFilterDuplexProb(unittest.TestCase):
    # setup the file paths, some of which are valid and others invalid
    def setUp(self):
        self.filtered_probe_filename = "tests/files/valid_filtered.bed"
        self.blank = "tests/files/blank.txt"

        self.valid_sam_filename = "tests/files/valid_test.sam"
        self.valid_bed_filename = "tests/files/valid_test.bed"

        self.invalid_sam_filename = "tests/files/invalid_test.sam"
        self.invalid_bed_filename = "tests/files/invalid_test.bed"

        self.valid_temp = 42
        self.valid_filter_prob = 0.2

    # remove the file that we create
    def tearDown(self):
        new_file = f"{self.filtered_probe_filename}_pDup_filtered.bed"
        if os.path.exists(new_file):
            os.remove(new_file)

    # non existent sam file
    def test_missing_sam(self):
        non_existent_sam = "nan.sam"
        with self.assertRaises(FileNotFoundError):
            filter_duplex_prob(non_existent_sam, self.valid_bed_filename, self.valid_temp, self.valid_filter_prob)

    # non existent bed file
    def test_missing_bed(self):
        non_existent_bed = "nan.bed"
        with self.assertRaises(FileNotFoundError):
            filter_duplex_prob(self.valid_sam_filename, non_existent_bed, self.valid_temp, self.valid_filter_prob)
    
    # wrong sam file type
    def test_wrong_sam_file_type(self):
        with self.assertRaises(ValueError):
            filter_duplex_prob(self.blank, self.valid_bed_filename, self.valid_temp, self.valid_filter_prob)

    # wrong bed file type
    def test_wrong_bed_file_type(self):
        with self.assertRaises(ValueError):
            filter_duplex_prob(self.valid_sam_filename, self.blank, self.valid_temp, self.valid_filter_prob)

    # bad bed file
    def test_calc_duplex_prob_wrong_format_bed(self):
        with self.assertRaises(ValueError):
            filter_duplex_prob(self.valid_sam_filename, self.invalid_bed_filename, self.valid_temp, self.valid_filter_prob)
        
    # test our output
    def test_filtered_probe_file(self):
        with open(self.filtered_probe_filename, 'r') as file:
            lines = file.readlines()

        header_parts = lines[0].split(' ')
        self.assertTrue(header_parts[0].isnumeric())
        number_of_probes = int(header_parts[0])

        # check each line's format
        for line in lines[1:]:
            probe_info = line.strip().split('\t')
            probe_info = [info.strip() for info in probe_info]

            # check that each probe came with 3 pieces of info
            self.assertEqual(len(probe_info), 3) 
            probe_number, sequence, probabilities_str = probe_info

            # make sure each part of the probe info is the correct type
            self.assertTrue(probe_number.isnumeric())
            self.assertIsInstance(sequence, str)
            probabilities = eval(probabilities_str)  
            self.assertIsInstance(probabilities, list)
            self.assertTrue(all(isinstance(p, float) for p in probabilities))

        # check total number of probes
        self.assertEqual(len(lines) - 1, number_of_probes)  

if __name__ == '__main__':
    unittest.main()
