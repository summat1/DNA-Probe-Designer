import unittest
from unittest.mock import patch, mock_open
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

current_path = os.path.dirname(os.path.abspath(__file__))
relative_path = os.path.normpath(os.path.join(current_path, '../DNA-Probe-Designer'))
sys.path.append(relative_path)
from duplex_prob import calc_duplex_prob, filter_duplex_prob, plot_duplex_prob

# test the calc_duplex_prob function using test input files
# mostly checking the sam file here since its used only in this function
class TestCalcDuplexProb(unittest.TestCase):
    # valid files should create the correct output
    # invalid files should toss errors
    def setUp(self):
        self.valid_sam_filename = "tests/files/valid_test.sam"
        self.valid_bed_filename = "tests/files/valid_test.bed"

        self.invalid_sam_filename = "tests/files/invalid_test.sam" 

        self.temps = np.array([32, 37, 42, 47, 52, 57])
        self.valid_temp = 37
        self.invalid_temp = 100

    # if it executes successfully we should return an array with probs between 0 and 1
    def test_calc_duplex_prob_valid(self):
        dupe_probs = calc_duplex_prob(self.valid_sam_filename, self.valid_bed_filename, self.valid_temp)
        self.assertTrue(all(0 <= p <= 1 for p in dupe_probs))
        self.assertIsInstance(dupe_probs, np.ndarray)

    # test with invalid temperature input ensuring correct error message
    def test_calc_duplex_prob_invalid_temp(self):
        with self.assertRaises(ValueError) as error:
            calc_duplex_prob(self.valid_sam_filename, self.valid_bed_filename, self.invalid_temp)
        self.assertEqual(str(error.exception), f"Invalid temperature value: {self.invalid_temp}. Valid values are {self.temps}")

    # bad sam file
    def test_calc_duplex_prob_wrong_format_sam(self):
        with self.assertRaises(ValueError) as error:
            calc_duplex_prob(self.invalid_sam_filename, self.valid_bed_filename, self.valid_temp)
        self.assertEqual(str(error.exception), "Probe sequence must contain only G, A, T, C in SAM file.")

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
        with self.assertRaises(FileNotFoundError) as error:
            filter_duplex_prob(non_existent_sam, self.valid_bed_filename, self.valid_temp, self.valid_filter_prob)
        self.assertEqual(str(error.exception), f"The file {non_existent_sam} does not exist.")

    # non existent bed file
    def test_missing_bed(self):
        non_existent_bed = "nan.bed"
        with self.assertRaises(FileNotFoundError) as error:
            filter_duplex_prob(self.valid_sam_filename, non_existent_bed, self.valid_temp, self.valid_filter_prob)
        self.assertEqual(str(error.exception), f"The file {non_existent_bed} does not exist.")

    # wrong sam file type
    def test_wrong_sam_file_type(self):
        with self.assertRaises(ValueError) as error:
            filter_duplex_prob(self.blank, self.valid_bed_filename, self.valid_temp, self.valid_filter_prob)
        self.assertEqual(str(error.exception), f"The file {self.blank} is not a SAM file.")

    # wrong bed file type
    def test_wrong_bed_file_type(self):
        with self.assertRaises(ValueError) as error:
            filter_duplex_prob(self.valid_sam_filename, self.blank, self.valid_temp, self.valid_filter_prob)
        self.assertEqual(str(error.exception), f"The file {self.blank} is not a BED file.")

    # bad bed file
    def test_calc_duplex_prob_wrong_format_bed(self):
        with self.assertRaises(ValueError) as error:
            filter_duplex_prob(self.valid_sam_filename, self.invalid_bed_filename, self.valid_temp, self.valid_filter_prob)
        self.assertEqual(str(error.exception), "Probe sequence must contain only G, A, T, C in BED file.")
        
    # test our output
    def test_filtered_probe_file(self):
        with open(self.filtered_probe_filename, 'r') as file:
            lines = file.readlines()

        header_parts = lines[0].split(' ')
        self.assertTrue(header_parts[0].isnumeric())
        number_of_probes = int(header_parts[0])

        # check our filtered file format
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

# test the plotting function
class testPlotDuplexProb(unittest.TestCase):
    # provide a filtered test file input
    def setUp(self):
        self.filtered_file = "tests/files/valid_filtered.bed"
    
    # ensure that the figures are plotting and showing with 'all' param
    @patch('matplotlib.pyplot.show')
    @patch('matplotlib.pyplot.plot')
    @patch('matplotlib.pyplot.figure')
    def test_plot_duplex_prob_all(self, mock_figure, mock_plot, mock_show):
        plot_duplex_prob(self.filtered_file, 'all')
        mock_figure.assert_called()
        mock_plot.assert_called()
        mock_show.assert_called()

    # ensure that the figures are plotting and showing when we select one probe
    @patch('matplotlib.pyplot.show')
    @patch('matplotlib.pyplot.plot')
    @patch('matplotlib.pyplot.figure')
    def test_plot_duplex_prob_single(self, mock_figure, mock_plot, mock_show):
        plot_duplex_prob(self.filtered_file, 1)
        mock_figure.assert_called()
        mock_show.assert_called()

    # check for value error with invalid param
    def test_plot_duplex_prob_invalid_input(self):
        with self.assertRaises(ValueError):
            plot_duplex_prob(self.filtered_file, 'invalid')

if __name__ == '__main__':
    unittest.main()
