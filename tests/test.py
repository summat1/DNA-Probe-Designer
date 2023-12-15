import unittest
import numpy as np
from unittest.mock import patch, mock_open
from ..duplex_prob import calc_duplex_prob
from secondary_structure import filter_secondary_structure

class TestDuplexProb(unittest.TestCase):
    def setUpTestCase(self):
        self.valid_sam_filename = "valid_test.sam"
        self.valid_bed_filename = "valid_test.bed"
        self.temp = 37

    def test_calc_duplex_prob_valid(self):
        result = calc_duplex_prob(self.valid_sam_filename, self.valid_bed_filename, self.temp)

        # Check that the result is an np.ndarray
        self.assertIsInstance(result, np.ndarray)

        # Check the length of the result to ensure it matches expected output size
        # The expected length should be based on the actual content of your .sam and .bed files
        expected_length = 0
        self.assertEqual(len(result), expected_length)


    def test_calc_duplex_prob_invalid_temp(self):
        # Test with invalid temperature input
        invalid_temp = 100
        with self.assertRaises(ValueError):  # Assuming ValueError is raised for invalid temperature
            calc_duplex_prob(self.valid_sam_filename, self.valid_bed_filename, invalid_temp)

    # Additional test cases...

class TestSecondaryStructure(unittest.TestCase):
    def setUp(self):
        # Setup any common data or configurations for the tests
        self.valid_bed_filename = "valid.bed"
        self.filter_MFE = -10.0

    def test_filter_secondary_structure_valid(self):
        # Test with valid input
        mock_bed_content = "mocked BED content"
        with patch('builtins.open', mock_open(read_data=mock_bed_content)) as mock_file:
            filter_secondary_structure(self.valid_bed_filename, self.filter_MFE)
            mock_file.assert_called_with(self.valid_bed_filename)
            # Add assertions to verify the output file content

    def test_filter_secondary_structure_invalid_MFE(self):
        # Test with invalid MFE value
        invalid_MFE = "not a float"
        with self.assertRaises(TypeError):  # Assuming TypeError is raised for invalid MFE
            filter_secondary_structure(self.valid_bed_filename, invalid_MFE)

    # Additional test cases...

if __name__ == '__main__':
    unittest.main()
