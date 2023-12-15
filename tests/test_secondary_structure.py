import unittest
import numpy as np
import sys
import os

current_path = os.path.dirname(os.path.abspath(__file__))
relative_path = os.path.normpath(os.path.join(current_path, '../DNA-Probe-Designer'))
sys.path.append(relative_path)
from secondary_structure import filter_secondary_structure

# test secondary structure check
class TestSecondaryStructure(unittest.TestCase):
    # we should pass in a filtered probes file
    def setUp(self):
        self.valid_bed_filename = "tests/files/valid_filtered.bed"
        self.filter_MFE = 0

    # remove what we create
    def tearDown(self):
        new_file = f"{self.valid_bed_filename}_pDup_MFE_filtered.bed"
        if os.path.exists(new_file):
            os.remove(new_file)
    
    def test_filter_secondary_structure_valid(self):
        filter_secondary_structure(self.valid_bed_filename, self.filter_MFE)

    def test_filter_secondary_structure_invalid_MFE(self):
        invalid_MFE = "not a float"
        with self.assertRaises(TypeError):
            filter_secondary_structure(self.valid_bed_filename, invalid_MFE)

if __name__ == '__main__':
    unittest.main()
