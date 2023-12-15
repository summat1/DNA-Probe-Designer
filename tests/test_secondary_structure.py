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
    # we need to pass in a filtered probes file
    def setUp(self):
        self.valid_bed_filename = "tests/files/valid_filtered.bed"
        self.output_file = f"{self.valid_bed_filename}_pDup_MFE_filtered.bed"
        self.filter_MFE = 0

    # remove what we create
    def tearDown(self):
        new_file = f"{self.valid_bed_filename}_pDup_MFE_filtered.bed"
        if os.path.exists(new_file):
            os.remove(new_file)
    
    # parse the output file that we create to make sure it follows
    # certain rules of types and structure
    def test_filter_secondary_structure_valid(self):
        filter_secondary_structure(self.valid_bed_filename, self.filter_MFE)
        with open(self.output_file, 'r') as file:
            next(file)
            for line in file:
                parts = line.split('\t')
                if len(parts) < 3:
                    self.fail("Filtered probe file schema is incorrect.")
                
                # check if first part is an int
                probe_num = parts[0]
                try:
                    probe_num = int(probe_num)
                except ValueError:
                    self.fail(f"The probe number is not an integer: {probe_num}")

                # check if second part is a seq, only GATC
                sequence = parts[1]
                if not set(sequence.strip()).issubset(set("GATC")):
                    self.fail(f"The sequence contains invalid characters: {sequence}")

                # ensure the scores is a list of floats
                scores = parts[2].strip('[]').split(',')
                for score in scores:
                    with self.subTest(element=score):
                        # ran into issue with last element in the list
                        score = score.strip().strip('[]')
                        try:
                            score = float(score)
                        except ValueError:
                            self.fail(f"Score '{score}' is not a float.")

    # pass in a non-float for MFE, should raise a type error
    def test_filter_secondary_structure_invalid_MFE(self):
        invalid_MFE = "not a float"
        with self.assertRaises(TypeError) as error:
            filter_secondary_structure(self.valid_bed_filename, invalid_MFE)
        self.assertEqual(str(error.exception), (f"MFE must be a float, unacceptable value: {invalid_MFE}"))

if __name__ == '__main__':
    unittest.main()
