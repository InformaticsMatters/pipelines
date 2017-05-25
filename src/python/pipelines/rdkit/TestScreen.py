import sys
import unittest

import screen


class TestScreen(unittest.TestCase):

    def test_upper(self):
        self.assertEqual('foo'.upper(), 'FOO')

    def test_from_file(self):
        sys.argv.append('--input')
        sys.argv.append('data/dhfr_3d.sdf.gz')
        sys.argv.append('--smiles')
        sys.argv.append('C1N=C(C2=CC=CC=C2)C2=CC=CC=C2C2=C1C=NC(NC1=CC=CC=C1)=N2')
        screen.main()

if __name__ == '__main__':
    unittest.main()