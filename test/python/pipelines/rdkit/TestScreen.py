import unittest

#from .context import pipelines

class TestScreen(unittest.TestCase):

    def test_upper(self):
        self.assertEqual('foo'.upper(), 'FOO')

#    def test_from_file(self):
#        args = {'input':'data/dhfr_3d.sdf.gz', 'simmin': 0.5, 'smiles': 'C1N=C(C2=CC=CC=C2)C2=CC=CC=C2C2=C1C=NC(NC1=CC=CC=C1)=N2'}
#        count = screen.execute(args)
#        self.assertEqual(count, 1)

if __name__ == '__main__':
    unittest.main()