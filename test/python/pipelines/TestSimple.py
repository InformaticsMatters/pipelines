import unittest

#from .context import pipelines

class TestSimple(unittest.TestCase):

    def test_upper(self):
        self.assertEqual('foo'.upper(), 'FOO')

#    def test_simple(self):
#        sys.argv.append('--name')
#        sys.argv.append('Tim')
#        result = pipelines.simple.main()
#        self.assertTrue(result)


if __name__ == '__main__':
    unittest.main()