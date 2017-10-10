import unittest
from context import esri

class TestInit(unittest.TestCase):

    def test_init_size(self):
        # test default size of new AsciiGrid
        grid = esri.AsciiGrid()
        self.assertEqual(grid.ncols, 10)
        self.assertEqual(grid.nrows, 10)

if __name__ == '__main__':
    unittest.main()