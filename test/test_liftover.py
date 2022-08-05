import unittest

from igv_compress.liftover import Chain, load_liftover
from igv_compress.feature import parse, Feature
import pathlib


class LiftoverTest(unittest.TestCase):


    def test_chain(self):
        '''
        chain 1000 chr1 190687455 + 4784516 190687455 chr1 1650847 + 0 1650847 1
        1725 14679 0
        3914 2468 0
        1535
        '''

        alignments = [(1725, 14679, 0), (3914, 2468, 0), (1535,)]

        chain = Chain("chr1", 190687455, 4784516, 190687455, "chr1", 1650847, 0, 1650847, 1)
        chain.set_alignments(alignments)

        # Feature overlapping end of block
        boundary = 4784516 + 1725
        f = Feature("chr1", boundary - 10, boundary + 10)
        mapped = chain.map(f)
        self.assertEqual(1, len(mapped))
        self.assertEqual(1725 - 10, mapped[0].start)
        self.assertEqual(10, mapped[0].size())

        # Feature within a block
        f = Feature("chr1", 4784516, 4784516 + 10)
        mapped = chain.map(f)
        self.assertEqual(1, len(mapped))
        self.assertEqual(0, mapped[0].start)
        self.assertEqual(10, mapped[0].size())

        # Feature overlapping beginning of block
        f = Feature("chr1", 4784516 - 10, 4784516 + 10)
        mapped = chain.map(f)
        self.assertEqual(1, len(mapped))
        self.assertEqual(0, mapped[0].start)
        self.assertEqual(10, mapped[0].size())

        # Feature overlapping 2 blocks
        boundary = 4784516 + 1725 + 14679
        f = Feature("chr1", 4784516 - 10, boundary + 10)
        mapped = chain.map(f)
        self.assertEqual(2, len(mapped))

        mapped.sort(key = lambda f: f.start)
        self.assertEqual(0, mapped[0].start)
        self.assertEqual(1725, mapped[0].size())

        self.assertEqual(1725, mapped[1].start)
        self.assertEqual(10, mapped[1].size())


        f = Feature("chr1", 10, 20)
        mapped = chain.map(f)
        self.assertEqual(0, len(mapped))

    def test_load(self):

        bedfile = str((pathlib.Path(__file__).parent / "data/K27Ac.chain").resolve())
        liftover = load_liftover(bedfile)

        # Feature overlapping end of block
        boundary = 4784516 + 1725
        f = Feature("chr1", boundary - 10, boundary + 10)
        mapped = liftover.map(f)
        self.assertEqual(1, len(mapped))
        self.assertEqual(1725 - 10, mapped[0].start)
        self.assertEqual(10, mapped[0].size())

        # Feature within a block
        f = Feature("chr1", 4784516, 4784516 + 10)
        mapped = liftover.map(f)
        self.assertEqual(1, len(mapped))
        self.assertEqual(0, mapped[0].start)
        self.assertEqual(10, mapped[0].size())

        # Feature overlapping beginning of block
        f = Feature("chr1", 4784516 - 10, 4784516 + 10)
        mapped = liftover.map(f)
        self.assertEqual(1, len(mapped))
        self.assertEqual(0, mapped[0].start)
        self.assertEqual(10, mapped[0].size())

        # Feature overlapping 2 blocks
        boundary = 4784516 + 1725 + 14679
        f = Feature("chr1", 4784516 - 10, boundary + 10)
        mapped = liftover.map(f)
        self.assertEqual(2, len(mapped))

        mapped.sort(key = lambda f: f.start)
        self.assertEqual(0, mapped[0].start)
        self.assertEqual(1725, mapped[0].size())

        self.assertEqual(1725, mapped[1].start)
        self.assertEqual(10, mapped[1].size())


        f = Feature("chr1", 10, 20)
        mapped = liftover.map(f)
        self.assertEqual(0, len(mapped))
