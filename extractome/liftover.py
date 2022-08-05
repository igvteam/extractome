from intervaltree import IntervalTree
from extractome.stream import getstream


def load_liftover(chains_file):
    chains = []
    chain = None
    data = []

    with  getstream(chains_file) as f:
        for line in f:
            line = line.strip()
            if len(line) == 0:
                # end of chain
                if chain is not None:
                    chain.set_alignments(data)
                    chains.append(chain)
                    chain = None
                    data = []
            elif line.startswith("chain"):
                t = line.split()
                # chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
                chain = Chain(t[2], int(t[3]), int(t[5]), int(t[6]), t[7], int(t[8]), int(t[10]), int(t[11]), t[12])
            else:
                t = line.split()
                for i in range(len(t)):
                    t[i] = int(t[i])
                data.append(t)

    # Last data
    if chain is not None and len(data) > 0:
        chain.set_alignments(data)
        chains.append(chain)

    return Liftover(chains)


class Liftover:

    def __init__(self, chains):

        self.chain_dict = {}

        for c in chains:
            self.chain_dict[c.tName] = c

    def map(self, feature):

        if feature.chr in self.chain_dict:
            return self.chain_dict[feature.chr].map(feature)


# chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
class Chain:

    def __init__(self, tName, tSize, tStart, tEnd, qName, qSize, qStart, qEnd, id):
        self.tName = tName
        self.tSize = tSize
        self.tStart = tStart
        self.tEnd = tEnd
        self.qName = qName
        self.qSize = qSize
        self.qStart = qStart
        self.qEnd = qEnd
        self.id = id

    def set_alignments(self, data):
        self.tree = self.build_tree(data)

    def build_tree(self, data):

        tree = IntervalTree()

        tStart = self.tStart
        qStart = self.qStart
        for a in data:
            '''
            a is a tuple with 3 elements, with the exception of the last tuple which has only the first
            size -- the size of the ungapped alignment
            dt -- the difference between the end of this block and the beginning of the next block (reference/target sequence)
            dq -- the difference between the end of this block and the beginning of the next block (query sequence)
            '''
            tEnd = tStart + a[0]
            qEnd = qStart + a[0]
            tree[tStart:tEnd] = (qStart, qEnd)

            if len(a) == 3:
                tStart += (a[0] + a[1])
                qStart += a[0]

        return tree

    def map(self, f):

        mapped = []

        intervals = self.tree[f.start:f.end]
        if intervals is not None:

            # A feature might span multiple intervals
            for i in intervals:
                ds = f.start - i.begin

                mf = f.clone()
                mf.chr = self.qName
                mf.start = max(i.data[0], i.data[0] + ds)
                mf.end = min(i.data[1], i.data[0] + ds + f.size())
                mapped.append(mf)

        return mapped
