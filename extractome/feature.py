from extractome.stream import getstream

class Feature:

    def __init__(self, chr, start, end, name='', score=1000, color=None):
        self.chr = chr
        self.start = start
        self.end = end
        self.name = name
        self.score = score
        self.strand = '+'
        self.color = color

    def clone(self):
        return Feature(self.chr, self.start, self.end, self.name, self.score, self.strand, self.color)

    def size(self):
        return self.end - self.start

    def setcolor(self, color):
        self.color = color

    def tostring(self):
        if self.color is not None:
            return f"{self.chr}\t{self.start}\t{self.end}\t{self.name}\t{self.score}\t{self.strand}\t{self.start}\t{self.end}\t{self.color}"
        else:
            return f"{self.chr}\t{self.start}\t{self.end}\t{self.name}\t{self.score}\t{self.strand}"

def parse(path, format=None):
    '''
    Parse a feature file and return an array of feature objects.  Supported formats are bed, gff, and gtf.
    :param path: Path to feature file, which can be local or url
    :param format: File format, bed | gtf | gff
    :param split_bool: Boolean specifying whether view is multi locus or not
    :return: List of feature objects {chr, start, end, text, name}
    '''

    f = None
    try:
        f = getstream(path)
        if not format:
            format = infer_format(path)
        if format == 'bed':
            return parse_bed(f)
        elif format == 'gff'  or format == 'gff3' or format == 'gtf':
            return parse_gff(f)
        elif format == 'tab':
            return parse_tab(f)
        elif format == 'bedpe':
            return parse_bedpe(f)
        elif format == 'refgene':
            return parse_refgene(f)
        else:
            raise Exception("Unknown file format: " + path)
    finally:
        if f:
            f.close()


def parse_bed(f):
    features = []
    for line in f:
        if not (line.startswith('#') or line.startswith('track') or line.startswith('browser')):
            tokens = line.rstrip('\n').rstrip('\r').split('\t')
            if len(tokens) >= 3:
                chr = tokens[0]
                start = int(tokens[1])
                end = int(tokens[2])
                name = tokens[3] if len(tokens) > 3 else ''
                score = tokens[4] if len(tokens) > 4 else None
                features.append(Feature(chr, start, end, name, score))
    return features


def parse_refgene(f):
    features = []
    for line in f:
        if not (line.startswith('#') or line.startswith('track') or line.startswith('browser')):
            tokens = line.rstrip('\n').rstrip('\r').split('\t')
            if len(tokens) >= 3:
                chr = tokens[2]
                start = int(tokens[4])
                end = int(tokens[5])
                name = tokens[12] if len(tokens) > 12 else ''
                features.append(Feature(chr, start, end, line, name))
    return features


def parse_bedpe(f):
    features = []
    for line in f:
        if not (line.startswith('#') or line.startswith('track') or line.startswith('browser')):
            tokens = line.rstrip('\n').rstrip('\r').split('\t')
            if len(tokens) >= 6:
                chr = tokens[0]
                start = int(tokens[1])
                end = int(tokens[2])
                chr2 = tokens[3]
                start2 = int(tokens[4])
                end2 = int(tokens[5])
                name = tokens[6] if len(tokens) > 6 else ''
                features.append(PEFeature(chr, start, end, line, name, chr2, start2, end2))
    return features


def parse_gff(f):
    features = []
    for line in f:
        if not (line.startswith('#') or line.startswith('track') or line.startswith('browser')):
            tokens = line.rstrip('\n').rstrip('\r').split('\t')
            # if we encounter a blank or malformed line (no start and end coords), skip it
            if len(tokens) < 5:
                continue
            chr = tokens[0]
            start = int(tokens[3]) - 1
            end = int(tokens[4])
            name = ''
            features.append(Feature(chr, start, end, line, name))

    return features


def parse_tab(f):
    rows = []
    for line in f:
        if not (line.startswith('#') or line.startswith('track') or line.startswith('browser')):
            tokens = line.rstrip('\n').rstrip('\r').split('\t')
            if len(tokens) > 2:
                rows.append(tokens)
    return rows


def infer_format(filename):
    '''
    Infer the genomic file format from the filename.  First known formats are checked.  Next presenece of
    the magic string "refgene" in the filename is checked for UCSC refgene files.  This is a legacy
    IGV convention.  The order is important, a recognized extension wins.
    NOTE: Formats are for output data uris.  CRAM format is converted to BAM before output.
    :param filename:
    :return:
    '''
    filename = filename.lower()
    if (filename.endswith(".gz")):
        filename = filename[:-3]

    if filename.endswith(".bam"):
        return "bam"
    if filename.endswith(".cram"):
        return "cram"
    elif filename.endswith(".vcf"):
        return "vcf"
    elif filename.endswith(".bed"):
        return "bed"
    elif filename.endswith(".bedpe"):
        return "bedpe"
    elif filename.endswith(".gff") or filename.endswith(".gff3"):
        return "gff"
    elif filename.endswith(".gtf"):
        return "gtf"
    elif filename.find("refgene"):
        return "refgene"
    else:
        idx = filename.rfind(".")
        if idx > 0:
            return filename[idx + 1:]
        else:
            return None
        idx = filename.rfind(".")
        if idx > 0:
            return filename[idx + 1:]
        else:
            return None
