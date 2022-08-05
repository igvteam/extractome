import os
import argparse

from extractome.fasta import FastaReader
from extractome.genome import get_genome
from extractome.feature import parse
from extractome.liftover import load_liftover

'''
This is the main function for the application.
'''
def extract_genome(args):

    # Create output directory
    if not os.path.exists(args.output):
        os.mkdir(args.output)

    # Read the region data (bed file)
    region_dict = build_region_dict(parse(args.regions, 'bed'))
    chrlist = list(region_dict.keys())
    chrlist.sort()

    # Check for optional genome argument.  If supplied an igv.js genome json definition is used in lieu of a fasta file
    if args.genome is not None:
        genome = get_genome(args.genome)
        if args.fasta is None:
            args.fasta = genome["fastaURL"]

    '''
    Create fasta
    '''
    fasta_reader = FastaReader(args.fasta)
    # output_file = os.path.join(args.output, f"{args.name}.fa")
    # with open(output_file, "w") as o:
    #     for chr in chrlist:
    #         o.write(f">{chr}\n")
    #         regions = region_dict[chr]
    #         for r in regions:
    #             seq = fasta_reader.slice({"chr": chr, "start": r.start + 1, "end": r.end})
    #             o.write(seq)
    #         o.write("\n")

    ''' 
    Create .chain file
    '''
    chromsizes = {}
    for chr in chrlist:
        regions = region_dict[chr]
        size = 0;
        for r in regions:
            size += r.end - r.start
        chromsizes[chr] = size

    chains_file = os.path.join(args.output, f"{args.name}.chain")
    with open(chains_file, "w") as o:
        id = 0
        for chr in chrlist:
            id += 1
            regions = region_dict[chr]
            qstart = 0
            qsize = chromsizes[chr]
            tstart = regions[0].start
            tsize = fasta_reader.size(chr) - tstart
            o.write(f"chain 1000 {chr} {tsize} + {tstart} {tsize} {chr} {qsize} + {qstart} {chromsizes[chr]} {id}\n")

            for i in range(len(regions) - 1):
                # size dt dq
                r = regions[i]
                rnext = regions[i + 1]
                o.write(f"{r.size()} {rnext.start - r.end} {0}\n")
            o.write(f"{rnext.size()}\n\n")

    '''
    Create bed file for marking regions -- alternating colors
    '''

    liftover = load_liftover(chains_file)
    color1 = '100,200,100'
    color2 = '100,100,200'
    color = color1

    output_file = os.path.join(args.output, f"{args.name}.regions.bed")
    with open(output_file, "w") as o:
        for chr in chrlist:
            regions = region_dict[chr]
            for r in regions:
                liftovers = liftover.map(r)
                for f in liftovers:
                    f.setcolor(color)
                    o.write(f.tostring() + "\n")
                    color = color1 if color is color2 else color2





def build_region_dict(region_list):
    dict = {}
    for region in region_list:
        if region.chr not in dict:
            dict[region.chr] = []
        dict[region.chr].append(region)

    for rlist in dict.values():
        rlist.sort(key=lambda r: r.start)

    return dict


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("regions", help="bed file defining regions, required")
    parser.add_argument("--fasta", default=None, help="reference fasta file, required if --genome is not specified")
    parser.add_argument("--genome", help="igv.js genome id (e.g. hg38)")
    parser.add_argument("--name", help="xome name", default="Xome")
    parser.add_argument("--output", help="output directory name", default="output")
    args = parser.parse_args()
    extract_genome(args)


if __name__ == "__main__":
    main()
