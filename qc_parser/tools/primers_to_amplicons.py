#!/usr/bin/env python
'''
Convert the nCoV primer scheme to a unique amplicon BED file.
'''
from qc_parser.models import primers as pr

def primers_to_amplicons(args):
    primers = pr.read_bed_file(args.primers)
    primer_pairs = pr.create_primer_pairs(primers=primers)
    amplicon_ranges = pr.create_amplicons(primer_pairs=primer_pairs,
                                        offset=args.offset,
                                        type=args.bed_type,
                                        prefix=args.primer_prefix)

    with open(args.output, 'w') as file_o:
        for line in amplicon_ranges:
            file_o.write('\t'.join(line))
            file_o.write('\n')
    file_o.close()
