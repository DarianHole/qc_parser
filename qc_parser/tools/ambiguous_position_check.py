'''
Parse an alleles.tsv file to identify possible problematic positions
'''
import qc_parser.models
from collections import defaultdict

def ambiguous_positions_check(args):

    alleles = qc_parser.models.Alleles(file=args.alleles)

    ambiguous_counts = defaultdict(int)
    ambiguous_alleles = defaultdict(dict)
    for sample, records in alleles.data.items():
        for position in records:
            aa = records[position]['alt']
            if aa != 'N' and qc_parser.models.is_variant_iupac(aa):
                p = int(position)
                ambiguous_counts[p] += 1
                ambiguous_alleles[p][aa] = 1

    print("position\tcount\talleles")
    for position in sorted(ambiguous_counts.keys()):
        count = ambiguous_counts[position]
        if count >= args.min_count:
            print("%d\t%d\t%s" % (int(position), count, ",".join(ambiguous_alleles[position])))
