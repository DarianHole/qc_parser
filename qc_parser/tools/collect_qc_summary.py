#!/usr/bin/env python
'''
A script for aggregating sample QC summary files into a single file.
'''
from qc_parser.models.qc import collect_qc_summary_data, write_qc_summary_header

def collect_qc_summary(args):
    summary_data = collect_qc_summary_data(path=args.path)
    summary_data.sort()

    write_qc_summary_header()
    for summary_line in summary_data:
        print(summary_line)
