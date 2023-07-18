#!/usr/bin/env python
'''
A Python package for summarizing QC data from the ncov-tools pipeline.
'''
import sys
import qc_parser.models.qc as qc
import qc_parser.models
import csv

def get_qc(args):
    qc_line = dict()
    qc_line.update({'sample' : args.sample})

    try:
        meta = qc_parser.models.Meta(file=args.meta)
        meta.import_metadata()
        qc_line.update(meta.data[args.sample])
    except:
        qc_line.update({'qpcr_ct' : 'NA', 'collection_date' : 'NA',
                        'num_months' : 'NA', 'num_weeks' : 'NA'})

    if args.platform == 'illumina':
        if str(args.variants).endswith('.variants.tsv'):
            vars = qc_parser.models.Variants(file=args.variants)
            qc_line.update(vars.get_total_variants())
        elif str(args.variants).endswith('.vcf') or str(args.variants).endswith('.vcf.gz'):
            vars = qc_parser.models.Vcf(file=args.variants)
            qc_line.update(vars.get_variant_counts())
        else:
            sys.exit('Must be a valid variant.tsv or .vcf file for the Illumina platform')
    elif args.platform == 'oxford-nanopore':
        if str(args.variants).endswith('.vcf') or str(args.variants).endswith('.vcf.gz'):
            vars = qc_parser.models.Vcf(file=args.variants)
            qc_line.update(vars.get_variant_counts())
        else:
            sys.exit('Must be a valid VCF file for the Oxford-Nanopore platform')


    alleles = qc_parser.models.Alleles(file=args.alleles)
    qc_line.update(alleles.get_variant_counts(sample=args.sample))

    cons = qc_parser.models.Consensus(file=args.consensus)
    qc_line.update(cons.count_iupac_in_fasta())
    qc_line.update(cons.get_genome_completeness())

    coverage = qc_parser.models.PerBaseCoverage(file=args.coverage)
    qc_line.update(coverage.get_coverage_stats())

    # Add the lineage from the Pangolin report
    try:
        lineage = qc_parser.models.Lineage(file=args.lineage, pangolin_ver=args.pangolin_version)
        lineage.create_lineage_dictionary()
        qc_line.update({"lineage" : lineage.lineage_dict[args.sample]["lineage"]})
        qc_line.update({"lineage_notes" : lineage.lineage_dict[args.sample]["notes"]})
        qc_line.update({"scorpio_call" : lineage.lineage_dict[args.sample]["scorpio_call"]})
    except:
        qc_line.update({"lineage" : "none"})
        qc_line.update({"lineage_notes" : "none"})
        qc_line.update({"scorpio_call" : "none"})

    # Add the watch list mutations
    try:
        watchlist = qc_parser.models.WatchList(file=args.mutations)
        qc_line.update({"mutations" : watchlist.get_mutation_string(sample=args.sample)})
    except:
        qc_line.update({"mutations" : "none"})

    # Get a list of consequences from the SNPEff variant annotations
    frameshift_indels = False
    try:
        annotations = qc_parser.models.Snpeff(file=args.aa_table)
        annotations.get_list_of_consequences()
        if annotations.has_frameshift():
            frameshift_indels = True
    except:
        pass

    # Produce warning flags
    qc_flags = list()
    if qc_line['genome_completeness'] < 0.5:
        qc_flags.append("INCOMPLETE_GENOME")
    elif qc_line['genome_completeness'] < 0.9:
        qc_flags.append("PARTIAL_GENOME")

    #num_indel_non_triplet = qc_line['num_variants_indel'] - qc_line['num_variants_indel_triplet']
    #if num_indel_non_triplet > 0:
    if frameshift_indels:
        qc_flags.append("POSSIBLE_FRAMESHIFT_INDELS")

    if qc_line['num_consensus_iupac'] > 5:
        qc_flags.append("EXCESS_AMBIGUITY")

    # the mixture report is currently generated for illuina runs, ont is
    # not supported at this time
    if args.mixture:
        if args.platform == 'illumina':
            mixture = set()
            with open(args.mixture, 'r') as mfh:
                reader = csv.DictReader(mfh, delimiter='\t')
                for record in reader:
                    mixture.add(record['sample_a'])
            if args.sample in mixture:
                qc_flags.append("POSSIBLE_MIXTURE")

    # Calculate number of variants per week, while accounting for incompleteness
    if qc_line['num_weeks'] != 'NA':

        if qc_line['genome_completeness'] > 0.1:
            scaled_variants = qc_line['num_variants_snvs'] / qc_line['genome_completeness']

            # very conservative upper limit on the number of acceptable variants
            # samples that fail this check should be manually reviewed incorporating other
            # evidence (frameshift indels, not failed outright)
            # Removing EXCESS_VARIANTS flag due to discovery of lineage carrying
            # high number of mutations -- December 17, 2020.
            #variant_threshold = qc_line['num_weeks'] * 0.75 + 15
            #if scaled_variants > variant_threshold:
            #    qc_flags.append("EXCESS_VARIANTS")
            qc_line['scaled_variants_snvs'] = "%.2f" % (scaled_variants)
        else:
            qc_line['scaled_variants_snvs'] = "NA"
    else:
        qc_line['scaled_variants_snvs'] = "NA"

    qc_flag_str = "PASS"
    if len(qc_flags) > 0:
        qc_flag_str = ",".join(qc_flags)

    qc_line.update({'qc_pass' : qc_flag_str})
    qc_line.update({'run_name' : args.run_name})

    qc.write_qc_summary_header()
    qc.write_qc_summary(summary=qc_line)
