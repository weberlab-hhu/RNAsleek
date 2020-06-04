#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re
import copy
import argparse


FAILPASS = [
    "basic_statistics",
    "per_base_sequence_quality",
    "per_sequence_quality_scores",
    "per_base_sequence_content",
    "per_sequence_gc_content",
    "per_base_n_content",
    "sequence_length_distribution",
    "sequence_duplication_levels",
    "overrepresented_sequences",
    "adapter_content",
    "kmer_content"
]


IN_SUMMARY = FAILPASS + ['%GC', 'total_deduplicated_percentage', 'avg_sequence_length']


def skip_indeces(df, to_skip):
    to_keep = []
    for i in range(df.shape[0]):
        if df.index[i] not in to_skip:
            to_keep.append(i)
        else:
            print('skipping {}'.format(df.index[i]))
    print(to_keep)
    return df.iloc[to_keep]


def confirm_sort(dfs):
    for df in dfs[1:]:
        assert all([x == y for x, y in zip(dfs[0].index, df.index)]), "0: {}, other: {}".format(
            len(dfs[0].index), len(df.index)
        )


def fix_sort(df_ref, df):
    df = df.reindex(index=df_ref.index)
    return df


def weighted_averages(df, weights, weightable):
    weights = df[weights]
    scaled = weights * df.loc[:, weightable].T
    summed = np.sum(scaled, axis=1)
    tot = np.sum(weights)
    return summed / tot


def collapse_fastqc2samples(fastqc, run2sample):
    num_sum = {float('nan'): float('nan'), 'fail': 0., 'warn': 0.5, 'pass': 1.}

    fastqc = copy.deepcopy(fastqc)
    fastqc.loc[:, "Run"] = [re.sub('_.*', '', x) for x in fastqc.loc[:, 'Sample']]
    fastqc.loc[:, 'SampleSample'] = [run2sample[x] for x in fastqc.loc[:, "Run"]]
    x = copy.deepcopy(fastqc.loc[:, FAILPASS])
    for column in x.columns:
        x.loc[:, column] = [num_sum[var] for var in x.loc[:, column]]

    fastqc.loc[:, 'N_Bases'] = fastqc.loc[:, 'Total Sequences'] * fastqc.loc[:, 'avg_sequence_length']
    fastqc[FAILPASS] = x
    grouped = fastqc.groupby(fastqc.loc[:, 'SampleSample'])
    collapsed = grouped.apply(weighted_averages, weights='N_Bases', weightable=IN_SUMMARY)
    # percents to fractions
    percents = ['%GC', 'total_deduplicated_percentage']
    collapsed.loc[:, percents] = collapsed.loc[:, percents] / 100
    # also, summed reads, percent to fraction as applicable
    return collapsed


def collapse_trimmomatic2samples(trimmomatic, run2sample):
    trimmomatic = copy.deepcopy(trimmomatic)
    trimmomatic.loc[:, "Run"] = [re.sub('_1', '', x) for x in trimmomatic.loc[:, 'Sample']]
    trimmomatic.loc[:, 'SampleSample'] = [run2sample[x] for x in trimmomatic.loc[:, "Run"]]
    # sum up raw values
    raw_counts = ['input_read_pairs', 'input_reads', 'surviving', 'forward_only_surviving', 'reverse_only_surviving',
                  'dropped']
    grouped = trimmomatic.loc[:, raw_counts].groupby(trimmomatic.loc[:, "SampleSample"])
    summed = grouped.apply(np.sum)
    # collapse SE and PE to fragments
    summed.loc[:, "fragment_counts"] = np.sum(summed.loc[:, ['input_read_pairs', 'input_reads']], axis=1)
    percents = (summed.loc[:, ['surviving',
                               'forward_only_surviving',
                               'reverse_only_surviving',
                               'dropped']].T / summed.loc[:, 'fragment_counts']).T

    percents.columns = ['surviving', 'F_surviving', 'R_surviving', 'dropped']
    return summed.loc[:, 'fragment_counts'], percents


def median_diff_coverage(coverage):
    # 101 "positions", so ignore middle. [::-1] reverses
    x = coverage.iloc[51:].values - coverage.iloc[:50][::-1].values
    return np.quantile(x, 0.5)


def main(srafile, pdfout, basedir, skip=None):

    sraruninfo = pd.read_csv(srafile)

    run2sample = {}
    for item in sraruninfo.loc[:, ['Run', 'Sample']].groupby(sraruninfo.loc[:, 'Sample']):
        for run in item[1].iloc[:, 0]:
            run2sample[run] = item[0]

    # import all QC output stats
    fastqc_trimmed = pd.read_csv(basedir + '/multiqc/multiqc_data/multiqc_fastqc.txt', delimiter='\t')
    fastqc_untrimmed = pd.read_csv(basedir + '/multiqc_untrimmed/multiqc_data/multiqc_fastqc.txt', delimiter='\t')
    trimmomatic = pd.read_csv(basedir + '/multiqc/multiqc_data/multiqc_trimmomatic.txt', delimiter='\t')
    hisat = pd.read_csv(basedir + '/multiqc/multiqc_data/multiqc_bowtie2.txt', sep='\t')
    picard = pd.read_csv(basedir + '/multiqc/multiqc_data/multiqc_picard_RnaSeqMetrics.txt', sep='\t')

    cov = pd.read_csv(basedir + '/coverage.tsv', sep='\t')

    # columns from fastqc that can be usefully displayed
    # note ignoring "per_tile_sequence_quality", bc not in all fastq files and not so meaningful here

    # clean, tidy, summarize sample-centric files as necessary
    hisat["SampleSample"] = [re.sub('.hisat.err', '', x) for x in hisat["Sample"]]

    hisat = hisat.set_index('SampleSample')

    picard = picard.set_index('Sample')

    bases = picard.loc[:, ['CODING_BASES', 'UTR_BASES', 'INTRONIC_BASES', 'INTERGENIC_BASES']]
    bases = bases.apply(lambda x: np.sum(x[:2]) / np.sum(x), axis=1)

    right_strand = picard['PCT_CORRECT_STRAND_READS'] / 100

    qc_untrimmed = collapse_fastqc2samples(fastqc_untrimmed, run2sample)
    qc_trimmed = collapse_fastqc2samples(fastqc_trimmed, run2sample)
    fragcounts, trimmo = collapse_trimmomatic2samples(trimmomatic, run2sample)
    unmapped = 1 - hisat.loc[:, "overall_alignment_rate"] / 100

    qc_merge = qc_untrimmed.copy().loc[:, IN_SUMMARY[:-1]]
    qc_merge.columns = ['basic_statistics', 'base_quality',
                        'seq_quality', 'base_ATCG',
                        'seq_gc_content', 'base_n_content',
                        'length', 'duplication',
                        'overrepr._seqs', 'adapters',
                        'kmer_content', 'f_GC', 'f_deduplicated']
    cleaning = qc_trimmed.copy()
    cleaning = cleaning.loc[:, ["sequence_length_distribution", "adapter_content"]]
    cleaning.columns = ["TRIMMED_length", "TRIMMED_adapters"]

    qc_merge = qc_merge.merge(cleaning, left_index=True, right_index=True)
    qc_merge = qc_merge.merge(trimmo, left_index=True, right_index=True)

    cov = cov.reindex(columns=qc_untrimmed.index)

    coverage = cov.apply(median_diff_coverage, axis=0)
    cov = cov.T  # flip back before plotting, but then it can be handled like the others

    plottables = {"qc_merge": qc_merge,
                  "fragcounts": fragcounts,
                  "qc_trimmed": qc_trimmed,
                  "right_strand": right_strand,
                  "unmapped": unmapped,
                  "bases": bases,
                  "coverage": coverage,
                  "cov": cov}
    # drop indeces as requested
    if skip is not None:
        with open(skip) as f:
            to_skip = f.readlines()
        to_skip = [x.rstrip() for x in to_skip]
        # filter every table
        for key in plottables:
            plottables[key] = skip_indeces(plottables[key], to_skip)

    # check and fix sorting
    try:
        confirm_sort(list(plottables.values()))
    except AssertionError:
        for key in plottables:
            if not key == "qc_merge":
                plottables[key] = fix_sort(plottables['qc_merge'], plottables[key])
        confirm_sort(list(plottables.values()))

    qc_merge = plottables['qc_merge']  # shortcut for most used
    # AND PLOTTING :-)

    grey = "0.85"
    bar_colors = ["#4d0000"] * 4 + ["#03396c"] * 4
    # so 40 samples have length 15
    labspace = 1.55
    figsize = (14, labspace + (15 - labspace) * qc_merge.shape[0] / 40)

    # qc untrimmed
    fig, (axqc, axfrags, axlen, axstrand,
          axmap, axexpect,  ax3p, axcov) = plt.subplots(nrows=1,
                                                        ncols=8,
                                                        figsize=figsize,
                                                        sharey=True,
                                                        gridspec_kw={"width_ratios": [2, 0.5, 0.5, 0.5, 1, 1, 1, 1]})

    # fastqc pre/post trim, trimmomatic
    axqc.set_xticks(list(range(qc_merge.shape[1])))
    axqc.set_yticks(list(range(qc_merge.shape[0])))
    axqc.set_xticklabels(qc_merge.columns, rotation=90)
    axqc.set_yticklabels(qc_merge.index)
    axqc.imshow(qc_merge, aspect="auto")
    axqc.margins(x=0)

    current_lim = list(axqc.get_ylim())
    current_lim[1] -= 1  # extra space for coverage

    axqc.set_ylim(current_lim)
    for ybar in range(0, qc_merge.shape[0], 4):
        axqc.axhline(y=ybar - 0.5, color="black")

    # million read pairs (fragments)
    plt.setp(axfrags.get_yticklabels(), visible=False)

    axfrags.set_xlabel("Frag. * 10 ^ 6")
    axfrags.barh(width=plottables['fragcounts'] / 10**6, y=qc_merge.index, color=bar_colors)
    axfrags.axvline(x=5, color=grey)

    # sequence length
    plt.setp(axlen.get_yticklabels(), visible=False)
    axlen.set_xlabel("length")
    axlen.barh(width=plottables["qc_trimmed"].loc[:, 'avg_sequence_length'], y=qc_merge.index, color=bar_colors)
    axlen.axvline(x=100, color=grey)

    # stranded check
    plt.setp(axstrand.get_yticklabels(), visible=False)
    axstrand.set_xlabel("right strand")
    axstrand.set_xlim((0., 1.))
    axstrand.barh(width=plottables["right_strand"], y=qc_merge.index, color=bar_colors)
    axstrand.axvline(x=0.5, color=grey)


    # unmapped rate
    plt.setp(axmap.get_yticklabels(), visible=False)
    axmap.axvline(x=0.1, color=grey)
    axmap.set_xlabel("unmapped")
    axmap.barh(width=plottables["unmapped"], y=qc_merge.index, color=bar_colors)

    # mapping where transcripts are expected
    plt.setp(axexpect.get_yticklabels(), visible=False)
    axexpect.axvline(x=0.1, color=grey)
    axexpect.set_xlabel("not CDS|UTR")
    axexpect.barh(width=1 - plottables["bases"], y=qc_merge.index, color=bar_colors)

    # 3' bias
    plt.setp(ax3p.get_yticklabels(), visible=False)
    ax3p.axvline(x=0, color=grey)
    ax3p.set_xlabel("median 3' diff")
    ax3p.barh(width=plottables["coverage"], y=qc_merge.index, color=bar_colors)
    ax3p.set_xticklabels(list(ax3p.get_xticks()), rotation=90)

    # coverage curvegs
    plt.setp(axcov.get_yticklabels(), visible=False)
    axcov.set_xlabel("coverage")
    axcov.axvline(x=50, color=grey)
    cov = plottables["cov"].T
    for i in range(cov.shape[1]):
        axcov.plot(list(range(101)), cov.iloc[:, i] * -1 + i + 0.1, color=bar_colors[i % 8])

    # fraction of figsize reserved for top and bottom margin
    total_margins = labspace / (labspace + (15 - labspace) * cov.shape[1] / 40)
    fig.subplots_adjust(bottom=total_margins * 0.95,
                        top=1 - total_margins * 0.05)

    plt.savefig(pdfout, format="pdf")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('directory', help='project working directory')
    parser.add_argument('run_file', help='file with SRR numbers, one per line')
    parser.add_argument('-o', '--pdf_out', required=True, help="output filename for plots")
    parser.add_argument('-s', '--skip_from_file',
                        help="file of samples IDs to skip (useful when a few samples skew the plots too much")
    args = parser.parse_args()
    main(args.run_file, args.pdf_out, args.directory, args.skip_from_file)





