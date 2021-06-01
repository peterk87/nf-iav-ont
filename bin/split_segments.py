#!/usr/bin/env python
# coding: utf-8
import os
import re

import click
import pandas as pd
from Bio import SeqIO


def segment_seqs_over_first_quartile(df_md, seg, threshold='25%'):
    q1 = df_md.loc[seg,'seq_length'].describe()[threshold]
    dfseg = df_md.loc[seg,:]
    dfseg = dfseg[dfseg.seq_length >= q1]
    return {acc:seg for acc in dfseg.accession}


@click.command()
@click.option('-s', '--influenza-sequences',
              required=True, help='NCBI Influenza DB sequences (i.e. influenza.fna from FTP site')
@click.option('-m', '--influenza-metadata',
              required=True,
              help='NCBI Influenza DB metadata (i.e. genomeset.dat from FTP site)')
@click.option('-o', '--outdir',
              default='./',
              help='Output destination for influenza genome segment sequence FASTA files')
@click.option('-N', '--exclude-seqs-sequential-ns',
              type=int,
              default=10,
              help='Exclude sequences with X or more sequential "N" characters (default: 10)')
@click.option('--seq-length-threshold-stat',
              type=click.Choice(['mean', 'min', '25%', '50%', '75%']),
              default='25%',
              help=('Pandas describe method statistic to use as min segment'
                    ' sequence threshold. (default "25%" or sequences must be '
                    'at least as long as the sequence at the 25th percentile '
                    'for each segment)'))
def split(influenza_sequences, influenza_metadata, outdir, exclude_seqs_sequential_ns, seq_length_threshold_stat):
    """Split influenza segment sequences into separate files.

    Exclude sequences that have too many consecutive Ns (default >= 10) and are shorter than the 25th percentile sequence for each segment (by default)
    """
    md_cols = [('accession', str),
               ('host', 'category'),
               ('segment', 'category'),
               ('subtype', 'category'),
               ('country', 'category'),
               ('date', 'category'),
               ('seq_length', 'uint16'),
               ('virus_name', 'category'),
               ('age', 'category'),
               ('gender', 'category'),
               ('group_id', 'category'),]
    df_md = pd.read_csv(influenza_metadata, sep='\t', names=[name for name, _ in md_cols], dtype={name:t for name,t in md_cols})
    df_md.set_index('segment', inplace=True)

    acc_regex = re.compile(r'gi\|\d+\|gb\|(\w+)\|\w+')
    regex_Ns = re.compile(r'N{' + str(exclude_seqs_sequential_ns) + r',}')
    accs_to_write = {}
    for seg in df_md.index.unique():
        accs_to_write.update(segment_seqs_over_first_quartile(df_md, seg))
    seg_to_fh = {seg: open(os.path.join(outdir, f'{seg}.fasta'), 'w') for seg in df_md.index.unique()}

    try:
        count = 0
        for rec in SeqIO.parse(influenza_sequences, 'fasta'):
            count += 1
            if count % 10000 == 0:
                print(f'Parsed {count}')
            m = acc_regex.match(rec.id)
            if m:
                acc = m.group(1)
                seg = accs_to_write.get(acc, None)
                fh = seg_to_fh.get(seg, None)
                if fh:
                    m = regex_Ns.search(str(rec.seq))
                    if m:
                        print(f'Too many Ns in a row in {rec.description}. N count={sum([1 for nt in rec.seq if nt == "N"])}')
                        continue
                    fh.write(f'>{rec.description}\n{rec.seq}\n')
    except Exception as e:
        raise e
    finally:
        for fh in seg_to_fh.values():
            fh.close()
    print('Done!')


if __name__ == '__main__':
    split()
