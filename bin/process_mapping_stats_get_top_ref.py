#!/usr/bin/env python
# coding: utf-8
import click
import pandas as pd
from Bio import SeqIO


def get_seq(fludb_fasta, seq_id):
    for rec in SeqIO.parse(fludb_fasta, 'fasta'):
        if rec.id == seq_id:
            return rec


@click.command()
@click.option('-s', '--influenza-sequences',
              required=True, help='NCBI Influenza DB sequences (i.e. influenza.fna from FTP site')
@click.option('-m', '--influenza-metadata',
              required=True,
              help='NCBI Influenza DB metadata (i.e. genomeset.dat from FTP site)')
@click.option('-d', '--depths-table', required=True,
              help='Samtools depth output table with header containing fields: ["genome", "position", "depth"].')
@click.option('-i', '--idxstats-table', required=True,
              help='Samtools idxstats output')
@click.option('-o', '--output-ref-fasta',
              default='top_ref.fasta',
              help='Top mapping reference sequence FASTA output path.')
@click.option('-c', '--output-cov-stats',
              default='mapping_coverage_stats_summary.tsv',
              help='Reference sequence mapping summary stats table.')
def process(influenza_sequences, influenza_metadata, depths_table, idxstats_table, output_ref_fasta, output_cov_stats):
    """Process mapping stats and get the top mapping reference sequence.

    The key parameter for determining the top reference sequence is the number 
    of mapped positions ('mapped_positions'). The reference genome with the most
    mapped positions will be output as the top reference for secondary mapping.
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
    df_idxstats = pd.read_csv(idxstats_table, sep='\t',
                              names='genome length mapped unmapped'.split())
    df_idxstats['accession'] = df_idxstats.genome.str.extract(r'gi\|\w+\|gb\|(\w+)\|.*')
    df_idxstats_md = pd.merge(df_idxstats, df_md, on='accession', how='left')
    df_depths = pd.read_csv(depths_table, sep='\t')
    df_depths['accession'] = df_depths.genome.str.extract(r'gi\|\w+\|gb\|(\w+)\|.*')
    df_depths['is_mapped'] = df_depths['depth'] > 0
    grouped_by_acc = df_depths.drop(columns=['position', 'genome']).groupby(by='accession')
    df_cov_summary = pd.DataFrame(
        dict(mapped_positions=grouped_by_acc.is_mapped.sum(),
             mean_coverage=grouped_by_acc.depth.mean(),
             median_coverage=grouped_by_acc.depth.median(), 
             max_coverage=grouped_by_acc.depth.max()))
    df_cov_summary.reset_index(inplace=True)
    df_cov_summary.query('mapped_positions > 0')
    df_mapping_cov_summary = pd.merge(df_cov_summary, df_idxstats_md.query('mapped > 0'), on='accession', how='right')
    df_mapping_cov_summary['p_mapped_positions'] = df_mapping_cov_summary['mapped_positions'] / df_mapping_cov_summary['length']
    df_mapping_cov_summary.sort_values('mapped_positions', ascending=False, inplace=True)
    df_mapping_cov_summary['H_type'] = df_mapping_cov_summary.subtype.str.extract(r'H(\d+)')
    df_mapping_cov_summary['N_type'] = df_mapping_cov_summary.subtype.str.extract(r'N(\d+)')
    top_seq_id = list(df_mapping_cov_summary.genome)[0]
    print(top_seq_id)
    print(df_mapping_cov_summary)

    rec = get_seq(influenza_sequences, top_seq_id)
    print(rec)

    df_mapping_cov_summary.to_csv(output_cov_stats, sep='\t', index=False)
    SeqIO.write([rec], output_ref_fasta, 'fasta')


if __name__ == '__main__':
    process()
