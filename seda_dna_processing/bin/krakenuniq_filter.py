#!/usr/bin/env python3

# -*- coding: utf-8 -*-


from sys import argv

import pandas as pd

import argparse


parser = argparse.ArgumentParser(

prog='krakenuniq_filter.py',

description='Filters KrakenUniq report by Species',

epilog='Siobhan is cool!')



parser.add_argument('--n_unique_kmers', help='Minimum number of unique k-mers', default=1000)

parser.add_argument('--n_tax_reads', help='Minimum number of reads assigned to taxonomy', default=100)

parser.add_argument('--ratio', help='Minimum ratio kmers:taxReads', default=1)

parser.add_argument('--rank', help='Rank to filter', default="species")

parser.add_argument('--only_reads', help='Do not filter by ratio', action='store_true')

parser.add_argument('--only_ratio', help='Do not filter by reads', action='store_true')

parser.add_argument('--krakenuniq_report', help='Path to krakenuniq report', required=True, default=None)



args = parser.parse_args()


krakenuniq_output = args.krakenuniq_report

n_unique_kmers = int(args.n_unique_kmers)

n_tax_reads = int(args.n_tax_reads)

ratio = float(args.ratio)

only_reads = args.only_reads

only_ratio = args.only_ratio

rank = args.rank


# Read and filter KrakenUniq output

kraken_output_df = pd.read_csv(krakenuniq_output, delimiter="\t", comment="#")


# Filter by n_unique_kmers

kraken_output_df["ratio"] = kraken_output_df["kmers"] / kraken_output_df["taxReads"]

kraken_output_df = kraken_output_df[kraken_output_df["rank"] == rank]

print(f"Data set with only species: {kraken_output_df.shape}")


if not only_ratio:

	kraken_output_df = kraken_output_df[kraken_output_df["kmers"] > n_unique_kmers]

	print(f"Data set dimensions after breadth of coverage filter: {kraken_output_df.shape}")

# Filter by n_reads classified as taxa

kraken_output_df = kraken_output_df[kraken_output_df["taxReads"] > n_tax_reads]

print(f"Data set dimensions after depth of coverage filter: {kraken_output_df.shape}")


if not only_reads:

	kraken_output_df = kraken_output_df[kraken_output_df["taxReads"] > 0]

	kraken_output_df = kraken_output_df[kraken_output_df["ratio"] > ratio]

print(f"Data set dimensions after ratio filter: {kraken_output_df.shape}")


kraken_output_df = kraken_output_df.sort_values(by="%", ascending=False)

kraken_output_df.to_csv(f"{krakenuniq_output}.filtered", sep="\t", index=False)