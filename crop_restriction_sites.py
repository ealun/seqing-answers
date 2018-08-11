#!/usr/bin/env python
import argparse
from collections import defaultdict
import gzip
import itertools
import multiprocessing
import os

from Bio import SeqIO
from joblib import Parallel, delayed


_SUFFIX = 'fastq.gz'


def make_file_names(name, directory):
    left = os.path.join(directory, '{}_1.{}'.format(name, _SUFFIX))
    left_out = os.path.join(directory, '{}_1_resitefound.{}'.format(name, _SUFFIX))
    right = os.path.join(directory, '{}_2.{}'.format(name, _SUFFIX))
    right_out = os.path.join(directory, '{}_2_resitefound.{}'.format(name, _SUFFIX))
    return left, left_out, right, right_out


def crop_sites(args):
    filenames = os.listdir(args.directory)
    base_names = set([filename.split('_')[0] for filename in filenames])
    num_cores = multiprocessing.cpu_count()
    print('Number of cores is: ', num_cores)

    def process_files(name):
        left, left_out, right, right_out = make_file_names(name, args.directory)
        counts = defaultdict(int)
        with gzip.open(left_out, 'wt') as left_outfile, gzip.open(right_out, 'wt') as right_outfile, gzip.open(left, 'rt') as left_infile, gzip.open(right, 'rt') as right_infile:
            print('Processing sequences for ', name)
            for left_seq, right_seq in itertools.izip(SeqIO.parse(left_infile, 'fastq'), SeqIO.parse(right_infile, 'fastq')):
                left_match = left_seq.seq.find(args.restriction_site)
                right_match = right_seq.seq.find(args.restriction_site)

                left_seq = left_seq if left_match < 0 else left_seq[:left_match + 4]
                right_seq = right_seq if right_match < 0 else right_seq[:right_match + 4]

                if left_match < 0 and right_match < 0:
                    counts['none'] += 1
                    continue

                match_type = 'both'
                if left_match > -1 and not right_match > -1:
                    match_tye = 'left'
                elif right_match > -1 and not left_match > -1:
                    match_type = 'right'

                SeqIO.write(left_seq, left_outfile, 'fastq')
                SeqIO.write(right_seq, right_outfile, 'fastq')

                counts[match_type] += 1
        print('Completed processing sequences for ', name)
        return name, dict(counts)
    results = Parallel(n_jobs=num_cores)(delayed(process_files)(name) for name in base_names)
    for result in results:
        print(result)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--restriction-site',
                        default='GATCGATC')
    parser.add_argument('--directory',
                        default='data')
    parser.add_argument('--min-length',
                        default=14)
    return parser.parse_args()


if __name__ == '__main__':
    crop_sites(parse_args())
