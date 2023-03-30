#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import os
import zsx_some_tools as st
from optparse import OptionParser
import sys


def extract_motif(sequence_, motif_, extract_1=12, extract_2=False, chr_=''):
    n_ = len(motif_[0])
    if not extract_2:
        extract_2 = extract_1

    start = extract_1
    end = int(len(sequence_)) - (extract_2 + n_ - 1)

    result_ = []
    for i_ in range(start, end):
        mot = sequence_[i_: i_ + n_]
        if mot in motif_:
            seq_out = sequence_[i_ - extract_1: i_ + extract_2 + n_]
            seq_left = seq_out[: extract_1]
            seq_mid = seq_out[extract_1: extract_1 + n_]
            seq_right = seq_out[extract_1 + n_:]

            result_ += [[seq_left, seq_mid, seq_right, chr_]]

    return result_


def main():
    usage = "usage: python %prog -i pathin -o pathout"
    description = "motif extract -i -o options are needed, -m -l -r options are optional"
    op = OptionParser(version="%prog 2.0", description=description, usage=usage, add_help_option=False)
    op.add_option("-h", "--help", action="help",
                  help="Show this help message and exit.")

    op.add_option("-i", "--inputpath", dest="pathin", type="str",
                  help="The input path-dir.")

    op.add_option("-o", "--outputpath", dest="pathout", type="str",
                  help="The output path-dir.")

    op.add_option("-m", "--motif", dest="motif", type="str", default='WRC',
                  help="motif name such as WRC.")

    op.add_option("-l", "--left_num", dest="extrac_left_num", type=int, default=6,
                  help="The number of nt to extract from left of the motif, default is 6.")

    op.add_option("-r", "--right_num", dest="extrac_right_num", type=int, default=None,
                  help="The number of nt to extract from right of the motif, default is None, "
                       "which means this parameter will be set the same as -l.")

    (options, args) = op.parse_args()

    if not options.pathin:
        op.print_help()
        sys.exit(1)
    if not options.pathout:
        op.print_help()
        sys.exit(1)

    pathin = options.pathin
    pathout = options.pathout
    motif_name = options.motif
    left_num = options.extrac_left_num
    right_num = options.extrac_right_num

    pathin = st.path_diagnosis(pathin)
    pathout = st.path_diagnosis(pathout)
    st.mkdir(pathout)

    motif = st.design_motif(motif_name)

    if right_num is None:
        right_num = left_num

    file_list = st.listdir(pathin, target='.fa|.seq|.fasta', sep='|', logic_and=False)
    for file in file_list:
        path = pathin + file
        file_name = file.split('.')[0] + '_motif_extract.txt'

        seq = ''
        result = []
        i = 0
        for line in open(path):
            if '>' in line:
                if i > 0:
                    result += extract_motif(seq, motif, extract_1=left_num, extract_2=right_num, chr_=name)

                name = line.strip('>').strip('\n')
                if '_' in name:
                    ok = False
                else:
                    ok = True
                seq = ''
                i += 1
            else:
                if not ok:
                    continue
                seq += line.strip('\n').upper()

        result += extract_motif(seq, motif, extract_1=left_num, extract_2=right_num, chr_=name)
        result = pd.DataFrame(result)
        result.columns = ['left', 'motif', 'right', 'chr']
        result.to_csv(pathout + file_name, sep='\t', index=False, header=False)


if __name__ == "__main__":
    main()
    