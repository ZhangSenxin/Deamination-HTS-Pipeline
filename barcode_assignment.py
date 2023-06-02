#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import zsx_some_tools as st
from scipy import spatial
import sys
import argparse
import warnings

warnings.filterwarnings("ignore")


def get_index(data, da0, Class0, barcode):
    ind_all = pd.DataFrame(None)
    barcode_len = barcode.shape[0]
    for i in range(len(Class0)):
        l = len(data.loc[da0.iloc[Class0[i], :].index])
        if l > barcode_len:
            raise ValueError('The number of barcodes is less than the number of a certain class. '
                             'Please increase the combination of barcodes.')

        ind = pd.DataFrame(da0.iloc[Class0[i], :].index)
        index = pd.concat([ind, barcode.iloc[:l, :]], axis=1)
        ind_all = ind_all.append(index)

    return ind_all


def main():
    usage = "usage: python %(prog)s -i input_path -o output_path -m meta"
    description = "-i -o -m option is needed"
    parser = argparse.ArgumentParser(prog="%prog 3.0", description=description, usage=usage, add_help=False)
    parser.add_argument("-h", "--help", action="help",
                        help="Show this help message and exit.")

    parser.add_argument("-i", "--input_path", dest="input_path", type=str,
                        help="input dir-path within \'inbarcode_class.pkl\'")

    parser.add_argument("-o", "--output_path", dest="output_path", type=str,
                        help="path to save the result table, .txt, .csv .xlsx and so on")

    parser.add_argument("-m", "--meta_file", dest="meta", type=str,
                        help="Input reference file .txt")

    args = parser.parse_args()
    if not args.input_path:
        parser.print_help()
        sys.exit(1)
    if not args.output_path:
        parser.print_help()
        sys.exit(1)
    if not args.meta:
        parser.print_help()
        sys.exit(1)

    path = args.input_path
    save_path = args.output_path
    meta_path = args.meta

    path = st.path_diagnosis(path)

    data, Class = st.read_file(path + 'inbarcode_class.pkl')
    barcode = st.read_file(meta_path)

    index_all = get_index(data, data, Class, barcode)
    index_all.columns = ['1', 'bar1', 'bar2']
    index_all = index_all.sort_values(by='1')
    index_all.index = index_all.iloc[:, 0]
    index_all = index_all.iloc[:, 1:]
    data.iloc[:, [2, 3]] = index_all

    st.write_file(save_path, data)


if __name__ == "__main__":
    main()
