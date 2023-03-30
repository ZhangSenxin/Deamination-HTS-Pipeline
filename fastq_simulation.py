#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import zsx_some_tools as st
import sys
import argparse


def main():
    usage = "usage: python %(prog)s -i input_path -o output_path -m meta"
    description = "-i -o -m option is needed"
    parser = argparse.ArgumentParser(prog="%(prog)s 1.0", description=description, usage=usage, add_help=False)
    parser.add_argument("-h", "--help", action="help",
                        help="Show this help message and exit.")

    parser.add_argument("-i", "--input_path", dest="input_path", type=str,
                        help="input file path")

    parser.add_argument("-o", "--output_path", dest="output_path", type=str,
                        help="path to save the result table, .txt, .csv .xlsx and so on")

    parser.add_argument("-r", "--repeat", dest="repeat", type=int, default=10,
                        help="repeat number")

    args = parser.parse_args()
    if not args.input_path:
        parser.print_help()
        sys.exit(1)
    if not args.output_path:
        parser.print_help()
        sys.exit(1)

    path = args.input_path
    save_path = args.output_path
    repeat = args.repeat

    data = st.read_file(path, index_col=0)

    with open(save_path + 'R1.fastq', "w+") as f1, open(save_path + 'R2.fastq', "w+") as f2:
        for i in range(data.shape[0]):
            a = '@' + data.iloc[i, 0]
            b = data.iloc[i, 2] + data.iloc[i, 1] + data.iloc[i, 3]
            b2 = st.complementary(b)
            c = '+'
            d = 'J' * len(b)

            for re in range(repeat):
                f1.write(a + '\n')
                f1.write(b + '\n')
                f1.write(c + '\n')
                f1.write(d + '\n')

                f2.write(a + '\n')
                f2.write(b2 + '\n')
                f2.write(c + '\n')
                f2.write(d + '\n')


if __name__ == "__main__":
    main()
