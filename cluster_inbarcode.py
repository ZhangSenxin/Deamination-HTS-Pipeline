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


# take C and T as the same
def ba_ba(chain):
    ba = []
    for i in range(len(chain)):
        if chain[i] == "A":
            ba += [1]
        elif chain[i] == "T":
            ba += [2]
        elif chain[i] == "G":
            ba += [3]
        elif chain[i] == "C":
            ba += [2]
    return ba


def xi_xi(chain1, chain2):
    chain = abs(chain1 - chain2)
    xb = sum([bool(i) for i in chain])
    return xb


def dis_matrix(dis, set_seq):
    a = 0
    dd = np.zeros([len(set_seq), len(set_seq)])
    for i in range(len(set_seq)):
        for j in range(i + 1, len(set_seq)):
            dd[i, j] = dis[a]
            dd[j, i] = dd[i, j]
            a += 1
    dd = pd.DataFrame(dd)
    return dd


def distance1(da):
    set_seq = list(da.iloc[:, 1])
    set_seq_a = [ba_ba(set_seq[i][:9]) for i in range(len(set_seq))]
    set_seq_a = pd.DataFrame(set_seq_a)
    set_seq_a = np.array(set_seq_a)

    dis = spatial.distance.pdist(set_seq_a, xi_xi)
    dd = dis_matrix(dis, set_seq)

    return dd


def distance2(da):
    set_seq = list(da.iloc[:, 1])
    print(set_seq)
    set_seq_a = [ba_ba(set_seq[i][-9:]) for i in range(len(set_seq))]
    set_seq_a = pd.DataFrame(set_seq_a)
    set_seq_a = np.array(set_seq_a)
    print(set_seq_a)

    dis = spatial.distance.pdist(set_seq_a, xi_xi)
    dd = dis_matrix(dis, set_seq)

    return dd


def draw_class(dd_, threshold=2):
    set_all = set(list(range(len(dd_))))
    set_in = set()
    set_out = set_all - set_in
    Class = []
    while len(set_out) > 0:
        set_all = set_out
        set_in = set()
        set_out = set_all - set_in
        a = list(set_out)[0]
        for i in list(dd_.loc[dd_.iloc[:, a] < threshold].index):
            set_in.add(i)
        t = 0
        while len(set_in) > t:
            t = len(set_in)
            for i in list(set_in):
                for j in list(dd_.loc[dd_.iloc[:, i] < threshold].index):
                    set_in.add(j)
        set_out = set_all - set_in
        Class += [list(set_in)]

    return Class


def main():
    usage = "usage: python %(prog)s -i input_path -o output_path"
    description = "-i -o option is needed"
    parser = argparse.ArgumentParser(prog="%prog 1.0", description=description, usage=usage, add_help=False)
    parser.add_argument("-h", "--help", action="help",
                  help="Show this help message and exit.")

    parser.add_argument("-i", "--input_path", dest="input_path", type=str,
                  help="input sequence file path.")

    parser.add_argument("-o", "--output_path", dest="output_path", type=str,
                  help="output dir-path.")

    args = parser.parse_args()
    if not args.input_path:
        parser.print_help()
        sys.exit(1)
    if not args.output_path:
        parser.print_help()
        sys.exit(1)

    path = args.input_path
    save_path = args.output_path

    save_path = st.path_diagnosis(save_path)
    st.mkdir(save_path)

    data = st.read_file(path)
    data = pd.concat([data, pd.DataFrame(None, columns=['bar1', 'bar2', 'nono'])], axis=1)
    data.loc[:, ['bar1', 'bar2', 'nono']] = 0

    dd1 = distance1(data)
    dd2 = distance2(data)

    # Using matrix corresponding element multiplication and addition,
    # lock the distance between (2,0) and (0,2), and apply a penalty of 3 t0 excluding them.
    # Requires about 10GB of RAM
    dd = dd1 + dd2
    ddd = dd1 + dd2
    dddd = dd1 * dd2

    ddd[dd != 2] = 0
    dddd[dddd == 0] = -1
    dddd[dddd != -1] = 0
    gg = ddd * dddd
    dd[gg != 0] = 3

    Class = draw_class(dd, threshold=2)
    st.write_file(save_path + 'inbarcode_class.pkl', [data, Class])
    max_num = max([len(Class[u]) for u in range(len(Class))])
    print(max_num)


if __name__ == "__main__":
    main()
