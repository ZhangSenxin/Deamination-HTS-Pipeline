#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import gzip
import time
import sys
import os
from fastq_alignment import *
from optparse import OptionParser
import some_tools as st

nt = ['A', 'T', 'G', 'C']
one_hot = {}
auto_index = {}
for i in range(len(nt)):
    aa = [0]*len(nt)
    aa[i] = 1
    one_hot[nt[i]] = aa
    auto_index[nt[i]] = i


# one-hot coding
def hot(chain):
    ba = []
    for ab in chain:
        ba += one_hot[ab]
    return ba


# id coding
def ba_ba(chain):
    ba = []
    for seq in range(len(chain)):
        if chain[seq] == "A":
            ba += [1]
        elif chain[seq] == "T":
            ba += [2]
        elif chain[seq] == "G":
            ba += [3]
        elif chain[seq] == "C":
            ba += [4]
    return ba


# one-hot coding
def mut(chain):
    mutation_sum = []
    for mu in chain:
        res = [0, 0, 0, 0]
        if mu == 1:
            res[0] = 1
        elif mu == 2:
            res[1] = 1
        elif mu == 3:
            res[2] = 1
        elif mu == 4:
            res[3] = 1  
        mutation_sum += [res]
    return np.array(mutation_sum)


def main():
    usage = "usage: python %prog -m meta -f file -s save"
    description = "split reads according to reference, -m -f -s option is needed"
    op = OptionParser(version="%prog 1.0", description=description, usage=usage, add_help_option=False)
    op.add_option("-h", "--help", action="help",
                  help="Show this help message and exit.")

    op.add_option("-m", "--meta_file", dest="meta", type="str",
                  help="Input reference file .txt")

    op.add_option("-f", "--file_path", dest="file", type="str",
                  help="data file dir-path.")

    op.add_option("-s", "--save_path", dest="save", type="str",
                  help="save dir-path.")

    op.add_option("-c", "--chain", dest="chain", type="int", default=1,
                  help="+/- chain type, 1 is + chain, 0 is - chain")

    (options, args) = op.parse_args()
    if not options.meta:
        op.print_help()
        sys.exit(1)
    if not options.file:
        op.print_help()
        sys.exit(1)
    if not options.save:
        op.print_help()
        sys.exit(1)

    meta_path = options.meta
    path = options.file
    out_path = options.save
    chain_type = options.chain

    path = st.path_diagnosis(path)
    out_path = st.path_diagnosis(out_path)

    reference = pd.read_csv(meta_path, sep='\t', index_col=0)
    # double strand: 6 barcode(3**6, C==T) and 9 in-barcode
    ref_barcode = np.array([hot(chain[:15]+chain[-15:]) for chain in reference.iloc[:, 1]])
    dir_list = os.listdir(path)

    st.mkdir(out_path)

    file_num = 0
    t = []  # save the time each step used
    name_file_list = []

    # 进入不同文件夹提取序列信息
    for direction in dir_list:
        if direction[:3] != 'cut':
            continue
        path1 = path + direction + '/'
        file_list = os.listdir(path1)
        for file in file_list:
            if 'R2' in file:
                continue
            t0 = time.time()

            file_name = file.rsplit('_', 1)[0]
            i = 0
            out = []
            with gzip.open(path1 + file, 'r') as f1:
                try:
                    con = f1.readlines()
                except:
                    return print(file, 'error')
                for line in con:
                    if i % 4 == 1:  # 只提取序列信息
                        line = line.decode('utf-8').strip('\n')
                        if 'N' in line:  # 含有N的极少，直接跳过
                            i += 1
                            continue
                        out += [line]
                    i += 1
            if len(out) == 0:  # 空文件
                result = np.zeros([50, 9])
                result = pd.DataFrame(result)
                result.columns = ['Pos', 'Base', 'Reads', 'Subs', 'Y', 'A', 'T', 'G', 'C']
                result.to_csv(out_path + 'NF_' + file_name + '.txt', index=False, sep='\t')
                name_file_list += [['Not Found', file_name]]
                continue

            # 计算序列类别数
            s_out = list(set(out))
            countlist = [[seq, out.count(seq)] for seq in s_out]
            countlist = pd.DataFrame(countlist)
            countlist = countlist.sort_values(by=1, ascending=False)

            # find reference
            s_out0 = np.array([hot(s_out[test][:15] + s_out[test][-15:]) for test in range(len(s_out))])
            K = np.array([np.sum(abs(ref_barcode - s_out0[test]), axis=1) for test in range(len(s_out0))])
            index = np.argmin(np.sum(K, axis=0))
            ref = reference.iloc[index, 1][6:-6]
            ref_name = reference.iloc[index, 0]
            ref_name = ref_name.replace('*', '_')

            # find mutation
            result = np.zeros([50, 9])
            result = pd.DataFrame(result)
            result.columns = ['Pos', 'Base', 'Reads', 'Subs', 'Y', 'A', 'T', 'G', 'C']
            result.iloc[:, 0] = range(1, 51)
            result.iloc[:, 1] = list(ref)
            result.iloc[:, 2] = len(out)

            # 对比mutation并计算
            chain1 = np.array(ba_ba(ref))
            for seq in range(countlist.shape[0]):
                counts = countlist.iloc[seq, 1]
                chain2 = np.array(ba_ba(countlist.iloc[seq, 0][6:-6]))
                if len(chain2) != 50:
                    continue
                chain = chain1 - chain2
                if sum(abs(chain)) == 0:
                    continue
                chain[chain != 0] = chain2[chain != 0]
                result.iloc[:, -4:] = result.iloc[:, -4:] + mut(chain) * counts
            result.loc[:, 'Subs'] = np.sum(result.iloc[:, -4:], axis=1)
            result.loc[:, 'Y'] = np.round(result.loc[:, 'Subs'] / result.loc[:, 'Reads'], 5)

            name_file_list += [[ref_name, file_name]]

            if not chain_type:
                result = result.loc[::-1]
            result.to_csv(out_path + ref_name + '.txt', index=False, sep='\t')

            t1 = time.time()         # a naive way to estimate time needed
            file_num += 1
            t += [t1 - t0]
            T = np.mean(t) * (len(reference)-file_num)
            S = int(T % 60)
            M = int(T/60)
            txt = str(file_num) + '/' + str(len(reference)) + ' ' + str(M) + 'min' + str(S) + 's'
            print(txt, end='\r')

    name_file_list = pd.DataFrame(name_file_list)
    name_file_list.to_csv(out_path + 'name_file_list.txt', index=False, header=False, sep='\t')


# 运行
main()
