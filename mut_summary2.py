#!/usr/bin/env python
# coding: utf-8

import os
import pandas as pd
import numpy as np
import time
import sys
from optparse import OptionParser
import zsx_some_tools as st


def WRC_mutation(data, name):
    data = data.reset_index(drop=True)
    seq = list(data.loc[:, 'Base'])
    summary_array = []
    for j in range(int(len(seq) - 2)):
        mer3 = seq[j: j + 3]
        if (mer3[0] == 'A' or mer3[0] == 'T') and (mer3[1] == 'A' or mer3[1] == 'G') and mer3[2] == 'C':
            summary_array += [[name] + list(data.loc[j + 2, ['Pos', 'Reads', 'A', 'T', 'G']])]

    return summary_array


def SYC_mutation(data, name):
    data = data.reset_index(drop=True)
    seq = list(data.loc[:, 'Base'])
    summary_array = []
    for j in range(int(len(seq) - 2)):
        mer3 = seq[j: j + 3]
        if (mer3[0] == 'G' or mer3[0] == 'C') and (mer3[1] == 'C' or mer3[1] == 'T') and mer3[2] == 'C':
            summary_array += [[name] + list(data.loc[j + 2, ['Pos', 'Reads', 'A', 'T', 'G']])]

    return summary_array


def AGC_mutation(data, name):
    data = data.reset_index(drop=True)
    seq = list(data.loc[:, 'Base'])
    summary_array = []
    for j in range(int(len(seq) - 2)):
        if seq[j: j + 3] == list('AGC'):
            summary_array += [[name] + list(data.loc[j + 2, ['Pos', 'Reads', 'A', 'T', 'G']])]

    return summary_array


def AGCT_mutation(data, name):
    data = data.reset_index(drop=True)
    seq = list(data.loc[:, 'Base'])
    summary_array = []
    for j in range(int(len(seq) - 3)):
        if seq[j: j + 4] == list('AGCT'):
            summary_array += [[name] + list(data.loc[j + 2, ['Pos', 'Reads', 'A', 'T', 'G']])]

    return summary_array


def WGCW_mutation(data, name):
    data = data.reset_index(drop=True)
    seq = list(data.loc[:, 'Base'])
    summary_array = []
    for j in range(int(len(seq) - 3)):
        mer3 = seq[j: j + 4]
        if (mer3[0] == 'A' or mer3[0] == 'T') and mer3[1] == 'G' \
                and mer3[2] == 'C' and (mer3[3] == 'A' or mer3[3] == 'T'):
            summary_array += [[name] + list(data.loc[j + 2, ['Pos', 'Reads', 'A', 'T', 'G']])]

    return summary_array


def C_mutation(data, name):
    data = data.reset_index(drop=True)
    seq = list(data.loc[:, 'Base'])
    summary_array = []
    for j in range(len(seq)):
        if seq[j] == 'C':
            summary_array += [[name] + list(data.loc[j, ['Pos', 'Reads', 'A', 'T', 'G']])]

    return summary_array


def A_mutation(data, name):
    data = data.reset_index(drop=True)
    seq = list(data.loc[:, 'Base'])
    summary_array = []
    for j in range(len(seq)):
        if seq[j] == 'A':
            summary_array += [[name] + list(data.loc[j, ['Pos', 'Reads', 'A', 'T', 'G']])]

    return summary_array


def T_mutation(data, name):
    data = data.reset_index(drop=True)
    seq = list(data.loc[:, 'Base'])
    summary_array = []
    for j in range(len(seq)):
        if seq[j] == 'T':
            summary_array += [[name] + list(data.loc[j, ['Pos', 'Reads', 'A', 'T', 'G']])]

    return summary_array


def G_mutation(data, name):
    data = data.reset_index(drop=True)
    seq = list(data.loc[:, 'Base'])
    summary_array = []
    for j in range(len(seq)):
        if seq[j] == 'G':
            summary_array += [[name] + list(data.loc[j, ['Pos', 'Reads', 'A', 'T', 'G']])]

    return summary_array


def main():
    usage = "usage: python %prog -f file -s save"
    description = "find target pos from meta, -f -s option is needed"
    op = OptionParser(version="%prog 1.0", description=description, usage=usage, add_help_option=False)
    op.add_option("-h", "--help", action="help",
                  help="Show this help message and exit.")

    op.add_option("-m", "--meta_file", dest="meta", type="str",
                  help="Input reference file .txt")

    op.add_option("-f", "--file_path", dest="file", type="str",
                  help="data file dir-path.")

    op.add_option("-s", "--save_path", dest="save", type="str",
                  help="save dir-path.")

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
    save_path = options.save

    path = st.path_diagnosis(path)
    save_path = st.path_diagnosis(save_path)
    st.mkdir(save_path)

    label = pd.read_csv(meta_path, sep='\t')
    label = label.fillna('Na')
    file_list = os.listdir(path)

    # 建立变量，用于存储后续信息
    ATGC = ['WRC', 'SYC', 'AGC', 'AGCT', 'WGCW', 'C', 'A', 'T', 'G']
    for atgc in ATGC:
        region = ['FR1', 'FR2', 'FR3', 'CDR1', 'CDR2', 'CDR3', 'L', 'offtarget', 'exon']
        for reg in region:
            code = reg + ' = []'
            exec(code, None, globals())

        # 提取exon 和 offtarget信息
        i = 0
        for file in file_list:
            if 'name' in file:
                i += 1
                continue
            data = pd.read_csv(path + file, sep='\t')
            if 'exon' in file:
                code = 'exon += ' + atgc + '_mutation(data, file.strip(\'.txt\'))'
                exec(code, {'data': data, 'file': file}, globals())
            elif 'IGHV' in file:
                i += 1
                continue
            else:
                code = 'offtarget += ' + atgc + '_mutation(data, file.strip(\'.txt\'))'
                exec(code, {'data': data, 'file': file}, globals())
            i += 1
            print(i, '/', len(file_list), ' ' * 10, end='\r')

        # 提取IGHV信息，根据label文件，提取对应的region('FR1', 'FR2', 'FR3', 'CDR1', 'CDR2', 'CDR3', 'L')
        mismatch_file = []
        for i in range(label.shape[0]):
            file = label.iloc[i, 1] + '.txt'
            file = file.replace('*', '_')
            isExists = os.path.exists(path + file)
            if isExists:
                data = pd.read_csv(path + file, sep='\t')
            else:
                mismatch_file += [file.strip('.txt')]
                continue
            lab = label.iloc[i, 3].split('_')

            if len(lab) == 1:  # 仅在一个区域
                code = lab[0] + ' += ' + atgc + '_mutation(data, file.strip(\'.txt\'))'
                exec(code, {'data': data, 'file': file}, globals())

            else:  # 多个区域的第一个区域
                region_1 = lab[1]
                range_1 = int(lab[0].split('-')[1]) - int(lab[0].split('-')[0]) + 1
                da = data.iloc[:range_1]
                code = region_1 + ' += ' + atgc + '_mutation(da, file.strip(\'.txt\'))'
                exec(code, {'da': da, 'file': file}, globals())

                lab2 = label.iloc[i, 4].split('_')
                if len(lab2) != 1:  # 多个区域的第二个区域
                    region_2 = lab2[1]
                    range_2 = int(lab2[0].split('-')[1]) - int(lab2[0].split('-')[0]) + range_1 + 1
                    da = data.iloc[range_1:range_2]
                    code = region_2 + ' += ' + atgc + '_mutation(da, file.strip(\'.txt\'))'
                    exec(code, {'da': da, 'file': file}, globals())

                    lab3 = label.iloc[i, 5].split('_')
                    if len(lab3) != 1:  # 多个区域的第三个区域
                        region_3 = lab3[1]
                        range_3 = int(lab3[0].split('-')[1]) - int(lab3[0].split('-')[0]) + range_2 + 1
                        da = data.iloc[range_2:range_3]
                        code = region_3 + ' += ' + atgc + '_mutation(da, file.strip(\'.txt\'))'
                        exec(code, {'da': da, 'file': file}, globals())
            print(i, '/', label.shape[0], ' ' * 10, end='\r')

        # 信息赋值给各类result
        for reg in region:
            code1 = reg + ' = pd.DataFrame(' + reg + ')'
            exec(code1, {'reg': reg}, globals())

            seq_list = ['A', 'T', 'G', 'C']
            if atgc in ['WRC', 'SYC', 'AGC', 'AGCT', 'WGCW']:
                seq_list.remove('C')
            else:
                seq_list.remove(atgc)
            code2 = reg + '.columns = [\'ID\', \'Pos\', \'Reads\'] + ' + str([reg + '_' + atgc.split('_')[0] +
                                                                              '>' + nt + '_Num' for nt in seq_list])
            try:
                exec(code2, {'reg': reg, 'seq_list': seq_list, 'atgc': atgc}, globals())
            except ValueError:
                print(reg)
                exec('print(' + reg + ')', {'reg': reg}, globals())

            code3 = 'col = ' + reg + '.columns'
            exec(code3, {'reg': reg}, globals())

            for nt in range(3, 6):
                try:
                    code4 = reg + '.loc[:, col[nt].rsplit(\'_\', 1)[0] + \'_freq\'] = ' \
                            + reg + '.loc[:, col[nt]] / ' + reg + '.loc[:, \'Reads\']'
                    exec(code4, {'reg': reg, 'col': col, 'nt': nt}, globals())

                except IndexError:
                    code = reg + '= pd.DataFrame([\'\'])'
                    exec(code, {'reg': reg}, globals())
                    break

            code = reg + '.to_csv(save_path + \'' + atgc + '_' + reg + '_result.txt\', sep=\'\t\')'
            exec(code, {'save_path': save_path, 'atgc': atgc, 'reg': reg}, globals())

        print(atgc, 'result finished')

    mismatch_file = pd.DataFrame(mismatch_file)
    mismatch_file.to_csv(save_path + 'mismatch_file.txt', index=False, sep='\t')


# 运行
main()
