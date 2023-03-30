#!/usr/bin/env python
# coding: utf-8

import os
import pandas as pd
import numpy as np
import zsx_some_tools as st
import sys
import argparse


# mot is a list of motif(s), even if only one motif, it also should be bagging in a list, like ['AGCT']
def find_motif(chain, mot: list, pos=0):
    chain_1 = np.zeros(len(chain))
    length = len(mot)
    for nt in range(len(chain) - length + 1):
        mer = ''.join(chain[nt: nt + length])
        if mer in mot:
            chain_1[nt + pos] = 1

    return chain_1


def find_Pos(chain):
    chain_1 = np.zeros(len(chain))
    chain_2 = np.zeros(len(chain))
    chain_3 = np.zeros(len(chain))
    chain_4 = np.zeros(len(chain))
    chain_5 = np.zeros(len(chain))
    chain_6 = np.zeros(len(chain))
    chain_7 = np.zeros(len(chain))
    chain_8 = np.zeros(len(chain))

    for nt in range(len(chain)):
        if chain[nt] == 'A':
            chain_1[nt] = 1
        elif chain[nt] == 'T':
            chain_2[nt] = 1
        elif chain[nt] == 'G':
            chain_3[nt] = 1
        elif chain[nt] == 'C':
            chain_4[nt] = 1

    for nt in range(int(len(chain) - 2)):
        mer3 = chain[nt: nt + 3]
        if (mer3[0] == 'A' or mer3[0] == 'T') and (mer3[1] == 'A' or mer3[1] == 'G') and mer3[2] == 'C':
            chain_5[nt + 2] = 1

    for nt in range(int(len(chain) - 2)):
        if chain[nt: nt + 3] == 'AGC':
            chain_6[nt + 2] = 1

    for nt in range(int(len(chain) - 3)):
        if chain[nt: nt + 4] == 'AGCT':
            chain_7[nt + 2] = 1

    for nt in range(int(len(chain) - 3)):
        mer3 = chain[nt: nt + 4]
        if (mer3[0] == 'A' or mer3[0] == 'T') and mer3[1] == 'G' \
                and mer3[2] == 'C' and (mer3[3] == 'A' or mer3[3] == 'T'):
            chain_8[nt + 2] = 1

    return chain_1, chain_2, chain_3, chain_4, chain_5, chain_6, chain_7, chain_8


def main():
    usage = "usage: python %prog -m reference_path -r region_path -f file_path -s save_path -i mismatch_path"
    description = "mapping to full-length reference, -m -r -f -s -i options are needed"
    parser = argparse.ArgumentParser(prog="%prog 1.0", description=description, usage=usage, add_help=False)
    parser.add_argument("-h", "--help", action="help",
                  help="Show this help message and exit.")

    parser.add_argument("-m", "--reference_path", dest="reference_path", type=str,
                  help="Input reference file .txt")

    parser.add_argument("-r", "--region_path", dest="region_path", type=str,
                  help="Input region file .txt")

    parser.add_argument("-f", "--file_path", dest="file_path", type=str,
                  help="data file dir-path.")

    parser.add_argument("-s", "--save_path", dest="save_path", type=str,
                  help="save dir-path.")

    parser.add_argument("-i", "--mismatch_path", dest="mismatch_path", type=str,
                  help="mismatch file .txt")

    args = parser.parse_args()
    if not args.reference_path:
        parser.print_help()
        sys.exit(1)
    if not args.region_path:
        parser.print_help()
        sys.exit(1)
    if not args.file_path:
        parser.print_help()
        sys.exit(1)
    if not args.save_path:
        parser.print_help()
        sys.exit(1)
    if not args.mismatch_path:
        parser.print_help()
        sys.exit(1)

    reference_path = args.reference_path
    region_path = args.region_path
    path = args.file_path
    out_path = args.save_path
    mis_path = args.mismatch_path
    st.mkdir(out_path)

    path = st.path_diagnosis(path)
    out_path = st.path_diagnosis(out_path)

    full_length_reference = pd.read_csv(reference_path)
    IGHV_ref = pd.read_csv(region_path, sep='\t')
    IGHV_ref.iloc[:, 0] = [IGHV_ref.iloc[i, 1].rsplit('_', 1)[1] for i in range(IGHV_ref.shape[0])]
    IGHV_ref.columns = ['Range', 'name', 'ref', 'region1', 'region2', 'region3']
    IGHV_ref = IGHV_ref.fillna('Na')

    # mask_nucleotide_num
    mask = 6
    # start_end = [0, 78, 169, 218, 282, 322, 452]
    range_list = ['1-50', '26-75', '51-100', '76-125', '101-150', '126-175', '151-200',
                  '176-225', '201-250', '226-275', '251-300', '276-325', '301-350', '326-375']

    mismatch_file = []
    for line in open(mis_path):
        mismatch_file += [line.strip('\n') + '.txt']

    C_Pos = []     # 记录C的位置信息
    WRC_Pos = []   # 记录WRC的位置信息
    AGC_Pos = []  # 记录AGC的位置信息
    AGCT_Pos = []  # 记录AGCT的位置信息
    WGCW_Pos = []  # 记录WGCW的位置信息
    A_Pos = []
    T_Pos = []
    G_Pos = []
    Num = []
    A_Mut = []
    T_Mut = []
    G_Mut = []
    C_Mut = []
    range_Num = []

    for i in range(full_length_reference.shape[0]):
        sequence = full_length_reference.iloc[i, 1]
        chain1, chain2, chain3, chain4, chain5, chain6, chain7, chain8 = find_Pos(sequence)
        A_Pos += [chain1]
        T_Pos += [chain2]
        G_Pos += [chain3]
        C_Pos += [chain4]
        WRC_Pos += [chain5]
        AGC_Pos += [chain6]
        AGCT_Pos += [chain7]
        WGCW_Pos += [chain8]

        Num += [np.zeros(len(sequence))]
        A_Mut += [np.zeros(len(sequence))]
        T_Mut += [np.zeros(len(sequence))]
        G_Mut += [np.zeros(len(sequence))]
        C_Mut += [np.zeros(len(sequence))]
        ran_Num = []
        for ran in range_list:
            da = IGHV_ref.loc[IGHV_ref.iloc[:, 0] == ran]
            bottom = int(ran.split('-')[0])-1
            top = int(ran.split('-')[1])
            seq = sequence[bottom: top]
            if sum(da.iloc[:, 2] == seq) == 0:
                continue
            file = da.loc[da.iloc[:, 2] == seq].iloc[0, 1] + '.txt'
            file = file.replace('*', '_')
            if file in mismatch_file:
                continue
            seq_data = pd.read_csv(path + file, sep='\t')
            ran_Num += [seq_data.loc[seq_data.index[0], 'Reads']]
            Num[i][bottom + mask: top - mask] += np.array(seq_data.loc[:, 'Reads'])[mask: -mask]
            A_Mut[i][bottom + mask: top - mask] += np.array(seq_data.loc[:, 'A'])[mask: -mask]
            T_Mut[i][bottom + mask: top - mask] += np.array(seq_data.loc[:, 'T'])[mask: -mask]
            G_Mut[i][bottom + mask: top - mask] += np.array(seq_data.loc[:, 'G'])[mask: -mask]
            C_Mut[i][bottom + mask: top - mask] += np.array(seq_data.loc[:, 'C'])[mask: -mask]

        range_Num += [ran_Num]
        print(i+1, '/', full_length_reference.shape[0], end='\r')

    C_Posdf = pd.DataFrame(C_Pos).fillna('')
    WRC_Posdf = pd.DataFrame(WRC_Pos).fillna('')
    AGC_Posdf = pd.DataFrame(AGC_Pos).fillna('')
    AGCT_Posdf = pd.DataFrame(AGCT_Pos).fillna('')
    WGCW_Posdf = pd.DataFrame(WGCW_Pos).fillna('')
    A_Posdf = pd.DataFrame(A_Pos).fillna('')
    T_Posdf = pd.DataFrame(T_Pos).fillna('')
    G_Posdf = pd.DataFrame(G_Pos).fillna('')
    Numdf = pd.DataFrame(Num).fillna('')
    A_Mutdf = pd.DataFrame(A_Mut).fillna('')
    T_Mutdf = pd.DataFrame(T_Mut).fillna('')
    G_Mutdf = pd.DataFrame(G_Mut).fillna('')
    C_Mutdf = pd.DataFrame(C_Mut).fillna('')
    range_Numdf = pd.DataFrame(range_Num).fillna('')

    C_Posdf.to_csv(out_path + 'C_Pos.txt', sep='\t', index=False)
    WRC_Posdf.to_csv(out_path + 'WRC_Pos.txt', sep='\t', index=False)
    AGC_Posdf.to_csv(out_path + 'AGC_Pos.txt', sep='\t', index=False)
    AGCT_Posdf.to_csv(out_path + 'AGCT_Pos.txt', sep='\t', index=False)
    WGCW_Posdf.to_csv(out_path + 'WGCW_Pos.txt', sep='\t', index=False)
    A_Posdf.to_csv(out_path + 'A_Pos.txt', sep='\t', index=False)
    T_Posdf.to_csv(out_path + 'T_Pos.txt', sep='\t', index=False)
    G_Posdf.to_csv(out_path + 'G_Pos.txt', sep='\t', index=False)
    Numdf.to_csv(out_path + 'Num.txt', sep='\t', index=False)
    A_Mutdf.to_csv(out_path + 'A_Mut.txt', sep='\t', index=False)
    T_Mutdf.to_csv(out_path + 'T_Mut.txt', sep='\t', index=False)
    G_Mutdf.to_csv(out_path + 'G_Mut.txt', sep='\t', index=False)
    C_Mutdf.to_csv(out_path + 'C_Mut.txt', sep='\t', index=False)
    range_Numdf.to_csv(out_path + 'range_Num.txt', sep='\t', index=False)


if __name__ == "__main__":
    main()
