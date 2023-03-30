#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import zsx_some_tools as st
from collections import defaultdict
from openpyxl import Workbook
from openpyxl import load_workbook
import os
import sys
import argparse
import warnings
warnings.filterwarnings("ignore")


def FR():
    return 'FR'


def read_txt(path, name, time=0):
    read_list = []
    path = path + str(time) + 'min/'

    for line in open(path + name):
        line = line.strip('\n').strip('\t').split('\t')
        line = [lin.split('.0')[0] for lin in line]
        if '' in line:
            read_list += [[]]
        else:
            read_list += [np.array(list(map(int, line)))]

    return read_list[1:]


# density
def density_summary(V_reg, Pos):
    dens_summary = []
    for g in range(V_reg.shape[0]):
        gene = V_reg.index[g]
        region = [0] + list(V_reg.loc[gene])
        info = []
        for r in range(7):
            reg = [int(region[r]), int(region[r + 1])]
            pos = Pos[g][reg[0]: reg[1]]
            if np.sum(pos) == 0:
                info += [-1]
                continue
            info += [np.mean(pos)]
        dens_summary += [info]

    dens_summary = pd.DataFrame(dens_summary)
    dens_summary.index = V_reg.index
    dens_summary.columns = V_reg.columns
    dens_summary = dens_summary.fillna(-1)

    return dens_summary


# freq
def length_summary(V_reg, ref, Num_info):
    len_summary = []
    for g in range(V_reg.shape[0]):
        gene = V_reg.index[g]
        sequence = ref.loc[gene, ref.columns[0]]
        region = [0] + list(V_reg.loc[gene])
        info = []
        for r in range(7):
            seq = sequence[int(region[r]): int(region[r + 1])]
            num = Num_info[g]
            if np.sum(num) == 0:
                info += [-1]
                continue
            info += [len(seq)]
        len_summary += [info]

    len_summary = pd.DataFrame(len_summary)
    len_summary.index = V_reg.index
    len_summary.columns = V_reg.columns
    len_summary = len_summary.fillna(-1)

    return len_summary


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


def spe_produce(summary, spe_dict):
    spe_info = np.array([spe_dict[i.rsplit('_', 1)[0]] for i in summary.index])
    spe_set = list(set(spe_info))

    infomation = []
    for spe in spe_set:
        spe_data = summary.loc[spe_info == spe]
        info = []
        for r in range(7):
            spe_da = spe_data.iloc[:, [r]]
            spe_da = spe_da.loc[spe_da.iloc[:, 0] != -1]
            if spe_da.shape[0] == 0:
                info += [-1]
                continue
            info += [np.mean(spe_da)[0]]
        infomation += [info]
    infomation = pd.DataFrame(infomation)
    infomation.index = spe_set
    infomation.columns = summary.columns

    return infomation


def HM_produce(summary, spe):
    spe_info = np.array([i.rsplit('_', 1)[0] for i in summary.index])
    infomation = summary.loc[spe_info == spe]

    return infomation


def align_process(list_like_, seq_):
    out_list_ = []
    i_ = 0
    for s_ in seq_:
        if s_ == '-':
            out_list_ += [-1]
        else:
            out_list_ += [list_like_[i_]]
            i_ += 1

    return np.array(out_list_)


def align_process2(list_like_, seq_):
    out_list_ = []
    i_ = 0
    for s_ in seq_:
        if s_ == '-':
            out_list_ += [0]
        else:
            out_list_ += [list_like_[i_]]
            i_ += 1

    return np.array(out_list_)


def small_try(da_, reg_name):
    try:
        return list(da_.loc[[reg_name], 1])
    except IndexError:
        return []


def main():
    usage = "usage: python %prog -m reference_path -r region_path -f file_path -s save_path -p species_path"
    description = "data prepare, -m -r -f -s -p options are needed"
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

    parser.add_argument("-p", "--species_path", dest="species_path", type=str,
                        help="species table .txt")

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
    if not args.species_path:
        parser.print_help()
        sys.exit(1)

    reference_path = args.reference_path
    region_path = args.region_path
    path = args.file_path
    save_path = args.save_path
    species_path = args.species_path

    path = st.path_diagnosis(path)
    save_path = st.path_diagnosis(save_path)
    st.mkdir(save_path)

    full_len = pd.read_csv(reference_path, index_col=0)
    V_region = pd.read_csv(region_path, index_col=0)

    variables = ['C_Pos', 'WRC_Pos', 'AGC_Pos', 'AAC_Pos', 'TGC_Pos', 'TAC_Pos', 'AGCT_Pos', 'WGCW_Pos',
                 'A_Pos', 'T_Pos', 'G_Pos', 'Num', 'A_Mut', 'T_Mut', 'G_Mut', 'C_Mut', 'range_Num']

    full_len_class_use = full_len.reset_index()
    spe_info = st.read_file(species_path)
    spedict = {}
    for i in range(spe_info.shape[0]):
        spedict[spe_info.iloc[i, 2]] = spe_info.iloc[i, 3]
    species = np.array([spedict[i.rsplit('_', 1)[0]] for i in full_len_class_use.iloc[:, 0]])
    spe_all = list(set(species))

    time = 20

    # read info data
    for var in variables:
        code = var + ' = []'
        exec(code, {'var': var}, globals())
        mid_data = read_txt(path, var + '.txt', time=time)
        code2 = var + ' = mid_data'
        exec(code2, {'mid_data': mid_data}, globals())

    all_freq = [[] for _ in range(12)]
    for spe in spe_all:
        if spe != 'Mouse':
            index = list(full_len_class_use.loc[species == spe].index)
        else:
            # 去除mouse中带有S的IgV
            da = full_len_class_use.loc[species == 'Mouse']
            index = list(da.index[[True if 'S' not in i else False for i in da.loc[:, 'index']]])

        inter = [[] for _ in range(12)]
        for ind in index:
            name = full_len.index[ind]
            # seq = full_len.loc[name].iloc[0]

            num = Num[ind]
            WRC_pos = WRC_Pos[ind] * (num > 1000)
            AGC_pos = AGC_Pos[ind] * (num > 1000)
            AAC_pos = AAC_Pos[ind] * (num > 1000)
            TGC_pos = TGC_Pos[ind] * (num > 1000)
            TAC_pos = TAC_Pos[ind] * (num > 1000)
            C_pos = C_Pos[ind] * (num > 1000)
            T_mut = T_Mut[ind]
            C_freq = T_mut * C_pos / num

            region_ = V_region.loc[name][:7]
            CDRs = list(range(region_[1], region_[2])) + list(range(region_[3], region_[4])) + list(
                range(region_[5], region_[6]))
            reg_dict_ = defaultdict(FR)
            for i_ in CDRs:
                reg_dict_[i_] = 'CDR'
            for i_ in range(region_[0]):
                reg_dict_[i_] = 'L'

            WRC_p = np.array(range(len(WRC_pos)))[WRC_pos == 1]
            AGC_p = np.array(range(len(AGC_pos)))[AGC_pos == 1]
            AAC_p = np.array(range(len(AAC_pos)))[AAC_pos == 1]
            TGC_p = np.array(range(len(TGC_pos)))[TGC_pos == 1]
            TAC_p = np.array(range(len(TAC_pos)))[TAC_pos == 1]
            C_p = np.array(range(len(C_pos)))[C_pos == 1]

            WRC_reg = [reg_dict_[i] for i in WRC_p]
            AGC_reg = [reg_dict_[i] for i in AGC_p]
            AAC_reg = [reg_dict_[i] for i in AAC_p]
            TGC_reg = [reg_dict_[i] for i in TGC_p]
            TAC_reg = [reg_dict_[i] for i in TAC_p]
            C_reg = [reg_dict_[i] for i in C_p]

            WRC_mut = pd.DataFrame([WRC_reg, C_freq[WRC_pos == 1]]).T
            WRC_mut = WRC_mut.set_index(0)
            AGC_mut = pd.DataFrame([AGC_reg, C_freq[AGC_pos == 1]]).T
            AGC_mut = AGC_mut.set_index(0)
            AAC_mut = pd.DataFrame([AAC_reg, C_freq[AAC_pos == 1]]).T
            AAC_mut = AAC_mut.set_index(0)
            TGC_mut = pd.DataFrame([TGC_reg, C_freq[TGC_pos == 1]]).T
            TGC_mut = TGC_mut.set_index(0)
            TAC_mut = pd.DataFrame([TAC_reg, C_freq[TAC_pos == 1]]).T
            TAC_mut = TAC_mut.set_index(0)
            C_mut = pd.DataFrame([C_reg, C_freq[C_pos == 1]]).T
            C_mut = C_mut.set_index(0)

            mutation_info = [small_try(WRC_mut, 'FR'), small_try(WRC_mut, 'CDR')]
            mutation_info += [small_try(AGC_mut, 'FR'), small_try(AGC_mut, 'CDR')]
            mutation_info += [small_try(AAC_mut, 'FR'), small_try(AAC_mut, 'CDR')]
            mutation_info += [small_try(TGC_mut, 'FR'), small_try(TGC_mut, 'CDR')]
            mutation_info += [small_try(TAC_mut, 'FR'), small_try(TAC_mut, 'CDR')]
            mutation_info += [small_try(C_mut, 'FR'), small_try(C_mut, 'CDR')]

            for i_ in range(12):
                inter[i_] += mutation_info[i_]

        for i in range(12):
            all_freq[i] += [inter[i]]

    for i, name in enumerate(['WRC_FR_pool', 'WRC_CDR_pool',
                              'AGC_FR_pool', 'AGC_CDR_pool',
                              'AAC_FR_pool', 'AAC_CDR_pool',
                              'TGC_FR_pool', 'TGC_CDR_pool',
                              'TAC_FR_pool', 'TAC_CDR_pool',
                              'C_FR_pool', 'C_CDR_pool']):
        result = pd.DataFrame(all_freq[i], index=spe_27).T
        result.to_csv(save_path + name + '.txt', sep='\t', index=False)


if __name__ == "__main__":
    main()
