#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import os
import zsx_some_tools as st
from optparse import OptionParser
import sys
import warnings
warnings.filterwarnings("ignore")


def read_txt(path_, name_, time_=0):
    read_list = []
    path_ = path_ + str(time_) + 'min/'

    for line in open(path_ + name_):
        line = line.strip('\n').strip('\t').split('\t')
        line = [lin.split('.0')[0] for lin in line]
        if '' in line:
            read_list += [[]]
        else:
            read_list += [np.array(list(map(int, line)))]

    return read_list[1:]


def main():
    usage = "usage: python %prog -m meta -i pathin -o pathout"
    description = "combine V, -m -i -o option is needed."
    op = OptionParser(version="%prog 1.3", description=description, usage=usage, add_help_option=False)
    op.add_option("-h", "--help", action="help",
                  help="Show this help message and exit.")

    op.add_option("-m", "--meta", dest="meta", type="str", default=0,
                  help="meta file.")

    op.add_option("-i", "--inputpath", dest="pathin", type="str",
                  help="The input dir-path.")

    op.add_option("-o", "--outputpath", dest="pathout", type="str",
                  help="The output dir-path.")

    (options, args) = op.parse_args()
    if not options.meta:
        op.print_help()
        sys.exit(1)
    if not options.pathin:
        op.print_help()
        sys.exit(1)
    if not options.pathout:
        op.print_help()
        sys.exit(1)

    meta_path = options.meta
    path = options.pathin
    save_path = options.pathout

    path = st.path_diagnosis(path)
    save_path = st.path_diagnosis(save_path)
    st.mkdir(save_path)

    all_len = pd.read_csv(meta_path, index_col=0)
    variables = ['C_Pos', 'WRC_Pos', 'AGC_Pos', 'AGCT_Pos', 'WGCW_Pos', 'A_Pos', 'T_Pos', 'G_Pos',
                 'Num', 'A_Mut', 'T_Mut', 'G_Mut', 'C_Mut', 'range_Num']

    folders = st.listdir(path)
    for folder in folders:
        time = int(folder.strip('min'))
        save_path_use = save_path + folder + '/'
        st.mkdir(save_path_use)

        # read info data
        for var in variables:
            code = var + ' = []'
            exec(code, {'var': var}, globals())
            mid_data = read_txt(path, var + '.txt', time_=time)
            code2 = var + ' = mid_data'
            exec(code2, {'mid_data': mid_data}, globals())

        for i, name in enumerate(all_len.index):
            seq = all_len.loc[name, '1']
            Base = pd.DataFrame(list(seq))

            A_mut = A_Mut[i]
            T_mut = T_Mut[i]
            G_mut = G_Mut[i]
            C_mut = C_Mut[i]
            num = Num[i]

            info = pd.DataFrame([A_mut, T_mut, G_mut, C_mut]).T
            subs = np.sum(info, axis=1)
            Y = pd.DataFrame(subs / num)
            num = pd.DataFrame(num)
            subs = pd.DataFrame(subs)
            pos = pd.DataFrame(list(range(1, len(seq) + 1)))

            infomation = pd.concat([pos, Base, num, subs, Y, info], axis=1)
            infomation.columns = ['Pos', 'Base', 'Reads', 'Subs', 'Y', 'A', 'T', 'G', 'C']
            infomation.to_csv(save_path_use + name + '.txt', index=False, sep='\t')


if __name__ == "__main__":
    main()
