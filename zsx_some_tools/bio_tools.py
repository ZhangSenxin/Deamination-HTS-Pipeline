# -*- coding: UTF-8 -*-
# @Author: Zhang Senxin

"""
Functions for bio-analysis
"""

from collections import defaultdict
from functools import reduce
import gzip
import os
import pandas as pd
import numpy as np


# 全排列函数
"""
Permutations list.
:x: two-layer nested list with str as fundamental element.
:code: joiner string, default is ''
:return: str list
"""
fn = lambda x, code='': reduce(lambda y, z: [str(i) + code + str(j) for i in y for j in z], x)


# PDB file related infomation
class PDBfilereader(object):
    from collections import defaultdict as __defaultdict
    from collections import namedtuple as __namedtuple

    __columns = ['ATOM', 'SPACE2', 'serial', 'SPACE1', 'name', 'altLoc', 'resName', 'SPACE1', 'chainID', 'resSeq',
                 'iCode', 'SPACE3', 'x', 'y', 'z', 'occupancy', 'tempFactor', 'SPACE6', 'segID', 'element', 'charge']

    __split = [0, 4, 6, 11, 12, 16, 17, 20, 21, 22, 26, 27,
               30, 38, 46, 54, 60, 66, 72, 76, 78, 80]

    __type = [str, str, int, str, str, str, str, str, str, int,
              str, str, float, float, float, float, float, str, str, str, str]

    __round = [False, False, False, False, False, False, False, False, False, False,
               False, False, 3, 3, 3, 2, 2, False, False, False, False]

    __direction = ['left', 'left', 'right', 'left', 'middle', 'left', 'left', 'left', 'left', 'right',
                    'left', 'left', 'right', 'right', 'right', 'right', 'right', 'left', 'left', 'right', 'left']

    __info = __namedtuple('pdb_info', ['split_info', 'type_info', 'round_info', 'direction_info'])
    __info_dict = __defaultdict(tuple)
    for __i, __col in enumerate(__columns):
        __info_use = __info(__split[__i: __i + 2], __type[__i], __round[__i], __direction[__i])
        __info_dict[__col] = __info_use

    def __init__(self, silence=False):
        if not silence:
            print('pdb_columns = pdb_file_reader.columns\n'
                  'pdb_split = pdb_file_reader.split\n'
                  'pdb_type = pdb_file_reader.typing\n'
                  'pdb_round = pdb_file_reader.rounding\n'
                  'pdb_direction = pdb_file_reader.direction\n'
                  'pdb_info_dict = pdb_file_reader.info_dict')

    @property
    def columns(self):
        return self.__columns

    @property
    def split(self):
        return self.__split

    @property
    def typing(self):
        return self.__type

    @property
    def rounding(self):
        return self.__round

    @property
    def direction(self):
        return self.__direction

    @property
    def info_dict(self):
        return self.__info_dict


# amino acid related infomation
class Aminoacid(object):
    from collections import defaultdict as __defaultdict

    __seqnum = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '*']
    __aa_dict = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"}

    __nt_dict = __defaultdict(list)
    for __nt in __aa_dict.keys():
        __nt_dict[__aa_dict[__nt]] += [__nt]

    # 三元组的命名为被排除者的后继字母
    __degeneracy = {'W': ['A', 'T'], 'S': ['G', 'C'], 'K': ['G', 'T'], 'M': ['A', 'C'], 'R': ['A', 'G'],
                    'Y': ['C', 'T'],
                    'B': ['T', 'C', 'G'], 'D': ['T', 'A', 'G'], 'H': ['T', 'C', 'A'], 'V': ['A', 'C', 'G'],
                    'N': ['A', 'C', 'T', 'G']}

    def __init__(self, silence=False):
        if not silence:
            print('aa_dict = amino_acid_info.aa_dict\n'
                  'nt_dict = amino_acid_info.nt_dict\n'
                  'degeneracy = amino_acid_info.degeneracy\n'
                  'seqnum = amino_acid_info.seqnum')

    @property
    def seqnum(self):
        return self.__seqnum

    @property
    def aa_dict(self):
        return self.__aa_dict

    @property
    def nt_dict(self):
        return self.__nt_dict

    @property
    def degeneracy(self):
        return self.__degeneracy


# PDB file read write
def read_pdb_file(path_):
    """
    Read pdb file.
    :param path_: pdb file path
    :return: pdb file
    """
    def read_pdb_line(line__, pdb_columns_, pdb_info_dict_):
        line_info_ = []
        for col__ in pdb_columns_:
            split_info_ = pdb_info_dict_[col__].split_info
            type_info_ = pdb_info_dict_[col__].type_info

            try:
                info_ = type_info_(line__[split_info_[0]: split_info_[1]].strip(' '))
            except ValueError:
                info_ = ''

            line_info_ += [info_]

        return line_info_

    pdb_file_reader = PDBfilereader(silence=True)
    pdb_columns = pdb_file_reader.columns
    pdb_info_dict = pdb_file_reader.info_dict

    data_ = []
    with open(path_, 'r+') as pdb_file_:
        for line_ in pdb_file_:
            line_ = line_.strip('\n')
            data_ += [read_pdb_line(line_, pdb_columns, pdb_info_dict)]
    data_ = pd.DataFrame(data_, columns=pdb_columns)

    return data_


def write_pdb_file(path_, data_):
    """
    Save df to pdb file (Format requirements are strict)
    :param path_: save path
    :param data_: pdb df
    :return: None
    """
    from .basic_tools import exactly_round
    def write_pdb_block(string__, value__, col__, i__, pdb_info_dict_):
        split_info = pdb_info_dict_[col__].split_info
        round_info = pdb_info_dict_[col__].round_info
        direction_info = pdb_info_dict_[col__].direction_info

        if round_info:
            try:
                value__ = exactly_round(value__, round_info)
            except ValueError:
                value__ = ''

        value__ = str(value__)
        length_exp = split_info[1] - split_info[0]
        length_true = len(value__)
        if length_true == length_exp:
            string__ += value__
        elif length_true > length_exp:
            raise ValueError('Value in row \'' + str(i__) + '\' and in col \'' +
                             col__ + '\' (\'' + value__ + '\') is too long to be set in a PDB file.')
        else:
            diff = length_exp - length_true
            if direction_info == 'right':
                value__ = ' ' * diff + value__
            elif direction_info == 'left':
                value__ = value__ + ' ' * diff
            elif direction_info == 'middle':
                value__ = ' ' + value__ + ' ' * (diff - 1)
            string__ += value__

        return string__

    pdb_file_reader = PDBfilereader(silence=True)
    pdb_info_dict = pdb_file_reader.info_dict

    with open(path_, 'w+') as pdb_file_:
        for i_ in range(data_.shape[0]):
            string_ = ''
            for j_, col_ in enumerate(data_.columns):
                value_ = data_.iloc[i_, j_]
                string_ = write_pdb_block(string_, value_, col_, i_, pdb_info_dict)

            pdb_file_.write(string_ + '\n')


# keep for those old code
def amino_acid():
    """
    Replaced by class Aminoacid()
    :return: amino acid related infomation
    """
    seqnum = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '*']
    aa_dict = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"}

    nt_dict = defaultdict(list)
    for nt in aa_dict.keys():
        nt_dict[aa_dict[nt]] += [nt]

    # 三元组的命名为被排除者的后继字母
    other_name = {'W': ['A', 'T'], 'S': ['G', 'C'], 'K': ['G', 'T'], 'M': ['A', 'C'], 'R': ['A', 'G'], 'Y': ['C', 'T'],
                  'B': ['T', 'C', 'G'], 'D': ['T', 'A', 'G'], 'H': ['T', 'C', 'A'], 'V': ['A', 'C', 'G']}

    return seqnum, aa_dict, nt_dict, other_name


def design_motif(motif_name_):
    """
    Use nucleotide degeneracy name to design motifs. For example, 'WRC' motif equal to ['AAC', 'AGC', 'TAC', 'TGC']
    :param motif_name_: str, degeneracy name
    :return: list of nucleotide motif
    """
    amino_acid_info_ = Aminoacid(silence=True)
    degeneracy_ = amino_acid_info_.degeneracy

    motif_list_ = []
    for alphabet_ in motif_name_:
        if alphabet_ in degeneracy_.keys():
            motif_list_ += [degeneracy_[alphabet_]]
        else:
            motif_list_ += [[alphabet_]]

    return fn(motif_list_)


# dna to amino acid
def translate(dna_reference, aa_dictionary, silence=False):
    """
    Translate nucleotide (DNA) sequence to amino acid (protein) sequence
    :param dna_reference: DNA cequence
    :param aa_dictionary: dict, keys are dna codons, vlaues are amino acid names.
    :param silence: print warning message or not
    :return: amino acid sequence
    """
    length = len(dna_reference)
    if length % 3 != 0:
        if silence:
            return None
        else:
            return print('DNA length can not be divided by 3.')

    aa_seq = ''
    for i in range(length // 3):
        aa_seq += aa_dictionary[dna_reference[3 * i: 3 * (i + 1)]]

    return aa_seq


def find_motif(sequence_, motif_, mot_p=2):
    """
    To find motif positions in such sequence
    :param sequence_: str, using sequence
    :param motif_: list, motif that designed by design_motif()
    :param mot_p: position start with the i_th nucleotide of the motif
    :return: list of positions that exist the motif
    """
    result_ = []
    mot_l = len(motif_[0])
    for i_ in range(len(sequence_) - len(motif_[0]) + 1):
        mot_ = sequence_[i_: i_ + mot_l]
        if mot_ in motif_:
            result_ += [i_ + mot_p]

    return result_


def count_motif(sequence_, motif_):
    """
    To calculate the number of motif in such sequence
    :param sequence_: str, using sequence
    :param motif_: list, motif that designed by design_motif()
    :return: int, motif number
    """
    num_ = 0
    mot_l = len(motif_[0])
    for i_ in range(len(sequence_) - len(motif_[0]) + 1):
        mot_ = sequence_[i_: i_ + mot_l]
        if mot_ in motif_:
            num_ += 1

    return num_


def extract_motif(sequence_, motif_, extract_1=6, extract_2=False):
    """

    :param sequence_:
    :param motif_:
    :param n_:
    :param extract_1:
    :param extract_2:
    :return:
    """
    start = extract_1
    end = int(len(sequence_)) - (extract_2 + n_ - 1)
    n_ = len(motif_[0])

    result_ = []
    for i_ in range(start, end):
        mot = sequence_[i_: i_ + n_]
        if mot in motif_:
            seq_out = sequence_[i_ - extract_1: i_ + extract_2 + n_]
            seq_left = seq_out[: extract_1]
            seq_mid = seq_out[extract_1: extract_1 + n_]
            seq_end = seq_out[extract_1 + n_:]

            result_ += [[seq_left, seq_mid, seq_end]]

    return result_


def read_fasta(path_, split=None):
    """
    Read fasta file as a df
    :param path_: fasta file path
    :param split: id info will split by the given separation and keep the front part
    :return: df, index is id, and only one column is sequence
    """

    i_ = 0
    result_ = []
    fasta2 = ''
    for line_ in open(path_):
        line_ = line_.strip("\n").strip("\r")
        if '>' in line_:
            if i_:
                result_ += [[fasta1, fasta2]]
                fasta2 = ''
            fasta1 = line_.split('>', 1)[-1].split(split)[0]   # need more parameter
            i_ = 1
        else:
            fasta2 += line_

    result_ += [[fasta1, fasta2]]
    result_ = pd.DataFrame(result_, columns=['ID', 'seq'])
    result_ = result_.set_index('ID')

    return result_


def write_fasta(path_, data__):
    """
    Save fasta file
    :param path_: save path
    :param data_: df
    :return: None
    """
    if data__.shape[1] == 1:
        data_ = data__.reset_index()
    elif data__.shape[1] > 2:
        data_ = data__.iloc[:, :2]
    else:
        data_ = data__
    with open(path_, 'w+') as file_:
        for i in range(data_.shape[0]):
            file_.write('>' + str(data_.iloc[i, 0]) + '\t')
            file_.write(str(data_.iloc[i, 1]) + '\t')


def read_fastq(path_, split=None):
    """
    Read fastq file as a df
    :param path_: fastq file path
    :param split: id info will split by the given separation and keep the front part
    :return: df, index is id, and three columns for sequence, addition info, and sequencing quality
    """
    fastq1 = []
    fastq2 = []
    fastq3 = []
    fastq4 = []

    for i_, line_ in enumerate(open(path_)):
        line_ = line_.strip("\n").strip("\r")
        if i_ % 4 == 0:
            line_ = line_.split(sep='@')[-1].split(split)[0]   # more parameter
            fastq1 += [line_]
        if i_ % 4 == 1:
            fastq2 += [line_]
        if i_ % 4 == 2:
            fastq3 += [line_]
        if i_ % 4 == 3:
            fastq4 += [line_]

    fastq1 = pd.DataFrame(fastq1)
    fastq2 = pd.DataFrame(fastq2)
    fastq3 = pd.DataFrame(fastq3)
    fastq4 = pd.DataFrame(fastq4)

    fastq = pd.concat([fastq1, fastq2, fastq3, fastq4], axis=1)
    fastq.columns = ['ID', 'seq', '+', 'qua']
    fastq = fastq.set_index('ID')

    return fastq


def write_fastq(path_, data__):
    """
    Save fastq file
    :param path_: save path
    :param data_: df
    :return: None
    """
    if data__.shape[1] == 3:
        data_ = data__.reset_index()
    elif data__.shape[1] > 4:
        data_ = data__.iloc[:, :4]
    else:
        data_ = data__
    with open(path_, 'w+') as file_:
        for i in range(data_.shape[0]):
            file_.write('@' + str(data_.iloc[i, 0]) + '\t')
            file_.write(str(data_.iloc[i, 1]) + '\t')
            file_.write(str(data_.iloc[i, 2]) + '\t')
            file_.write(str(data_.iloc[i, 3]) + '\t')


def read_fastq_gz(path, split='.', num_limit=0, decode='utf-8'):
    """
    Read fastq.gz file as a df
    :param path_: fastq file path
    :param split: id info will split by the given separation and keep the front part
    :param num_limit: stop reading at which line
    :param decode: decode parameter
    :return: df, index is id, and three columns for sequence, addition info, and sequencing quality
    """
    fastq1 = []
    fastq2 = []
    fastq3 = []
    fastq4 = []
    i = 0
    if num_limit == 0:
        num_limit = 99999999
    with gzip.open(path, 'r') as file_:
        for line_ in file_:
            line_ = line_.decode(decode).strip('\n').strip("\r")
            if ('\x00' in line_) and i > 0:
                continue
            if i % 4 == 0:
                line_ = line_.split(sep='@')[-1].split(sep=split)[0]   # more parameter
                fastq1 += [line_]
            elif i % 4 == 1:
                fastq2 += [line_]
            elif i % 4 == 2:
                fastq3 += [line_]
            else:
                fastq4 += [line_]
            i += 1
            if i == num_limit * 4:
                break

    fastq1 = pd.DataFrame(fastq1)
    fastq2 = pd.DataFrame(fastq2)
    fastq3 = pd.DataFrame(fastq3)
    fastq4 = pd.DataFrame(fastq4)

    fastq = pd.concat([fastq1, fastq2, fastq3, fastq4], axis=1)
    fastq.columns = ['ID', 'seq', '+', 'qua']
    fastq = fastq.set_index('ID')

    return fastq


def write_fastq_gz(path_, data_):
    """
    Save fasta.gz file
    :param path_: save path
    :param data_: df
    :return: None
    """
    with gzip.open(path_, 'w+') as file_:
        for i in range(data_.shape[0]):
            file_.write(('@' + str(data_.iloc[i, 0]) + '\t').encode())
            file_.write((str(data_.iloc[i, 1]) + '\t').encode())
            file_.write((str(data_.iloc[i, 2]) + '\t').encode())
            file_.write((str(data_.iloc[i, 3]) + '\t').encode())


def read_fastq_gz_onlyseq(path_, split='.', num_limit=False, decode=None):
    """
    Read fastq.gz file as a df only contain id and sequence
    :param path_: fastq file path
    :param split: id info will split by the given separation and keep the front part
    :param num_limit: stop reading at which line.
    :param decode: decode parameter
    :return:
    """
    fastq1 = []
    fastq2 = []
    i = 0

    with gzip.open(path_, 'r') as f1:
        # con = f1.readlines()
        for line in f1:
            line = line.decode(decode).strip('\n').strip("\r")
            if ('\x00' in line) and i > 0:
                continue
            if i % 4 == 0:
                line = line.split(sep='@')[-1].split(sep=split)[0]  # more parameter
                fastq1 += [line]
            elif i % 4 == 1:
                fastq2 += [line]

            i += 1
            if i == num_limit * 4:
                break

    fastq1 = pd.DataFrame(fastq1)
    fastq2 = pd.DataFrame(fastq2)

    fastq = pd.concat([fastq1, fastq2], axis=1)
    fastq.columns = ['ID', 'seq']
    fastq = fastq.set_index('ID')

    return fastq


def complementary(sequence, reverse=True):
    """
    Reverse and complementary the sequence
    :param sequence: sequence
    :param reverse: reverse or not
    :return: processed sequence
    """
    if reverse:
        sequence = sequence[::-1]
    cha = ""
    for c in range(len(sequence)):
        if sequence[c] == "A":
            cha += "T"
        elif sequence[c] == "T":
            cha += "A"
        elif sequence[c] == "C":
            cha += "G"
        else:
            cha += "C"

    return cha
