import glob
import re
import sys
import unittest
from collections import Counter, defaultdict
from itertools import combinations, product

import pysam
from xopen import xopen

from celescope.tools import utils
from celescope.tools.__init__ import PATTERN_DICT
from celescope.__init__ import ROOT_PATH, HELP_DICT
from celescope.tools.step import Step, s_common
import json
import pandas as pd

MIN_T = 10

bp=dict(B="C/G/T",
D="A/G/T",
H="A/C/T",
K="G/T",
M="A/C",
N="A/C/G/T",
R="A/G",
S="C/G",
V="A/C/G",
W="A/T",
Y="C/T",)

BP={key:val.split("/") for key,val in bp.items()}

def get_probe_mismatch(seq):
    """
    this is probe function
    """
    seq_set=set()
    mismatch_dict={}
    seq = seq.upper()
    seq_list = [[x] for x in list(seq)]
    _locs = {}
    for key,_ in BP.items():
        loc = [each.start() for each in re.finditer(key, seq)]
        for x in loc:
            _locs[x] = key
    for key,val in _locs.items():
        seq_list[key] = BP[val]  
    for pos in product(*seq_list):
        seq_set.add(''.join(pos))
    for mismatch_seq in seq_set:
        mismatch_dict[mismatch_seq]=seq
    return mismatch_dict

def get_probe_all_mismatch(seq_list):
    """
    this is probe function
    """
    mismatch_dict={}
    correct_set=set()
    for seq in seq_list:
        seq = seq.upper()
        seq=seq.strip()
        if seq.startswith('>'):
            continue
        correct_set.add(seq)
        mismatch_dict.update(get_probe_mismatch(seq))
    return correct_set,mismatch_dict

def findall_mismatch(seq, n_mismatch=1, bases='ACGTN'):
    """
    choose locations where there's going to be a mismatch using combinations
    and then construct all satisfying lists using product

    Return:
    all mismatch <= n_mismatch set. 

    >>> answer = set(["TCG", "AAG", "ACC", "ATG", "ACT", "ACN", "GCG", "ANG", "ACA", "ACG", "CCG", "AGG", "NCG"])
    >>> seq_set = findall_mismatch("ACG")
    >>> seq_set == answer
    True
    """
    seq_set = set()
    seq_len = len(seq)
    if n_mismatch > seq_len:
        n_mismatch = seq_len
    for locs in combinations(range(seq_len), n_mismatch):
        seq_locs = [[base] for base in seq]
        for loc in locs:
            seq_locs[loc] = list(bases)
        for poss in product(*seq_locs):
            seq_set.add(''.join(poss))
    return seq_set

def hamming_fac(x,y,mismatch=3):
    """
    Calculate the Hamming distance
    """
    dis = 0
    for index,i in enumerate(x):
        if i != y[index]:
            dis += 1
            if dis > mismatch:
                return dis
    return dis

@utils.add_log
def creat_probe_file(probefile,sample):
    probe_list,_ = utils.read_one_col(probefile)
    probe_list[1::2] = [x.upper() for x in probe_list[1::2]]
    probe_dict = dict(zip(probe_list[1::2],probe_list[::2]))
    _,probe_mismatch_dict = get_probe_all_mismatch(probe_list)
    _probe = pd.DataFrame([probe_mismatch_dict])
    df = _probe.T
    df.columns = ['probe']
    df['seq'] = df.index
    df.reset_index(inplace=True)
    df.drop(labels='index',axis=1,inplace=True)
    df.probe=df.probe.map(probe_dict)
    s = df.groupby('probe').cumcount().add(1).astype(str)
    df['probe'] += ('_' + s).replace('_1', '')
    with open(f"./{sample}_probe.fasta","w") as f:
        for i in range(df.shape[0]):
            f.write('%s\n%s\n' % (df.probe[i],df.seq[i]))
    return f"./{sample}_probe.fasta"

def read_fasta(fasta_file, equal=False):
    # seq must have equal length
    dicts = {}
    LENGTH = 0
    with open(fasta_file, "rt") as f:
        while True:
            line1 = f.readline()
            line2 = f.readline()
            if not line2:
                break
            name = line1.strip(">").strip()
            seq = line2.strip()
            seq_length = len(seq)
            if equal:
                if LENGTH == 0:
                    LENGTH = seq_length
                else:
                    if LENGTH != seq_length:
                        raise Exception(f"{fasta_file} have different seq length")
            dicts[name] = seq
    if equal:
        return dicts, LENGTH
    return dicts

def merge_dicts(dict1, dict2):
    merged_dict = {}
    n = 0 
    for key in set(dict1.keys()).union(dict2.keys()):
        ##print(n)
        ##print(key)
        if key in dict1 and key in dict2:
            tmp_dict = {}
            list_of_dicts = [dict1.get(key), dict2.get(key)]
            #print(list_of_dicts)
            for single_dict in list_of_dicts:
                #print(single_dict)
                for k, v in single_dict.items():
                    if k in tmp_dict:
                        tmp_dict[k] = int(tmp_dict[k]) + int(v)
                    else:
                        tmp_dict[k] = v
            
            merged_dict[key] = tmp_dict
            #print(merged_dict)
            #break
        elif key in dict1:
            merged_dict[key] = dict1[key]
        else:
            merged_dict[key] = dict2[key]
        ##print(merged_dict)
        ##print("--")
        ##n = n+1
        ##if n>10:
        ##    break
    return merged_dict

# check
#dict1 = {"t1": {"t11":1},"t2":{'t21': 1}}
#dict2 = {"t2":{'t21': 1, "a":1},"t3":{"t31":3}}
#merged_dict = merge_dicts(dict1, dict2)
#merged_dict


def merge_dict_probe(data, dict_update, probe_names, probe):
    probe_sub = [item for item in probe_names if probe in item]
    #print(probe_sub)
    dict_tmp = {}
    for item in probe_sub:
        #print(item)
        dict_tmp = merge_dicts(dict_tmp, data[item])
    dict_update[probe] = dict_tmp


class Chemistry():
    """
    Auto detect chemistry from R1-read
    """

    def __init__(self, fq1):
        '''
        'scopeV2.0.1': 'C8L16C8L16C8L1U8T18'
        'scopeV2.1.1': 'C8L16C8L16C8L1U12T18'
        'scopeV2.2.1': 'C8L16C8L16C8L1U12T18' with 4 types of linkers
        'scopeV3.0.1': 'C9L16C9L16C9L1U12T18' with 4 types of linkers
        '''
        self.fq1 = fq1
        self.fq1_list = fq1.split(',')
        self.n_read = 10000

        self.pattern_dict_v2, * \
            _, self.linker_1_v2_set_list, self.linker_1_v2_mismatch_list = Barcode.parse_chemistry('scopeV2.1.1')
        self.pattern_dict_v2, * \
            _, self.linker_4_v2_set_list, self.linker_4_v2_mismatch_list = Barcode.parse_chemistry('scopeV2.2.1')
        self.pattern_dict_v3, *_, self.linker_v3_set_list, self.linker_v3_mismatch_list = Barcode.parse_chemistry('scopeV3.0.1')
        self.pattern_dict_flv, *_, self.linker_flv_set_list, self.linker_flv_mismatch_list = Barcode.parse_chemistry('flv')
        self.pattern_dict_flv_rna, *_, self.linker_flv_rna_set_list, self.linker_flv_rna_mismatch_list = Barcode.parse_chemistry('flv_rna')


    @utils.add_log
    def check_chemistry(self):
        """check chemistry in the fq1_list"""
        chemistry_list = []
        for fastq1 in self.fq1_list:
            print(fastq1)
            chemistry = self.get_chemistry(fastq1)
            chemistry_list.append(chemistry)
        if len(set(chemistry_list)) != 1:
            Chemistry.check_chemistry.logger.warning('multiple chemistry found!' + str(chemistry_list))
        return chemistry_list


    def seq_chemistry(self, seq):
        """
        Returns: chemistry or None

        >>> runner = Chemistry("fake_fq1_string")
        >>> seq = "TCGACTGTCATCCACGTGCTTGAGATTCTAGGATTCAGCATGCGGCTACGTGCACGAGACATATCAATGGGTTTTCTTGTTGCTTTTTTTTTTTTTTTTTTTTTTTT"
        >>> runner.seq_chemistry(seq)
        'scopeV3.0.1'

        >>> seq = "GTCGTAGAATCCACGTGCTTGAGACTCAATGATCAGCATGCGGCTACGGCGATTAACGTTGAATGTTTTTTTTTTTTTTTTTTTTT"
        >>> runner.seq_chemistry(seq)
        'scopeV2.0.1'

        >>> seq = "NCAGATTC" + "ATCCACGTGCTTGAGA" + "GTACGCAA" + "TCAGCATGCGGCTACG" + "CTGAGCCA" + "C" + "TCCGAAGCCCAT" + "TTTTTTTTTTTTTTTTTTTTTTTTTTATTGC"
        >>> runner.seq_chemistry(seq)
        'scopeV2.1.1'

        >>> seq = "NCAGATTC" + "TCGGTGACAGCCATAT" + "GTACGCAA" + "CGTAGTCAGAAGCTGA" + "CTGAGCCA" + "C" + "TCCGAAGCCCAT" + "TTTTTTTTTTTTTTTTTTTTTTTTTTATTGC"
        >>> runner.seq_chemistry(seq)
        'scopeV2.2.1'

        """

        # check flv_rna first. otherwise it may be considered as scopeV2.1.1 and scopeV2.2.1
        linker_flv_rna = Barcode.get_seq_str(seq, self.pattern_dict_flv_rna["L"])
        bool_valid, _, _ = Barcode.check_seq_mismatch(
            [linker_flv_rna], self.linker_flv_rna_set_list, self.linker_flv_rna_mismatch_list)
        if bool_valid:
            return "flv_rna"


        linker_v2 = Barcode.get_seq_str(seq, self.pattern_dict_v2["L"])
        bool_valid, _, _ = Barcode.check_seq_mismatch(
            [linker_v2], self.linker_1_v2_set_list, self.linker_1_v2_mismatch_list)
        if bool_valid:
            if seq[65:69] == "TTTT":
                return "scopeV2.0.1"
            else:
                return "scopeV2.1.1"

        linker_v3 = Barcode.get_seq_str(seq, self.pattern_dict_v3["L"])
        bool_valid, _, _ = Barcode.check_seq_mismatch(
            [linker_v3], self.linker_v3_set_list, self.linker_v3_mismatch_list)
        if bool_valid:
            return "scopeV3.0.1"

        linker_v2 = Barcode.get_seq_str(seq, self.pattern_dict_v2["L"])
        bool_valid, _, _ = Barcode.check_seq_mismatch(
            [linker_v2], self.linker_4_v2_set_list, self.linker_4_v2_mismatch_list)
        if bool_valid:
            return "scopeV2.2.1"

        linker_flv = Barcode.get_seq_str(seq, self.pattern_dict_flv["L"])
        bool_valid, _, _ = Barcode.check_seq_mismatch(
            [linker_flv], self.linker_flv_set_list, self.linker_flv_mismatch_list)
        if bool_valid:
            return "flv"

        return

    @utils.add_log
    def get_chemistry(self, fq1):
        results = defaultdict(int)

        with pysam.FastxFile(fq1) as fh:
            for _ in range(self.n_read):
                entry = fh.__next__()
                seq = entry.sequence
                chemistry = self.seq_chemistry(seq)
                if chemistry:
                    results[chemistry] += 1
        # if it is 0, then no other linker types
        if results["scopeV2.2.1"] != 0:
            results["scopeV2.2.1"] += results["scopeV2.1.1"]
        sorted_counts = sorted(results.items(), key=lambda x: x[1], reverse=True)
        self.get_chemistry.logger.info(sorted_counts)

        chemistry, read_counts = sorted_counts[0][0], sorted_counts[0][1]
        percent = float(read_counts) / self.n_read
        if percent < 0.5:
            self.get_chemistry.logger.warning("Valid chemistry read counts percent < 0.5")
        if percent < 0.1:
            self.get_chemistry.logger.error("Valid chemistry read counts percent < 0.1")
            raise Exception(
                'Auto chemistry detection failed! ' + HELP_DICT['chemistry']
            )
        Chemistry.get_chemistry.logger.info(f'chemistry: {chemistry}')

        return chemistry


class Barcode(Step):
    """
    ## Features

    - Demultiplex barcodes.
    - Filter invalid R1 reads, which includes:
        - Reads without linker: the mismatch between linkers and all linkers in the whitelist is greater than 2.  
        - Reads without correct barcode: the mismatch between barcodes and all barcodes in the whitelist is greater than 1.  
        - Reads without polyT: the number of T bases in the defined polyT region is less than 10.
        - Low quality reads: low sequencing quality in barcode and UMI regions.

    ## Output

    - `01.barcode/{sample}_2.fq(.gz)` Demultiplexed R2 reads. Barcode and UMI are contained in the read name. The format of 
    the read name is `{barcode}_{UMI}_{read ID}`.
    """

    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        self.fq1_list = args.fq1.split(",")
        self.fq2_list = args.fq2.split(",")
        self.fq_number = len(self.fq1_list)
        if self.fq_number != len(self.fq2_list):
            raise Exception('fastq1 and fastq2 do not have same file number!')
        if args.chemistry == 'auto':
            ch = Chemistry(args.fq1)
            self.chemistry_list = ch.check_chemistry()
        else:
            self.chemistry_list = [args.chemistry] * self.fq_number
        self.barcode_corrected_num = 0
        self.linker_corrected_num = 0
        self.total_num = 0
        self.clean_num = 0
        self.no_polyT_num = 0
        self.lowQual_num = 0
        self.no_linker_num = 0
        self.no_barcode_num = 0
        self.barcode_qual_Counter = Counter()
        self.umi_qual_Counter = Counter()
        self.pattern = args.pattern
        self.linker = args.linker
        self.whitelist = args.whitelist
        self.lowNum = args.lowNum
        self.lowQual = args.lowQual
        self.filterNoPolyT = args.filterNoPolyT
        self.allowNoLinker = args.allowNoLinker
        self.nopolyT = args.nopolyT  # true == output nopolyT reads
        self.noLinker = args.noLinker
        self.output_R1 = args.output_R1
        self.bool_flv = False
        self.probefile=args.probe_file
        self.probefilemode=args.probe_file_mode
        self.ampfile=args.amp_file

        self._assay = self.get_slot_key(
            slot='metrics',
            step_name='sample',
            key='Assay',
        )
        self._assay = self._assay.split(' ')[-1]
        # flv_trust4, flv_CR
        if 'flv' in self._assay:
            self.bool_flv = True
            self.barcode_read_Counter = Counter()
            if self._assay == 'flv_trust4':
                if args.match_dir == 'None':
                    raise FileNotFoundError('Match_dir required when running flv_trust4')
                self.match_barcodes = set(utils.get_barcode_from_match_dir(args.match_dir)[0]) # barcode set of flv_rna.
                self.match_num = 0 # record read number match with flv_rna.
                self.match_cbs = set() # record barcode number match with flv_rna.

        # out file
        if args.gzip:
            self.suffix = ".gz"
        else:
            self.suffix = ""
        self.out_fq2 = f'{self.out_prefix}_2.fq{self.suffix}'
        self.out_fq1 = f'{self.out_prefix}_1.fq{self.suffix}'
        if self.nopolyT:
            self.nopolyT_1 = f'{self.out_prefix}_noPolyT_1.fq'
            self.nopolyT_2 = f'{self.out_prefix}_noPolyT_2.fq'
        if self.noLinker:
            self.noLinker_1 = f'{self.out_prefix}_noLinker_1.fq'
            self.noLinker_2 = f'{self.out_prefix}_noLinker_2.fq'

        self.open_files()

    @staticmethod
    def get_seq_str_no_exception(seq, sub_pattern_dict):
        """get subseq with intervals in arr and concatenate"""
        return ''.join([seq[item[0]: item[1]] for item in sub_pattern_dict])

    @staticmethod
    def get_seq_str(seq, sub_pattern_dict):
        """
        Get subseq with intervals in arr and concatenate

        Args:
            seq: str
            sub_pattern_dict: [[0, 8], [24, 32], [48, 56]]

        Returns:
            str

        Raise:
            IndexError: if sequence length is not enough

        >>> sub_pattern_dict = [[0, 8]]
        >>> seq = "A" * 7
        >>> Barcode.get_seq_str(seq, sub_pattern_dict)
        Traceback (most recent call last):
        ...
        IndexError: sequence length is not enough in R1 read: AAAAAAA
        >>> seq = "A" * 8
        >>> Barcode.get_seq_str(seq, sub_pattern_dict)
        'AAAAAAAA'
        """
        seq_len = len(seq)
        ans = []
        for item in sub_pattern_dict:
            start, end = item[0], item[1]
            if end > seq_len:
                raise IndexError(f"sequence length is not enough in R1 read: {seq}")
            else:
                ans.append(seq[start:end])
        return ''.join(ans)

    @staticmethod
    def get_seq_list(seq, pattern_dict, abbr):
        """
        >>> pattern_dict = Barcode.parse_pattern("C2L3C2")
        >>> seq = "AAGGGTT"
        >>> Barcode.get_seq_list(seq, pattern_dict, "C")
        ['AA', 'TT']
        """
        
        return [seq[item[0]: item[1]] for item in pattern_dict[abbr]]

    @staticmethod
    @utils.add_log
    def parse_pattern(pattern):
        """
        >>> pattern_dict = Barcode.parse_pattern("C8L16C8L16C8L1U12T18")
        >>> pattern_dict['C']
        [[0, 8], [24, 32], [48, 56]]
        >>> pattern_dict['L']
        [[8, 24], [32, 48], [56, 57]]
        """
        pattern_dict = defaultdict(list)
        p = re.compile(r'([CLUNT])(\d+)')
        tmp = p.findall(pattern)
        if not tmp:
            Barcode.parse_pattern.logger.error(f'Invalid pattern: {pattern}')
            sys.exit()
        start = 0
        for item in tmp:
            end = start + int(item[1])
            pattern_dict[item[0]].append([start, end])
            start = end
        return pattern_dict

    @staticmethod
    def get_abbr_len(pattern_dict, abbr):
        """
        >>> pattern_dict = Barcode.parse_pattern("C8L16C8L16C8L1U12T18")
        >>> Barcode.get_abbr_len(pattern_dict, 'C')
        24
        >>> Barcode.get_abbr_len(pattern_dict, 'L')
        33
        """
        length = 0
        for item in pattern_dict[abbr]:
            length += item[1] - item[0]

        return length
        

    @staticmethod
    def get_scope_bc(chemistry):
        """Return (linker file path, whitelist file path)"""
        try:
            linker_f = glob.glob(f'{ROOT_PATH}/data/chemistry/{chemistry}/linker*')[0]
            whitelist_f = f'{ROOT_PATH}/data/chemistry/{chemistry}/bclist'
        except IndexError:
            return None, None
        return linker_f, whitelist_f

    @staticmethod
    def ord2chr(q, offset=33):
        return chr(int(q) + offset)

    @staticmethod
    def qual_int(char, offset=33):
        return ord(char) - offset

    @staticmethod
    def low_qual(quals, minQ, num):
        # print(ord('/')-33)           14
        return True if len([q for q in quals if Barcode.qual_int(q) < minQ]) > num else False

    @staticmethod
    def findall_mismatch(seq, n_mismatch=1, bases='ACGTN'):
        """
        choose locations where there's going to be a mismatch using combinations
        and then construct all satisfying lists using product

        Return:
        all mismatch <= n_mismatch set. 

        >>> answer = set(["TCG", "AAG", "ACC", "ATG", "ACT", "ACN", "GCG", "ANG", "ACA", "ACG", "CCG", "AGG", "NCG"])
        >>> seq_set = Barcode.findall_mismatch("ACG")
        >>> seq_set == answer
        True
        """
        seq_set = set()
        seq_len = len(seq)
        if n_mismatch > seq_len:
            n_mismatch = seq_len
        for locs in combinations(range(seq_len), n_mismatch):
            seq_locs = [[base] for base in seq]
            for loc in locs:
                seq_locs[loc] = list(bases)
            for poss in product(*seq_locs):
                seq_set.add(''.join(poss))
        return seq_set

    @staticmethod
    @utils.add_log
    def get_mismatch_dict(seq_list, n_mismatch=1):
        """
        Return:
        mismatch dict. Key: mismatch seq, value: seq in seq_list

        >>> seq_list = ["AACGTGAT", "AAACATCG"]
        >>> mismatch_dict = Barcode.get_mismatch_dict(seq_list)
        >>> mismatch_dict["AACGTGAA"] == "AACGTGAT"
        True
        """
        mismatch_dict = {}

        for seq in seq_list:
            seq = seq.strip()
            if seq == '':
                continue
            for mismatch_seq in Barcode.findall_mismatch(seq, n_mismatch):
                mismatch_dict[mismatch_seq] = seq

        return mismatch_dict

    @staticmethod
    def check_seq_mismatch(seq_list, correct_set_list, mismatch_dict_list):
        '''
        Return bool_valid, bool_corrected, corrected_seq

        >>> seq_list = ['ATA', 'AAT', 'ATA']
        >>> correct_set_list = [{'AAA'},{'AAA'},{'AAA'}]
        >>> mismatch_dict_list = [Barcode.get_mismatch_dict(['AAA'])] * 3

        >>> Barcode.check_seq_mismatch(seq_list, correct_set_list, mismatch_dict_list)
        (True, True, 'AAAAAAAAA')

        >>> seq_list = ['AAA', 'AAA', 'AAA']
        >>> Barcode.check_seq_mismatch(seq_list, correct_set_list, mismatch_dict_list)
        (True, False, 'AAAAAAAAA')
        '''
        bool_valid = True
        bool_corrected = False
        corrected_seq = ''
        for index, seq in enumerate(seq_list):
            if seq not in correct_set_list[index]:
                if seq not in mismatch_dict_list[index]:
                    bool_valid = False
                    return bool_valid, bool_corrected, corrected_seq
                else:
                    bool_corrected = True
                    corrected_seq += mismatch_dict_list[index][seq]
            else:
                corrected_seq += seq
        return bool_valid, bool_corrected, corrected_seq

    @staticmethod
    def parse_whitelist_file(files: list, n_pattern: int, n_mismatch: int):
        """
        files: file paths
        n_pattern: number of sections in pattern
        n_mismatch: allowed number of mismatch bases
        Returns:
            white_set_list
            mismatch_list
        """
        n_files = len(files)
        if n_files == 1 and n_pattern > 1:
            files = [files[0]] * n_pattern
        elif n_files != n_pattern:
            sys.exit(f'number of whitelist files({n_files} files:{files}) != n_pattern({n_pattern})')
        
        white_set_list, mismatch_list = [], []
        for f in files:
            barcodes, _ = utils.read_one_col(f)
            white_set_list.append(set(barcodes))
            barcode_mismatch_dict = Barcode.get_mismatch_dict(barcodes, n_mismatch)
            mismatch_list.append(barcode_mismatch_dict)

        return white_set_list, mismatch_list

    @staticmethod
    def parse_chemistry(chemistry):
        """
        Returns: pattern_dict, barcode_set_list, barcode_mismatch_list, linker_set_list, linker_mismatch_list
        """
        pattern = PATTERN_DICT[chemistry]
        pattern_dict = Barcode.parse_pattern(pattern)
        linker_file, whitelist_file = Barcode.get_scope_bc(chemistry)

        barcode_set_list, barcode_mismatch_list = Barcode.parse_whitelist_file([whitelist_file], n_pattern=len(pattern_dict['C']), n_mismatch=1)
        linker_set_list, linker_mismatch_list = Barcode.parse_whitelist_file([linker_file],n_pattern=1, n_mismatch=2)

        return pattern_dict, barcode_set_list, barcode_mismatch_list, linker_set_list, linker_mismatch_list

    @staticmethod
    def check_polyT(seq, pattern_dict, min_polyT_count=MIN_T):
        """
        Return:
            True if polyT is found
        """
        seq_polyT = Barcode.get_seq_str(seq, pattern_dict['T'])
        n_polyT_found = seq_polyT.count('T')
        if n_polyT_found >= min_polyT_count:
            return True
        return False

    def open_files(self):
        if self.output_R1 or self.bool_flv:
            self.fh_fq1 = xopen(self.out_fq1, 'w')
        self.fh_fq2 = xopen(self.out_fq2, 'w')

        if self.nopolyT:
            self.fh_nopolyT_fq1 = xopen(self.nopolyT_1, 'w')
            self.fh_nopolyT_fq2 = xopen(self.nopolyT_2, 'w')

        if self.noLinker:
            self.fh_nolinker_fq1 = xopen(self.noLinker_1, 'w')
            self.fh_nolinker_fq2 = xopen(self.noLinker_2, 'w')

    def close_files(self):
        if self.output_R1 or self.bool_flv:
            self.fh_fq1.close()
        self.fh_fq2.close()

        if self.nopolyT:
            self.fh_nopolyT_fq1.close()
            self.fh_nopolyT_fq2.close()
        
        if self.noLinker:
            self.fh_nolinker_fq1.close()
            self.fh_nolinker_fq2.close()

    @utils.add_log
    def add_step_metrics(self):

        self.add_metric(
            name='Raw Reads',
            value=self.total_num,
            help_info='total reads from FASTQ files'
        )
        self.add_metric(
            name='Valid Reads',
            value=self.clean_num,
            total=self.total_num,
            help_info='reads pass filtering(filtered: reads without poly T, reads without linker, reads without correct barcode or low quality reads)'
        )

        BarcodesQ30 = sum([self.barcode_qual_Counter[k] for k in self.barcode_qual_Counter if k >= self.ord2chr(
            30)]) / float(sum(self.barcode_qual_Counter.values())) * 100
        BarcodesQ30 = round(BarcodesQ30, 2)
        BarcodesQ30_display = f'{BarcodesQ30}%'
        self.add_metric(
            name='Q30 of Barcodes',
            value=BarcodesQ30,
            display=BarcodesQ30_display,
            help_info='percent of barcode base pairs with quality scores over Q30',
        )

        UMIsQ30 = sum([self.umi_qual_Counter[k] for k in self.umi_qual_Counter if k >= self.ord2chr(
            30)]) / float(sum(self.umi_qual_Counter.values())) * 100
        UMIsQ30 = round(UMIsQ30, 2)
        UMIsQ30_display = f'{UMIsQ30}%'
        self.add_metric(
            name='Q30 of UMIs',
            value=UMIsQ30,
            display=UMIsQ30_display,
            help_info='percent of UMI base pairs with quality scores over Q30',
        )

        self.add_metric(
            name='No PolyT Reads',
            value=self.no_polyT_num,
            total=self.total_num,
            show=False
        )

        self.add_metric(
            name='Low Quality Reads',
            value=self.lowQual_num,
            total=self.total_num,
            show=False,
        )

        self.add_metric(
            name='No Linker Reads',
            value=self.no_linker_num,
            total=self.total_num,
            show=False,
        )

        self.add_metric(
            name='No Barcode Reads',
            value=self.no_barcode_num,
            total=self.total_num,
            show=False,
        )

        self.add_metric(
            name='Corrected Linker Reads',
            value=self.linker_corrected_num,
            total=self.total_num,
            show=False,
        )

        self.add_metric(
            name='Corrected Barcode Reads',
            value=self.barcode_corrected_num,
            total=self.total_num,
            show=False,
        )

        if self.clean_num == 0:
            raise Exception('no valid reads found! please check the --chemistry parameter.' + HELP_DICT['chemistry'])
        
        if self._assay == 'flv_trust4':
            self.add_metric(
                name='Valid Matched Reads',
                value=self.match_num,
                total=self.total_num,
                help_info='reads match with flv_rna cell barcodes'
            )

            self.add_metric(
                name='Matched Barcodes',
                value=len(self.match_cbs),
                help_info='barcodes match with flv_rna'
            )


    @utils.add_log
    def run(self):
        """
        Extract barcode and UMI from R1. Filter reads with 
            - invalid polyT
            - low quality in barcode and UMI
            - invalid inlinker
            - invalid barcode
            
        for every sample
            get chemistry
            get linker_mismatch_dict and barcode_mismatch_dict
            for every read in read1
                filter
                write valid R2 read to file
        """
        is_probe = False
        
        for i in range(self.fq_number):

            chemistry = self.chemistry_list[i]
            lowNum = int(self.lowNum)
            lowQual = int(self.lowQual)
            if chemistry == 'scopeV1':
                lowNum = min(0, lowNum)
                lowQual = max(10, lowQual)
                self.run.logger.info(f'scopeV1: lowNum={lowNum}, lowQual={lowQual} ')
            # get linker and whitelist
            bc_pattern = PATTERN_DICT[chemistry]
            if (bc_pattern):
                linker_file, whitelist_file = self.get_scope_bc(chemistry)
                whitelist_files = [whitelist_file]
            else:
                bc_pattern = self.pattern
                linker_file = self.linker
                whitelist_file = self.whitelist
                whitelist_files = whitelist_file.split(',')
            if not bc_pattern:
                raise Exception("invalid bc_pattern!")

            pattern_dict = self.parse_pattern(bc_pattern)

            bool_T = True if 'T' in pattern_dict else False
            bool_L = True if 'L' in pattern_dict else False
            bool_whitelist = (whitelist_file is not None) and whitelist_file != "None"
            C_len = sum([item[1] - item[0] for item in pattern_dict['C']])

            #QC
            bool_probe=False
            if self.probefile:
                probefile=creat_probe_file(self.probefile,self.sample)
                bool_probe = True
                count_dic = utils.genDict(dim=3)
                valid_count_dic = utils.genDict(dim=2)
                probe_dic = read_fasta(probefile)
                reads_without_probe = 0

            bool_amp=False
            ampfile=self.ampfile
            if ampfile and ampfile != 'None':
                bool_amp=True
                count_amp_dic=utils.genDict(dim=3)
                valid_count_amp_dic=utils.genDict(dim=2)
                amp_dic=read_fasta(ampfile)
                reads_without_amp = 0           
            #

            if bool_probe:
                is_probe = True
                prior_probe_mismatch_dict = {}
                for probe_name in probe_dic:
                    prior_seq = probe_dic[probe_name][:9]
                    for mismatch_seq in findall_mismatch(prior_seq,n_mismatch=1):
                        prior_probe_mismatch_dict[mismatch_seq] = probe_name
                re_dict = {prior_probe:re.compile(prior_probe) for prior_probe in prior_probe_mismatch_dict}

            with pysam.FastxFile(self.fq1_list[i], persist=False) as fq1, \
                    pysam.FastxFile(self.fq2_list[i], persist=False) as fq2:
                for entry1, entry2 in zip(fq1, fq2):
                    header1, seq1, qual1 = entry1.name, entry1.sequence, entry1.quality
                    header2, seq2, qual2 = entry2.name, entry2.sequence, entry2.quality
                    self.total_num += 1

                    # polyT filter
                    if bool_T and self.filterNoPolyT:
                        if not Barcode.check_polyT(seq1, pattern_dict):
                            self.no_polyT_num += 1
                            if self.nopolyT:
                                self.fh_nopolyT_fq1.write(
                                    '@%s\n%s\n+\n%s\n' % (header1, seq1, qual1))
                                self.fh_nopolyT_fq2.write(
                                    '@%s\n%s\n+\n%s\n' % (header2, seq2, qual2))
                            continue

                    # lowQual filter
                    C_U_quals_ascii = Barcode.get_seq_str(
                        qual1, pattern_dict['C'] + pattern_dict['U'])
                    if lowQual > 0 and Barcode.low_qual(C_U_quals_ascii, lowQual, lowNum):
                        self.lowQual_num += 1
                        continue

                    # linker filter
                    if bool_L and (not self.allowNoLinker):
                        seq_str = Barcode.get_seq_str(seq1, pattern_dict['L'])
                        bool_valid, bool_corrected, _ = Barcode.check_seq_mismatch(
                            [seq_str], linker_set_list, linker_mismatch_list)
                        if not bool_valid:
                            self.no_linker_num += 1
                            if self.noLinker:
                                self.fh_nolinker_fq1.write(
                                    f'@{header1}\n{seq1}\n{seq_str}\n{qual1}\n')
                                self.fh_nolinker_fq2.write(
                                    '@%s\n%s\n+\n%s\n' % (header2, seq2, qual2))
                            continue
                        elif bool_corrected:
                            self.linker_corrected_num += 1

                    # barcode filter
                    seq_list = self.get_seq_list(seq1, pattern_dict, 'C')
                    if self.bool_flv:
                        seq_list = [utils.reverse_complement(seq) for seq in seq_list[::-1]]
                    if bool_whitelist:
                        bool_valid, bool_corrected, corrected_seq = Barcode.check_seq_mismatch(
                            seq_list, barcode_set_list, barcode_mismatch_list)

                        if not bool_valid:
                            self.no_barcode_num += 1
                            continue
                        elif bool_corrected:
                            self.barcode_corrected_num += 1
                        cb = corrected_seq
                    else:
                        cb = "".join(seq_list)

                    self.clean_num += 1
                    self.barcode_qual_Counter.update(C_U_quals_ascii[:C_len])
                    self.umi_qual_Counter.update(C_U_quals_ascii[C_len:])

                    umi = Barcode.get_seq_str(seq1, pattern_dict['U'])

                    #QC
                    read_name_probe = 'None'
                    ## Extract partial reads statistics
                    if self.total_num > 1000000:
                        bool_probe = False
                    
                    if bool_probe:
                        # valid count
                        valid_count_dic[cb][umi] += 1
                        # output probe UMi and read count
                        """
                        # 耗时长需改进
                        find_probe = False
                        for probe_name in probe_dic:
                            probe_seq = probe_dic[probe_name]
                            probe_seq = probe_seq.upper()
                            
                            # 设置mismatch
                            for _i in range(70,130-(len(probe_seq))):
                                target_seq = seq1[_i:(_i+(len(probe_seq)))]                                                              
                                if hamming(target_seq,probe_seq) < 3:
                                    count_dic[probe_name][cb][umi] += 1
                                    find_probe = True
                                    break
                            if find_probe == True:
                                break
                        if not find_probe:
                            reads_without_probe += 1
                        
                        """
                        # factor
                        #re_dict = {prior_probe:re.compile(prior_probe) for prior_probe in prior_probe_mismatch_dict}
                        start_pos_dict = {prior_probe:[int(each.start())+70 for each in re_dict[prior_probe].finditer(seq1[70:130])] for prior_probe in re_dict}
                        start_pos_dict = {k:v for k,v in start_pos_dict.items() if len(v) != 0} 
                        probe_name_old = ''
                        find_probe = ''
                        for prior_probe in start_pos_dict:
                            find_probe = False
                            probe_name = prior_probe_mismatch_dict[prior_probe]
                            probe_seq = probe_dic[probe_name] 
                            probe_len = len(probe_seq)
                            for pos in start_pos_dict[prior_probe]:
                                targe_seq = seq1[pos:pos+probe_len]
                                if probe_name_old == probe_name:
                                    break
                                elif hamming_fac(targe_seq,probe_seq) < 4:
                                    count_dic[probe_name][cb][umi] += 1
                                    find_probe = True
                                    probe_name_old = probe_name
                                    break                       
                            if find_probe:
                                break
                        if not find_probe:
                            reads_without_probe += 1
                        

                    if bool_amp:
                        valid_count_amp_dic[cb][umi] += 1

                        find_amp=False
                        for amp_name in amp_dic:
                            amp_seq=amp_dic[amp_name]
                            if seq2[0:23].find(amp_seq.upper()) != -1:
                                count_amp_dic[amp_name][cb][umi] += 1
                                read_name_amp=amp_name
                                find_amp=True
                                break
                        if not find_amp:
                            reads_without_amp += 1
                    #


                    if self.bool_flv:
                        qual1 = 'F' * len(cb + umi)
                        self.barcode_read_Counter.update(cb)
                        if self._assay == 'flv_trust4' and cb in self.match_barcodes:
                            self.match_num += 1
                            self.match_cbs.add(cb)
                            if self.barcode_read_Counter[cb] <= 80000:
                                self.fh_fq2.write(f'@{cb}_{umi}_{self.total_num}\n{seq2}\n+\n{qual2}\n')
                                self.fh_fq1.write(f'@{cb}_{umi}_{self.total_num}\n{cb}{umi}\n+\n{qual1}\n')
                        elif self._assay == 'flv_CR':
                            self.fh_fq2.write(f'@{cb}_{umi}_{self.total_num}\n{seq2}\n+\n{qual2}\n')
                            self.fh_fq1.write(f'@{cb}_{umi}_{self.total_num}\n{cb}{umi}\n+\n{qual1}\n')

                    else:
                        self.fh_fq2.write(f'@{cb}_{umi}_{self.total_num}\n{seq2}\n+\n{qual2}\n')
                        if self.output_R1:
                            self.fh_fq1.write(f'@{cb}_{umi}_{self.total_num}\n{seq1}\n+\n{qual1}\n')                   
            
                self.run.logger.info(self.fq1_list[i] + ' finished.')

        self.close_files()
        self.add_step_metrics()

        #QC
        if is_probe:
            bool_probe = True
        if bool_probe:
            # total probe summary
            total_umi = 0
            total_valid_read = 0
            for cb in valid_count_dic:
                total_umi += len(valid_count_dic[cb])
                total_valid_read += sum(valid_count_dic[cb].values())
            Barcode.run.logger.info("total umi:"+str(total_umi))
            Barcode.run.logger.info("total valid read:"+str(total_valid_read))
            Barcode.run.logger.info("reads without probe:"+str(reads_without_probe))

            # probe summary
            count_list = []
            for probe_name in probe_dic:
                UMI_count = 0
                read_count = 0
                if probe_name in count_dic:
                    for cb in count_dic[probe_name]:
                        UMI_count += len(count_dic[probe_name][cb])
                        read_count += sum(count_dic[probe_name][cb].values())
            
            probe_names = sorted(count_dic.keys())

            if self.probefilemode  == "16S":
                dict_update = {}
                merge_dict_probe(count_dic, dict_update, probe_names, "16S-V9-probe")
                merge_dict_probe(count_dic, dict_update, probe_names, "P479-probe")
                merge_dict_probe(count_dic, dict_update, probe_names, "1061R-probe")
                count_dic = dict_update

            #json_file = self.outdir + '/' + self.sample + '_probe_dic.json'
            #with open(json_file, 'w', encoding='utf-8') as file:
            #    json.dump(count_dic, file, ensure_ascii=False, indent=4)                         

            count_probe_list2 = []
            for probe_name in count_dic.keys():
                UMI_count = 0
                read_count = 0
                #if amp_name in count_amp_dic:
                for cb in count_dic[probe_name]:
                    UMI_count += len(count_dic[probe_name][cb])
                    read_count += sum(count_dic[probe_name][cb].values())
                #break
                count_probe_list2.append({"probe_name": probe_name, "UMI_count": UMI_count, "read_count": read_count})
            df_probe_count = pd.DataFrame(count_probe_list2, columns=["probe_name", "read_count", "UMI_count"])

            def format_percent(x):
                x = str(round(x*100, 2))+"%"
                return x

            df_probe_count["read_fraction"] = (df_probe_count["read_count"]/total_valid_read).apply(format_percent)
            df_probe_count["UMI_fraction"] = (df_probe_count["UMI_count"]/total_umi).apply(format_percent)
            df_probe_count.sort_values(by="UMI_count", inplace=True, ascending=False)

            df_probe_count_file = self.outdir + '/' + self.sample + '_probe_count_stat.tsv'
            df_probe_count.to_csv(df_probe_count_file, sep="\t", index=False)


#            for probe_name in count_dic.keys():   
#                count_list.append(
#                    {"probe_name": probe_name, "UMI_count": UMI_count, "read_count": read_count})

#            df_count = pd.DataFrame(count_list, columns=[
#                                "probe_name", "read_count", "UMI_count"])
#            def format_percent(x):
#                x = str(round(x*100, 2))+"%"
#                return x
#            df_count["read_fraction"] = (
#                df_count["read_count"]/total_valid_read).apply(format_percent)
#            df_count["UMI_fraction"] = (
#                df_count["UMI_count"]/total_umi).apply(format_percent)
#            df_count.sort_values(by="UMI_count", inplace=True, ascending=False)
#            df_count_file = self.outdir + '/' + self.sample + '_probe_count.tsv'
#            df_count.to_csv(df_count_file, sep="\t", index=False)
        
        if bool_amp:
            total_amp_umi = 0
            total_valid_amp_read = 0

            for cb in valid_count_amp_dic:
                total_amp_umi += len(valid_count_amp_dic[cb])
                total_valid_amp_read += sum(valid_count_amp_dic[cb].values())
            Barcode.run.logger.info("total amp umi:"+str(total_amp_umi))
            Barcode.run.logger.info("total valid amp read:"+str(total_valid_amp_read))
            Barcode.run.logger.info("reads without amp:"+str(reads_without_amp))


            # probe summary
            count_amp_list = []
            for amp_name in amp_dic:
                UMI_count = 0
                read_count = 0
                if amp_name in count_amp_dic:
                    for cb in count_amp_dic[amp_name]:
                        UMI_count += len(count_amp_dic[amp_name][cb])
                        read_count += sum(count_amp_dic[amp_name][cb].values())
                count_amp_list.append(
                    {"amp_name": amp_name, "UMI_count": UMI_count, "read_count": read_count})
            df_amp_count = pd.DataFrame(count_amp_list, columns=[
                                "amp_name", "read_count", "UMI_count"])
            def format_percent(x):
                x = str(round(x*100, 2))+"%"
                return x
            df_amp_count["read_fraction"] = (
                df_amp_count["read_count"]/total_valid_amp_read).apply(format_percent)
            df_amp_count["UMI_fraction"] = (
                df_amp_count["UMI_count"]/total_amp_umi).apply(format_percent)
            df_amp_count.sort_values(by="UMI_count", inplace=True, ascending=False)
            df_amp_count_file = self.outdir + '/' + self.sample + '_amp_count.tsv'
            df_amp_count.to_csv(df_amp_count_file, sep="\t", index=False)   

            #amp_dic_json = self.outdir + '/' + self.sample + '_amp_dic.json'
            #with open(amp_dic_json, 'w', encoding='utf-8') as file:
            #    json.dump(amp_dic, file, ensure_ascii=False, indent=4)
            
            #count_amp_dic_json = self.outdir + '/' + self.sample + '_count_amp_dic.json'
            #with open(count_amp_dic_json, 'w', encoding='utf-8') as file:
            #    json.dump(count_amp_dic, file, ensure_ascii=False, indent=4)

            probe_names = sorted(count_amp_dic.keys())

            if self.probefilemode  == "16S":
                dict_update = {}
                merge_dict_probe(count_amp_dic, dict_update, probe_names, "16S-V2-R2")
                merge_dict_probe(count_amp_dic, dict_update, probe_names, "16S-V8-R2")
                merge_dict_probe(count_amp_dic, dict_update, probe_names, "515f-R2")
                merge_dict_probe(count_amp_dic, dict_update, probe_names, "P799-R2")
                merge_dict_probe(count_amp_dic, dict_update, probe_names, "V6-F4-R")

                count_amp_dic = dict_update
                count_amp_dic.keys()

            count_amp_list2 = []
            for amp_name in count_amp_dic.keys():
                UMI_count = 0
                read_count = 0
                #if amp_name in count_amp_dic:
                for cb in count_amp_dic[amp_name]:
                    UMI_count += len(count_amp_dic[amp_name][cb])
                    read_count += sum(count_amp_dic[amp_name][cb].values())
                #break
                count_amp_list2.append({"amp_name": amp_name, "UMI_count": UMI_count, "read_count": read_count})
            df_amp_count = pd.DataFrame(count_amp_list2, columns=["amp_name", "read_count", "UMI_count"])

            def format_percent(x):
                x = str(round(x*100, 2))+"%"
                return x

            df_amp_count["read_fraction"] = (df_amp_count["read_count"]/total_valid_amp_read).apply(format_percent)
            df_amp_count["UMI_fraction"] = (df_amp_count["UMI_count"]/total_amp_umi).apply(format_percent)
            df_amp_count.sort_values(by="UMI_count", inplace=True, ascending=False)

            df_amp_count_file = self.outdir + '/' + self.sample + '_amp_count_stat.tsv'
            df_amp_count.to_csv(df_amp_count_file, sep="\t", index=False)  
        #



@utils.add_log
def barcode(args):
    with Barcode(args, display_title='Demultiplexing') as runner:
        runner.run()


def get_opts_barcode(parser, sub_program=True):
    parser.add_argument(
        '--chemistry',
        help='Predefined (pattern, barcode whitelist, linker whitelist) combinations. ' + HELP_DICT['chemistry'],
        choices=list(PATTERN_DICT.keys()),
        default='auto'
    )
    parser.add_argument(
        '--pattern',
        help="""The pattern of R1 reads, e.g. `C8L16C8L16C8L1U12T18`. The number after the letter represents the number 
        of bases.  
- `C`: cell barcode  
- `L`: linker(common sequences)  
- `U`: UMI    
- `T`: poly T""",
    )
    parser.add_argument(
        '--whitelist',
        help='Cell barcode whitelist file path, one cell barcode per line.'
    )
    parser.add_argument(
        '--linker',
        help='Linker whitelist file path, one linker per line.'
    )
    parser.add_argument(
        '--lowQual',
        help='Default 0. Bases in cell barcode and UMI whose phred value are lower than \
lowQual will be regarded as low-quality bases.',
        type=int,
        default=0
    )
    parser.add_argument(
        '--lowNum',
        help='The maximum allowed lowQual bases in cell barcode and UMI.',
        type=int,
        default=2
    )
    parser.add_argument(
        '--nopolyT',
        help='Outputs R1 reads without polyT.',
        action='store_true',
    )
    parser.add_argument(
        '--noLinker',
        help='Outputs R1 reads without correct linker.',
        action='store_true',
    )
    parser.add_argument(
        '--filterNoPolyT',
        help="Filter reads without PolyT.",
        action='store_true'
    )
    parser.add_argument(
        '--allowNoLinker',
        help="Allow valid reads without correct linker.",
        action='store_true'
    )
    parser.add_argument(
        '--gzip',
        help="Output gzipped fastq files.",
        action='store_true'
    )
    parser.add_argument(
        '--output_R1',
        help="Output valid R1 reads.",
        action='store_true'
    )

    parser.add_argument(
        '--probe_file', help="probe fasta file"  # js
    )
    parser.add_argument(
        '--probe_file_mode', help="probe fasta file mode"  # js
    )
    parser.add_argument(
        '--amp_file', help="amplification fasta file"
    )

    if sub_program:
        parser.add_argument('--fq1', help='R1 fastq file. Multiple files are separated by comma.', required=True)
        parser.add_argument('--fq2', help='R2 fastq file. Multiple files are separated by comma.', required=True)
        parser.add_argument('--match_dir', help='Matched scRNA-seq directory, required for flv_trust4')
        parser = s_common(parser)

    return parser


if __name__ == '__main__':
    unittest.main()
