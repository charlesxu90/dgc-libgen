"""
Adapted from https://github.com/bioverse/dynamcc/blob/master/src/DYNAMCC_0.py

Copyright (c) 2016, Andrea Halweg-Edwards, Gur Pines, Assaf Pines, James Winkler
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of DYNAMCC nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""
import os
import argparse
import itertools
import multiprocessing
from loguru import logger
from collections import defaultdict
from .CodonWorker import CodonWorker
from . import util
from loguru import logger

NAT_AAs = ['A', 'R', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'N', 'X']
RULES = {'R': ['A', 'G'],
             'Y': ['C', 'T'],
             'M': ['A', 'C'],
             'K': ['G', 'T'],
             'S': ['C', 'G'],
             'W': ['A', 'T'],
             'H': ['A', 'C', 'T'],
             'B': ['C', 'G', 'T'],
             'V': ['A', 'C', 'G'],
             'D': ['A', 'G', 'T'],
             'N': ['A', 'C', 'G', 'T'],
             'A': ['A'],
             'C': ['C'],
             'G': ['G'],
             'T': ['T']
             }

def CalcCombinations(filtered_dict):
    """Calculates the total combinations of codons possible. This is n! for
    all combinations of codons in the filtered_dict.

    Parameters
    ----------
    filtered_dict : dict
        Dictionary formatted in the same way as EditUsageDict(). Dictionary of
        lists of dictionaries for codon usage. Technically, any
        dictionary that has single letter amino acid symbols as keys would
        work

    Returns
    -------
    total : int
        This is the total number of combinations of codons that are possible
        after the user has selected which codon(s) to remove.

    Examples
    --------
    >>> CalcCombinations(filtered_dict)
    """
    total = 1
    for key in filtered_dict:
        total *= len(filtered_dict[key])
    return total


def BuildCodonCount(new_dict, key_order):
    """Iterate through a codon usage dictionary, and build a list that contains
    the sum of number of codons for each key (amino acid) at the indices

    Parameters
    ----------
    new_dict : dict
        This is the dictionary returned by either RemoveCodonByRank or
        ByUsage. The dictionary formatted the same as EditUsageDict().
        Dictionary of lists of dictionaries except that the codons with
        usage frequency below user specified rank or usage have been removed

    Returns
    -------
    codon_count : list
        A list containing the number of codons for each amino acid (this is
        a list of integers)
    """
    codon_count = []
    for key in key_order:
        codon_count.append(len(new_dict[key]))
    return codon_count


def FindMinimumThreshold(filtered_dict):
    """Iterates through filtered_dict and builds a list of all the codons with
    the hightest usage (because the list of codons is ordered by usage,
    this script just grabs the first index of each value's list). Then the
    script returns the lowest number in the list. This number will be the
    minimum threshold (i.e. users should select a cut-off below this number
    otherwise they will be omitting codons that they did not intend to omit).

    Parameters
    ----------
    filtered_dict : dict
        Dictionary formatted in the same way as EditUsageDict(). Dictionary of
        lists of dictionaries for codon usage. Technically, any
        dictionary that has single letter amino acid symbols as keys would
        work

    Returns
    -------
    float(min(usage_list)) : float
        Among all the codons remaing (after the user has thrown out particular
        residues), this number represents the lowest usage frequency of all
        codons with the highest usage frequency. In other words, given a list
        of the highest usage frequency for each codon, return the minimum
        number.

    Examples
    --------
    >>> min_threshold = FindMinimumThreshold(filtered_dict)
    """
    usage_list = []
    for key in filtered_dict:
        highest_usage = filtered_dict[key][0]
        for codon in highest_usage:
            usage_list.append(highest_usage[codon])
    return float(min(usage_list))


def RemoveLowCodons(threshold, filtered_dict):
    """Given the user input for the codon usage threshold (i.e. the user
    would like to remove all codons with a usage frequency below the threshold)
    the script builds a new dictionary (in the same format as the input
    dictionary) that does not contain the codons with usage frequency below
    threshold.

    Parameters
    ----------
    threshold : float
        This is the codon usage frequency input by the user that specifies
        which codons to remove based on usage frequency.
    filtered_dict : dict
        Dictionary formatted in the same way as EditUsageDict(). Dictionary of
        lists of dictionaries for codon usage. Technically, any
        dictionary that has single letter amino acid symbols as keys would
        work

    Returns
    -------
    new_dict : dict
        Dictionary formatted the same as input except that the codons with
        usage frequency below user specified threshold have been removed

    Examples
    --------
    """
    new_dict = {}
    for key1 in filtered_dict:
        new_dict[key1] = [x for x in filtered_dict[key1] if x[1] > threshold]
        if (len(new_dict[key1]) == 0):
            raise SyntaxError('No codons avaliable for key: %s' % key1)
    return new_dict


def RemoveCodonByRank(rank, sorted_dict):
    """This script builds a new dictionary with the same format as the input
    dictionary except that is will not contain codons below the 'rank'
    specified by the user. In this case, rank is an integer that corresponds
    to 1 + the index of a list of codons ordered by usage frequency. In other
    words, given a list of codon, the codon with the highest usage frequency
    will be given rank 1 and the codon with the lowest usage frequency will be
    given the highest number (depending on how many codons code for the
    particular amino acid).

    Parameters
    ----------
    rank : int
        the user specified codon rank that has been chosen as a cutoff value
        (i.e. codons with a rank below the specified rank will not be
        included)
    sorted_dict : dict
        Dictionary formatted in the same way as EditUsageDict(). Dictionary of
        lists of dictionaries for codon usage. Technically, any
        dictionary that has single letter amino acid symbols as keys would
        work

    Returns
    -------
    new_dict : dict
        Dictionary formatted the same as input except that the codons with
        usage frequency below user specified rank have been removed

    Examples
    --------
    >>> rank = int(raw_input("Set codon rank threshold: "))
    >>> RemoveCodonByRank(rank, filtered_dict)
    """
    new_dict = {}
    for key in sorted_dict:
        if len(sorted_dict[key]) > rank:
            new_dict[key] = sorted_dict[key][0:rank]
        else:
            new_dict[key] = sorted_dict[key]

    return new_dict


def SetRedundancy():
    """
    Parameters
    ----------
    none

    Returns
    -------
    redundancy : int
        Return the integer entered by the user
    """
    redundancy = int(input("Set redundancy" +
                           " (0 for no redundancy at all): "))
    return redundancy


def CalcTotCombinations(combinations, codon_count, new_dict, redundancy):
    """
    Parameters
    ----------
    combinations : int
    codon_count : int
    new_dict : dict
    redundancy : int

    Returns
    -------
    tot_combinations : int
    """
    tot_combinations = (combinations *
                        (codon_count - len(new_dict)) ** redundancy)
    return tot_combinations


def codon_exploder(codon_count):

    temp_list = []

    for item in codon_count:
        temp_list.append(list(range(0, item)))

    combinations = itertools.product(*temp_list)
    return combinations


def start_multiprocessing(new_dict, rules_dict, selection, codon_count, redundancy, processes=3):

    if (selection != 'R' and selection != 'U'):
        logger.info('Unknown ranking method selected.')
        return

    worker_array = []
    output_queue = multiprocessing.Queue()
    for thread in range(0, processes):
        input_queue = multiprocessing.Queue()
        w = CodonWorker(input_queue, output_queue, new_dict, rules_dict, selection, redundancy)
        worker_array.append(w)

    for w in worker_array:
        w.start()

    cyclical_adder = itertools.cycle(list(range(0, processes)))

    combinations = codon_exploder(codon_count)

    for combo in combinations:
        worker_array[next(cyclical_adder)].push(combo)

    for worker in worker_array:
        worker.push(None)

    output = []
    while (len(output) != len(worker_array)):
        output.append(output_queue.get(block=True))

    BestReduceSize = 22
    BestRatio = 0

    best_result = None
    for result in output:
        # print('result: ', result)

        reduced_list = result['BestReducedList']
        total_ratio = result['Ratio']
        reduced_size = result['ReduceSize']

        if selection == 'R' and reduced_size < BestReduceSize or (reduced_size == BestReduceSize and total_ratio < BestRatio):
            best_result = result

        elif selection == 'U' and reduced_size < BestReduceSize or (reduced_size == BestReduceSize and total_ratio > BestRatio):
            best_result = result

        BestReduceSize = best_result['ReduceSize']
        BestRatio = best_result['Ratio']

    # for key in best_result:
    #     print(key, best_result[key])

    return best_result


def expand_codons(best_reduced_list, rules):
    # logger.debug(f"best_reduced_list: {best_reduced_list}")
    exploded_codons = {}
    codon_list = []
    for codon in best_reduced_list:
        exploded_codons[codon] = list(codon)
        codon_list.append(list(codon))
    # logger.debug(f"exploded_codons: {exploded_codons}")
    # logger.debug(f"codon_list: {codon_list}")

    exploded_codons_copy1 = {}
    for key in exploded_codons:
        exploded_codons_copy1[key] = []
    for codon in exploded_codons:
        for j in range(len(exploded_codons[codon])):
            exploded_codons_copy1[codon].append(rules[exploded_codons[codon][j]])
    # logger.debug(f"exploded_codons_copy1: {exploded_codons_copy1}")

    exploded_codons_copy2 = {}
    for key in exploded_codons:
        exploded_codons_copy2[key] = []
    for codon in exploded_codons_copy1:
        combos = list(itertools.product(*exploded_codons_copy1[codon]))
        for combo in combos:
            exploded_codons_copy2[codon].append(combo)
    # logger.debug(f"exploded_codons_copy2: {exploded_codons_copy2}")

    exploded_codons = {}
    for key in exploded_codons_copy2:
        exploded_codons[key] = []
        for value in exploded_codons_copy2[key]:
            joined_codon = ''.join(list(value))
            exploded_codons[key].append(joined_codon)
    return exploded_codons

def load_ecoli_codons():
    # get the script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    sorted_dict = util.BuildUsageDict(os.path.join(script_dir, 'ecoli_codons.txt'))
    return sorted_dict

def get_aa_top1_codon(aa, sorted_dict,):
    try:
        return sorted_dict[aa][0][0], None
    except KeyError:
        logger.info(f"Invalid amino acid: {aa} with {sorted_dict[aa]}")
        return None, None


def get_dg_codon_dict(aas, keep_or_remove='keep', compression='rank', rank=2):
    remove_aa = aas if keep_or_remove is False else list(set(NAT_AAs) - set(aas))
    script_dir = os.path.dirname(os.path.abspath(__file__))
    sorted_dict = util.BuildUsageDict(os.path.join(script_dir, 'ecoli_codons.txt'))
    rules_dict, inverse_dict = util.BuildRulesDict(os.path.join(script_dir, 'rules.txt'))
    filtered_dict = util.EditUsageDict(remove_aa, sorted_dict)

    if compression == 'rank':
        new_dict = RemoveCodonByRank(rank, filtered_dict)
    else:
        new_dict = RemoveLowCodons(0.1, filtered_dict)

    codon_order = list(new_dict.keys())
    codon_count = BuildCodonCount(new_dict, codon_order)

    rank_method = 'R' if compression == 'rank' else 'U'
    best_result = start_multiprocessing(new_dict, rules_dict, rank_method, codon_count, redundancy=0, processes=3)

    exploded_codons = expand_codons(best_result['BestReducedList'], RULES)
    
    return exploded_codons, sorted_dict



def main(args):
    aas = args.input_aa.split(',')
    assert len(aas) > 0, 'No amino acids provided'
    assert all(aa in NAT_AAs for aa in aas), 'Invalid amino acid provided'

    assert args.keep_or_remove in ['keep', 'remove'], 'Invalid keep_or_remove argument'

    assert args.compression in ['rank', 'usage'], 'Invalid compression argument'
    assert args.rank > 0, 'Invalid rank argument'

    exploded_codons, sorted_dict = get_dg_codon_dict(aas, args.keep_or_remove, args.compression, args.rank)
    logger.info(f"Exploded codons: {exploded_codons}")

    codon_dict = util.BuildCodonDict(sorted_dict)
    
    for dg_codon in exploded_codons: 
        logger.info(f"DG codon: {dg_codon}")
        for codon in exploded_codons[dg_codon]:
            logger.info(f"codon: {codon}, rank: {codon_dict[codon]}")


def get_args():
    parser = argparse.ArgumentParser(description='Optimize protein property via reinforcement learning')
    parser.add_argument('-i', '--input_aa', type=str, help='Input amino acid list')
    parser.add_argument('-k', '--keep_or_remove', type=str, default='keep', help='Keep or remove, keep by defult')
    parser.add_argument('-c', '--compression', type=str, default='rank',  help='Compression method: rank or usage, rank by default')
    parser.add_argument('-r', '--rank', type=int, default=2, help='Top n codons to usage, defult is n=2')

    return parser.parse_args()


if __name__ == '__main__':
    args = get_args()
    main(args)
