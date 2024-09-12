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

from time import time
import multiprocessing
from collections import defaultdict
from .Recursive import Recursive


class CodonWorker(multiprocessing.Process):

    def __init__(self, input_queue, output_queue, new_dict, rules_dict, method, redundancy):

        super(CodonWorker, self).__init__()

        self.input_queue = input_queue
        self.output_queue = output_queue
        self.position_cache = defaultdict(dict)

        self.method = method
        self.new_dict = new_dict
        self.rules_dict = rules_dict
        self.redundancy = redundancy

    def push(self, combination):

        self.input_queue.put(combination)

    def run(self):

        new_dict = self.new_dict
        rules_dict = self.rules_dict
        redundancy = self.redundancy

        BestReduceSize = 20

        if(self.method == 'R'):
            BestRatio = 1000000
        elif(self.method == 'U'):
            BestRatio = 0

        BestList = []
        BestReducedList = []
        t = 0

        inverse_rule_dict = {}

        for key in rules_dict:
            inverse_rule_dict[rules_dict[key]] = key

        combo = []

        while(True):

            combo = self.input_queue.get(block=True)

            if(combo == None):
                break

            codons, ratios = CodonWorker.CreateListFromIndex(combo, new_dict)

            t += 1

            if redundancy != 0:
                pass
            else:
                recursive = Recursive(codons, rules_dict, inverse_rule_dict)
                reduced_list = recursive.Reduce(None)

                total_usage_frequency = 0
                for frequency in ratios:
                    total_usage_frequency += frequency

                total_ranking_sum = sum(combo)

                if(self.method == 'R' and (len(reduced_list) < BestReduceSize or (len(reduced_list) == BestReduceSize and total_ranking_sum < BestRatio))):
                    BestList = [codons, ratios]
                    #BestList = [codons, ratios, ranks]
                    BestReduceSize = len(reduced_list)
                    BestRatio = total_ranking_sum
                    BestReducedList = reduced_list

                elif self.method == 'U' and (len(reduced_list) < BestReduceSize or (len(reduced_list) == BestReduceSize and total_usage_frequency > BestRatio)):
                    BestList = [codons, ratios]
                    BestReduceSize = len(reduced_list)
                    BestRatio = total_usage_frequency
                    BestReducedList = reduced_list

        information_dict = {
            "BestReducedList":BestReducedList,
            "BestList": BestList,
            "ReduceSize": BestReduceSize,
            "Ratio": BestRatio
        }
        self.output_queue.put(information_dict)

    @staticmethod
    def CreateListFromIndex(empty_list, new_dict):
        codons = []
        ratios = []
        #ranks = []

        counter = 0

        for aa in new_dict:
            #print aa,empty_list[i],empty_list,i
            (codon, ratio) = new_dict[aa][empty_list[counter]]
            codons.append(codon)
            ratios.append(ratio)
            #ranks.append(counter + 1) ## need the index of ...
            counter+=1

        return codons, ratios
