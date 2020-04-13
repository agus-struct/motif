# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt


test_protein_sequence = "MSDSRSDDRHPSVTPEQLATLSHEFRTPLNGVLGMARLLENTKLTAEQRSYVTALRESGDHLLSLVNDVLHFARLGAAAIELSLAPVDIEGLLRQVAELMSPRAHEKGIEIAWAVSSPLPTILADEGRLRQILLNFAGNAVKFTEAGGVLLTASAIDGGRVRFSVADTGPGVAPDARARIFEAFVQTDVTHATQLGGAGLGLAIVSRLSAAMGGAVGVGGELGQGAEFWFEAPFATAAAPLRAAPLEGRNVAIASPNAIVRAATARQIEAAGGRAYAAVDIASALAGAPADAVLLIDAALSGPRGALKPPAGRRSVVLLTPEQRDRIDRLKAAGFSGYLIKPLRAASLVAQVLQAVTADGVAEDEPAHDDRIAGAVASGARVLLAEDNPINALLARTLLEREGCIVDRVADGEQAIAAASAGVYDLILMDLRMPGLTGIEAARALRAKGVATPIAALTADAFDEDRRTCLAAGMDDFLVKPLTQEALRDALKRWTTGGVSGGWTKPATRAKVAG"
test_sequences = ['FITKS', 'YLLKD', 'YLSKG', 'FLSKR', 'IVLKQ', 'YLLKE', 'IINKD', 'AISKT',
                  'FIFKD', 'FVSKR', 'FVSKK', 'FVSKC', 'YLVKP', 'YLLKP', 'YILKP', 'FILKP',
                  'FIHKP', 'FVTKP', 'HFAKP', 'YIMKP', 'YIPKP', 'FIEKP', 'TVDKP', 'FLQKP',
                  'FIAKP', 'FLSKP', 'FVEKP', 'FLTKP', 'YLPKP', 'FLTKP', 'FLTKP', 'YLVKP',
                  'YLIKP', 'YVVKP', 'CLFKP', 'CLFKP', 'FLSKP', 'VIVKP', 'ILAKP', 'VLSKP',
                  'YLAKP', 'VLLKP', 'FLTKP', 'YITKP', 'QISKP', 'HLTKP', 'YLSKP', 'YLPKP',
                  'YLVKP', 'YVVKP', 'YVTKP', 'FIAKP']


class motif():
    
    def __init__(self, motif_sequences, *args, **kwargs):
        
        self.__sequences = motif_sequences
        self.__counts    = self.__counts_func()
        self.__rcounts   = self.__rcounts_func()
        self.__pwm       = self.__pwm_func(*args, **kwargs)
        self.__consensus, self.__mv = self.__consensus_func()

    
    def sequences(self):
        return self.__sequences
    
    def consensus(self):
        return self.__consensus
    
    def consensusval(self):
        return self.__mv
    
    def length(self):
        return len(self.__consensus)
    
    def counts(self):
        return self.__counts
    
    def rcounts(self):
        return self.__rcounts
    
    def pwm(self):
        return self.__pwm
    
    def scan(self, protein_sequence, plot=False, toList=False):
        
        result = self.__scan_func(protein_sequence)
            
        if plot is True:
            y = [s for s, p, f in result]
            x = range(len(y))
            #w = len(self.__pwm)
            #xticks = list(protein_sequence[:-w])
            avg = sum(y)/len(y)
            mv = self.__mv
            plt.figure(figsize=(14,4))
            plt.axhline(avg,              color="black",  linestyle="dotted")
            plt.axhline(avg+(mv-avg)/3,   color="red",    linestyle="dotted")
            plt.axhline(avg+(mv-avg)/3*2, color="orange", linestyle="dotted")
            plt.axhline(mv,               color="green",  linestyle="dotted")
            plt.step(x, y, where="post")
            plt.xlim(min(x), max(x)+1)
            plt.ylabel("Log Odds")
            plt.xlabel("Protein sequence position")
            
            #plt.xticks(x, xticks, fontsize=8)
            #plt.grid()
            plt.show()
            
        return result if toList is False else list(result)
    
    def find(self, protein_sequence, score_threshold=float("-inf"), return_best=True):
        
        scoreboard = []
        for s, p, f in self.__scan_func(protein_sequence):
            if s > score_threshold:
                scoreboard.append([s, p, f])
                
        scoreboard = sorted(scoreboard, key=lambda x: x[0])
        
        if len(scoreboard) == 0:
            return None
        elif return_best:
            spf = scoreboard[-1]
            return spf[0], spf[1], spf[2]
        else:
            return scoreboard
        

    def __str__(self):
        return self.__consensus
    
    def __counts_func(self, alphabet=['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']):
        
        sequences = self.__sequences
        counts_list = []
        col_range = range(min([len(s) for s in sequences]))
        sequences = np.asarray([list(s.upper()) for s in sequences])    # Split strings into uppercase-character lists
                                                                        # and convert into 2D numpy array
        for col_index in col_range:
            d = dict.fromkeys(alphabet, 0)
            uniques, counts = np.unique(sequences[:, col_index], return_counts=True)
            for i, u in enumerate(uniques):
                d[u] = float(counts[i])
                
            counts_list.append(d)
        return counts_list

    def __rcounts_func(self):
        
        relative_counts_list = []
        for d in self.__counts:
            totcounts = sum(d.values())
            e = {}
            for u in d.keys():
                e[u] = d[u]/totcounts
            relative_counts_list.append(e)
        return relative_counts_list
    
    def __pwm_func(self, pseudocounts=0.5):
        
        # TODO: add background instead of simply dividing by 1/totletters (i.e. multiplying by totletters)
        
        pseudocounts = max(1e-3, pseudocounts) # Safety
        pwm_list = []
        for d in self.__counts:
            totletters = len(d.keys())
            totcounts  = sum([val+pseudocounts for val in d.values()])
            e = {}
            for u in d.keys():
                e[u] = np.log2((d[u]+pseudocounts)/totcounts * totletters)
            pwm_list.append(e)
        
        
        return pwm_list
    
    def __consensus_func(self):
        cons = ""
        maxval = 0
        for pos in self.__pwm:
            aa, val = "", float("-inf")
            for letter, weight in pos.items():
                if weight > val:
                    aa = letter
                    val = weight
            cons += aa
            maxval += val
        return cons, maxval
    
    def __scan_func(self, sequence):
        w = len(self.__pwm)
        for position in range(0, len(sequence)-w):
            fragment = sequence[position:position+w]
            score_sum = 0
            for position_fragment, letter in enumerate(fragment):
                score_sum += self.__pwm[position_fragment][letter]
            average_score = np.average(score_sum)
            yield average_score, position, fragment


            
def parse_fasta(fasta_filepath): 
    fasta_dict = {}
    name   = ""
    seq    = ""
    with open(fasta_filepath, "r") as fh:
        for line in fh:
            line = line.strip().replace("-", "")
            if line.startswith(">"):
                # Add previous to db
                if name != "": fasta_dict[name.split(">")[-1]] = seq.replace(" ", "").replace("*", "")
                # Start saving current entry
                name = line
                seq = ""
            else:
                seq = "".join((seq, line))
                
    # For last sequence, after iteration finished
    fasta_dict[name.split(">")[-1]] = seq.replace(" ", "").replace("*", "")
    
    return fasta_dict