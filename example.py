# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 2020

@author: E. Agustoni (GitHub: agus-struct)
"""

from motif import motif

test_protein_sequence = "MSDSRSDDRHPSVTPEQLATLSHEFRTPLNGVLGMARLLENTKLTAEQRSYVTALRESGDHLLSLVNDVLHFARLGAAAIELSLAPVDIEGLLRQVAELMSPRAHEKGIEIAWAVSSPLPTILADEGRLRQILLNFAGNAVKFTEAGGVLLTASAIDGGRVRFSVADTGPGVAPDARARIFEAFVQTDVTHATQLGGAGLGLAIVSRLSAAMGGAVGVGGELGQGAEFWFEAPFATAAAPLRAAPLEGRNVAIASPNAIVRAATARQIEAAGGRAYAAVDIASALAGAPADAVLLIDAALSGPRGALKPPAGRRSVVLLTPEQRDRIDRLKAAGFSGYLIKPLRAASLVAQVLQAVTADGVAEDEPAHDDRIAGAVASGARVLLAEDNPINALLARTLLEREGCIVDRVADGEQAIAAASAGVYDLILMDLRMPGLTGIEAARALRAKGVATPIAALTADAFDEDRRTCLAAGMDDFLVKPLTQEALRDALKRWTTGGVSGGWTKPATRAKVAG"

test_sequences = ['FITKS', 'YLLKD', 'YLSKG', 'FLSKR', 'IVLKQ', 'YLLKE', 'IINKD', 'AISKT',
                  'FIFKD', 'FVSKR', 'FVSKK', 'FVSKC', 'YLVKP', 'YLLKP', 'YILKP', 'FILKP',
                  'FIHKP', 'FVTKP', 'HFAKP', 'YIMKP', 'YIPKP', 'FIEKP', 'TVDKP', 'FLQKP',
                  'FIAKP', 'FLSKP', 'FVEKP', 'FLTKP', 'YLPKP', 'FLTKP', 'FLTKP', 'YLVKP',
                  'YLIKP', 'YVVKP', 'CLFKP', 'CLFKP', 'FLSKP', 'VIVKP', 'ILAKP', 'VLSKP',
                  'YLAKP', 'VLLKP', 'FLTKP', 'YITKP', 'QISKP', 'HLTKP', 'YLSKP', 'YLPKP',
                  'YLVKP', 'YVVKP', 'YVTKP', 'FIAKP']

test_motif = motif(test_sequences)

scan_output = test_motif.scan(test_protein_sequence,
                              plot=True)

sorted_scan_output = sorted(scan_output,
                            key=lambda x: x[0],
                            reverse=True)
    
print("n\tscore\tpos\tfragment")
for i, j in enumerate(sorted_scan_output[:5]):
    print("%i\t%.2f\t%i\t%s"%(i+1, j[0], j[1], j[2]))
