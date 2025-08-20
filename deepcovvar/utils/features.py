#!/usr/bin/python3

"""
Title: Script for protein feature generation
Author: Naveen Duhan
"""

from .feature_data import *

class FEATURE:

    def __init__(self):
        self.AA = AA
        self._conj = conj
        self.gr = gr

    def CKSAAP(self, seq, gap=5, order='alphabetically'):
        """Composition of k-spaced amino acid pairs"""
        if gap < 0:
            print('Error: the gap should be equal or greater than zero' + '\n\n')
            return 0

        if len(seq) < gap + 2:
            print('Error: all the sequence length should be larger than the (gap value) + 2 = ' + str(gap + 2) + '\n\n')
            return 0
        
        AA = myAAorder[order]
        encodings = []
        aaPairs = []
        for aa1 in AA:
            for aa2 in AA:
                aaPairs.append(aa1 + aa2)
        
        for g in range(gap + 1):
            for aa in aaPairs:
                count = 0
                for i in range(len(seq) - g - 1):
                    if seq[i] == aa[0] and seq[i + g + 1] == aa[1]:
                        count += 1
                encodings.append(count)
        return encodings

    def hybrid(self, seq, method1, method2):
        """Hybrid features combining two methods"""
        features1 = method1(seq)
        features2 = method2(seq)
        return features1 + features2

    def hybrid2(self, seq, method1, method2, method3):
        """Hybrid features combining three methods"""
        features1 = method1(seq)
        features2 = method2(seq)
        features3 = method3(seq)
        return features1 + features2 + features3

    def hybrid3(self, seq, method1, method2, method3, method4):
        """Hybrid features combining four methods"""
        features1 = method1(seq)
        features2 = method2(seq)
        features3 = method3(seq)
        features4 = method4(seq)
        return features1 + features2 + features3 + features4

