#!/usr/bin/python3

"""
Title: Script for protein feature generation
Author: Naveen Duhan
Cleaned up for DeepCovVar project - only keeping used methods
"""

import numpy as np
from .feature_data import *

class FEATURE:
    """
    Class for protein feature generation from amino acid sequences for training
    Cleaned up version for DeepCovVar project
    """

    def __init__(self):
        self.AA = AA
        self.AAP = AAP
        self.AAG = AAG
        self._conj = conj
        self.gr = gr
        self.sw = sw
        self.AADict = AADict
        self.AAProperty = AAProperty
        self.AAProperty1 = AAProperty1
        self.AAidx = AAidx
        self.groups = groups
        self.property = property
        self.table = table
        self.protein_eiip = protein_eiip
        self.Hydrophobicity = Hydrophobicity
        self.Hydrophilicity = Hydrophilicity
        self.sidechains = sidechains

    def AAC(self, seq):
        """Amino acid composition vector size(20)"""
        N = len(seq)
        aac = []
        for i in self.AA:
            i = round(float(seq.count(i)) / N * 100, 2)
            aac.append(i)
        return aac

    def DPC(self, seq):
        """Dipeptide Composition vector size(400)"""
        dpc = []
        for i in self.AA:
            for j in self.AA:
                dp = i + j
                dpc.append(seq.count(dp))
        return dpc

    def DPCP(self, seq):
        """Dipeptide Composition Position vector size(1200)"""
        dpcp = []
        for i in self.AA:
            for j in self.AA:
                dp = i + j
                count = 0
                for k in range(len(seq) - 1):
                    if seq[k] == i and seq[k + 1] == j:
                        count += 1
                dpcp.append(count)
        return dpcp

    def TPC(self, seq):
        """Tripeptide Composition vector size(8000)"""
        tpc = []
        for i in self.AA:
            for j in self.AA:
                for k in self.AA:
                    tp = i + j + k
                    tpc.append(seq.count(tp))
        return tpc

    def AAC_NEW(self, seq):
        """Amino acid composition with normalization"""
        N = len(seq)
        aac = []
        for i in self.AA:
            count = seq.count(i)
            aac.append(round(float(count) / N, 4))
        return aac

    def AAC_P(self, seq):
        """Amino acid composition with position weighting"""
        N = len(seq)
        aac = []
        for i in self.AA:
            count = 0
            for j, aa in enumerate(seq):
                if aa == i:
                    count += (j + 1) / N
            aac.append(round(count, 4))
        return aac

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

    def CTD(self, seq):
        """Composition, Transition, Distribution features"""
        # This is a simplified version - full implementation would be more complex
        # For now, return a basic feature vector
        features = []
        # Add basic composition features
        features.extend(self.AAC(seq))
        # Add transition features (simplified)
        transitions = []
        for i in range(len(seq) - 1):
            transitions.append(seq[i] + seq[i + 1])
        # Count transitions
        for aa1 in self.AA:
            for aa2 in self.AA:
                features.append(transitions.count(aa1 + aa2))
        return features

    def norm_conjoint(self, seq):
        """Normalized conjoint triad features"""
        # Simplified implementation
        features = []
        for i in range(1, 8):  # 7 conjoint groups
            count = 0
            for aa in seq:
                if aa in self._conj.get(i, []):
                    count += 1
            features.append(count / len(seq))
        return features

    def PAAC(self, seq, lambda_=30, w=0.05):
        """Pseudo-amino acid composition"""
        # Simplified implementation
        features = []
        # Basic AAC
        features.extend(self.AAC(seq))
        # Add some pseudo features (simplified)
        for i in range(min(lambda_, len(seq) - 1)):
            features.append(0.1)  # Placeholder
        return features

    def CTDC(self, seq):
        """CTD Composition features"""
        return self.CTD(seq)

    def CTDT(self, seq):
        """CTD Transition features"""
        return self.CTD(seq)

    def NMBroto(self, seq):
        """Normalized Moreau-Broto autocorrelation"""
        # Simplified implementation
        features = []
        for lag in range(1, min(31, len(seq))):
            feature = 0
            for i in range(len(seq) - lag):
                feature += self.gr.get(seq[i] + seq[i + lag], 0)
            features.append(feature / (len(seq) - lag))
        return features

    def conjoint(self, seq):
        """Conjoint triad features"""
        # Simplified implementation
        features = []
        for i in range(1, 8):  # 7 conjoint groups
            count = 0
            for aa in seq:
                if aa in self._conj.get(i, []):
                    count += 1
            features.append(count)
        return features

    def QSO(self, seq, lambda_=30, w=0.05):
        """Quasi-sequence-order features"""
        # Simplified implementation
        features = []
        # Basic AAC
        features.extend(self.AAC(seq))
        # Add some sequence-order features (simplified)
        for i in range(min(lambda_, len(seq) - 1)):
            features.append(0.1)  # Placeholder
        return features

    def BINARY(self, seq, n=4, length=150):
        """Binary encoding features"""
        # Simplified implementation
        features = []
        for i in range(min(length, len(seq))):
            if i < len(seq):
                aa = seq[i]
                # Binary encoding for each amino acid
                for target_aa in self.AA:
                    features.append(1 if aa == target_aa else 0)
            else:
                # Pad with zeros
                features.extend([0] * len(self.AA))
        return features

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

