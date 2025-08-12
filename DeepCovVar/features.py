#!/usr/bin/python3

"""
Title: Script for protein feature generation
Author: Naveen Duhan

"""

# import statement

import copy, math,re
import numpy as np
from .feature_data import *

# Feature class

class FEATURE:
    """
    Class for protein feature generation from amino acid sequences for training
    """


    def __init__(self):
        self.AA=AA
        self.AAP=AAP
        self.AAG=AAG
        self._conj=conj
        self.gr=gr
        self.sw=sw
        self.AADict=AADict
        self.AAProperty=AAProperty
        self.AAProperty1=AAProperty1
        self.AAidx=AAidx
        self.groups=groups
        self.property=property
        self.table=table
        self.protein_eiip= protein_eiip
        self.Hydrophobicity = Hydrophobicity
        self.Hydrophilicity = Hydrophilicity
        self.sidechains = sidechains
        

    # Amino acid composition vector size(20)
    def AAC(self, seq):
        N = len(seq)
        aac = []
        for i in self.AA:
            i = round(float(seq.count(i)) / N * 100, 2)
            aac.append(i)
        return aac
        # Amino acid composition N-terminal residue

    def AAC_nt(self, seq, N):
        aac_nt = []
        sub_str = seq[0:N]
        for i in self.AA:
            count = 0
            for j in sub_str:
                tp = j
                if tp == i:
                    count += 1
                nt = round(float(count / N) * 100, 2)
                aac_nt.append(nt)
        return aac_nt

        # Amino acid composition C-terminal residue

    def AAC_ct(self, seq, N):
        aac_ct = []
        sub_str = seq[-N:]
        for i in self.AA:
            count = 0
            for j in sub_str:
                tp = j
                if tp == i:
                    count += 1
                nt = round(float(count / N) * 100, 2)
                aac_ct.append(nt)
        return aac_ct
    # Dipeptide Composition vector size(400)
    
    def CKSAAP(self,seq, gap=5, order='alphabetically'):
        if gap < 0:
            print('Error: the gap should be equal or greater than zero' + '\n\n')
            return 0

        if len(seq) < gap+2:
            print('Error: all the sequence length should be larger than the (gap value) + 2 = ' + str(gap+2) + '\n\n')
            return 0
        AA = myAAorder[order] 
        encodings = []
        aaPairs = []
        for aa1 in AA:
            for aa2 in AA:
                aaPairs.append(aa1 + aa2)
        # header = ['#']
        # for g in range(gap+1):
        #     for aa in aaPairs:
        #         header.append(aa + '.gap' + str(g))
        # encodings.append(header)
        
        code =[]
        for g in range(gap+1):
            myDict = {}
            for pair in aaPairs:
                myDict[pair] = 0
            sum = 0
            for index1 in range(len(seq)):
                index2 = index1 + g + 1
                if index1 < len(seq) and index2 < len(seq) and seq[index1] in AA and seq[index2] in AA:
                    myDict[seq[index1] + seq[index2]] = myDict[seq[index1] + seq[index2]] + 1
                    sum = sum + 1
            for pair in aaPairs:
                encodings.append(myDict[pair] / sum)
        
        return encodings


    def DPC(self, seq):
        N = len(seq)
        dpc = []
        for i in self.AA:
            for j in self.AA:
                dp = i + j
                dp = round(float(seq.count(dp)) / (N - 1) * 100, 2)
                dpc.append(dp)
        return dpc
    
    def DPCP(self, seq):
        N = len(seq)
        dpc = []
        for i in self.AA:
            for j in self.AA:
                dp = i + j
                dpk = round(float(seq.count(dp)) / (N - 1) * 100, 2)
                dpc.append(dpk)
                dpc.append(seq.count(dp))
                dpc.append(self.sw[dp])
        return dpc

    # Tripeptide composition vector size(8000)
    def TPC(self, seq):
        N = len(seq)
        tpc = []
        for i in self.AA:
            for j in self.AA:
                for k in self.AA:
                    tripep = i + j + k

                    tripep = round(float(seq.count(tripep)) / (N - 1) * 100, 2)

                    tpc.append(tripep)

        return tpc
    
    def AAC_NEW(self, seq):
        N = len(seq)
        aac = []
        for i in AA:
            j = round(float(seq.count(i)) / N * 100, 2)
            data = table[i]
            aac.append(seq.count(i))
            aac.append(j)

            for k in self.table[i]:
                aac.append(k)
            aac.append(self.protein_eiip[i])
        return aac
    
    def AAC_P(self, seq):
        N = len(seq)
        aac = []
        for i in AA:
            j = round(float(seq.count(i)) / N * 100, 2)
            aac.append(seq.count(i))
            aac.append(j)

            for k in self.table[i]:
                aac.append(k)
            aac.append(self.Hydrophobicity[i])
            aac.append(self.Hydrophilicity[i])
            aac.append(self.sidechains[i])
            aac.append(protein_eiip[i])
        return aac



    # Composition Transition Distribution (CTD) vector size(168)
    def str2n(self, seq, property):
        prot = copy.deepcopy(seq)
        for k, m in property.items():
            for index in m:
                prot = str.replace(prot, index, k)
        return prot

    def CalculateComposition(self, protseq, property):
        tprotseq = self.str2n(protseq, property)
        result = []
        N = len(tprotseq)
        a = round(float(tprotseq.count('1')) / N, 5)
        b = round(float(tprotseq.count('2')) / N, 5)
        c = round(float(tprotseq.count('3')) / N, 5)
        result.append(a)
        result.append(b)
        result.append(c)
        return result

    def CalculateTransition(self, protseq, property):
        tprotseq = self.str2n(protseq, property)

        result = []
        N = len(tprotseq)

        a = round(float(tprotseq.count('12') + tprotseq.count('21')) / (N - 1), 5)
        b = round(float(tprotseq.count('13') + tprotseq.count('31')) / (N - 1), 5)
        c = round(float(tprotseq.count('23') + tprotseq.count('32')) / (N - 1), 5)
        result.append(a)
        result.append(b)
        result.append(c)
        return result

    def CalculateDistribution(self, protseq, property):

        tprotseq = self.str2n(protseq, property)
        # Result={}
        Result = []
        Num = len(tprotseq)
        temp = ('1', '2', '3')
        for i in temp:
            num = tprotseq.count(i)
            ink = 1
            indexk = 0
            cds = []
            while ink <= num:
                indexk = str.find(tprotseq, i, indexk) + 1
                cds.append(indexk)
                ink = ink + 1

            if cds == []:
                a = 0
                b = 0
                c = 0
                d = 0
                e = 0
                Result.append(a)
                Result.append(b)
                Result.append(c)
                Result.append(d)
                Result.append(e)
            else:

                a = round(float(cds[0]) / Num * 100, 5)
                b = round(float(cds[int(math.floor(num * 0.25)) - 1]) / Num * 100, 5)
                c = round(float(cds[int(math.floor(num * 0.5)) - 1]) / Num * 100, 5)
                d = round(float(cds[int(math.floor(num * 0.75)) - 1]) / Num * 100, 5)
                e = round(float(cds[-1]) / Num * 100, 5)
                Result.append(a)
                Result.append(b)
                Result.append(c)
                Result.append(d)
                Result.append(e)
        return Result
    def CTD(self,seq, desc='ctd'):
        result=[]
        get_C = True
        get_T = True
        get_D = True
        desc = desc.lower()
        if 'c' not in desc:
            get_C = False
        if 't' not in desc:
            get_T = False
        if 'd' not in desc:
            get_D = False
        for i in range(len(self.AAG)):
            property = self.AAP[i]
            # named = self.AAG[i]
            if get_C:
                result.extend(self.CalculateComposition(seq, property))
            if get_T:
                result.extend(self.CalculateTransition(seq, property))
            if get_D:
                result.extend(self.CalculateDistribution(seq, property))
        return result
    # Conjoint Triad
    def s2n(self,seq):
        conj = {}
        for i in self._conj:
            for j in self._conj[i]:
                conj[j] = i

        res = seq
        for i in conj:
            res = res.replace(i, str(conj[i]))
        return res
    def conjoint(self,seq):
        result = []

        proteinnum = self.s2n(seq)
        for i in range(1, 8):
            for j in range(1, 8):
                for k in range(1, 8):
                    temp = str(i) + str(j) + str(k)
                    count = proteinnum.count(temp)
                    result.append(count)

        return result
    # normalize conjoint vector size (343)
    def norm_conjoint(self,seq):
        data=self.conjoint(seq)
        fmax = max(data)
        fmin = min(data)
        result=[]
        for l in range(0, len(data)):
            val = float((data[l] - fmin) / fmax)

            result.append(val)
        return result

# Quasi order 1 and 2 based on Grantham chemical distance  vector size(168)
    # Sequence-order-coupling number
    def cn(self, seq, rank, matrix):
        length = len(seq)
        _cn = []
        for i in range(length - rank):

            a = seq[i]
            b = seq[i + rank]
            # print(a)
            # print(b)
            # print(self.gr[a + b])
            _cn.append(math.pow(matrix[a + b], 2))


        return _cn

    def cnu(self, seq, max, matrix):
        length = len(seq)
        _cn = []
        for i in range(max):
            _cn.append(self.cn(seq, i + 1, matrix))
        # print(_cn)
        return _cn[0]

    def QSO(self, seq, max, w):
        data = []
        data_sw=[]
        for i in range(max):
            data.append(self.cnu(seq, i + 1, self.gr))
        for i in range(max):
            data_sw.append(self.cnu(seq, i + 1, self.sw))
        caa = self.AAC(seq)
        # print(caa)
        qs1=[]
        qs2=[]
        qs11 = []
        qs22 = []
        a_gr = 1 + (w * sum(data[0]))
        a_sw = 1 + (w * sum(data_sw[0]))

        for i, j in enumerate(self.AA):
            # print(j)
            qs1.append(round(caa[i] /a_sw , 10))
            qs11.append(round(caa[i] /a_gr , 10))
            # result['QSO:' + j] = round(caa[i] / a,10)
        for i in range(20, 20 + max):
            qs2.append(float(round(w * data[0][i - 20] / a_sw, 10)))
            qs22.append(float(round(w * data_sw[0][i - 20] / a_gr, 10)))
        result=qs1
        result.extend(qs2)
        result.extend(qs11)
        result.extend(qs22)
        return result

    def Rvalue(self,aa1, aa2, AADict, Matrix):
        return sum([(Matrix[i][AADict[aa1]] - Matrix[i][AADict[aa2]]) ** 2 for i in range(len(Matrix))]) / len(Matrix)

    def PAAC(self,seqStr, lambdaValue=30, w=0.05, **kw):




        theta = []
        code = []
        for n in range(1, lambdaValue + 1):
            theta.append(
                sum([self.Rvalue(seqStr[j], seqStr[j + n], self.AADict, self.AAProperty1) for j in
                     range(len(seqStr) - n)]) / (
                        len(seqStr) - n))
        myDict = {}
        for aa in AA:
            myDict[aa] = seqStr.count(aa)
        code = [myDict[aa] / (1 + w * sum(theta)) for aa in AA]
        code = code + [(w * j) / (1 + w * sum(theta)) for j in theta]

        return code

    def NMBroto(self,seq, nlag=30):
        pstd = np.std(self.AAidx, axis=1)
        pmean = np.average(self.AAidx, axis=1)

        for i in range(len(self.AAidx)):
            for j in range(len(self.AAidx[i])):
                AAidx[i][j] = (self.AAidx[i][j] - pmean[i]) / pstd[i]
        code = []
        N = len(seq)
        for prop in range(8):
            for n in range(1, nlag + 1):
                if len(seq) > nlag:
                    # if key is '-', then the value is 0
                    rn = sum([AAidx[prop][self.AADict.get(seq[j], 0)] * AAidx[prop][self.AADict.get(seq[j + n], 0)] for j in
                              range(len(seq) - n)]) / (N - n)
                else:
                    rn = 0
                code.append(rn)
        return code

    def Count(self,seq1, seq2):
        sum = 0
        for aa in seq1:
            sum = sum + seq2.count(aa)
        return sum

    def CTDC(self,seq):
        code = []
        for p in self.property:
            c1 = self.Count(group1[p], seq) / len(seq)
            c2 = self.Count(group2[p], seq) / len(seq)
            c3 = 1 - c1 - c2
            code = code + [c1, c2, c3]
        return code
    def CTDT(self,sequence):
        code = []
        aaPair = [sequence[j:j + 2] for j in range(len(sequence) - 1)]
        for p in self.property:
            c1221, c1331, c2332 = 0, 0, 0
            for pair in aaPair:
                if (pair[0] in group1[p] and pair[1] in group2[p]) or (pair[0] in group2[p] and pair[1] in group1[p]):
                    c1221 = c1221 + 1
                    continue
                if (pair[0] in group1[p] and pair[1] in group3[p]) or (pair[0] in group3[p] and pair[1] in group1[p]):
                    c1331 = c1331 + 1
                    continue
                if (pair[0] in group2[p] and pair[1] in group3[p]) or (pair[0] in group3[p] and pair[1] in group2[p]):
                    c2332 = c2332 + 1
            code = code + [c1221 / len(aaPair), c1331 / len(aaPair), c2332 / len(aaPair)]
        return code
    def BINARY(self,seq,sampleCount,sampleLen):
        gap = (len(seq) - sampleCount * sampleLen) / sampleCount  #
        res =''
        if gap < 0:
            # Sequence is too short
            res = seq + ('-' * (sampleCount * sampleLen - len(seq)))
        else:
            pos = 0
            for i in range(0, sampleCount):
                res += seq[pos:pos + sampleLen]
                pos += math.floor(gap)

        if (len(res) != sampleCount * sampleLen):  # Check id
            print('Err' + str(len(res)))
            exit()
        code=[]
        for a in seq:
            if a == '-':
                code.append([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
            for aa in self.AA:
                tag = 1 if a == aa else 0
                code.append(tag)
        return code
    def hybrid(self,seq,f1,f2):
        a = f1(seq)
        a.extend(f2(seq))
        return a

    def hybrid2(self, seq, f1, f2, f3):
        a = f1(seq)
        a.extend(f2(seq))
        a.extend(f3(seq))
        return a

    def hybrid3(self, seq, f1, f2, f3, f4):
        a = f1(seq)
        a.extend(f2(seq))
        a.extend(f3(seq))
        a.extend(f4(seq))
        return a

