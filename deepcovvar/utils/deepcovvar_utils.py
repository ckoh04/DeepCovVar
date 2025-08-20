#!/usr/bin/python
# Script for pre-processing data
# Author: Naveen Duhan
import numpy as np
from .features import *
import re

def preprocess(_seq, feat, classLabel, size):
    """
    Preprocess sequences for feature extraction.
    
    Args:
        _seq: List of BioPython sequence objects
        feat: Feature type (AAC, DPC, DPCP, TPC, etc.)
        classLabel: Class label for the sequences
        size: Expected feature size
        
    Returns:
        tuple: (features_array, class_labels, sequence_ids)
    """
    tmp = []
    f = FEATURE()

    totalSequences = 0
    _classLabels = []
    seqID = []
    
    for s in _seq:
        res = ''
        totalSequences += 1
        seqS = str(s.seq)
        seqStr = re.sub('[^ARNDCQEGHILKMFPSTWYV-]', '', ''.join(seqS).upper())
        res += seqStr
        
        if feat == "AAC":
            tmp.append(f.AAC(res))
            _classLabels.append(classLabel)
            seqID.append(s.id)
        elif feat == "DPC":
            tmp.append(f.DPC(res))
            _classLabels.append(classLabel)
            seqID.append(s.id)
        elif feat == "DPCP":
            tmp.append(f.DPCP(res))
            _classLabels.append(classLabel)
            seqID.append(s.id)
        elif feat == "TPC":
            tmp.append(f.TPC(res))
            _classLabels.append(classLabel)
            seqID.append(s.id)
        elif feat == "AACN":
            tmp.append(f.AAC_NEW(res))
            _classLabels.append(classLabel)
            seqID.append(s.id)
        elif feat == "AACP":
            tmp.append(f.AAC_P(res))
            _classLabels.append(classLabel)
            seqID.append(s.id)
        elif feat == "CKSAAP":
            tmp.append(f.CKSAAP(res))
            _classLabels.append(classLabel)
            seqID.append(s.id)
        elif feat == "CTD":
            tmp.append(f.CTD(res))
            _classLabels.append(classLabel)
            seqID.append(s.id)
        elif feat == "norm_conjoint":
            tmp.append(f.norm_conjoint(res))
            _classLabels.append(classLabel)
            seqID.append(s.id)
        elif feat == "PAAC":
            tmp.append(f.PAAC(res))
            _classLabels.append(classLabel)
            seqID.append(s.id)
        elif feat == "CTDC":
            tmp.append(f.CTDC(res))
            _classLabels.append(classLabel)
            seqID.append(s.id)
        elif feat == "CTDT":
            tmp.append(f.CTDT(res))
            _classLabels.append(classLabel)
            seqID.append(s.id)
        elif feat == "NMBroto":
            tmp.append(f.NMBroto(res))
            _classLabels.append(classLabel)
            seqID.append(s.id)
        elif feat == "conjoint":
            tmp.append(f.conjoint(res))
            _classLabels.append(classLabel)
            seqID.append(s.id)
        elif feat == "QSO":
            tmp.append(f.QSO(res, 30, 0.1))
            _classLabels.append(classLabel)
            seqID.append(s.id)
        elif feat == "BINARY":
            tmp.append(f.BINARY(res, 4, 150))
            _classLabels.append(classLabel)
            seqID.append(s.id)
        elif feat == "DPC_AAC":
            tmp.append(f.hybrid(res, f.DPC, f.AAC))
            _classLabels.append(classLabel)
            seqID.append(s.id)
        elif feat == "DPC_NM":
            tmp.append(f.hybrid(res, f.DPC, f.NMBroto))
            _classLabels.append(classLabel)
            seqID.append(s.id)
        elif feat == "DPC_NM_CONJOINT":
            tmp.append(f.hybrid2(res, f.DPC, f.conjoint, f.NMBroto))
            _classLabels.append(classLabel)
            seqID.append(s.id)
        elif feat == "TPC_DPCP":
            tmp.append(f.hybrid(res, f.DPCP, f.TPC))
            _classLabels.append(classLabel)
            seqID.append(s.id)
        elif feat == "DPC_NM_CKSAAP":
            tmp.append(f.hybrid2(res, f.DPC, f.NMBroto, f.CKSAAP))
            _classLabels.append(classLabel)
            seqID.append(s.id)
        elif feat == "DPC_NM_CONJOINT_CKSAAP":
            tmp.append(f.hybrid3(res, f.DPC, f.conjoint, f.NMBroto, f.CKSAAP))
            _classLabels.append(classLabel)
            seqID.append(s.id)
        elif feat == "DPCP_NM_CONJOINT_CTD":
            tmp.append(f.hybrid3(res, f.DPCP, f.conjoint, f.NMBroto, f.CTD))
            _classLabels.append(classLabel)
            seqID.append(s.id)
        elif feat == "CKSAAP_NM_CONJOINT_CTD":
            tmp.append(f.hybrid3(res, f.CKSAAP, f.conjoint, f.NMBroto, f.CTD))
            _classLabels.append(classLabel)
            seqID.append(s.id)
        elif feat == "CKSAAP_NM_CTD":
            tmp.append(f.hybrid2(res, f.CKSAAP, f.NMBroto, f.CTD))
            _classLabels.append(classLabel)
            seqID.append(s.id)
        elif feat == "NM_CONJOINT_CTD":
            tmp.append(f.hybrid2(res, f.conjoint, f.NMBroto, f.CTD))
            _classLabels.append(classLabel)
            seqID.append(s.id)
        elif feat == "CKSAAP_NM":
            tmp.append(f.hybrid(res, f.CKSAAP, f.NMBroto))
            _classLabels.append(classLabel)
            seqID.append(s.id)
    
    tmp = np.array(tmp)
    aa = np.zeros((len(tmp), int(size)))
    
    for i, s in enumerate(tmp):
        aa[i] = list(s)

    return aa, _classLabels, seqID


def preprocessdf(_seq, feat, classLabel, size):
    """
    Preprocess sequences from DataFrame for feature extraction.
    
    Args:
        _seq: List of sequence strings
        feat: Feature type (AAC, DPC, DPCP, TPC, etc.)
        classLabel: Class label for the sequences
        size: Expected feature size
        
    Returns:
        tuple: (features_array, class_labels)
    """
    tmp = []
    f = FEATURE()

    totalSequences = 0
    _classLabels = []
    
    for s in _seq:
        res = ''
        totalSequences += 1
        seqS = str(s)
        seqStr = re.sub('[^ARNDCQEGHILKMFPSTWYV-]', '', ''.join(seqS).upper())
        res += seqStr
        
        if feat == "AAC":
            tmp.append(f.AAC(res))
            _classLabels.append(classLabel)
        elif feat == "DPC":
            tmp.append(f.DPC(res))
            _classLabels.append(classLabel)
        elif feat == "DPCP":
            tmp.append(f.DPC(res))
            _classLabels.append(classLabel)
        elif feat == "TPC":
            tmp.append(f.TPC(res))
            _classLabels.append(classLabel)
        elif feat == "AACN":
            tmp.append(f.AAC_NEW(res))
            _classLabels.append(classLabel)
        elif feat == "AACP":
            tmp.append(f.AAC_P(res))
            _classLabels.append(classLabel)
        elif feat == "CTD":
            tmp.append(f.CTD(res))
            _classLabels.append(classLabel)
        elif feat == "norm_conjoint":
            tmp.append(f.norm_conjoint(res))
            _classLabels.append(classLabel)
        elif feat == "PAAC":
            tmp.append(f.PAAC(res))
            _classLabels.append(classLabel)
        elif feat == "CTDC":
            tmp.append(f.CTDC(res))
            _classLabels.append(classLabel)
        elif feat == "CTDT":
            tmp.append(f.CTDT(res))
            _classLabels.append(classLabel)
        elif feat == "NMBroto":
            tmp.append(f.NMBroto(res))
            _classLabels.append(classLabel)
        elif feat == "conjoint":
            tmp.append(f.conjoint(res))
            _classLabels.append(classLabel)
        elif feat == "QSO":
            tmp.append(f.QSO(res, 30, 0.1))
            _classLabels.append(classLabel)
        elif feat == "BINARY":
            tmp.append(f.BINARY(res, 4, 150))
            _classLabels.append(classLabel)
        elif feat == "hybrid":
            tmp.append(f.hybrid(res, f.DPC, f.conjoint))
            _classLabels.append(classLabel)
        elif feat == "hybrid1":
            tmp.append(f.hybrid(res, f.DPC, f.NMBroto))
            _classLabels.append(classLabel)
        elif feat == "hybrid2":
            tmp.append(f.hybrid2(res, f.DPC, f.conjoint, f.NMBroto))
            _classLabels.append(classLabel)
    
    tmp = np.array(tmp)
    aa = np.zeros((len(tmp), int(size)))
    
    for i, s in enumerate(tmp):
        aa[i] = list(s)

    return aa, _classLabels
