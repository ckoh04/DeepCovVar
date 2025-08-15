#!/usr/bin/python
# Script for pre-processing data
# Author: Naveen Duhan
import numpy as np
from .features import *
import subprocess
import re, os
import math
import pandas as pd
from Bio import SeqIO

# Utility Functions
def softmax(z):
    assert len(z.shape) == 2
    s = np.max(z, axis=1)
    s = s[:, np.newaxis] # necessary step to do broadcasting
    e_x = np.exp(z - s)
    div = np.sum(e_x, axis=1)
    div = div[:, np.newaxis] # necessary step to do broadcasting
    return e_x / div
def shuffleTwoArrays(a, b):
    # Generate the permutation index array.
    permutation = np.random.permutation(a.shape[0])
    # Shuffle the arrays by giving the permutation in the square brackets.
    shuffled_a = a[permutation]
    shuffled_b = b[permutation]
    return shuffled_a, shuffled_b

def preprocess(_seq,feat, classLabel,size):
    tmp = []
    f=FEATURE()

    totalSequences = 0
    _classLabels = []
    seqID = []
    for s in _seq:
        res = ''
        totalSequences += 1
        seqS = str(s.seq)
        seqStr=re.sub('[^ARNDCQEGHILKMFPSTWYV-]', '', ''.join(seqS).upper())
        res+=seqStr
        if feat=="AAC":
            tmp.append(f.AAC(res))
            _classLabels.append(classLabel)
            seqID.append(s.id)
        elif feat=="DPC":
            tmp.append(f.DPC(res))
            _classLabels.append(classLabel)
            seqID.append(s.id)
        elif feat=="DPCP":
            tmp.append(f.DPCP(res))
            _classLabels.append(classLabel)
            seqID.append(s.id)
        elif feat=="TPC":
            tmp.append(f.TPC(res))
            _classLabels.append(classLabel)
            seqID.append(s.id)
        elif feat=="AACN":
            tmp.append(f.AAC_NEW(res))
            _classLabels.append(classLabel)
            seqID.append(s.id)
        elif feat=="AACP":
            tmp.append(f.AAC_P(res))
            _classLabels.append(classLabel)
            seqID.append(s.id)
        elif feat=="CKSAAP":
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
            tmp.append(f.QSO(res,30,0.1))
            _classLabels.append(classLabel)
            seqID.append(s.id)
        elif feat == "BINARY":
            tmp.append(f.BINARY(res,4,150))
            _classLabels.append(classLabel)
            seqID.append(s.id)
        elif feat == "DPC_AAC":
            tmp = tmp
            tmp.append(f.hybrid(res, f.DPC, f.AAC))
            _classLabels.append(classLabel)
            seqID.append(s.id)
        elif feat == "DPC_NM":
            tmp = tmp
            tmp.append(f.hybrid(res, f.DPC, f.NMBroto))
            _classLabels.append(classLabel)
            seqID.append(s.id)
        elif feat == "DPC_NM_CONJOINT":
            tmp = tmp
            tmp.append(f.hybrid2(res, f.DPC, f.conjoint, f.NMBroto))
            _classLabels.append(classLabel)
            seqID.append(s.id)
        elif feat == "TPC_DPCP":
            tmp = tmp
            tmp.append(f.hybrid(res, f.DPCP, f.TPC))
            _classLabels.append(classLabel)
            seqID.append(s.id)
        elif feat == "DPC_NM_CKSAAP": # features  2400 + 640 = 3040
            tmp = tmp
            tmp.append(f.hybrid2(res, f.DPC, f.NMBroto, f.CKSAAP))
            _classLabels.append(classLabel)
            seqID.append(s.id)
        elif feat == "DPC_NM_CONJOINT_CKSAAP":# features  400 + 240 + 343 + 2400 = 3383
            tmp = tmp
            tmp.append(f.hybrid3(res, f.DPC, f.conjoint, f.NMBroto, f.CKSAAP))
            _classLabels.append(classLabel)
            seqID.append(s.id)
        elif feat == "DPCP_NM_CONJOINT_CTD": # features  1200 + 240 + 343 + 168 = 1951
            tmp = tmp
            tmp.append(f.hybrid3(res, f.DPCP, f.conjoint, f.NMBroto, f.CTD))
            _classLabels.append(classLabel)
            seqID.append(s.id)
        elif feat == "CKSAAP_NM_CONJOINT_CTD": # features  2400 + 240 + 343 + 168 = 3151
            tmp = tmp
            tmp.append(f.hybrid3(res, f.CKSAAP, f.conjoint, f.NMBroto, f.CTD))
            _classLabels.append(classLabel)
            seqID.append(s.id)
        elif feat == "CKSAAP_NM_CTD": # features  2400 + 240 + 168 = 2808
            tmp = tmp
            tmp.append(f.hybrid2(res, f.CKSAAP, f.NMBroto, f.CTD))
            _classLabels.append(classLabel)
            seqID.append(s.id)
        elif feat == "NM_CONJOINT_CTD": # features 240 + 343 + 168 = 751
            tmp = tmp
            tmp.append(f.hybrid2(res, f.conjoint, f.NMBroto, f.CTD))
            _classLabels.append(classLabel)
            seqID.append(s.id)
        elif feat == "CKSAAP_NM": # features  2400 + 240 = 2640
            tmp = tmp
            tmp.append(f.hybrid(res, f.CKSAAP, f.NMBroto))
            _classLabels.append(classLabel)
            seqID.append(s.id)                           
    tmp=np.array(tmp)

    aa=np.zeros((len(tmp),int(size)))
    for i,s in enumerate(tmp):
        aa[i]=list(s)

    return aa, _classLabels, seqID

def preprocessdf(_seq,feat, classLabel,size):
    tmp = []
    f=FEATURE()

    totalSequences = 0
    _classLabels = []
    seqID = []
    for s in _seq:
        res = ''
        totalSequences += 1
        seqS = str(s)
        seqStr=re.sub('[^ARNDCQEGHILKMFPSTWYV-]', '', ''.join(seqS).upper())
        res+=seqStr
        if feat=="AAC":
            tmp.append(f.AAC(res))
            _classLabels.append(classLabel)
           
        elif feat=="DPC":
            tmp.append(f.DPC(res))
            _classLabels.append(classLabel)
            
        elif feat=="DPCP":
            tmp.append(f.DPC(res))
            _classLabels.append(classLabel)
          
        elif feat=="TPC":
            tmp.append(f.TPC(res))
            _classLabels.append(classLabel)
            
        elif feat=="AACN":
            tmp.append(f.AAC_NEW(res))
            _classLabels.append(classLabel)
            
        elif feat=="AACP":
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
            tmp.append(f.QSO(res,30,0.1))
            _classLabels.append(classLabel)
            
        elif feat == "BINARY":
            tmp.append(f.BINARY(res,4,150))
            _classLabels.append(classLabel)
            
        elif feat == "hybrid":
            tmp = tmp
            tmp.append(f.hybrid(res, f.DPC, f.conjoint))
            _classLabels.append(classLabel)
           
        elif feat == "hybrid1":
            tmp = tmp
            tmp.append(f.hybrid(res, f.DPC, f.NMBroto))
            _classLabels.append(classLabel)
           
        elif feat == "hybrid2":
            tmp = tmp
            tmp.append(f.hybrid2(res, f.DPC, f.conjoint, f.NMBroto))
            _classLabels.append(classLabel)
           
    tmp=np.array(tmp)

    aa=np.zeros((len(tmp),int(size)))
    for i,s in enumerate(tmp):
        aa[i]=list(s)

    return aa, _classLabels

def result_stat(result):
    tn = result[0]
    print(tn)
    fp = result[1]#11
    print(fp)
    fn = result[2]#56
    print(fn)
    tp = result[3]
    print(tp)
    sens = tp / (tp + fn)
    spec = tn / (tn + fp)
    prec = tp / (tp + fp)
    npv = tn / (tn+fn)
    fpr =  fp/ (fp+tn)
    fdr = fp /(fp+tp)
    fnr = fn /(fn+tp)
    acc = (tp + tn) / (tp + tn + fp + fn)
    f1score = (2 * tp) / ((2 * tp) + fp + fn)
    # mcc = ((tp * tn) - (fp * fn)) / math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    ######################### MCC ##############################
    # Calculate numerator and denominator using numpy
    # Convert to numpy float64 to handle large numbers
    tp, tn, fp, fn = np.float64(tp), np.float64(tn), np.float64(fp), np.float64(fn)
    
    # Calculate numerator
    numerator = (tp * tn) - (fp * fn)
    
    # Calculate denominator
    denominator = (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)
    
    # Check if denominator is too small or zero, which might cause instability
    if denominator == 0:
        print("Warning: Denominator is zero. Returning MCC = 0.")
        return 0.0
    
    # Calculate MCC
    mcc = numerator / np.sqrt(denominator)
    print("MCC:",mcc) 

    return [tp, fn, tn, fp, sens, spec, prec,npv,fpr,fdr,fnr, acc, f1score, mcc]

def sbatch_slurm(job_name='pyseqRNA', jname='test',command='command',  dep=''):
    if dep != '': dep = '--dependency=afterok:{} --kill-on-invalid-dep=yes '.format(dep)
    # command="srun -l --gres=gpu:tesla:1 -N1 "+command
    sbatch_command = "sbatch -J {} -o {}.out -e {}.err --gpus=a100:1 --nodes=1 --exclusive --time=7-0  --partition gpu-preempt  --wrap='{}' {}".format(job_name, jname,jname, command, dep)
    print(sbatch_command)
    sbatch_response = subprocess.getoutput(sbatch_command)
    print(sbatch_response)

    job_id = sbatch_response.split(' ')[-1].strip()
    return job_id

# Function for Loading environment module

# if 'MODULEPATH' not in os.environ:
#          f = open(os.environ['MODULESHOME'] + "/init/.modulespath", "r")
#          path = []
#          for line in f.readlines():
#                  line = re.sub("#.*$", '', line)
#                  if line != '':
#                          path.append(line)
#          os.environ['MODULEPATH'] = ':'.join(path)

# if 'MODULEPATH' not in os.environ:
#          os.environ['LOADEDMODULES'] = ''

# def module(*args):
#          if type(args[0]) == type([]):
#                  args = args[0]
#          else:
#                  args = list(args)
#          (output, error) = subprocess.Popen(['/usr/bin/modulecmd','python'] +
#                          args, stdout=subprocess.PIPE).communicate()
#          exec (output)


def converter(infile=None, fclass=None):
    rec=[]
    seq=[]
    df=pd.DataFrame()
    for record in SeqIO.parse(infile, "fasta"):
        seq.append(str(record.seq).upper())
        rec.append(record.id)
    df = df.assign(name=rec)
    df = df.assign(sequence=seq)
    df['label']=fclass
    df=df[['label','name','sequence' ]]
    # df.to_csv(ofile,index=False,sep="\t")
    return df
