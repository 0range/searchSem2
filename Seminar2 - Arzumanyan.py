import numpy as np
import matplotlib.pyplot as plt
#%matplotlib inline
import ctypes
from time import time


def UnaryCode(n):
    return n + 1


def TruncatedBinaryCode(n, m):
    b = int(np.ceil(np.log(m + 0.) / np.log(2.0)))
    small = (n < (2 ** b - m))
    return (b - 1) * small + b * (1 - small)


def GolombCompress(data, m):
    res = np.sum(UnaryCode(data / m)) + np.sum(TruncatedBinaryCode(data % m, m))
    return res


def Generate():
    filenames = []
    for i in (0.001, 0.004, 0.007, 0.01, 0.04, 0.07, 0.1, 0.5, 0.7):
        data = np.random.geometric(i, 1000000)
        filename = str(int(1000 * i))
        filenames.append(filename)
        with open(filename, 'w') as outputfile:
            for el in data:
                outputfile.write(str(el) + '\n')
    return filenames


def SaveCompCoeffs(filenames):
    compcoeffs = np.zeros((len(filenames), 512), dtype=np.float64)
    thetime = time()
    for i in xrange(len(filenames)):
        data = []
        j = 0
        with open(filenames[i], 'r') as inputfile:
            for line in inputfile:
                j += 1
                if len(line) == 0:
                    break
                n = int(line.strip())
                data.append(n)
                #if (j+1) % 200000 == 0:
                #    print i, j
        print i, "data size is", len(data), time() - thetime
        original_size = len(data) * ctypes.sizeof(ctypes.c_int) * 8
        data = np.array(data)
        for m in xrange(512):
            coded_size = GolombCompress(data, m + 1)
            #if m % 100 == 0:
            #    print m, time() - thetime
            compcoeffs[i][m] = (original_size + 0.) / coded_size   
    return compcoeffs


def Vurhis(p):
    for m in xrange(512):
        if p ** (m + 1) + p ** (m + 2) <= 1 and 1 < p ** (m + 1) + p ** m:
            return m+1#, p, p ** (m + 1) + p ** (m + 2),  p ** (m + 1) + p ** m
    return 512


def main():
    thetime = time()
    filenames = Generate()
    print time() - thetime

    thetime = time()
    compcoeffs = SaveCompCoeffs(filenames)
    print time() - thetime

    # best coeff
    totalcoeffs = np.ravel(np.dot(np.ones((1,9), dtype=np.float64), compcoeffs))
    m_opt = np.argmax(totalcoeffs)
    print 'm=', m_opt + 1 , compcoeffs[:, m_opt], (totalcoeffs[m_opt]) / compcoeffs.shape[0]

    # longest list
    l = np.argmin(compcoeffs[:, m_opt])
    with open(filenames[l], 'r') as inputfile:
        data = []
        for line in inputfile:
            if len(line) == 0:
                 break
            n = int(line.strip())
            data.append(n)
    data = np.array(data)
    p = (data.shape[0] + 0.) / np.sum(data)
    p_g = (int(filenames[l]) + 0.) / 1000
    print 'p from data =', p, ', p from generator =', p_g

    # best coeff for longest
    m_opt_l = np.argmax(compcoeffs[l,:])
    print 'best m ', m_opt_l+1, compcoeffs[l,m_opt_l]

    # best coeff from theory
    print 'best m theoretical', Vurhis(p), compcoeffs[l, Vurhis(p) - 1]

