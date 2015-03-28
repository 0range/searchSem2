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
    for m in xrange(2 ** 15):
        if p ** (m + 1) + p ** (m + 2) <= 1 and 1 < p ** (m + 1) + p ** m:
            return m+1#, p, p ** (m + 1) + p ** (m + 2),  p ** (m + 1) + p ** m
    return 2 ** 15


def loadLists(filename):
    """Load index from a file, returns dict key = word number
    value = list of (document number, count of word in the document)"""
    index = []
    i = 0
    with open(filename, 'r') as index_file:
        for line in index_file:
            linelist = line.strip().split()[1:]
            current_list = []
            for j in range(0, len(linelist)):
                current_pair = linelist[j].strip().split(':')
                for k in xrange(int(current_pair[1])):
                    current_list.append(int(current_pair[0]))
            index.append(sorted(current_list))
    return index


def main():
    """thetime = time()
    filenames = Generate()
    print time() - thetime

    thetime = time()
    compcoeffs = SaveCompCoeffs(filenames)
    print time() - thetime

    totalcoeffs = np.ravel(np.dot(np.ones((1,9), dtype=np.float64), compcoeffs))
    m_opt = np.argmax(totalcoeffs)

    fig = plt.figure(figsize=(16,6))
    colors = ('b','g','r','c','m','y','w','b','g')
    styles = ('-','-','-','-','-','-','-','--','--')
    axes = fig.add_subplot(1, 1, 1, axisbg='black')
    for i in xrange(9):
        plt.plot(compcoeffs[i, :], c=colors[i], ls=styles[i], label = filenames[i])
    xlim = plt.xlim([0,512])
    plt.axvline(x=m_opt, ymin=0, ymax=14, color = 'white')
    plt.legend(framealpha=0.2)
    """

    # now the task - optimal m for whole index
    thetime = time()
    index = loadLists('index.txt')
    print time() - thetime

    # compressing
    thetime = time()
    compressed_index = []
    empties = []
    for i in xrange(len(index)):
        if len(index[i]) > 0:
            curr_line = np.zeros((len(index[i]),), dtype=int)
            curr_line[0] = index[i][0]
            for j in xrange(1, len(index[i])):
                curr_line[j] = index[i][j] - index[i][j-1]
            compressed_index.append(curr_line)
        else:
            empties.append(i)
    print time() - thetime

    thetime = time()
    full_sizes = np.zeros((len(compressed_index)), dtype=int)
    for i in xrange(len(compressed_index)):
        full_sizes[i] =  len(compressed_index[i]) * ctypes.sizeof(ctypes.c_int) * 8
    print time() - thetime

    thetime = time()
    comp_sum_coeffs_first = np.zeros((16))
    for i in xrange(len(compressed_index)):
        if i % 100000 == 0:
            print i, time() - thetime
        for m in xrange(16):
            coded_size = GolombCompress(compressed_index[i], 2 ** m)
            comp_sum_coeffs_first[m] += coded_size   
    print time() - thetime

    #plt.plot(comp_sum_coeffs_first)
    #lim = plt.ylim([0, 0.1 * 10**10])

    m_min_first = 2 ** (np.argmin(comp_sum_coeffs_first) - 1)
    m_max_first = 2 ** (np.argmin(comp_sum_coeffs_first) + 1)

    thetime = time()
    comp_sum_coeffs_second = np.zeros((12))
    for i in xrange(len(compressed_index)):
        if i % 100000 == 0:
            print i, time() - thetime
        for m in xrange(12):
            coded_size = GolombCompress(compressed_index[i], m_min_first + 2 ** m)
            comp_sum_coeffs_second[m] += coded_size   
    print time() - thetime

    m_min_second = m_min_first + 2 ** 7
    m_max_second = 2 ** 11

    thetime = time()
    comp_sum_coeffs_third = np.zeros((16))
    for i in xrange(len(compressed_index)):
        if i % 100000 == 0:
            print i, time() - thetime
        for m in xrange(16):
            coded_size = GolombCompress(compressed_index[i], m_min_second + m * (m_max_second - m_min_second) / 15)
            comp_sum_coeffs_third[m] += coded_size   
    print time() - thetime

    m_min_third = m_min_second + (np.argmin(comp_sum_coeffs_third) - 1) * (m_max_second - m_min_second) / 15
    m_max_third = m_min_second + (np.argmin(comp_sum_coeffs_third) + 1) * (m_max_second - m_min_second) / 15

    thetime = time()
    comp_sum_coeffs_last = np.zeros(16)
    for i in xrange(len(compressed_index)):
        if i % 100000 == 0:
            print i, time() - thetime
        for m in xrange(16):
            coded_size = GolombCompress(compressed_index[i], m_min_third + m * (m_max_third - m_min_third) / 15)
            comp_sum_coeffs_last[m] += coded_size 
    print time() - thetime

    

