from sequence import *

import math
import numpy as np
from scipy import stats
import random
import string
import csv

import pylab as plt

from collections import defaultdict
from itertools import chain


def getGapIndices(sequence):
    """
    Return indices where gaps are
    """
    return [i for i in range(len(sequence)) if sequence[i] == '-']

def generateRandomStringID(length = 5, alpha = None, seed  = None):
    """
    @param length: length of string
    @return: string
    """
    if seed:
        random.seed(seed)
    if not alpha:
        alpha = string.ascii_uppercase + string.digits
    return  ''.join(random.choice(alpha) for x in range(length))


def kde_pvalue(observation_value, background_list):
    """ Un-normalised Kernel Density Estimate P-value"""
    kde = stats.gaussian_kde(background_list)
    print observation_value
    return kde.integrate_box_1d(low=min(background_list),high=float(observation_value))

def plotKernalDistribution(background_list, observation_values = None, title = None):
    kde = stats.gaussian_kde(background_list)
    if observation_values:
        max_val = max([max([x[0] for x in observation_values]), max(background_list)])
    else:
        max_val = max(background_list)
    min_val = min(background_list)
    x = np.linspace(min_val, max_val, 500)
    fig, ax1 = plt.subplots(figsize = (1,1))
    plt.xlabel("Value")
    plt.title("Kernal Distribution")
    y = kde(x)
    if observation_values:
        for observed_val in observation_values:
            ax1.set_ylim(max(y))
            ax1.axvline(x = observed_val[0], ls = '--', color = 'r', linewidth = 3)
            p_val = 1 - kde_pvalue(observed_val[0], background_list)
            plt.text(observed_val[0], max(y)/2.0, observed_val[1] + \
                            " observed:" + str(p_val)[:7], rotation = 'vertical',
                            size = 16)
    if title:
        plt.title(title)
    ax1.set_ylim(0,max(y) + max(y)/10)
    ax1.plot(x, y, 'r', lw = 3) #Distribution function
    ax1.hist(background_list, normed=1,alpha=.3) # Histogram
    plt.show()

def fromKeys(iterable, value):
    """
    defaultdict.fromkeys does not work they way I hoped.
    This just makes sure all values in iterable are in the
    dict.
    """
    d = defaultdict(float)
    for k in iterable:
        d[k] = value
    return d

def mean(X):
    sum_val = 0
    for x in X:
        sum_val += x
    return sum_val/len(X)

def meanvar(X):
    """ The mean and variance of the sample. """
    mu = mean(X)
    dev = 0
    for x in X:
        dev += (x - mu) * (x - mu)
    return mu, dev / len(X)

def getZScore(X, sample):
    (mu, var) = meanvar(X)
    return (sample - mu) / math.sqrt(var)
    
def getZScores(X):
    (mu, var) = meanvar(X)
    Y = [((x - mu) / math.sqrt(var)) for x in X]
    return Y

def readDelimitedFile(filename, d_sym = '\t', skip_headings = False):
    """
    @param filename: Delimited file to read
    @skip_headings: Set to True to skip first heading line
    @d_sym: delimiter symbol used, default is tab
    @return: Columns -> dict()
    """
    columns = []
    with open(filename,'r') as fh:
        lines = csv.reader(fh, delimiter = d_sym)
        for line in lines:
            columns.append(tuple(c for c in line))
        if skip_headings: columns = columns[1:]
    return columns


def kl(p, q):
    """Kullback-Leibler divergence D(P || Q) for discrete distributions
 
    Parameters
    ----------
    p, q : array-like, dtype=float, shape=n
        Discrete probability distributions.
    """
    p = np.asarray(p, dtype = np.float)
    q = np.asarray(q, dtype = np.float)
 
    return np.sum(np.where(p != 0, p * np.log(p / q), 0))

def nestedDictToMatrix(nested_dict):
    """
    Returns a matrix from a nested dictionary i.e. {Row : {Columns}}.
    Assumes each row will have the same number of column values
    """
    rows = sorted(nested_dict.keys())
    columns = sorted(nested_dict[rows[0]].keys())
    matrix = np.zeros((len(rows),len(columns)))
    for row in range(len(rows)):
        for column in range(len(columns)):
            matrix[row][column] = nested_dict[rows[row]][columns[column]]
    return matrix

def prettyPrintMatrixSimple(matrix):
    string = ''
    for row in matrix:
        string += '\t'.join(map(str, ["%.3f" % i for i in row]))
        string += '\n'
    return string

#Broken/inflexible
def prettyPrintMatrix(matrix,  column_labels, row_labels, name = None, type_float= True):
    """
    Returns a "pretty printed" string of a matrix
    """
    string = ''
    if name: string += "%s\n" % name
    string += "\t"
    if column_labels:
        # string += "\t\t".join(["\t".join(map(str, r))
                                        # for r in column_labels])
        string += "\t\t".join(["\t\t".join(map(str, [r for r in column_labels]))])
    string += "\n"
    r_cnt = 0
    while r_cnt < len(matrix):
        if row_labels: string += "%s\t" % row_labels[r_cnt]
        for v in matrix[r_cnt]:
            if type_float: string += "%0.2f\t" % v
            else: string += "%i\t\t" % v
        string += "\n"
        r_cnt += 1
    # string += "\n".join(["\t".join(map(str, r)) for r in matrix])
    return string

def groupByAttribute(object_list, attribute, function = None):
    """
    Group objects based on the value of a specified attribute and return
    dictionary, using attribute value as key and object as value.
    If function specified, will apply function to attribute, effecting key.
    """
    if function:
        out = dict()
        for obj in object_list:
            temp = function(obj.__dict__[attribute])
            if temp in out:
                out[temp].append(obj)
            else: out[temp] = [obj]
        # out = dict((function(obj.__dict__[attribute]), obj) for obj in object_list)
    else: out = dict((obj.__dict__[attribute], obj) for obj in object_list)
    return out

def whereMutated(template_seq, mutated_seq):
    """
    Returns a list of the positions where a mutation
    has occurred and the amino acids observed at that position.
    """
    return dict((i, mutated_seq[i]) for i in
    range(len(template_seq)) if template_seq[i] != mutated_seq[i])

def phylipNameFix(filename):
    pass

def getTreeSiteDist(tree, sites):
    """
    @param tree: PhyloTree with alignment attached
    @param sites: sites of interest
    @return: return matrix of distribution of symbols
    """
    matrix = np.zeros(len(tree.alignment.alphabet), len(tree.alignment.alphabet))


def getFrequencyMatrices(sequences, mut_positions, sampled_aas, counts = True):
    """
    Returns the observation frequency matrices for a collection of sequences.
    mut_positions -> [mutation locations in sequence]*
    sampled_aas -> [The amino acids sampled during the generation of the sequences.]
    *don't correct for python indexing, this is done for you
    """
    out = {}
    num_aas = len(sampled_aas)
    if counts:
        datatype = int
    else:
        datatype = float
    for x in xrange(len(mut_positions)):
        y = x + 1
        mp1 = mut_positions[x]-1
        while y < len(mut_positions):
            matrix = np.zeros((num_aas, num_aas),
                        dtype = datatype)
            mp2 = mut_positions[y]-1
            curname = "%i,%i" % (mp1 + 1, mp2 + 1)
            for i in xrange(num_aas):
                for j in xrange(num_aas):
                    for s in sequences:
                        if s[mp1] == sampled_aas[i] and s[mp2] == sampled_aas[j]:
                            matrix[i][j] += 1
            y += 1
            out[curname] = matrix
    return out

def calcDeltaDeltaGibbs(e_value, R, temperature):
    """
    e.g. temp = 0.304, R = 1.987, eval = 102 -> 102
    @param e_value: activity value
    @param R: gas constant R
    @param temperature: absolute temperature
    @return: delta delta G (change in gibbs free energy)
    """
    if e_value <= 0: e_value = 1
    try:
        return -R * temperature * log(e_value)
    except:
        print "ERROR", e_value, R, temperature


def get_partial_matrix(amatrix, exclude_idx):
    """
    Return a sub matrix from the full matrix, excluding specified row and column
    :param exclude_idx: row/column index to exclude
    :return: partial matrix
    """
    N = len(amatrix)
    rows = [[i] for i in xrange(N) if i != exclude_idx]
    return amatrix[[rows], chain(*rows)]

if __name__ == "__main__":
    pass