#!/usr/bin/python3
import sys, os, re, getopt, random, argparse
import subprocess 
import statistics 
import scipy, scipy.stats
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


#######################
# Degree Distribution
#######################
def plot_distribution(data):
    print("##########################")
    print(bcolors.HEADER + "Degree Distribution" + bcolors.ENDC)
    print("##########################")

    [N, n, m] = [2**45, 2**7, 2**44]

    print("Generate hypergeometric distribution...")
    rv = scipy.stats.hypergeom.rvs(N, n, m, size=10000000)
    plt.hist(rv, 50, label="Ideal", alpha=0.75)

    print("Generate plot...")
    plt.hist(data, 50, label="Experiments", alpha=0.75)
    plt.grid(True)
    plt.title("Hypergeometric distribution test")
    plt.xlabel("Hypergeometrix random values")
    plt.ylabel("#Occurences")
    plt.legend(loc=0)
    plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)
    pp.savefig()
    plt.close()


#######################
# Main
#######################

parser = argparse.ArgumentParser(description='Evaluate hypergeometric distribution.')
parser.add_argument('-f', type=str, nargs='?')
parser.add_argument('-i', type=str, nargs='?')
args = parser.parse_args()

in_file = args.i
out_file = args.f

# Init plots
pp = PdfPages(out_file + '.pdf')

# Network Generation
print("##########################")
print(bcolors.HEADER + "Hypergeometric Evaluation" + bcolors.ENDC)
print("##########################")

f = open(in_file)
lines = [int(x) for x in f.readlines()]

# Degree Distribution
plot_distribution(lines)

# Finalize
pp.close()
