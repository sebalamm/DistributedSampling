#!/usr/bin/python3
import sys, os, re, getopt, numpy, random, argparse, itertools
import scipy, scipy.stats
import matplotlib.pyplot as plt
import subprocess 
import statistics 
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
# Plot data
#######################
def plot_scaling(data):
    print("##########################")
    print(bcolors.HEADER + "Generate Plot" + bcolors.ENDC)
    print("##########################")
    print("Generate plot...")

    algorithm = list(set([t[0] for t in data]))
    runs = sorted(list(set([int(t[1]) for t in data if int(t[1]) >= limit])))

    ax = plt.subplot(1,1,1)
    for alg in sorted(algorithm):
        means = [float(t[2])/int(t[1]) * 1e9 for t in data if t[0] == alg and int(t[1]) >= limit]
        stdev = [float(t[3])/int(t[1]) * 1e9 for t in data if t[0] == alg and int(t[1]) >= limit]
        ax.errorbar(runs, means, stdev, marker=next(marker), label=alg)

    handles, labels = ax.get_legend_handles_labels()
    new_handles = [handles[0], handles[5], handles[3], handles[1], handles[2], handles[4]]
    new_labels = [labels[0], labels[5], labels[3], labels[1], labels[2], labels[4]]
    plt.legend(new_handles, new_labels, loc=2)
    plt.grid(True)
    plt.title(r"Running time for $N=2^{30}$ averaged over $N/n$ repetitions")
    plt.xscale("log", basex=2)
    # plt.xlabel(r"$\log_{10}$(Sample size)")
    plt.xlabel("Sample size $n$")
    plt.ylabel("Time per sample (ns)")
    # plt.ylabel("Cycles per sample")
    plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)
    plt.tight_layout()
    pp.savefig()
    plt.close()


#######################
# Main
#######################

parser = argparse.ArgumentParser(description='Evaluate different sampling methods.')
parser.add_argument('-i', type=str, nargs='?')
parser.add_argument('-f', type=str, nargs='?')
parser.add_argument('-l', type=int, nargs='?')
args = parser.parse_args()

input_file = args.i
output_file = args.f
limit = args.l
print(limit)

# Init plots
marker = itertools.cycle((',', '+', '.', 'o', '*')) 
pp = PdfPages(output_file + '.pdf')

# Title
print("##########################")
print(bcolors.HEADER + "Sampling (Direct Comparison)" + bcolors.ENDC)
print("##########################")

running_times = []

with open(input_file) as infile:
    for line in infile:
        if (line.strip().split()[0] != "#"):
            running_times.append(line.strip().split())

# Create plot
plot_scaling(running_times)

# Finalize
pp.close()
