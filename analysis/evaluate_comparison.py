#!/usr/bin/python3
import sys, os, re, getopt, numpy, random, argparse
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

    cycles = 1200000000

    runs = [int(t[1]) for t in data if t[0] == "D"]
    d_means = [float(t[2])/int(t[1]) * cycles for t in data if t[0] == "D"]
    d_stdev = [float(t[3]) for t in data if t[0] == "D"]
    plt.errorbar(runs, d_means, marker="^", label="D")

    h_means = [float(t[2])/int(t[1]) * cycles for t in data if t[0] == "H"]
    h_stdev = [float(t[3]) for t in data if t[0] == "H"]
    plt.errorbar(runs, h_means, marker="^", label="H")

    r_means = [float(t[2])/int(t[1]) * cycles for t in data if t[0] == "R"]
    r_stdev = [float(t[3]) for t in data if t[0] == "R"]
    plt.errorbar(runs, r_means, marker="^", label="R")

    plt.grid(True)
    plt.title(r"Running time for $N=2^{" + str(52) + "}$ averaged over $" + str(10) + r"$ repetitions")
    plt.xscale("log", basex=2)
    # plt.xlabel(r"$\log_{10}$(Sample size)")
    plt.xlabel("Sample size")
    # plt.ylabel("Time per sample(s)")
    plt.ylabel("Cycles per sample")
    plt.legend(loc=0)
    plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)
    pp.savefig()
    plt.close()


#######################
# Main
#######################

parser = argparse.ArgumentParser(description='Evaluate different sampling methods.')
parser.add_argument('-i', type=str, nargs='?')
parser.add_argument('-f', type=str, nargs='?')
args = parser.parse_args()

input_file = args.i
output_file = args.f

# Init plots
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
