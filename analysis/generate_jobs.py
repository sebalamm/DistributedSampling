#!/usr/bin/python
import sys, math, os, re, getopt, random, argparse

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
# Sampling procedures
#######################
def sample(n, N, k, P, iterations, method):
    output = "job_" + str(n) + "_" + str(N) + "_" + str(P)
    out = open(output, 'w')
    out.write("# @ job_name = sample_" + str(n) + "_" + str(N) + "_" + str(P) + "\n") 
    out.write("# @ commment = 'test sample weak-scaling'" + "\n")
    out.write("# @ error = $(job_name).$(job_id).out" + "\n")
    out.write("# @ output = $(job_name).$(job_id).out" + "\n")
    out.write("# @ environment = COPY_ALL" + "\n")
    out.write("# @ wall_clock_limit = 00:05:00" + "\n")
    out.write("# @ notification = error" + "\n")
    out.write("# @ notify_user = seba.lamm@gmail.com" + "\n")
    out.write("# @ job_type = blueqene" + "\n")
    out.write("# @ bg_size = 32" + "\n")
    out.write("# @ queue" + "\n")
    out.write("\n") 
    out.write("runjob --np " + str(P) + " --ranks-per-node 1 : ~/code/build/run_method" + method + " -n " + str(n) + " -N " + str(N) + " -k " + str(k) + " -i " + str(iterations) + " -output " + output)
    out.close()

#######################
# Main
#######################

parser = argparse.ArgumentParser(description='Generate jobs for parallel sampling.')
parser.add_argument('-P', type=int, nargs='+')
parser.add_argument('-n', type=int, nargs='+')
parser.add_argument('-N', type=int, nargs='?')
parser.add_argument('-k', type=int, nargs='?')
parser.add_argument('-i', type=int, nargs='?')
args = parser.parse_args()

i = args.i
P = args.P
n = args.n
N = args.N
k = args.k

# Title
print("##########################")
print(bcolors.HEADER + "Sampling (Generate Jobs)" + bcolors.ENDC)
print("##########################")

for sample_size in n:
	for p in P:
	    total_n = int(math.log(p, 2)) + sample_size
	    iterations = int((2**i)/(2**total_n))
	    sample(total_n, N, k, p, iterations, "P")


