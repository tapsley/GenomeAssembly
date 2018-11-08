from random import randint

import sys
import time

from lib import fasta_to_string
from lib import generateContigs
from lib import construct_graph
from lib import calc_degrees
from lib import find_start_node
from lib import find_eulerian_path
from lib import window


class assembly():
    def __init__(self, kmers):
        self.kmers = kmers

        self.contig_len = -1
        self.genome = ""

        # Contigs
        #self.contigs = generateContigs(self.kmers)
        self.contigs = ["ACGT"]
        self.contig_len = len(self.contigs)
        self.n50 = calc_n50(self.contigs)

        # Genome
        self.genome = assemble_genome(self.kmers)
        self.genome_len = len(self.genome)


def calc_n50(contigs):
    contigs.sort(key=len)
    # Even
    if len(contigs) % 2 == 0:
        left = int((len(contigs)/2))
        right = left + 1
        return int((len(contigs[left]) + len(contigs[right]) / 2))
    else: # Odd
        median = int(len(contigs) / 2)
        return len(contigs[median])

def fasta_to_kmer_counts(fasta_lines, k):
    kmer_counts = {}

    for line in fasta_lines:
        line_kmers = list(window(line, k))
        for kmer in line_kmers:
            if kmer not in kmer_counts:
                kmer_counts[kmer] = 1
            else:
                kmer_counts[kmer] += 1

    return kmer_counts

def build_kmer_list(fasta_lines,k , error_cutoff):

    kmer_counts = fasta_to_kmer_counts(fasta_lines, k)

    kmers = []
    for kmer, count in kmer_counts.items():
        if count > error_cutoff:
            kmers.append(kmer)

    return kmers

def kmer_count_average(fasta_lines):
    kmer_counts = fasta_to_kmer_counts(fasta_lines)

    sum = 0
    for _, count in kmer_counts.items():
        sum += count

    return sum / len(kmer_counts)

def random_search(fasta_lines, count=100):
    #average_kmer_count = kmer_count_average(fasta_lines)
    err_min = 0
    err_max = 0

    k_min = 5
    k_max = 30

    best_n50 = -1
    best_k = 0
    best_err = 0
    best_genome = ""

    for i in range(count):
        local_n50 = 0
        err = randint(err_min, err_max)
        k = randint(k_min, k_max)
        kmers = build_kmer_list(fasta_lines, k , err)
        result = assembly(kmers)

        if best_n50 < result.n50:
            best_n50 = result.n50
            best_k = k
            best_err = err
            best_genome = result.genome

    return best_n50, best_k, best_err, best_genome



def assemble_genome(kmers):
    graph = construct_graph(kmers)
    degrees = calc_degrees(graph)
    start_node = find_start_node(degrees)
    eulerian_path = find_eulerian_path(start_node, graph, degrees, [])

    assembled_genome = eulerian_path[0][:-1]
    for node in eulerian_path:
        assembled_genome = assembled_genome + node[-1:]

    return assembled_genome

def main():
    fasta_lines = fasta_to_string("files/synthetic.noerror.small.fasta")
    #fasta_lines = fasta_to_string("files/example.data.fasta")
    fasta_lines = fasta_to_string("files/real.error.large.fasta")


    start = time.time()
    n50, k, err, genome = random_search(fasta_lines, 1)
    print("n50: ", n50)
    print("k:  ", k)
    print("err: ", err)
    print("genome: ", genome)

    print("")
    print(time.time() - start)





sys.setrecursionlimit(1000000000)
main()

# Todo: calc n50
# Todo: remove errors by removing infrequent kmers
# Todo: refactor into functions with params
# Todo: automate testing and data collecting, use a random algorithm to find local max n50
