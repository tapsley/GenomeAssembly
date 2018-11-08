from random import randint

import sys
import time
import numpy as np

from lib import fasta_to_string
from lib import generateContigs
from lib import construct_graph
from lib import calc_degrees
from lib import find_start_node
from lib import find_eulerian_path
from lib import window

from contigGraph import generate_contigs_graph


class assembly():
    def __init__(self, kmers, k, percent):
        self.kmers = kmers
        self.k = k
        self.percent = percent

        self.contig_len = -1
        self.genome = ""

        # Contigs
        self.contigs = generate_contigs_graph(self.kmers)
        self.contig_len = len(self.contigs)
        self.n50 = calc_n50(self.contigs)

        # Genome
        self.genome = assemble_genome(self.kmers)
        self.genome_len = len(self.genome)

    def print(self):
        print("n50: ", self.n50)
        print("k:  ", self.k)
        print("percent: ", self.percent)
        print("genome: ", self.genome)
        return


def calc_n50(contigs):
    lengths = [len(contig) for contig in contigs]
    median = np.ma.median(lengths)
    return median


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

# kmers must be sorted!
def remove_bottom_percent(kmers, percent):
    start = int(len(kmers) * percent)


    v = [t[0] for t in kmers[start:]]
    return v


def random_search(fasta_lines, find_best_error=False):
    k_min = 15 # Much smaller than this will cause stack overflow issues! Linux works better than Windows
    k_max = int(len(fasta_lines[0]) * .95) # Max size is 80% of line length

    max_percent = .1

    best_result = None
    step = .001 # Remove bottom 1% more each iteration

    for k in range(k_min, k_max):
        percent = 0

        kmer_counts = fasta_to_kmer_counts(fasta_lines, k)
        keys = list(kmer_counts.keys())
        values = list(kmer_counts.values())
        kmer_count_list = list(zip(keys, values))

        kmer_count_list.sort(key=lambda x: x[1])

        while True:
            kmers = remove_bottom_percent(kmer_count_list, percent)

            result = assembly(kmers, k, percent)
            print(k, result.n50, result.percent, result.contig_len)


            if best_result == None or\
                    best_result.n50 < result.n50 or\
                    best_result.n50 == result.n50 and best_result.genome_len < result.genome_len:

                best_result = result

                percent += step
                if percent > max_percent:
                    break
            else:
                if best_result == None:
                    best_result = result
                break

            if not find_best_error:
                break

    return best_result



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
    #fasta_lines = fasta_to_string("files/synthetic.noerror.small.fasta")
    #fasta_lines = fasta_to_string("files/example.data.fasta")
    fasta_lines = fasta_to_string("files/real.error.large.fasta")
    #fasta_lines = fasta_to_string("files/real.error.small.fasta")



    start = time.time()
    result = random_search(fasta_lines, find_best_error=True)

    print()
    result.print()

    print("")
    print(time.time() - start)




if __name__ == "__main__":
    sys.setrecursionlimit(10000000)
    main()

# Todo: remove errors by removing infrequent kmers
# Todo: automate testing and data collecting, use a random algorithm to find local max n50
