from random import randint

import sys
import time
import csv
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
        #self.contigs = generateContigs(self.kmers)
        self.contig_len = len(self.contigs)
        self.n50, self.contig_max = calc_n50_and_max(self.contigs)

        # Genome
        self.genome = assemble_genome(self.kmers)
        self.genome_len = len(self.genome)

    def print(self):
        print("n50: ", self.n50)
        print("k:  ", self.k)
        print("percent: ", self.percent)
        print("genome: ", self.genome)
        return


def calc_n50_and_max(contigs):
    lengths = [len(contig) for contig in contigs]
    median = np.ma.median(lengths)
    return median, np.max(lengths)


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
    results = []
    k_min = 13 #15 # Much smaller than this will cause stack overflow issues! Linux works better than Windows
    k_max = int(len(fasta_lines[0]) * .95) # Max size is 80% of line length

    #if k_max > 50:
     #   k_max = 50

    max_percent = .15

    best_result = None
    step = .025 # Remove bottom 1% more each iteration

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
            results.append(result)
            print(k, result.n50, result.contig_len, result.genome_len, percent, sep="\t")


            if best_result == None or\
                    best_result.n50 < result.n50 or\
                    best_result.n50 == result.n50 and best_result.contig_max < result.contig_max:

                best_result = result

                percent += step
                if percent > max_percent:
                    break
            else:
                if best_result is None:
                    best_result = result
                break

            if not find_best_error:
                break

    return best_result, results



def assemble_genome(kmers):
    graph = construct_graph(kmers)
    degrees = calc_degrees(graph)
    start_node = find_start_node(degrees)
    eulerian_path = find_eulerian_path(start_node, graph, degrees, [])

    assembled_genome = eulerian_path[0][:-1]
    for node in eulerian_path:
        assembled_genome = assembled_genome + node[-1:]

    return assembled_genome

def write_results_to_csv(file_name, results, best_result):
    file = open(file_name, "w", newline='')
    csv_writer = csv.writer(file)

    header = ["K", "N50", "Contig Count", "Longest Contig Length", "Contigs", "Bottom % Removed", "Genome Length", "Genome"]
    rows = [["*The first result is the best result!"], header]

    row = []
    row.append(best_result.k)
    row.append(best_result.n50)
    row.append(best_result.contig_len)
    row.append(best_result.contig_max)
    row.append("\n".join(best_result.contigs))
    row.append(best_result.percent)
    row.append(best_result.genome_len)
    row.append(best_result.genome.strip())
    rows.append(row)

    rows.append(["-------","-------","-------","-------","-------","-------","-------","-------"])

    for r in results:
        row = []
        row.append(r.k)
        row.append(r.n50)
        row.append(r.contig_len)
        row.append(r.contig_max)
        row.append("\n".join(r.contigs))
        row.append(r.percent)
        row.append(r.genome_len)
        row.append(r.genome.strip())

        rows.append(row)

    csv_writer.writerows(rows)
    file.close()
    pass

def main():



    #fasta_lines = fasta_to_string("files/synthetic.example.noerror.small.fasta")
    #fasta_lines = fasta_to_string("files/synthetic.noerror.small.fasta")
    #fasta_lines = fasta_to_string("files/synthetic.noerror.large.fasta")
    #fasta_lines = fasta_to_string("files/example.data.fasta")
    #fasta_lines = fasta_to_string("files/real.error.large.fasta")
    #fasta_lines = fasta_to_string("files/real.error.small.fasta")

    tests = [
        ("files/example.data.fasta", "results/example.data.csv", False),
        ("files/synthetic.example.noerror.small.fasta", "results/synthetic.example.noerror.small.csv", False),
        ("files/synthetic.noerror.small.fasta", "results/synthetic.noerror.small.csv", False),
        ("files/synthetic.noerror.large.fasta", "results/synthetic.noerror.large.csv", False),
        ("files/real.error.small.fasta", "results/real.error.small.csv", True),
        ("files/real.error.large.fasta", "results/real.error.large.csv", True)
    ]


    test = tests[5]

    fasta_lines = fasta_to_string(test[0])
    best_result, results = random_search(fasta_lines, find_best_error=test[2])


    print("Writing to file:", test[1])
    write_results_to_csv(test[1], results, best_result)
    print("Done!")



if __name__ == "__main__":
    sys.setrecursionlimit(10000000)
    main()

