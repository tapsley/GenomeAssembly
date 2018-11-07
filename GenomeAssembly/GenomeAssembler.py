from lib import fasta_to_string
from lib import generateContigs
from lib import construct_graph
from lib import calc_degrees
from lib import find_start_node
from lib import find_eulerian_path
from lib import window

def n50(data):
    return -1

def build_kmer_list(fata_lines,k,  error_cutoff):
    return []

def random_search():
    pass

def assemble_genome(kmers):
    return ""

def main():
    dna = fasta_to_string("files/synthetic.noerror.small.fasta")
    solutionFile = open("files/example.solution.fasta", 'r')
    solution = solutionFile.readline()
    solutionFile.close()
    #print "\n".join(dna)



    kmers = set()
    k = 20
    for d in dna:
        kmers.update(list(window(d, k)))


    contigs = generateContigs(list(kmers))

    print("Contigs")
    for kmer in contigs:
        print(kmer)


    print()
    print("Assembled Genomre")
    graph = construct_graph(kmers)
    degrees = calc_degrees(graph)
    start_node = find_start_node(degrees)
    eulerian_path = find_eulerian_path(start_node, graph, degrees, [])

    assembled_genome = eulerian_path[0][:-1]
    for node in eulerian_path:
        assembled_genome = assembled_genome + node[-1:]

    print(assembled_genome)
    print(len(assembled_genome))
    #print(solution)

main()

# Todo: calc n50
# Todo: remove errors by removing infrequent kmers
# Todo: refactor into functions with params
# Todo: automate testing and data collecting, use a random algorithm to find local max n50
