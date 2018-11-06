from lib import fasta_to_string
from lib import generateContigs
from lib import construct_graph
from lib import calc_degrees
from lib import find_start_node
from lib import find_eulerian_path



def main():
    dna = fasta_to_string("files/example.data.fasta")
    solutionFile = open("files/example.solution.fasta", 'r')
    solution = solutionFile.readline()
    solutionFile.close()
    #print "\n".join(dna)

    contigs = generateContigs(dna)

    #print "Contigs"
    #for kmer in contigs:
    #    print kmer


    graph = construct_graph(dna)
    degrees = calc_degrees(graph)
    start_node = find_start_node(degrees)
    eulerian_path = find_eulerian_path(start_node, graph, degrees, [])

    assembled_genome = eulerian_path[0][:-1]
    for node in eulerian_path:
        assembled_genome = assembled_genome + node[-1:]

    print assembled_genome
    print solution

main()