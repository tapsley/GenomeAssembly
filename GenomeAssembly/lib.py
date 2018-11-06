"""
Given a string containing the file path to a fasta file
Return a string of genes
"""


def fasta_to_string(file_name):
    lines = []
    with open(file_name, 'r') as f:
        count = 0
        for line in f:
            count += 1
            if count % 2 == 0:  # get every other line
                line = line.strip()
                lines.append(line)
    f.close()
    return lines


# Given an array of dna sequences, returns an array of contigs

def generateContigs(dna):
    final_kmers = []
    pos = len(dna[0]) - 1
    oneInOneOut = []
    for i, kmer1 in enumerate(dna):
        # Look at each kmer and find the ones with only 1 in and 1 out
        prefix1 = kmer1[0:pos]
        suffix1 = kmer1[len(kmer1) - pos:len(kmer1)]
        ins = 0
        outs = 0
        for j, kmer2 in enumerate(dna):
            prefix2 = kmer2[0:pos]
            suffix2 = kmer2[len(kmer2) - pos:len(kmer2)]
            if prefix1 == prefix2:
                ins += 1
            if prefix1 == suffix2:
                outs += 1
        if ins == 1 and outs == 1:
            oneInOneOut.append(kmer1)
        else:
            final_kmers.append(kmer1)

    i = 0
    while len(oneInOneOut) > 0:
        # This is where we glue them together
        kmer = oneInOneOut[i]
        prefix = kmer[0:pos]
        for x in final_kmers:
            if prefix in x:
                final_kmers.append(x + kmer[pos:len(kmer)])
                final_kmers.remove(x)
                oneInOneOut.remove(kmer)
                i = -1
                break
        i += 1

    # for i, kmer in enumerate(sorted(final_kmers)):
    #    print kmer

    return final_kmers


# Given a start node, a deBruijn graph, the degrees, and an empty path array
# Returns an Eulerian Path

def find_eulerian_path(node, graph, degrees, path):
    path += [node] # Same as cycle.append(node)

    if node not in degrees or degrees[node][1] == 0:
        return path

    while len(graph[node]) > 0:
        temp_node = graph[node][0]
        graph[node].remove(temp_node)

        sub_path = find_eulerian_path( temp_node, graph, degrees, [])

        path = path[:1] + sub_path + path[1:]

    return path

# Given the degrees of a graph
# returns the start node

def find_start_node(degrees):
    start_node = '0'

    for node, degree in degrees.items():
        if degree[0] < degree[1]:
            start_node = node

    return start_node


# Given a deBruijn graph
# returns the degrees

def calc_degrees(graph):
    degrees = {}

    for node, neighbors in graph.items():
        degrees[node] = (0, len(neighbors))

    for _, neighbors in graph.items():
        for node in neighbors:
            if node in degrees:
                degrees[node] = (degrees[node][0] + 1, degrees[node][1])
    return degrees


def getPrefix(kmer):
    return kmer[:-1]

def getSuffix(kmer):
    return kmer[1:]

# Given a list of dna sequences
# returns the edges in the form of an adjacency list

def GetEdges(patterns):
    adjList = []
    for pattern in patterns:
        prefix = getPrefix(pattern)
        suffix = getSuffix(pattern)
        adjList.append([prefix, suffix])
    return adjList


# Given a list of edges
# returns a deBruijn graph

def construct_graph(dna):
    data = GetEdges(dna)
    graph = {}
    for thing in data:
        if thing[0] in graph:
            graph[thing[0]]= graph[thing[0]] + thing[1:]
        else:
            graph[thing[0]] = list(thing[1:])

    return graph
