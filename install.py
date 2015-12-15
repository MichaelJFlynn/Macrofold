import linecache

def loadUNAfold(unafoldPath):
    f = open('include/loadData.h', 'w')
    hybridSS = open(unafoldPath + '/src/hybrid-ss.c', 'r')

    f.write("void loadData() {\n")

    f.write("}\n")

    f.close()

loadUNAfold("/Users/michaelflynn/unafold-scratch")
