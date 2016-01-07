# for reading rawdata/loop.txt and outputting :
# -data/hairpin.csv
# -data/bulge.csv 
# -data/internal.csv
def convertHairpinBulgeInternal():
    # open files
    loopFile = open("rawdata/loop.txt", 'r')
    hairpinFile = open("data/hairpin.csv", 'w')
    bulgeFile = open("data/bulge.csv", 'w')
    internalFile = open("data/internal.csv", 'w')

    ## write headers
    hairpinFile.write("Length,Energy\n")
    bulgeFile.write("Length,Energy\n")
    internalFile.write("Length,Energy\n")

    skips = 4 # hardcoded number of lines to skip
    for line in loopFile:
        if(skips > 0):
            skips -= 1
            continue
        tokens = line.split()
        if tokens[1] is ".":
            continue
        # size column
        for f in [hairpinFile, bulgeFile, internalFile]:
            f.write(tokens[0])
            f.write(",")
        
        # internal energy column
        if tokens[1] is ".":
            internalFile.write("Inf")
        else:
            internalFile.write(tokens[1])
        internalFile.write("\n")

        # bulge energy column
        if tokens[2] is ".":
            bulgeFile.write("Inf")
        else:
            bulgeFile.write(tokens[2])
        bulgeFile.write("\n")

        # hairpin energy column
        if tokens[3] is ".":
            hairpinFile.write("Inf")
        else:
            hairpinFile.write(tokens[3])
        hairpinFile.write("\n")
    
    loopFile.close()
    hairpinFile.close()
    bulgeFile.close()
    internalFile.close()


def intToBase(num):
    if num == 0: 
        return 'A'
    elif num == 1:
        return 'C'
    elif num == 2:
        return 'G'
    elif num == 3:
        return 'U'

def convertStack():
    stackFile = open("rawdata/stack.txt", 'r')
    outputFile = open("data/stack.csv", 'w')

    outputFile.write("5primeOuter,5primeInner,3primeInner,3primeOuter,Energy\n")
    skip = 26
    outer5primeIndex = 0
    inner5primeIndex = 0
    for line in stackFile:
        if skip > 0:
            skip -= 1
            continue
        tokens = line.split()
        for i in range(len(tokens)):
            energy = tokens[i]
            if energy is ".":
                continue
                ##energy = "Inf"
            inner3prime = intToBase(i % 4)
            outer3prime = intToBase(i // 4)
            outer5prime = intToBase(outer5primeIndex)
            inner5prime = intToBase(inner5primeIndex)
            outputFile.write("{0},{1},{2},{3},{4}\n".format(outer5prime, inner5prime, inner3prime, outer3prime, energy))
        inner5primeIndex += 1
        if inner5primeIndex >= 4:
            inner5primeIndex = 0
            outer5primeIndex += 1
            skip = 11
    stackFile.close()
    outputFile.close()

def convertDangles():
    dangleFile = open("rawdata/dangle.txt", 'r')
    dangle5 = open("data/dangle5.csv", 'w')
    dangle3 = open("data/dangle3.csv", 'w')


    dangle5.write("5primeBase,3primeBase,DangleBase,Energy\n")
    dangle3.write("5primeBase,3primeBase,DangleBase,Energy\n")
    skip = 10
    fivePrime = True
    index5prime = 0

    for line in dangleFile:
        if skip > 0:
            skip -= 1
            continue
        tokens = line.split()
        for i in range(len(tokens)):
            energy = tokens[i]
            if energy is ".":
                continue
                #energy = "Inf"
            dangleBase = intToBase(i % 4)
            base3prime = intToBase(i // 4)
            base5prime = intToBase(index5prime)
            if fivePrime:
                dangle5.write("{0},{1},{2},{3}\n".format(base5prime, base3prime, dangleBase, energy))
            else:
                dangle3.write("{0},{1},{2},{3}\n".format(base5prime, base3prime, dangleBase, energy))
            
        index5prime += 1
        skip = 10
        if index5prime >= 4:
            fivePrime = False
            index5prime = 0
        
    dangleFile.close()
    dangle5.close()
    dangle3.close()
    

def convertTerminalStack():
    # open files
    inputFile = open("rawdata/tstack.txt", 'r')
    outputFile = open("data/tstack.csv", 'w')

    outputFile.write("5primeOuter,5primeInner,3primeInner,3primeOuter,Energy\n")
    skip = 26
    outer5primeIndex = 0
    inner5primeIndex = 0
    for line in inputFile:
        if skip > 0:
            skip -= 1
            continue
        tokens = line.split()
        for i in range(len(tokens)):
            energy = tokens[i]
            if energy is ".":
                continue;
                #energy = "Inf"
            inner3prime = intToBase(i % 4)
            outer3prime = intToBase(i // 4)
            outer5prime = intToBase(outer5primeIndex)
            inner5prime = intToBase(inner5primeIndex)
            outputFile.write("{0},{1},{2},{3},{4}\n".format(outer5prime, inner5prime, inner3prime, outer3prime, energy))
        inner5primeIndex += 1
        if inner5primeIndex >= 4:
            inner5primeIndex = 0
            outer5primeIndex += 1
            skip = 11
    inputFile.close()
    outputFile.close()
    

convertHairpinBulgeInternal()
convertTerminalStack()
convertDangles()
convertStack()
