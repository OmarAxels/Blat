# búa til index fyrir reference, notið t.d. k=11 og sleppið algengum k-merum
# framkvæma seed og extend fyrir allar raðir í transcript, þegar staðsetning er fundin þá að framkvæma alignment
# finna spliced alignment fyrir þessi transcript

# Gagnasettið er stutt reference röð, 2Mb úr erfðamengi mannsins og það á að framvæma samröðun með röðunum sem eru í transcripts.fa (sömu og eru í kallisto). 
# Skilið samröðunum sem blat finnur í úttak. Þið getið athugað hvernig þetta lítur út með því að nota http://genome.ucsc.edu/cgi-bin/hgBlat til að tékka ykkur af.

import numpy as np
import math

# Reads query sequences and genome subsequence from files
def readSequencesFromFiles(subsequenceFile, transcriptsFile):
    dataFile = open(subsequenceFile).readlines()
    transcriptsFile = open(transcriptsFile).readlines()
    data = dataFile[1].strip()
    transcripts = []
    for i in range(1, len(transcriptsFile), 2):
        transcripts.append(transcriptsFile[i].strip())
    return data, transcripts

# Encodes A, C, G and T into 0, 1, 2 and 3 respectively
def letterToNumber(letter):
    if(letter == "A"):
        return 0
    if(letter == "C"):
        return 1
    if(letter == "G"):
        return 2
    if(letter == "T"):
        return 3
    else:
        raise RuntimeError("Letter was", letter,'but letter must be A, C, G or T')

# Hashes a kmer into a number between 0 and 4^k - 1
def hashit(kmer):
    value = 0
    for i in range(len(kmer)):
        number = letterToNumber(kmer[i])
        number *= pow(4, i)
        value += number
    return value

# Creates a list mapping hashed k-mers found in data to their index in data
def createIndexedList(k, data, startingPosInData, checkForCommonKmers, maxHitsPerKmer):
    hashLength = pow(4,k)
    indexedList = [ [] for i in range(hashLength) ]
    commonKmers = np.zeros(hashLength)
    for i in range(0, len(data), k):
        kmer = data[i:i+k]
        if(len(kmer) == k):
            kmerHashed = hashit(kmer)
            if(commonKmers[kmerHashed] == 0):
                indexedList[kmerHashed].append(i + startingPosInData)
                if(checkForCommonKmers):
                    if(len(indexedList[kmerHashed]) > maxHitsPerKmer):
                        commonKmers[kmerHashed] = 1
                        indexedList[kmerHashed] = []
    return indexedList


# Builds a list of “hits” where the query and the target match.
# Each hit contains [kmer, query position, database position, diagonal position]. 
def getHits(query, data, indexedList, k, alignmentStage):
    hitList = []
    for queryPos in range(len(query)-k):
        kmer = query[queryPos:queryPos+k] 
        kmerHashed = hashit(kmer)
        kmerIndexes = indexedList[kmerHashed]
        # If the kmer exists in data we mark it as a hit
        if(kmerIndexes != []):
            # If we are in the alignment stage we extend the k-mer by 1 untill we find 
            # a unique kmer in the homologous area or k has reached 25
            if(alignmentStage and len(kmerIndexes) > 1):
                limit = 25
                hit = extendKmer(query, data, queryPos, kmerIndexes, k, 25)
                hitList.append(hit)      
            else:
                for databasePos in kmerIndexes:
                    diagonalPos = databasePos - queryPos
                    hitList.append([kmer, queryPos, databasePos, diagonalPos])
    return hitList

# We extend the k-mer by 1 until we find 
# a unique kmer in the homologous area or k has reached the limit
def extendKmer(query, data, queryPos, hits, k, maxSize):
    queryStartPos = queryPos
    directionBackwards = True 
    j = 0
    while(len(hits) > 0 and j + k <= maxSize):
        if(directionBackwards):
            j += 1
        else:
            k += 1
        hitList = []
        for i in range(len(hits)):   
            #Expand each hit in both directions allowing no mismatches
            if(directionBackwards):
                databaseStartPos = hits[i] - 1
                if(queryStartPos - j >= 0 and databaseStartPos >= 0):
                    if(query[queryStartPos - j] != data[databaseStartPos]):
                        hitList.append(databaseStartPos + 1)
                        directionBackwards = False 
                    else:
                        hitList.append(databaseStartPos)
                else:
                    directionBackwards = False 
                    j -= 1

            if(not directionBackwards):
                databaseEndPos = databaseStartPos + k + j + 1            
                if(query[queryStartPos + k] == data[databaseEndPos]):
                    hitList.append(databaseStartPos)
            
        if(hitList == []):
            hits = [hits[0]]
            break
        else:
            hits = hitList
    
    hit = [query[queryStartPos-j:queryStartPos+k], queryStartPos-j, hits[0], hits[0] - queryStartPos + j]

    return hit


# Creates empty buckets
def createBuckets(size, data):
    bucketNumber = math.ceil(len(data)/size)
    buckets = []
    for i in range(bucketNumber):
        buckets.append([])
    return buckets

# Sorts buckets by the diagonal position of the hits
def sortBuckets(buckets):
    for bucket in buckets:
        bucket.sort(key=lambda x: x[3])
# Sorts proto-clumps by the database position of the hits
def sortProtoClumps(protoClumps):
    for protoClump in protoClumps:
        protoClump.sort(key=lambda x: x[2])
# Sorts hitList by the database position of the hits
def sortHitList(hitList):
    hitList.sort(key=lambda x: x[2])       

# Puts hits into buckets based on their database position
# and sorts each bucket by the diagonal position of the hits
def putHitsIntoBuckets(size, data, hitList):
    buckets = createBuckets(size, data)
    for hit in hitList:
        index = hit[2] // size
        buckets[index].append(hit)
    sortBuckets(buckets)
    return buckets

# Builds proto-clumps where each clump contains hits that share the same diagonal position
def buildProtoClumps(buckets):
    protoClumps = []
    protoClumpInUse = False
    for bucket in buckets:
        for i in range(len(bucket)-1):
            if(bucket[i][3] == bucket[i+1][3]):
                if(protoClumpInUse):
                    protoClump.append(bucket[i+1])
                else:
                    protoClump = [bucket[i], bucket[i+1]]
                    protoClumpInUse = True
            else:
                if(protoClumpInUse):
                    protoClumps.append(protoClump)
                    protoClumpInUse = False
    sortProtoClumps(protoClumps)
    return protoClumps

# Merges clumps that are within merge distance (300)
def mergeClumps(clumps, mergeDistance):
    for i in range(len(clumps) - 1):
        lastClumpIndex = len(clumps[i]) - 1 #index of last hit in clump i
        lastClump = clumps[i][lastClumpIndex]
        firstClump =  clumps[i+1][0]
        distanceBetweenClumps = firstClump[2] - lastClump[2]
        if(distanceBetweenClumps < 0):
            raise RuntimeError("Order of clumps not as expected")
        if(distanceBetweenClumps <= 300):
            clumps[i+1] = clumps[i] + clumps[i+1]
    return clumps

# Builds clumps from protoclumps so that hits within clumps are
# within window limit (W = 100) and each clump has at least minimumHits (2)
# Clumps within merge distance (300) are merged
def buildClumps(protoClumps, W, minimumHits):
    mergeDistance = 300
    clumps = []
    clumpInUse = False
    for protoClump in protoClumps:
        for i in range(len(protoClump) - 1):
            proximity = protoClump[i+1][2] - protoClump[i][2] 
            if(proximity < W):
                if(clumpInUse):
                    clump.append(protoClump[i+1])
                else:
                    clump = [protoClump[i], protoClump[i+1]]
                    clumpInUse = True
            # Add clump to list if next hit is further than
            # W away or if there are no more hits
            if(proximity > W or i == len(protoClump) - 2):
                if(clumpInUse):
                    clumps.append(clump)
                    clumpInUse = False
    clumps = mergeClumps(clumps, mergeDistance)
    return clumps

# Builds a homologous region from each clump and adds 500 additional bases
# on each side to form the final homologous region
def buildHomologousRegions(clumps, data):
    additionalBases = 500
    homologousRegions = []
    for clump in clumps:
        startPos = int(clump[0][2])
        if(startPos > additionalBases):
            startPos -= additionalBases
        else:
            startPos = 1     
        endPos = int(clump[len(clump)-1][2])
        if(endPos + additionalBases < len(data)):
            endPos += additionalBases
        else:
            endPos = len(data)   
        homologousRegions.append([startPos, endPos+1])
    return homologousRegions

# Removes hits that are empty
def removeEmptyHits(hitList):
    i = 0
    while(i < len(hitList)):
        if(i < len(hitList)):
            if(hitList[i] == []):
                del hitList[i]
            else:
                i += 1
        else:
            break
    return hitList


# Merges hits that overlap and thus reduces the size of the hitlist
def mergeHits(hitList):
    sortHitList(hitList)
    for i in range(len(hitList)-1):
        databasePos1 = hitList[i][2]
        databasePos2 = hitList[i+1][2]
        k1 = len(hitList[i][0])
        k2 = len(hitList[i+1][0])
        databaseEndPos1 = databasePos1 + k1
        databaseEndPos2 = databasePos2 + k2
        diagonalPos1 = hitList[i][3]
        diagonalPos2 = hitList[i+1][3]
        if(databaseEndPos1 >= databasePos2 and diagonalPos1 == diagonalPos2):
            if(databaseEndPos1 >= databaseEndPos2):
                hitList[i+1] = hitList[i]
                hitList[i] = []
            else:
                kmer = hitList[i][0][:databasePos2-databasePos1] + hitList[i+1][0]
                queryPos1 = hitList[i][1]
                hitList[i+1] = [kmer, queryPos1, databasePos1,  diagonalPos1]
                hitList[i] = []

    hitList = removeEmptyHits(hitList)
    return hitList

def spliceAlignments(alignments, query, data):
    splicedAlignments = []
    for i in range(len(alignments)-1):
        queryPos1 = int(alignments[i][1])
        queryEndPos1 = queryPos1 + len(alignments[i][0])
        queryPos2 = int(alignments[i+1][1])
        queryEndPos2 = queryPos2 + len(alignments[i+1][0])
        databasePos1 = alignments[i][2]
        databasePos2 = alignments[i+1][2]
        databaseEndPos1 = databasePos1 + len(alignments[i][0])
        databaseEndPos2 = databasePos2 + len(alignments[i+1][0])
        if(queryEndPos1 >= queryPos2 and queryPos1 < queryPos2):
            score = queryEndPos2 - queryPos1
            queryStart = queryPos1 + 1
            queryEnd = queryStart + score
            qsize = len(query)
            identity = "100%"
            chrom = "chr12"
            dataBaseStart = 53000001 + databasePos1
            databaseEnd = 53000001 + databaseEndPos2
            span = databaseEnd - dataBaseStart + 1
            splicedAlignment = [score, queryStart, queryEnd, qsize, identity, chrom, dataBaseStart, databaseEnd, span]
        else:
            score = queryEndPos1 - queryPos1
            queryStart = queryPos1 + 1
            queryEnd = queryStart + score
            qsize = len(query)
            identity = "100%"
            chrom = "chr12"
            dataBaseStart = 53000001 + databasePos1
            databaseEnd = 53000001 + databaseEndPos1
            span = databaseEnd - dataBaseStart + 1
        
        splicedAlignment = [score, queryStart, queryEnd, qsize, identity, chrom, dataBaseStart, databaseEnd, span]
        splicedAlignments.append(splicedAlignment)
    finalAlignments = []
    for finalAlignment in splicedAlignments:
        # Alignments with less than 20 score are removed
        if(finalAlignment[0] > 20):
            finalAlignments.append(finalAlignment)
    return finalAlignments   

def getAlignments(homologousRegions, query, data, k):
    alignments = []
    for region in homologousRegions:
        startPos = region[0]
        endPos = region[1]
        indexedList = createIndexedList(k, data[startPos:endPos], startPos, False, 0)
        hitList = getHits(query, data, indexedList, k, True)
        mergedHits = mergeHits(hitList)
        for hit in mergedHits:
            alignment = extendKmer(query, data, hit[1], [hit[2]], len(hit[0]), len(query))
            alignments.append(alignment)

    alignments = spliceAlignments(alignments, query, data)
    return alignments


def main():
    subsequenceFile, transcriptsFile = "subseq.fasta", "transcripts.fasta"
    data, transcripts = readSequencesFromFiles(subsequenceFile, transcriptsFile)
    k = 11
    maxHitsPerKmer = 10
    # Create an indexed list of k-mers from the data
    # Common k-mers that appear more than 10 times are removed
    indexedList = createIndexedList(k, data, 0, True, maxHitsPerKmer)

    # Run a query for all transcripts and print results
    for i in range(len(transcripts)):
        # Compare kmers from the query to the kmers in data and create a hit list
        query = transcripts[i]
        hitList = getHits(query, data, indexedList, k, False)
        # Split the hit list into buckets of 64K each, based on the database position
        # Hits within buckets have been sorted by the diagonal position
        bucketSize = 64000
        buckets = putHitsIntoBuckets(bucketSize, data, hitList)
        # Build proto-clumps from buckets
        # Hits within proto-clumps have been sorted by the database position
        protoClumps = buildProtoClumps(buckets)
        # Build clumps from protoclumps so that hits within clumps are
        # within window limit (W = 100) and each clump has at least minimumHits (2)
        # Clumps within merge distance (300) are merged
        W = 100
        minimumHits = 2
        clumps = buildClumps(protoClumps, W, minimumHits)
        # Build a homologous region from each clump and adds 500 additional bases
        # on each side to form the final homologous region
        homologousRegions = buildHomologousRegions(clumps, data)
        # From the homologous regions form alignments
        alignments = getAlignments(homologousRegions, query, data, k)
        if(alignments != []):
            print("Query   Score   Start  End     QSize Identity    Chrom      Start            End     Span")      
            for alignment in alignments:
                print(i+1, end = "       ")
                print(*alignment, sep = "     ")
            print()
        else:
            print("No alignment found for query", i+1, "\n")
            


if __name__ == "__main__":
    main()

