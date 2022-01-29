"""
Protein Sequencing Project
Name:
Roll Number:
"""

from asyncio import new_event_loop
from asyncore import read
from hashlib import new
from itertools import count
from json import JSONDecodeError
from re import U
from tkinter.filedialog import Open
import hw6_protein_tests as test

project = "Protein" # don't edit this

### WEEK 1 ###

'''
readFile(filename)
#1 [Check6-1]
Parameters: str
Returns: str
'''
def readFile(filename):
  new_file = open(filename,"r").read()
  file_staring =""
  for i in new_file.splitlines():
      file_staring+=i
  return file_staring


'''
dnaToRna(dna, startIndex)
#2 [Check6-1]
Parameters: str ; int
Returns: list of strs
'''
def dnaToRna(dna, startIndex):
    new_list=[]
    dna=dna.replace("T", "U")
    for i in range(startIndex,len(dna),3):
        new_list.append(dna[i:i+3])
    for j in new_list:
        if j=="UAA"or j=="UGA" or j=="UAG":
          new=new_list.index(j)
          return new_list[:new+1]   
    return new_list    


'''
makeCodonDictionary(filename)
#3 [Check6-1]
Parameters: str
Returns: dict mapping strs to strs
'''
def makeCodonDictionary(filename):
    import json
    new_dicts={}
    new_open = open(filename,"r")
    new_read = json.load(new_open)
    for i in new_read:
        for j in new_read[i]:
            j = j.replace("T","U")
            new_dicts[j]=i
    return new_dicts

'''
generateProtein(codons, codonD)
#4 [Check6-1]
Parameters: list of strs ; dict mapping strs to strs
Returns: list of strs
'''
def generateProtein(codons, codonD):
    new_list=[]
    for i in codons:
        for j in codonD:
            if i==j:
                new_list.append(codonD[j])
                if new_list[0] == 'Met':
                    new_list[0]='Start'
    
    return new_list


'''
synthesizeProteins(dnaFilename, codonFilename)
#5 [Check6-1]
Parameters: str ; str
Returns: 2D list of strs
'''
def synthesizeProteins(dnaFilename, codonFilename):
    new_list=[]
    unused=0
    first= readFile(dnaFilename)
    second = makeCodonDictionary(codonFilename)
    i=0

    while i < len(first):
        new= first[i:i+3]
        if new == 'ATG':
            newca=dnaToRna(first,i)
            newcall= generateProtein(newca,second)
            new_list.append(newcall)
            i+=3*len(newca)
        else:
            i=i+1
            unused=unused+1
    print(len(first),unused,len(new_list))            

    return new_list


def runWeek1():
    print("Human DNA")
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    print("Elephant DNA")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")


### WEEK 2 ###

'''
commonProteins(proteinList1, proteinList2)
#1 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs
Returns: 2D list of strs
'''
def commonProteins(proteinList1, proteinList2):
    new_list=[]
    x=proteinList1
    y=proteinList2
    for element1 in x:
        for element2 in y:
            if element1==element2 and element1 not in new_list:
                new_list.append(element2)
                print(new_list)
    return new_list
    


'''
combineProteins(proteinList)
#2 [Check6-2]
Parameters: 2D list of strs
Returns: list of strs
'''
def combineProteins(proteinList):
    new_list= [j for sub_list in proteinList for j in sub_list]

    return new_list


'''
aminoAcidDictionary(aaList)
#3 [Check6-2]
Parameters: list of strs
Returns: dict mapping strs to ints
'''
def aminoAcidDictionary(aaList):
    new_dict={}
    for i in aaList:
        if i not in new_dict:
            new_dict[i]=aaList.count(i)
    return new_dict


'''
findAminoAcidDifferences(proteinList1, proteinList2, cutoff)
#4 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs ; float
Returns: 2D list of values
'''
def findAminoAcidDifferences(proteinList1, proteinList2, cutoff):
    new_list1=combineProteins(proteinList1)
    new_list2=combineProteins(proteinList2)
    new_list3=aminoAcidDictionary(new_list1)
    new_list4=aminoAcidDictionary(new_list2)
    new_list=[]
    for i,j in new_list3.items():
        if i!="Start" and i!="Stop" and i not in new_list:
            new_list.append(i)
    for i,j in new_list4.items():
        if i!="Start" and i!="Stop" and i not in new_list:
            new_list.append(i)
    new_feq=[]
    for i in new_list:
        freq1=0
        freq2=0
        if i in new_list1:
            freq1=new_list3[i]/len(new_list1)
        if i in new_list2:
            freq2=new_list4[i]/len(new_list2)
        difference=freq1-freq2
        if difference>cutoff or difference<-cutoff:
            new_feq.append([i,freq1,freq2])
    return new_feq        
    


'''
displayTextResults(commonalities, differences)
#5 [Check6-2]
Parameters: 2D list of strs ; 2D list of values
Returns: None
'''
def displayTextResults(commonalities, differences):
    new_list=[]
    for i in commonalities:
        i.remove("Start")
        i.remove("Stop")
        if len(i)>1:
            add_symbol='-'.join(i)
            new_list.append([add_symbol])
        else:
            if i not in new_list:
                new_list.append(i)
    sort_list=sorted(new_list)
    for k in sort_list:
        for l in k:
            #print(l)
            for j in differences:
                print(j[0]+":"+str(round(j[1]*100,2))+"%"+ " in seq1 ,"+str(round(j[2]*100,2))+"%"+"in seq2")
                    
    return


def runWeek2():
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")

    commonalities = commonProteins(humanProteins, elephantProteins)
    differences = findAminoAcidDifferences(humanProteins, elephantProteins, 0.005)
    displayTextResults(commonalities, differences)


### WEEK 3 ###

'''
makeAminoAcidLabels(proteinList1, proteinList2)
#2 [Hw6]
Parameters: 2D list of strs ; 2D list of strs
Returns: list of strs
'''
def makeAminoAcidLabels(proteinList1, proteinList2):
    new_list=[]
    combine_list1=combineProteins(proteinList1)
    combine_list2=combineProteins(proteinList2)
    combine_list=combine_list1+combine_list2
    for i in combine_list:
        if i not in new_list:
            new_list.append(i)
    sorted_list=sorted(new_list)
    return sorted_list


'''
setupChartData(labels, proteinList)
#3 [Hw6]
Parameters: list of strs ; 2D list of strs
Returns: list of floats
'''
def setupChartData(labels, proteinList):
    new_list=[]
    combine_protine=combineProteins(proteinList)
    l2=aminoAcidDictionary(combine_protine)
    for j in labels:
        if j in l2:
            new_list.append(l2[j]/len(combine_protine))
        else:
            new_list.append(0)
    return new_list


'''
createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None)
#4 [Hw6] & #5 [Hw6]
Parameters: list of strs ; list of floats ; str ; list of floats ; str ; [optional] list of strs
Returns: None
'''
def createChart(xLabels, freqList1, label1, freqList2, label2, edgeList="black"):
    import matplotlib.pyplot as plt
    w = 0.35  # the width of the bars

    plt.bar(xLabels, freqList1, width=-w, align='edge', label=label1,edgecolor=edgeList)
    plt.bar(xLabels, freqList2, width= w, align='edge', label=label2,edgecolor=edgeList)

    plt.xticks(rotation="vertical")
    plt.legend()
    plt.title("create chart")

    plt.show()
    return


'''
makeEdgeList(labels, biggestDiffs)
#5 [Hw6]
Parameters: list of strs ; 2D list of values
Returns: list of strs
'''
def makeEdgeList(labels, biggestDiffs):
    return


'''
runFullProgram()
#6 [Hw6]
Parameters: no parameters
Returns: None
'''
def runFullProgram():
    return


### RUN CODE ###

# This code runs the test cases to check your work
if __name__ == "__main__":
    # print("\n" + "#"*15 + " WEEK 1 TESTS " +  "#" * 16 + "\n")
    # test.week1Tests()
    # print("\n" + "#"*15 + " WEEK 1 OUTPUT " + "#" * 15 + "\n")
    # runWeek1()
    # test.testReadFile()
    # test.testDnaToRna()
    # test.testMakeCodonDictionary()
    # test.testGenerateProtein()
    # test.testSynthesizeProteins()
    # test.testCommonProteins()
    # test.testCombineProteins()
    # test.testAminoAcidDictionary()
    # test.testFindAminoAcidDifferences()
    # runWeek2()
    # test.testMakeAminoAcidLabels()
    # test.testSetupChartData()
    test.testCreateChart()
    



    ## Uncomment these for Week 2 ##
    """
    print("\n" + "#"*15 + " WEEK 2 TESTS " +  "#" * 16 + "\n")
    test.week2Tests()
    print("\n" + "#"*15 + " WEEK 2 OUTPUT " + "#" * 15 + "\n")
    runWeek2()
    """

    ## Uncomment these for Week 3 ##
    """
    print("\n" + "#"*15 + " WEEK 3 TESTS " +  "#" * 16 + "\n")
    test.week3Tests()
    print("\n" + "#"*15 + " WEEK 3 OUTPUT " + "#" * 15 + "\n")
    runFullProgram()
    """
