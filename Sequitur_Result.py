
import re
import gc
import sequitur

def deleteWhiteSpacefromFasta(file1):
	file1 = re.sub(r'^.{5}', '', file1.read().lower()) #ignore the title
	file2 = re.sub(r'[^natcg]', '', file1) #delete unwanted chars
	return file2

def readFastaFile(fileName):
    fastaFile = open(fileName, 'r') # chrx.fa file extension here  #~250 MB
    fastaString = deleteWhiteSpacefromFasta(fastaFile)
    fastaFile.close()
    memoryCollected = gc.collect()
    print (memoryCollected, "collected… in readFastFile() ", fileName)
    return fastaString

sequitur.reset_avg()
sequitur.reset_graph()
sequitur.reset_numRules()

chrFileName = ['', "chr1.fa", "chr2.fa", "chr3.fa", "chr4.fa", "chr5.fa", "chr6.fa", "chr7.fa","chr8.fa", "chr9.fa", "chr10.fa", "chr11.fa", "chr12.fa", "chr13.fa", "chr14.fa", "chr15.fa", "chr16.fa","chr17.fa", "chr18.fa", "chr19.fa", "chr20.fa", "chr21.fa", "chr22.fa", "chrX.fa", "chrY.fa" ]
fastaString = ['']
for i in range(1,len(chrFileName)):
    fastaString.append(readFastaFile('/Users/dykim/Desktop/Ewha Bioinformatics/Sequitur/sequitur-master/hg19/' + chrFileName[i]))
f25 = open('ch_Sequitur_Output_25.txt','wt')
f30 = open('ch_Sequitur_Output_30.txt','wt')
for i in range(1,16):
    trainBEDFileName = "/Users/dykim/Desktop/Ewha Bioinformatics/정부과제/snpX200bp15states/All_200bp_DominantSub_TwoStates_"+str(i)+"_v5.bed"
    bedFile= open(trainBEDFileName,'r')
    bedFile.seek(0)

    count = 0

    f25.write('<chrom'+str(i)+'state>\n')
    f30.write('<chrom'+str(i)+'state>\n')

    while True:
            line = bedFile.readline()
            if not line:
                    break
            tokens = line.strip().split('\t')

            if int(tokens[3]) == 1:  
                    chrom = int(tokens[0])
                    start = int(tokens[1])
                    end = int(tokens[2])
                    dnaseq = fastaString[chrom][start:end]
                    rules = sequitur.run_sequitur(dnaseq)
                    graph = sequitur.num_graph()
                    count = count +1
                    ruleCount = sequitur.get_ruleCount()
                    if ruleCount>30:
                        f30.write('chorm'+str(chrom)+'\n')
                        f30.write(str(start)+'\t'+str(end)+'\n')
                        #f30.write('dnaseq\n')
                        #f30.write(rules+'\n\n')
                        f30.write('\n')
                     
                    elif ruleCount>25:
                        f25.write('chorm'+str(chrom)+'\n')
                        f25.write(str(start)+'\t'+str(end)+'\n')
                        #f25.write('dnaseq\n')
                        #f25.write(rules+'\n\n')   
                        f25.write('\n')
                    

    sequitur.reset_avg()
    sequitur.reset_graph()
    sequitur.reset_numRules()
    
f25.close()
f30.close()
