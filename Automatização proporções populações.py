import csv
import re
from scipy.stats import chi2

nucleotides = ['A', 'T', 'G', 'C']

symbolsResults = []
fileName = 'rs6826785.csv'
translationTable = str.maketrans({'\'': '', "\"":"", '\n':'', ',':''})
expectedQuiSquared = float(chi2.ppf(1-0.05, df=2))

with open(fileName, 'r') as file:
    title = file.readline()

    line = file.readline().translate(translationTable)
    while line:
        firstSemicolon = line.find(';')
        popSymbol = line[:firstSemicolon].replace(' ','')
        line = line[firstSemicolon+1:]

        activated = False
        p = int; q = int #Frequências alélicas em números absolutos
        pRate = float; qRate = float #Frequências alélicas em frequência

        observedTT = int; observedTC = int; observedCC = int #Números absolutos dos genótipos observados
        observedTTRate = float; observedTCRate = float; observedCCRate = float #Frequências dos genótipos observados
        n = int #Número total da população

        # for cell in line.split('"'):
        #     if re.match(r'[ATGC]:', cell): #Frequências alélicas
                
        rateMatches = re.findall(r'0.[0-9]{3}', line)
        if rateMatches:
            pRate = float(rateMatches[0])
            qRate = float(rateMatches[1])
            observedTTRate = float(rateMatches[2])
            observedCCRate = float(rateMatches[3])
            try:
                observedTCRate = float(rateMatches[4])
            except:
                observedTCRate = None

            valueMatches = re.findall(r'\(([0-9]{1,})\)', line)
            p = int(valueMatches[0])
            q = int(valueMatches[1])
            observedTT = int(valueMatches[2])
            observedCC = int(valueMatches[3])
            try:
                observedTC = int(valueMatches[4])
            except:
                observedTC = None

            n = (p + q)/2
            activated = True
         
        # elif re.match('[ATGC]|[ATGC]', cell): #Frequências genotípicas observadas
        if(activated):
            expectedTT = pRate**2 #Frequência genotipica esperada TT         
            expectedCC = qRate**2 #Frequência genotípica esperada TC

            expectedTTPop = expectedTT * n #Frequência genotipica esperada valores absolutos TT
            expectedCCPop = expectedCC * n #Frequência genotípica esperada valores absolutos TC

            varianceA = (observedTT - expectedTTPop) ** 2
            varianceB = (observedCC - expectedCCPop) ** 2

            if observedTCRate:
                expectedTC = 2*pRate*qRate #Frequência genotípica esperada CC
                expectedTCPop = expectedTC * n #Frequência genotípica esperada valores absolutos CC
                varianceC = (observedTC - expectedTCPop) ** 2
            else:
                varianceC = 0

            quiSquared = (varianceA + varianceB + varianceC) / n            
            nullHypothesisFailed = quiSquared > expectedQuiSquared

            symbolsResults.append([popSymbol, quiSquared, nullHypothesisFailed])
        line = file.readline().translate(translationTable)


with open(fileName+'_Results.csv', 'w', newline='') as file:
    writer = csv.writer(file, delimiter=';')
    for pop in symbolsResults:
        writer.writerow(pop)









                





        





