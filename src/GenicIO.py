# -------------------------------------------------------------------------------
# Name:        GenicIO.py
# Purpose: File processing for genomic files
#
# Author:      Tristan J. Hayeck
#
# Created:     2020
# Copyright:   (c) Tristan J. Hayeck 2020
# Licence:     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# -------------------------------------------------------------------------------
import sys, os
import numpy as np
import pandas as pd
import argparse
import pysam
import gzip
import datetime
ID_TYPE_COL = "idType"
ENSEBL = "ensembl"
CHROM_COL = "CHROM"
POS_COL = "POS"
MAF_COL = "MAF"
SUB_REGION_COL = "subRegion"
COUNTS_COL = "counts"
REF_COL = "REF"
ALT_COL = "ALT"
START_COL = "start"
END_COL = "end"


def checkFileExists(fileName):
    """
    Check if the file exists, otherwise raise an argparse ArgumentTypeError
    :param fileName:
    :return:
    """
    if os.path.isfile(fileName):
        return os.path.realpath(fileName)
    else:
        raise argparse.ArgumentTypeError("Error, %s does not exist" % fileName)

def extractMafTableFromVcf(populations, bedFileName, vcfFileName, outputFile, DEBUG):
    neutralIntervals = readBedFile(bedFileName)
    neutralMafTable = convertVcfToMafsTable(neutralIntervals, populations, vcfFileName, DEBUG)
    pd.DataFrame.to_csv(neutralMafTable, outputFile, sep="\t", index=False, header=True)


def readBedFile(bedFileName):
    bedCol=[CHROM_COL, START_COL, END_COL, SUB_REGION_COL]
    bedData = pd.read_csv(bedFileName, sep="\t", header=None,names=bedCol)
    return bedData


def getIntervalsFromBedFile(bedFileName,DEBUG=False):
    """
    Reads in the bed file and gets the intervals
    :param bedFileName: name of the bed file
    :return: intervals which should contain four column: the chromosome, start position, end position,
            and the sub-region or window
    """
    print("Running getIntervalsFromBedFile for in %s" % bedFileName)
    intervals = []
    with open(bedFileName,'r') as input:
        for line in input:
            intervals.append(line.strip().split('\t'))
    return intervals

def readFileAtLineNumber(lineNumber,fileName):
    with open(fileName) as inputFile:
        for lineIndex, line in enumerate(inputFile):
            if lineIndex == lineNumber:
                return line.replace('\n','').replace('\r','')
            elif lineIndex > lineNumber:
                print("Error for file %s line number %s is out of range %s" % (fileName,lineNumber,lineIndex))
                raise ValueError


def convertVcfToMafsTable(intervals, populations, tabixFileName, offset=0, ignoreIndels = True,
                      excludeLower = False,
                    CHROM_COL = "CHROM",POS_COL = "POS",ID_COL = "ID",REF_COL = "REF",ALT_COL = "ALT",
                      QUALITY_COL = "QUAL",FILTER_COL = "FILTER",INFO_COL = "INFO",ALL_MAF="AF",
                      DEBUG=False):
    """
    Get's the MAF values from a tabix sorted vcf file for the input populations over the input intervals.
    :param intervals: data frame with the chr:start-end to get MAF values
            start (int) start of the genomic region (0-based inclusive)
            stop (int)  end of the genomic region (0-based exclusive)
    :param populations: the populations with MAF values in the vcf to get
    :param tabixFileName: tabix sorted vcf file
    :param offset: Currently default to zero but, might be use case where this is helpful in the future
            maybe using this term to add to the start and end because sometimes this was an issue
    :param ignoreIndels: flag if we exclude indels from the analysis, default True
    :param excludeLower:
    :param DEBUG:
    :return:
    """
    NAN_VALUE = "."
    headerFields = [CHROM_COL, POS_COL, ID_COL, REF_COL, ALT_COL, QUALITY_COL, FILTER_COL, INFO_COL]
    columnNames = [CHROM_COL,FILTER_COL, POS_COL,REF_COL,ALT_COL]+ populations #SUB_REGION_COL
    totalSpan = (intervals[END_COL] - intervals[START_COL]).sum()
    # It is possiblye you could input a series of multi-allelic sites with indels included
    # and the total number of variants spanned is greater than totalSpan*3, but for almost all use
    # cases this is probably fine. For now leaving, in the future may modify
    crossPopulationMafTable = pd.DataFrame(columns=columnNames,index=range(totalSpan*3))
    crossPopulationMafTable[CHROM_COL] = crossPopulationMafTable[CHROM_COL].apply(str)
    crossPopulationMafTable[FILTER_COL] = crossPopulationMafTable[FILTER_COL].apply(str)
    crossPopulationMafTable[POS_COL] = -1
    crossPopulationMafTable[POS_COL] = crossPopulationMafTable[POS_COL].astype(int)
    crossPopulationMafTable[REF_COL] = crossPopulationMafTable[REF_COL].astype(str)
    crossPopulationMafTable[ALT_COL] = crossPopulationMafTable[ALT_COL].astype(str)
    #crossPopulationMafTable[SUB_REGION_COL] = crossPopulationMafTable[SUB_REGION_COL].apply(str)
    for popMaf in populations:
        crossPopulationMafTable[popMaf] = crossPopulationMafTable[popMaf].apply(np.double)
    regionTabix = pysam.Tabixfile(tabixFileName, 'r')

    lineCounter = 0
    excludeLowerCount = 0
    variantIndex = 0
    warningPos = 0
    for index, row in intervals.iterrows():
        #TODO this won't work on X,Y chrom
        #chromosome = str(row[CHROM_COL].astype(int))
        chromosome = str(np.int(row[CHROM_COL]))
        startPos= row[START_COL]
        endPos= row[END_COL]
        startPosOffset = np.int(startPos)+offset
        endPosOffset = np.int(endPos)+offset
        try:
            #--------------------#
            # start (int) start of the genomic region (0-based inclusive)
            # stop (int)  end of the genomic region (0-based exclusive)
            for tabixLine in regionTabix.fetch(chromosome, startPosOffset, endPosOffset):
                fields = tabixLine.rstrip().split('\t')
                infoColIndex = headerFields.index(INFO_COL)
                infoDct = dict([field.split("=", 1) for field in
                                fields[infoColIndex].split(";") if "=" in field])
                if np.int(fields[headerFields.index(POS_COL)]) <= startPosOffset or np.int(fields[headerFields.index(POS_COL)]) > endPosOffset :
                    if excludeLower and np.int(fields[headerFields.index(POS_COL)]) == startPosOffset:
                        excludeLowerCount+=1
                    else:
                        warningPos+=1
                        print("\n\nWarning, appears to be an issue with positions not between interval range!!")
                        print(row)
                        print("\n%s\n\n"%tabixLine)
                else:
                    curRef = fields[headerFields.index(REF_COL)]
                    #sometimes the REF maybe something like TCTG, so the alt will probably be an indel, so we may
                    #want to ignore it
                    if not ignoreIndels or len(curRef.replace(" ",""))<2:
                        #sometimes we have mulitallelic sites (multiple alternate alleles) and we should
                        # add a row for each
                        altAlleleCounter = 0
                        for curAlt in fields[headerFields.index(ALT_COL)].split(","):
                            if not ignoreIndels or len(curAlt.replace(" ", "")) < 2:
                                crossPopulationMafTable.at[variantIndex,CHROM_COL] = str(chromosome)
                                crossPopulationMafTable.at[variantIndex,FILTER_COL] = fields[headerFields.index(FILTER_COL)]
                                crossPopulationMafTable.at[variantIndex,POS_COL] = np.int(fields[headerFields.index(POS_COL)])
                                crossPopulationMafTable.at[variantIndex,REF_COL] = fields[headerFields.index(REF_COL)]
                                crossPopulationMafTable.at[variantIndex,ALT_COL] = curAlt
                                for popMaf in populations:
                                    missingMafForAPop = False
                                    if popMaf in infoDct.keys():
                                        curValue = infoDct[popMaf].split(",")[altAlleleCounter]
                                    else:
                                        missingMafForAPop=True
                                        curValue = 0
                                    if curValue == NAN_VALUE:
                                        missingMafForAPop = True
                                        crossPopulationMafTable.at[variantIndex,popMaf] = np.nan
                                    else:
                                        crossPopulationMafTable.at[variantIndex,popMaf] = np.double(curValue )
                                    if missingMafForAPop:
                                        print("\n--------------------")
                                        print("Warning, missing value for one or more populations Chrom = %s | Pos = %s |" %
                                        (fields[headerFields.index(CHROM_COL)],fields[headerFields.index(POS_COL)]))
                                        for pop in populations:
                                            if pop in infoDct.keys():
                                                print("%s = %s" %  (pop,infoDct[pop]))
                                            else:
                                                print("Missing value for Population = %s, setting to zero" % pop)
                                        print("\n--------------------\n")

                                altAlleleCounter += 1
                                variantIndex+=1
                lineCounter+=1
            # --------------------#
        except Exception as exception:
            print("Error at interval: ")
            print(row)
            print(tabixLine)
            print("for line %s: \n\n" % lineCounter)
            print(exception)
            print("----------------------")
    if excludeLower:
        print("For bed values that are exactly the start value, we exclude lower count = %s" % excludeLowerCount)
    crossPopulationMafTableCleaned = crossPopulationMafTable[crossPopulationMafTable[POS_COL] > -1]
    return crossPopulationMafTableCleaned




def getScoresOverIntervals(intervals, scoreFileName, header, outHeader=None, DEBUG=False, ):
    """
    Reads in the scores file and for each value in the interval returns a data frame with the position and
    corresponding scores
    :param intervals: data frame with the intervals, should have columns like a bed file
    :param header: column names for the scores file
    :param scoreFileName:
    :param outHeader: sometimes we don't want all the scores, so this allows us to specify the columns we want to reduce
    the size of the output. Default is None and then sets this to the full set of header values. Updated for more
    efficient process in getScorePercForAnnotations.py
    :return: scoresTable with the scores over the input intervals
    """
    if outHeader is None:
        outHeader = header
    regionTabix = pysam.Tabixfile(scoreFileName, 'r')
    intervalSpanCol = "intervalSpan"
    intervals[intervalSpanCol] = intervals[END_COL]-intervals[START_COL] #+1
    totPositionsSpanned = intervals[intervalSpanCol].sum()
    if intervals[intervals[intervalSpanCol]<1].shape[0]>0:
        print("Error, invalid intervals in getScoresOverIntervals: ")
        print(intervals[intervals[intervalSpanCol] < 1])
        raise ValueError
    scoresTable = pd.DataFrame(columns=outHeader,index=range(totPositionsSpanned))
    outputIndex = 0
    intervalIndicesFound = []
    for index, row in intervals.iterrows():
        chromosome = row[CHROM_COL]
        startPos = row[START_COL]
        endPos = row[END_COL]
        if DEBUG:
            print("%s from %s to %s  " % (index, startPos, endPos ))
        for tabixLine in regionTabix.fetch(chromosome,(startPos),(endPos)):
            intervalIndicesFound.append(index)
            if DEBUG:
                print(tabixLine)
            resultItems = tabixLine.rstrip().split('\t')
            colIndex = 0
            for col in header:
                if col in outHeader:
                    if col == CHROM_COL:
                        scoresTable.at[outputIndex, col] = np.int(str(resultItems[colIndex]).replace("chr",""))
                    elif col == POS_COL:
                        scoresTable.at[outputIndex,col] = np.int(resultItems[colIndex])
                        curPos= np.int(resultItems[colIndex])
                    else:
                        scoresTable.at[outputIndex, col] = np.double(resultItems[colIndex])
                colIndex+=1
            if int(curPos) < startPos or int(curPos)> endPos:
                print("Error, appears to be an issue in getLogRankScoresOverIntervals where tabix returned values not in range")
                print("%s from %s to %s vs current pos %s\n" % (chromosome, startPos, endPos,curPos))
                print(tabixLine)
                raise ValueError
            outputIndex+=1
    intervalIndicesNotFound =   set(intervals.index.tolist()) - set(intervalIndicesFound)
    if len(intervalIndicesNotFound) > 0:
        print("--------Warning Missing Intevals--------")
        print("indices: ")
        print(list(intervalIndicesNotFound))
        print("--------Warning Missing Intevals--------")
    return scoresTable