# -------------------------------------------------------------------------------
# Name:        PCIT.py
# Purpose: Genomic regions subject to purifying selection are more likely to carry disease causing mutations.
# Tests for genetic intolerance look for depletion of variation relative to expectation within a species, allowing
# species specific features to be identified. The Population Conditional Intolerance Test (PCIT) is a weighted and
# stratified log rank test to test for differences between the SFS within a given query window and the SFS estimated
# from intergenic neutral sequence across multiple ancestral populations. This implementation takes in a neutral
# SFS spectrum and tests across test region(s). 
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
import os, sys
import argparse
import numpy as np
import pandas as pd
import scipy
import scipy.stats as stats
import time
import datetime
import GenicIO
LOG_RANK_COL = "logrank"
LOG_RANK_PVAL_COL = "logrank_pValue"
SUB_REGION_COL = "subRegion"
CHROM_COL = "CHROM"
POS_COL = "POS"
ID_COL = "ID"
REF_COL = "REF"
ALT_COL = "ALT"
QUALITY_COL ="QUAL"
FILTER_COL = "FILTER"
INFO_COL = "INFO"
headerFields =[CHROM_COL, POS_COL, ID_COL, REF_COL, ALT_COL, QUALITY_COL, FILTER_COL, INFO_COL]
START_COL = "start"
END_COL = "end"
POP_COL = "population"
WEIGHT_COL = "weight"
P_VALUE_SUFFIX = "_negLog10P"
PASS_CRITERIA = "PASS"
LCR_STATUS = "LCR"
SEGDUP_STATUS = "SEGDUP"
STRATIFIED_LR_COL = "StratLogRank"
WEIGHT_STRAT_LR_COL = "WeightStratLogRank"


def runTestOverRegion(windowSize,allPopulations,weightedPops,totPosSpannedNeutral, PASS, vcfFileName, neutralMafTable, row, DEBUG=False):
    """
    This is esentially a wrapper function for calcLogRankAsScan which will run the analysis over the specified start
    and end positions.
    :param windowSize:
    :param allPopulations:
    :param weightedPops:
    :param totPosSpannedNeutral:
    :param PASS:
    :param vcfFileName:
    :param neutralMafTable:
    :param row:
    :param DEBUG:
    :return:
    """
    rowOverWindow = row.copy()
    rowOverWindow[START_COL] = np.maximum(0, row[START_COL] - np.int(np.ceil(windowSize / 2.)))
    rowOverWindow[END_COL] = row[END_COL] + np.int(np.ceil(windowSize / 2.))
    funcMafTable = GenicIO.convertVcfToMafsTable(rowOverWindow, allPopulations, vcfFileName,
                                                 DEBUG=DEBUG)
    if PASS:
        funcMafTablePASS = funcMafTable[funcMafTable[FILTER_COL] == PASS_CRITERIA]
    else:
        funcMafTablePASS = funcMafTable
    if DEBUG:
        print("funcMafTable table dimensions:")
        print(funcMafTable.shape)
        print("funcMafTable PASS criteria table dimensions:")
        print(funcMafTablePASS.shape)
        print("neutralMafTable")
        print(neutralMafTable.shape)
        print(neutralMafTable.head(5))
        print("\n-----------------------\n")
        print(rowOverWindow)
        print("-----------------------\n")
        print(funcMafTable)
        print("\n-----------------------\n")
        print("Warning, in debug, rasing erorr now to stop test, no actual issue")
    logRankPopTestsTable = LogRankTest.calcLogRankAsScan(totPosSpannedNeutral, windowSize, allPopulations,weightedPops,
                                                         neutralMafTable, funcMafTablePASS,
                                                         startPos=row[START_COL].values[0]+1,
                                                         endPos=row[END_COL].values[0]+1,
                                                         DEBUG=DEBUG)
    #TODO this won't work with X or Y chrom
    logRankPopTestsTable[CHROM_COL]=np.int(row[CHROM_COL].values[0])
    return logRankPopTestsTable


def runTestWithinSubPops(args):
    GenicIO.checkFileExists(args.vcfFileName)
    GenicIO.checkFileExists(args.vcfFileName + ".tbi")
    print("\nBoth the %s vcf and corresponding tbi file exist." % args.vcfFileName)
    GenicIO.checkFileExists(args.bedFileName)
    print("\nRunning on: %s" % args.bedFileName)
    intervals = GenicIO.readBedFile(args.bedFileName)
    neutralMafTable = pd.read_csv(args.neutralFile, sep="\t")
    weightedPops = pd.read_csv(args.weightedPopFile, sep="\t")
    pops = pd.read_csv(args.popFile, sep="\t")
    allPopulations = pops[POP_COL].tolist() + weightedPops[POP_COL].tolist()
    if args.lineIndex>0:
        row = intervals.iloc[[args.lineIndex - 1]]
        logRankPopTestsTable = runTestOverRegion(args.windowSize,allPopulations,weightedPops, args.totPosSpannedNeutral,args.PASS,
                                                 args.vcfFileName,neutralMafTable, row, args.DEBUG)
        outpuTableFileName = args.outputDir + "%s.txt" % str(row[SUB_REGION_COL].values[0])
    else:
        outpuTableFileName = args.outputDir
        firstIter = True
        basesSpanned = (intervals[END_COL]-intervals[START_COL]).sum()
        posIndex=0
        for rowIndex in intervals.index.values:
            row=intervals.iloc[[rowIndex]]
            logRankPopTestsTableRegion = runTestOverRegion(args.windowSize,allPopulations,weightedPops,args.totPosSpannedNeutral,args.PASS,
                                                           args.vcfFileName, neutralMafTable, row, args.DEBUG)
            if firstIter:
                logRankPopTestsTable = pd.DataFrame(columns=logRankPopTestsTableRegion.columns.to_list(),
                                                    index=range(basesSpanned))
                firstIter=False
            for i in range(logRankPopTestsTableRegion.shape[0]):
                logRankPopTestsTable.iloc[posIndex]=logRankPopTestsTableRegion.iloc[i]
                posIndex+=1
    logRankPopTestsTable[POS_COL]=logRankPopTestsTable[POS_COL].astype(np.int)
    outputCols = [CHROM_COL,POS_COL]+allPopulations + [pop + P_VALUE_SUFFIX for pop in allPopulations] +\
                 [STRATIFIED_LR_COL,STRATIFIED_LR_COL+ P_VALUE_SUFFIX,
                  WEIGHT_STRAT_LR_COL,WEIGHT_STRAT_LR_COL + P_VALUE_SUFFIX]
    if args.DEBUG:
        print("logRankPopTestsTable")
        print(logRankPopTestsTable.head())
    pd.DataFrame.to_csv(logRankPopTestsTable[outputCols], outpuTableFileName, sep="\t", index=False, header=True)



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--bedFileName", type=str,
                        help="The bed file with the chromosome, start, end, and subregions columns to run the scan over")
    parser.add_argument("-n", "--neutralFile", type=str,
                        default="Neutral10M.tsv",
                        help="The neutral file name with tab seperated MAFs, typically created by extractSubPopData")
    parser.add_argument("-l", "--lineIndex", type=int,
                        default=-1,
                        help="Line index for which region to read in over bed file, the idea is this allows for easy "
                             "parallelization if you want to split up bed regions for each job. "
                             "Default -1 will read in all regions.")
    parser.add_argument("-v", "--vcfFileName", type=str,
                        default="data/gnomad/gnomad.genomes.r2.1.1.sites.vcf.bgz",
                        help="Input vcf file to be analyzed (must have corresponding tabix tbi file)")
    parser.add_argument("-w", "--windowSize", type=int,
                        default=1000,
                        help="Size of window to search around position")
    parser.add_argument("-p", "--popFile", type=str,
                        default="/scr1/users/hayeckt/PCIT/input/PopFiles/Pop.txt",
                        help="File with the populations, seperate from the weighted set to be used in the combined analysis.")
    parser.add_argument("--weightedPopFile", type=str,
                        default="/scr1/users/hayeckt/PCIT/input/PopFiles/GroupsPerWeight.txt",
                        help="File with the populations and their weights")
    parser.add_argument("-t", "--totPosSpannedNeutral", type=int,
                        help="Total number of positions spanned by neutral variants. This is not the number of neutral "
                             "variants but the regions they represent because many sites are monomorphic, so the "
                             "corresponding vcf the neutral sample was taken from will only have items for "
                             "each variant found.")
    parser.add_argument('--dont-PASS', dest='PASS', action='store_false',
                        help="Variants must have quality PASS to be included")
    parser.add_argument("-o", "--outputDir", type=str,
                        default="output/",
                        help="output file name or the directory for the MAF tables if running over a single region")
    parser.add_argument('--debug', dest='DEBUG', action='store_true')
    parser.set_defaults(DEBUG=False)
    parser.set_defaults(PASS=True)
    start = datetime.datetime.now()
    args = parser.parse_args()
    print("Starting processing %s" % start)
    print(args)
    runTestWithinSubPops(args)
    done = datetime.datetime.now()
    elapsed = done - start
    duration=':'.join(str(elapsed).split(':')[1:])
    seconds = elapsed.total_seconds()
    print("The duration was: %s" % duration)
    print("(ie %s sec.)" % seconds)
    print("Finished processing %s" % done)


if __name__ == '__main__':
    main()










