# -------------------------------------------------------------------------------
# Name:        LogRankTest.py
# Purpose: Called from PCIT.py,  Population Conditional Intolerance Test (PCIT) is a weighted and
# stratified log rank test to test for differences between the SFS within a given query window and the SFS estimated
# from intergenic neutral sequence across multiple ancestral populations. This implementation takes in a neutral
# SFS spectrum and tests across test region(s). It is done typically as a scan (calcLogRankAsScan) where the data
# structures are updated over regions to increase efficiency by avoiding repopulating the data structure
# for running the test. The other main functions are collapseSurvivalData which formats the data for quick analysis
# and calcLogRankFromEventGroups
#
# Author:      Tristan J. Hayeck
#
# Created:     2020
# Copyright:   (c) Tristan J. Hayeck 2020
# Licence:     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# # IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# # FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# # AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# # LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# # OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# # SOFTWARE.
# -------------------------------------------------------------------------------
import os, sys
import numpy as np
import pandas as pd
import scipy
import scipy.stats as stats
import time
import datetime
import pysam

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
#populations MAF columns
AFR_MAF_COL = "AF_afr"
AMR_MAF_COL = "AF_amr"
ASJ_MAF_COL = "AF_asj"
EAS_MAF_COL = "AF_eas"
FIN_MAF_COL = "AF_fin"
NFE_MAF_COL = "AF_nfe"
OTH_MAF_COL = "AF_oth"
ALL_MAF = "AF"
allPopulations = [ALL_MAF, AFR_MAF_COL, AMR_MAF_COL, ASJ_MAF_COL, EAS_MAF_COL, FIN_MAF_COL, NFE_MAF_COL, OTH_MAF_COL]
GROUP_COL = "group"
FAILURE_COL = "Failure"
POP_1_COL = "POP1"
POP_2_COL = "POP2"
P_VALUE_SUFFIX = "_negLog10P"
#nulls in the gnomad file come up as "."
NULL_VALUE = "."
POP_COL = "population"
WEIGHT_COL = "weight"
#These are very closely related names 'weight' and 'weights', may want to change in the future
WEIGHTS_COL = "weights"
NUM_EVENT_AT_TIME_COL = "NumberEventsAtTime"
NUM_AT_RISK_COL = "NumberAtRisk"
START_AT_RISK_COL = "StartAtRisk"
GROUP_SUFFIX = "_G"
NUM_LR_SUFFIX = "_numeratorLR"
DENOM_LR_SUFFIX = "_denominatorLR"
STRATIFIED_LR_COL = "StratLogRank"
WEIGHT_STRAT_LR_COL = "WeightStratLogRank"



def calcWeightedLogRank(timesAndGroups, groupCol, gamma1=0, gamma2=0, DEBUG=False):
    """
    Note: not allowing for censoring, could probably add in pretty easily
    Calculates the log rank test (currently not weighted, will update this) where we assume failure time
    corresponds to the number of alleles, ie a singleton denotes failure at time 1, doubleton is at time 2.
    Another way of thinking about this is comparing time to fixation of the alleles.
    #
    d_0j: failures in group 0 at j-th failure time
    r_0j: number at risk in group 0 at j-th failure time
    d_1j: failures in group 0 at j-th failure time
    r_1j: number at risk in group 1 at j-th failure time
    d_j: failures in both groups at j-th failure time
    r_j: number at risk in both groups at j-th failure time
    o_j: observed failures for group 0, d_0j
    e_j: expected failures for group 0, E[d_0]
    v_j: variances of failures for group 0, E[d_0]
    #
    s_j: survival at j-th failure time, prod  s_j*(1-d_j/)
    #
    allowing for The G_w1_w2 family (Fleming and Harrington; 1991) of weights
    w_j = (S(j)^gamma1 (1-S(j))^gamma2, where w1 = 1 and w2 = 0 is the Wilcoxon test
    #
    main reference:
    Fleming T & Harrington D (1991). Counting Processes and Survival Analysis, Wiley, New York.
    quick references:
    https://cran.r-project.org/web/packages/survminer/vignettes/Specifiying_weights_in_log-rank_comparisons.html
    http://vassarstats.net/survival.html
    #
    Note: there is no censoring in this model, probably not to hard to adapt if wanted in the future.
    :param gamma1:
    :param gamma2:
    :param DEBUG:
    :return:
    """
    timesAndGroupsSorted = timesAndGroups.copy()
    timesAndGroupsSorted = timesAndGroupsSorted.sort_values(by=[FAILURE_COL ])
    numerator = 0
    denominator = 0
    s_j = 1.0
    r_0j = 1.0 * (timesAndGroupsSorted[timesAndGroupsSorted[groupCol] == 0]).shape[0]
    r_1j = 1.0 * (timesAndGroupsSorted[timesAndGroupsSorted[groupCol] == 1]).shape[0]
    numAtRisk = r_0j + r_1j
    tot_obs = 0
    tot_exp = 0
    if DEBUG:
        print("d_0j, r_0j, d_1j, r_1j, d_j, r_j,numAtRisk,s_j, w_j, o_j, e_j, v_j,numerator,denominator")
    uniqueTimes = timesAndGroupsSorted[FAILURE_COL ].unique()
    d_0j = 0
    d_1j = 0
    for curTime in uniqueTimes:
        prev_d_0j = d_0j
        prev_d_1j = d_1j
        # get all events at the time
        eventsAtTime = timesAndGroupsSorted[timesAndGroupsSorted[FAILURE_COL ] == curTime]
        # get the deaths in each group
        d_0j = eventsAtTime[eventsAtTime[groupCol] == 0].shape[0]
        d_1j = eventsAtTime[eventsAtTime[groupCol] == 1].shape[0]
        r_0j -= prev_d_0j
        r_1j -= prev_d_1j
        d_j = d_0j + d_1j
        r_j = r_0j + r_1j
        o_j = d_0j
        numAtRisk = r_j
        1.0 * numAtRisk - d_j
        if d_j != 0 and r_j != 1 and r_j != 0:
            w_j = 1.0 * np.power(s_j, gamma1) * np.power((1.0 - s_j), gamma2)
            e_j = r_0j * d_j / r_j
            v_j = (r_0j * r_1j * d_j * (r_j - d_j)) / (r_j ** 2 * (r_j - 1))
            numerator += (w_j * (o_j - e_j))
            denominator += (np.power(w_j, 2) * v_j)
            tot_obs += w_j * o_j
            tot_exp += w_j * e_j
            if DEBUG:
                print( "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s" % \
                      (d_0j, r_0j, d_1j, r_1j, d_j, r_j, numAtRisk, s_j, w_j, o_j, e_j, v_j, numerator, denominator))
        s_j = s_j * (1.0 - (1.0 * d_j / r_j))
    if DEBUG or denominator == 0:
        print("numerator = {%s}^2 and denominator = %s" % (numerator, denominator))
        print("tot obs = %s and tot exp = %s" % (tot_obs, tot_exp))
        if denominator == 0:
            print("Invalide denominator for log rank test")
            print("start and tail of table ")
            print(timesAndGroups.head(100))
            print(timesAndGroups.tail(100))
            print("weighed with gamma1 = %s, gamma2 = %s" % (gamma1, gamma2))
            raise ValueError
    logRank = numerator / np.sqrt(denominator)
    return logRank


def calcSparseWeightedLogRank(timesAndGroups, numZerosInGroup0,numZerosInGroup1, groupCol, gamma1=0, gamma2=0, DEBUG=False):
    """
    Seperately calculating the initial zero values (ie monomorphic sites) because including them in the
    times and groups involves sorting which if we can avoid it should speed up calculations a lot.

    Note: not allowing for censoring, could probably add in pretty easily
    Calculates the log rank test (currently not weighted, will update this) where we assume failure time
    corresponds to the number of alleles, ie a singleton denotes failure at time 1, doubleton is at time 2.
    Another way of thinking about this is comparing time to fixation of the alleles.
    #
    d_0j: failures in group 0 at j-th failure time
    r_0j: number at risk in group 0 at j-th failure time
    d_1j: failures in group 0 at j-th failure time
    r_1j: number at risk in group 1 at j-th failure time
    d_j: failures in both groups at j-th failure time
    r_j: number at risk in both groups at j-th failure time
    o_j: observed failures for group 0, d_0j
    e_j: expected failures for group 0, E[d_0]
    v_j: variances of failures for group 0, E[d_0]
    #
    s_j: survival at j-th failure time, prod  s_j*(1-d_j/)
    #
    allowing for The G_w1_w2 family (Fleming and Harrington; 1991) of weights
    w_j = (S(j)^gamma1 (1-S(j))^gamma2, where w1 = 1 and w2 = 0 is the Wilcoxon test
    #
    main reference:
    Fleming T & Harrington D (1991). Counting Processes and Survival Analysis, Wiley, New York.
    quick references:
    https://cran.r-project.org/web/packages/survminer/vignettes/Specifiying_weights_in_log-rank_comparisons.html
    http://vassarstats.net/survival.html
    #
    Note: there is no censoring in this model, probably not to hard to adapt if wanted in the future.
    :param gamma1:
    :param gamma2:
    :param DEBUG:
    :return:
    """
    timesAndGroupsSorted = timesAndGroups.copy()
    timesAndGroupsSorted = timesAndGroupsSorted.sort_values(by=[FAILURE_COL ])
    numerator = 0
    denominator = 0
    s_j = 1.0
    if DEBUG:
        print("numZerosInGroup0 = %s" %numZerosInGroup0)
        print("numZerosInGroup1 = %s" %numZerosInGroup1)
    r_0j = numZerosInGroup0 + 1.0 * (timesAndGroupsSorted[timesAndGroupsSorted[groupCol] == 0]).shape[0]
    r_1j = numZerosInGroup1 + 1.0 * (timesAndGroupsSorted[timesAndGroupsSorted[groupCol] == 1]).shape[0]
    numAtRisk = r_0j + r_1j
    r_j = r_0j + r_1j
    tot_obs = 0
    tot_exp = 0
    if DEBUG:
        print("d_0j, r_0j, d_1j, r_1j, d_j, r_j,numAtRisk,s_j, w_j, o_j, e_j, v_j,numerator,denominator")
    uniqueTimes = timesAndGroupsSorted[FAILURE_COL ].unique()
    d_0j = numZerosInGroup0
    d_1j = numZerosInGroup1
    d_j = d_0j + d_1j
    #effectively just adding first iteration here accounting for zero
    w_j = 1.0 * np.power(s_j, gamma1) * np.power((1.0 - s_j), gamma2)
    o_j = d_0j
    e_j = r_0j * d_j / r_j
    v_j = (r_0j * r_1j * d_j * (r_j - d_j)) / (r_j ** 2 * (r_j - 1))
    numerator += (w_j * (o_j - e_j))
    denominator += (np.power(w_j, 2) * v_j)
    tot_obs += w_j * o_j
    tot_exp += w_j * e_j
    for curTime in uniqueTimes:
        if curTime == 0:
            print("Error, inside of calcSparseWeightedLogRank should not be iterating over time %s " % curTime)
            print(timesAndGroupsSorted[timesAndGroupsSorted[FAILURE_COL] == curTime].head(10))
            print(timesAndGroupsSorted[timesAndGroupsSorted[FAILURE_COL] == curTime].tail(10))
            raise ValueError
        prev_d_0j = d_0j
        prev_d_1j = d_1j
        # get all events at the time
        eventsAtTime = timesAndGroupsSorted[timesAndGroupsSorted[FAILURE_COL ] == curTime]
        # get the deaths in each group
        d_0j = eventsAtTime[eventsAtTime[groupCol] == 0].shape[0]
        d_1j = eventsAtTime[eventsAtTime[groupCol] == 1].shape[0]
        r_0j -= prev_d_0j
        r_1j -= prev_d_1j
        d_j = d_0j + d_1j
        r_j = r_0j + r_1j
        o_j = d_0j
        numAtRisk = r_j
        1.0 * numAtRisk - d_j
        if d_j != 0 and r_j != 1 and r_j != 0:
            w_j = 1.0 * np.power(s_j, gamma1) * np.power((1.0 - s_j), gamma2)
            e_j = r_0j * d_j / r_j
            v_j = (r_0j * r_1j * d_j * (r_j - d_j)) / (r_j ** 2 * (r_j - 1))
            numerator += (w_j * (o_j - e_j))
            denominator += (np.power(w_j, 2) * v_j)
            tot_obs += w_j * o_j
            tot_exp += w_j * e_j
            if DEBUG:
                print("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s" % \
                      (d_0j, r_0j, d_1j, r_1j, d_j, r_j, numAtRisk, s_j, w_j, o_j, e_j, v_j, numerator, denominator))
        s_j = s_j * (1.0 - (1.0 * d_j / r_j))
    if DEBUG or denominator == 0:
        print("numerator = {%s}^2 and denominator = %s" % (numerator, denominator))
        print("tot obs = %s and tot exp = %s" % (tot_obs, tot_exp))
        if denominator == 0:
            print("Invalide denominator for log rank test")
            print("start and tail of table ")
            print(timesAndGroups.head(100))
            print(timesAndGroups.tail(100))
            print("weighed with gamma1 = %s, gamma2 = %s" % (gamma1, gamma2))
            raise ValueError
    logRank = numerator / np.sqrt(denominator)
    return logRank


def convertSFStoCounts(SFS):
    """
    Converts the SFS to the original counts. For example say the SFS = [2 4,3,2,2,1]
    then the counts are [0,0,1,1,1,1,2,2,2,3,3,4,4,5]
    :param SFS: Site Frequency Spectrum vector of length equal to the number of alleles
                with non-polymorphic, singletons, doubletons, ... ... number Individuals
                (or 2x number Individuals)
    """
    numSamples = np.sum(SFS)
    counts = np.zeros(numSamples, dtype=np.int64)
    curSample = 0
    for alleleIndex in range(SFS.shape[0]):
        numAlleles = SFS[alleleIndex]
        for nestedIndex in range(numAlleles):
            counts[curSample] = alleleIndex
            curSample += 1
    return counts


def collapseSurvivalData(events, groupSuffix, weights=1):
    """
    Reads in events and converts it to collapsed sorted survival data to later be processed to calculate the logrank
    statistic:
    For example
        reads in failure events: array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 5])
        group suffix "_G"
        weights = 1
                        weights_G0  NumberEventsAtTime_G0  StartAtRisk_G0  NumberAtRisk_G0
            Failure
            0                10                     10              25             25.0
            1                 9                      9               0             15.0
            2                 4                      4               0              6.0
            3                 1                      1               0              2.0
            5                 1                      1               0              1.0
    :param events:
    :param groupSuffix:
    :param weights:
    :return:
    """
    eventsDf = pd.DataFrame(events,columns=[FAILURE_COL])
    eventsDf[WEIGHTS_COL+groupSuffix]=weights
    eventsDf[NUM_EVENT_AT_TIME_COL+groupSuffix] = 1
    collapsedEvents = eventsDf.groupby(FAILURE_COL).sum()
    collapsedEvents.sort_index(inplace=True)
    collapsedEvents[START_AT_RISK_COL+groupSuffix] = 0
    collapsedEvents[START_AT_RISK_COL+groupSuffix].loc[0] = eventsDf[WEIGHTS_COL+groupSuffix].sum()
    collapsedEvents[NUM_AT_RISK_COL+groupSuffix] = collapsedEvents[START_AT_RISK_COL+groupSuffix].cumsum() - collapsedEvents[NUM_EVENT_AT_TIME_COL+groupSuffix].cumsum().shift(1).fillna(0)
    return collapsedEvents



def calcLogRank(eventsGroup0, eventsGroup1, DEBUG=False):
    """
    Calculates more effecient implementation of log rank test (one sided) and does not include censoring
    calls on collapseSurvivalData to get effective collapsing of the data
    :param eventsGroup0: events for first group (or ie allele counts or MAF)
    :param eventsGroup1: events for second group (or ie allele counts or MAF)
    :param DEBUG:
    :return:
    """
    eventsGroup0Collapsed = collapseSurvivalData(eventsGroup0, GROUP_SUFFIX + str(0))
    eventsGroup1Collapsed = collapseSurvivalData(eventsGroup1, GROUP_SUFFIX + str(1))
    collapsedData = eventsGroup0Collapsed.join(eventsGroup1Collapsed, how='outer')
    collapsedData = collapsedData.fillna(0)
    if DEBUG:
        print("Group 1: ")
        print(eventsGroup0Collapsed.head())
        print("Group 2: ")
        print(eventsGroup1Collapsed.head())
        print("Combined 1: ")
        print(collapsedData.head())
    eventGroups = collapsedData.filter(like=NUM_EVENT_AT_TIME_COL)
    #total in each group
    numInGroups_j = eventGroups.sum(0).values #we don't have censoring
    #total in each group at each time, recalculating just so we can combine with indices below
    numInGroupsAtTime_ij = (eventGroups.sum(0).values - eventGroups.cumsum(0).shift(1).fillna(0))
    #failures at time i for both groups, summing across 2 columns
    eventsAtTime_i = eventGroups.sum(1)
    #number at risk at each time
    numAtRiskAtTime_i = eventGroups.values.sum() - eventGroups.sum(1).cumsum().shift(1).fillna(0)
    #expected values
    expectedValue = numInGroupsAtTime_ij.mul(eventsAtTime_i / numAtRiskAtTime_i, axis='index').sum(0)
    # vector of observed minus expected, numerator for the statistic
    Z_j = numInGroups_j - expectedValue
    # compute covariance matrix
    factor = (((numAtRiskAtTime_i - eventsAtTime_i) / (numAtRiskAtTime_i - 1)).replace([np.inf, np.nan], 1)) * eventsAtTime_i / numAtRiskAtTime_i ** 2
    numInGroupsAtTime_ij['_'] = numAtRiskAtTime_i.values
    V_ = numInGroupsAtTime_ij.mul(np.sqrt(factor), axis='index').fillna(0)
    V = -np.dot(V_.T, V_)
    #since we have two groups
    ix = np.array([0,1]) #np.arange(n_groups)
    V[ix, ix] = V[ix, ix] - V[-1, ix]
    V = V[:-1, :-1]
    #we only want the one-sided p-value for two group test
    stat = Z_j[0]/np.sqrt(np.abs(V[0,0]))
    return stat



def testMethods():
    testFile = "WithinSubPop/MAF_GT_05_chr1_nonTeloCento_intergenic.masked.txt"
    testTable = pd.read_csv(testFile, sep='\t')
    testTableHead = testTable.head(100000)
    maf2 = testTableHead["AF_NFE"][0:80]
    countsAndGroup = pd.DataFrame(
        {GROUP_COL: np.zeros(testTableHead[ALL_MAF].shape[0]), FAILURE_COL: testTableHead[ALL_MAF]})
    countsAndGroup = countsAndGroup.append(
        pd.DataFrame({GROUP_COL: np.ones(maf2.shape[0]), FAILURE_COL: maf2}))
    test = calcWeightedLogRank(countsAndGroup, GROUP_COL, gamma1=0, gamma2=0, DEBUG=False)
    print("Original implementation: ")
    print(test)
    print("New implmentation: ")
    newTest = calcLogRank(testTableHead[ALL_MAF].to_numpy(), maf2.to_numpy())
    print(newTest)




def testCalcLogRank():
    """
    Output with debugging looks like:
    d_0j, r_0j, d_1j, r_1j, d_j, r_j,numAtRisk,s_j, w_j, o_j, e_j, v_j,numerator,denominator
    10,25.0,8,21.0,18,46.0,46.0,1.0,1.0,10,9.78260869565,2.77882797732,0.217391304348,2.77882797732
    9,15.0,6,13.0,15,28.0,28.0,0.608695652174,1.0,9,8.03571428571,1.79634353741,1.18167701863,4.57517151473
    4,6.0,3,7.0,7,13.0,13.0,0.282608695652,1.0,4,3.23076923077,0.869822485207,1.95090778786,5.44499399994
    1,2.0,1,4.0,2,6.0,6.0,0.130434782609,1.0,1,0.666666666667,0.355555555556,2.2842411212,5.80054955549
    0,1.0,1,3.0,1,4.0,4.0,0.0869565217391,1.0,0,0.25,0.1875,2.0342411212,5.98804955549
    1,1.0,0,2.0,1,3.0,3.0,0.0652173913043,1.0,1,0.333333333333,0.222222222222,2.70090778786,6.21027177772
    0,0.0,1,2.0,1,2.0,2.0,0.0434782608696,1.0,0,0.0,0.0,2.70090778786,6.21027177772
    numerator = {2.70090778786}^2 and denominator = 6.21027177772
    tot obs = 25.0 and tot exp = 22.2990922121
    1.08381324447
    1.17465114888
    Out[11]:
    1.0838132444657027
    :return:
    """
    funcSFS = np.array([10, 9, 4, 1, 0, 1])
    nonFuncSFS = np.array([8, 6, 3, 1, 1, 0, 0, 1, 1])
    numIndivs = 25
    funcCounts = convertSFStoCounts(funcSFS)
    nonFuncCounts = convertSFStoCounts(nonFuncSFS)
    countsAndGroup = pd.DataFrame({GROUP_COL: np.zeros(funcCounts.shape[0]), FAILURE_COL: funcCounts})
    countsAndGroup = countsAndGroup.append(
        pd.DataFrame({GROUP_COL: np.ones(nonFuncCounts.shape[0]), FAILURE_COL: nonFuncCounts}))
    test = calcWeightedLogRank(countsAndGroup, GROUP_COL, gamma1=0, gamma2=0, DEBUG=False)
    print("test result %s using calcWeightedLogRank" % test)
    funcCountsNonZero = funcCounts[np.nonzero(funcCounts)]
    nonFuncCountsNonZero = nonFuncCounts[np.nonzero(nonFuncCounts)]
    countsAndGroup = pd.DataFrame({GROUP_COL: np.zeros(funcCountsNonZero.shape[0]), FAILURE_COL: funcCountsNonZero})
    countsAndGroup = countsAndGroup.append(
        pd.DataFrame({GROUP_COL: np.ones(nonFuncCountsNonZero.shape[0]), FAILURE_COL: nonFuncCountsNonZero}))
    test = calcSparseWeightedLogRank(countsAndGroup,  np.sum(funcCounts==0),np.sum(nonFuncCounts==0),GROUP_COL,
                                     gamma1=0, gamma2=0, DEBUG=False)
    print("test result %s using calcSparseWeightedLogRank" % test)
    return test



def calculateNeutralVsFuncWithinPopLogRank(totPosSpannedNeutral, totPosSpannedFunc, populations,
                                           neutralMafTable, funcMafTable, addMonomorphicSites = True, DEBUG = True):
    logRankPopTests = {}
    if addMonomorphicSites:
        for pop in populations:
            neutralMafs = neutralMafTable[pop]
            funcMafs = funcMafTable[pop]
            neutralMafsNonZero = (neutralMafs.iloc[neutralMafs.nonzero()[0]]).dropna()
            funcMafsNonZero = (funcMafs.iloc[funcMafs.nonzero()[0]]).dropna()
            countsAndGroup = pd.DataFrame({GROUP_COL: np.zeros(neutralMafsNonZero.shape[0]),
                                           FAILURE_COL: neutralMafsNonZero})
            countsAndGroup = countsAndGroup.append(
                                pd.DataFrame({GROUP_COL: np.ones(funcMafsNonZero.shape[0]),
                                                FAILURE_COL: funcMafsNonZero}))
            if DEBUG:
                if np.sum(neutralMafsNonZero==0) > 0 or \
                    np.sum(funcMafsNonZero == 0) > 0 or \
                    np.sum(neutralMafsNonZero.isnull()) > 0 or \
                    np.sum(funcMafsNonZero.isnull())> 0:
                    print("neutralMafsNonZero number zeros: ")
                    print(np.sum(neutralMafsNonZero==0))
                    print("funcMafsNonZero number zeros: ")
                    print(np.sum(funcMafsNonZero == 0))
                    print("np.sum(neutralMafsNonZero.isnull())")
                    print(np.sum(neutralMafsNonZero.isnull()))
                    print("np.sum(funcMafsNonZero.isnull())")
                    print(np.sum(funcMafsNonZero.isnull()))
                print("neutralMafs.shape")
                print(neutralMafs.shape)
                print("funcMafs.shape")
                print(funcMafs.shape)
                print("neutralMafsNonZero.shape")
                print(neutralMafsNonZero.shape)
                print("funcMafsNonZero.shape")
                print(funcMafsNonZero.shape)
            numMonomorphicNeutral = totPosSpannedNeutral * 3 - neutralMafsNonZero.shape[0]
            numMonomorphicFunc = totPosSpannedFunc * 3 - funcMafsNonZero.shape[0]
            neutralMafsCleanedWithZeros = np.concatenate([neutralMafsNonZero, np.zeros(np.int(numMonomorphicNeutral))])
            funcMafsCleanedWithZeros = np.concatenate([funcMafsNonZero, np.zeros(np.int(numMonomorphicFunc))])
            logRankStat = calcLogRank(neutralMafsCleanedWithZeros,funcMafsCleanedWithZeros,DEBUG)
            logRankLogPValue = calcNegLog10PvalueOneSided(logRankStat)
            logRankPopTests[pop] = logRankStat
            logRankPopTests[pop+P_VALUE_SUFFIX] = logRankLogPValue
            if DEBUG:
                print("countsAndGroup.shape")
                print(countsAndGroup.shape)
                print("countsAndGroup.head()")
                print(countsAndGroup.head())
                print("%s = %s with -log_10(pValue) %s" % (pop,logRankStat,logRankLogPValue))
    else:
        print("For now only allowing the use of monomorphic sites because this appears to be the proper way to test")
        raise ValueError
    return logRankPopTests



def testStratLogRank():
    """
    This is a unit test comparint results of a stratified analysis using the R library
    survival. Created an arbitrary example with fake data where we have drug vs placebo
    (neutral vs functional) broken up into different regions.

    The R script, testStratLogRank.r runs the same test
    library(survival)
    eventTime=as.numeric(c(2,3,4,5,10,11,3,3,3,5,10,11)/100)
    status=rep(1,12)
    group=c(rep('placebo',3),rep('drug',3),rep('placebo',3),rep('drug',3))
    loc=c(rep('NYC',6),rep('Philly',6))
    eventGroups = data.frame(cbind(eventTime,status,group,loc))
    survfit <- survdiff(Surv(eventTime, status) ~ group+strata(loc))#, data = eventGroups
    survfit$chisq
    [1] 9.953437
    Could not find a way to directly test against weighted stratified.
    """
    NYC_COL = "NYC"
    PHILLY_COL = "Philly"
    placebo = NUM_EVENT_AT_TIME_COL + GROUP_SUFFIX + str(0)
    drug  = NUM_EVENT_AT_TIME_COL + GROUP_SUFFIX + str(1)
    nyc = pd.DataFrame({FAILURE_COL:[0.02, 0.03, 0.04, 0.05, 0.10, 0.11],
                        placebo:[1,1,1,0,0,0],drug:[0,0,0,1,1,1]})
    nyc = nyc.set_index(FAILURE_COL)
    philly=pd.DataFrame({FAILURE_COL:[0.03,0.05, 0.10, 0.11],
                        placebo:[3,0,0,0],drug:[0,1,1,1]})
    philly = philly.set_index(FAILURE_COL)
    nycNumerLR, nycDenomLR, nycLogRankStat = calcLogRankFromEventGroups(nyc)
    phillyNumerLR, phillyDenomLR, phillyLogRankStat = calcLogRankFromEventGroups(philly)
    weightedPops = pd.DataFrame({POP_COL:[NYC_COL,PHILLY_COL],WEIGHT_COL:[1,1]})
    strataValues = pd.DataFrame( {NYC_COL + NUM_LR_SUFFIX:[nycNumerLR,nycNumerLR],
                   NYC_COL + DENOM_LR_SUFFIX: [nycDenomLR,nycDenomLR],
                    PHILLY_COL + NUM_LR_SUFFIX:[phillyNumerLR,phillyNumerLR],
                    PHILLY_COL + DENOM_LR_SUFFIX: [phillyDenomLR,phillyDenomLR]},index=range(2))
    stratlogRankTest = strataValues.apply(lambda row:
                          calcStratStat(row, weightedPops, useWeights=False,
                                        DEBUG=False), axis=1)
    expectedResult=9.953437
    print("Stratified test %s expected %s diff of %s"%
          (np.round(stratlogRankTest[0]**2,6),expectedResult,
           (np.round(stratlogRankTest[0]**2,6)-expectedResult)) )



def calcLogRankFromEventGroups(eventGroups):
    """
    Calculates the log rank assuming we have two groups with input sorted like this:
                        NumberEventsAtTime_G0  NumberEventsAtTime_G1
        Failure
        0.000000              2949663.0                 2939.0
        0.000114                  640.0                    0.0
        0.000114                 1486.0                    1.0
        0.000115                 2162.0                    2.0
        0.000115                 2303.0                    0.0
    Note: this is a modified version of calcLogRank for easier optimization by just directly updating and
    passing in the eventGroups data structure instead of recalculating it every time.
    :param eventGroups: index with the failure time, then the event times for each of the two groups
    :return: one-sided log rank state
    """
    #total in each group
    numInGroups_j = eventGroups.sum(0).values #we don't have censoring
    #total in each group at each time, recalculating just so we can combine with indices below
    numInGroupsAtTime_ij = (eventGroups.sum(0).values - eventGroups.cumsum(0).shift(1).fillna(0))
    #failures at time i for both groups, summing across 2 columns
    eventsAtTime_i = eventGroups.sum(1)
    #number at risk at each time
    numAtRiskAtTime_i = eventGroups.values.sum() - eventGroups.sum(1).cumsum().shift(1).fillna(0)
    #expected values
    expectedValue = numInGroupsAtTime_ij.mul(eventsAtTime_i / numAtRiskAtTime_i, axis='index').sum(0)
    # vector of observed minus expected, numerator for the statistic
    Z_j = numInGroups_j - expectedValue
    # compute covariance matrix
    factor = (((numAtRiskAtTime_i - eventsAtTime_i) / (numAtRiskAtTime_i - 1)).replace([np.inf, np.nan], 1)) * eventsAtTime_i / numAtRiskAtTime_i ** 2
    numInGroupsAtTime_ij['_'] = numAtRiskAtTime_i.values
    V_ = numInGroupsAtTime_ij.mul(np.sqrt(factor), axis='index').fillna(0)
    V = -np.dot(V_.T, V_)
    #since we have two groups
    ix = np.array([0,1]) #np.arange(n_groups)
    V[ix, ix] = V[ix, ix] - V[-1, ix]
    V = V[:-1, :-1]
    #we only want the one-sided p-value for two group test
    stat = Z_j[0]/np.sqrt(np.abs(V[0,0]))
    return Z_j[0], V[0,0], stat



def calcNegLog10PvalueOneSided(stat):
    """
    Calculating a one-sided p-value for a normal distribution assuming that negative values are
    the extreme we care about. (ie -10 stat gives ~ 23.1 whereas 10~0).
    :param stat:
    :return:
    """
    return -1 * np.log10(stats.norm.cdf( stat))


def calcStratStat(row, weightedPops, useWeights, DEBUG=False):
    """
    Called on to calculate the stratified and weighted+stratified log rank statistics from the seperately
    calculated numerator and denominators within each group
    :param row:
    :param weightedPops:
    :param useWeights:
    :param DEBUG:
    :return:
    """
    numerator = 0
    denominator = 0
    if DEBUG:
        print("calcStratStat")
        print(row)
    for index, popRow in weightedPops.iterrows():
        weight = popRow[WEIGHT_COL] if useWeights else 1.0
        numerator += weight*row[popRow[POP_COL]+NUM_LR_SUFFIX]
        denominator += weight*weight*row[popRow[POP_COL]+DENOM_LR_SUFFIX]
        if DEBUG:
            print('weight = %s' % weight)
            print("running numerator / denominator = %s / %s ~ %s" %
                  (numerator , denominator, np.round(numerator/ denominator,3) ))
    if denominator<0:
        print("Error, invalid denominator for position %s denom = %s for stratified (weightin=%s) "
              % (row[POS_COL], denominator, useWeights))
        raise ValueError
    return numerator / np.sqrt(np.abs(denominator))


def calcLogRankAsScan(totPosSpannedNeutral, windowSize, populations, weightedPops,
                     neutralMafTable, funcMafTable, startPos, endPos,
                    DEBUG = False):
    """
    Main funciton that calculates the log rank statistics scanning through the funcMafTable across startPos to endPos
    using the input windowSize to define the test region.
    :param totPosSpannedNeutral: this corresonds to the number of positions as opposed to observed variants,
        for example if we look over 10 positions we may only see 3 SNPs which in the input vcf and
        corresponding MAF table this would only actually be 3 rows of items when really it was over 10 positoins.
    :param windowSize:
    :param populations:
    :param weightedPops:
    :param neutralMafTable:
    :param funcMafTable:
    :param startPos:
    :param endPos:
    :param DEBUG:
    :return:
    """
    outputCols = [POS_COL]+ populations + [pop + P_VALUE_SUFFIX for pop in populations] + \
                 [pop + NUM_LR_SUFFIX for pop in populations] + [pop + DENOM_LR_SUFFIX for pop in populations]
    logRankPopTests = pd.DataFrame(np.nan, index=range(startPos, endPos), columns=outputCols)
    for pop in populations:
        neutralMafs = neutralMafTable[pop]
        neutralMafsNonZero=neutralMafs[neutralMafs>0]
        numMonomorphicNeutral = totPosSpannedNeutral * 3 - neutralMafsNonZero.shape[0]
        neutralMafsCleanedWithZeros = np.concatenate([neutralMafsNonZero, np.zeros(np.int(numMonomorphicNeutral))])
        curStartPos = startPos - int(windowSize/2)
        curEndPos =   startPos + int(windowSize/2)
        funcMafTableOverWindow = funcMafTable[ (funcMafTable[POS_COL]>=curStartPos ) &
                                                           (funcMafTable[POS_COL]<=curEndPos)]
        funcMafs = funcMafTableOverWindow[pop]
        funcMafsNonZero = funcMafs[funcMafs>0]
        numMonomorphicFunc = windowSize * 3 - funcMafsNonZero.shape[0]
        if DEBUG:
            print("windowSize = %s; Functional numMonomorphicFunc = %s; and polymorphic = %s" %\
              (windowSize, numMonomorphicFunc,funcMafsNonZero.shape[0]))
        if numMonomorphicFunc < 0:
            print("\n\nWarning!!!")
            print("windowSize = %s; Functional numMonomorphicFunc = %s; and polymorphic = %s" % \
                  (windowSize, numMonomorphicFunc, funcMafsNonZero.shape[0]))
            print("This is maybe because we're looking over a small window, setting negative value to 0. \n\n")
            numMonomorphicFunc = np.max([numMonomorphicFunc,0])
            funcMafsCleanedWithZeros = np.concatenate([funcMafsNonZero, np.zeros(np.int(numMonomorphicFunc))])
        else:
            funcMafsCleanedWithZeros = np.concatenate([funcMafsNonZero, np.zeros(np.int(numMonomorphicFunc))])
        if DEBUG:
            if np.sum(neutralMafsNonZero==0) > 0 or \
                np.sum(funcMafsNonZero == 0) > 0 or \
                np.sum(neutralMafsNonZero.isnull()) > 0 or \
                np.sum(funcMafsNonZero.isnull())> 0:
                print("neutralMafsNonZero number zeros: ")
                print(np.sum(neutralMafsNonZero==0))
                print("funcMafsNonZero number zeros: ")
                print(np.sum(funcMafsNonZero == 0))
                print("np.sum(neutralMafsNonZero.isnull())")
                print(np.sum(neutralMafsNonZero.isnull()))
                print("np.sum(funcMafsNonZero.isnull())")
                print(np.sum(funcMafsNonZero.isnull()))
            print("neutralMafs.shape")
            print(neutralMafs.shape)
            print("funcMafs.shape")
            print(funcMafs.shape)
            print("neutralMafsNonZero.shape")
            print(neutralMafsNonZero.shape)
            print("funcMafsNonZero.shape")
            print(funcMafsNonZero.shape)
        #probably don't have to copy but just making sure we aren't modifying the starting data
        eventsGroup0 = neutralMafsCleanedWithZeros.copy()
        eventsGroup1 = funcMafsCleanedWithZeros.copy()
        eventsGroup0Collapsed = collapseSurvivalData(eventsGroup0, GROUP_SUFFIX + str(0))
        eventsGroup1Collapsed = collapseSurvivalData(eventsGroup1, GROUP_SUFFIX + str(1))
        collapsedData = eventsGroup0Collapsed.join(eventsGroup1Collapsed, how='outer')
        collapsedData = collapsedData.fillna(0)
        if DEBUG:
            print("Group 1: ")
            print(eventsGroup0Collapsed.head())
            print("Group 2: ")
            print(eventsGroup1Collapsed.head())
            print("Combined 1: ")
            print(collapsedData.head())
        eventGroups = collapsedData.filter(like=NUM_EVENT_AT_TIME_COL)
        numerLR, denomLR, logRankStat = calcLogRankFromEventGroups(eventGroups)
        logRankLogPValue = calcNegLog10PvalueOneSided(logRankStat)
        logRankPopTests.at[startPos, POS_COL] = startPos
        logRankPopTests.at[startPos,pop] = logRankStat
        logRankPopTests.at[startPos,pop+P_VALUE_SUFFIX] = logRankLogPValue
        logRankPopTests.at[startPos, pop+NUM_LR_SUFFIX] = numerLR
        logRankPopTests.at[startPos, pop + DENOM_LR_SUFFIX] = denomLR
        # now scan without having to totally repopulate data structures
        eventsGroup0Col = NUM_EVENT_AT_TIME_COL + GROUP_SUFFIX + '0'
        eventsGroup1Col = NUM_EVENT_AT_TIME_COL + GROUP_SUFFIX + '1'
        if DEBUG:
            print("before scan starts, eventGroups.head(50)")
            print(eventGroups.head(50))
            print("eventGroups[eventGroups[eventsGroup1Col]>0].shape")
            print(eventGroups[eventGroups[eventsGroup1Col] > 0].shape)
            print(eventGroups[eventGroups[eventsGroup1Col] > 0])
        prevLogRankStat = logRankStat
        for curPos in range(startPos+1, endPos):
            curStartPos = curPos - int(windowSize / 2)
            curEndPos = curPos + int(windowSize / 2)
            prevStartPos = curStartPos - 1
            prevStartPosMafs = funcMafTable[funcMafTable[POS_COL]==prevStartPos][pop].values
            curEndPosMafs = funcMafTable[funcMafTable[POS_COL]==curEndPos][pop].values
            needToShiftWindow = True
            if DEBUG:
                print("eventGroups[eventGroups[eventsGroup1Col]>0].shape")
                print(eventGroups[eventGroups[eventsGroup1Col] > 0].shape)
                print(eventGroups[eventGroups[eventsGroup1Col] > 0])
            #Put this in because sometimes in the highly polymorphic regions there would be issues and
            # we have more polymorphisms than 3*window size, trying to just re-read in the data every time this
            # happens, will slow things down should avoid breaking at multi-allelic sites
            highlyPolyMorphicRegion = (windowSize * 3 - funcMafsNonZero.shape[0]) < 0
            if curEndPosMafs.shape[0] > 1 or prevStartPosMafs.shape[0] > 1 or highlyPolyMorphicRegion:
                if DEBUG:
                    print("Multi-allilic at previous start of window %s or current window end %s " % (prevStartPos,curEndPos))
                    print("Previous start:")
                    print(funcMafTable[funcMafTable[POS_COL]==prevStartPos][pop])
                    print("current end:")
                    print(funcMafTable[funcMafTable[POS_COL]==curEndPos][pop])
                    print("Recalculating collapsed set of events ")
                funcMafTableOverWindow = funcMafTable[ (funcMafTable[POS_COL]>=curStartPos ) &
                                                                   (funcMafTable[POS_COL]<=curEndPos)]
                funcMafs = funcMafTableOverWindow[pop]
                funcMafsNonZero = funcMafs[funcMafs>0]
                numMonomorphicFunc = windowSize * 3 - funcMafsNonZero.shape[0]
                if numMonomorphicFunc < 0:
                    print("\n\nWarning!!!")
                    print("windowSize = %s; Functional numMonomorphicFunc = %s; and polymorphic = %s" % \
                          (windowSize, numMonomorphicFunc, funcMafsNonZero.shape[0]))
                    print("This is maybe because we're looking over a small window, setting negative value to 0. \n\n")
                    numMonomorphicFunc = np.max([numMonomorphicFunc, 0])
                    funcMafsCleanedWithZeros = np.concatenate([funcMafsNonZero, np.zeros(np.int(numMonomorphicFunc))])
                else:
                    funcMafsCleanedWithZeros = np.concatenate([funcMafsNonZero, np.zeros(np.int(numMonomorphicFunc))])
                eventsGroup0 = neutralMafsCleanedWithZeros.copy()
                eventsGroup1 = funcMafsCleanedWithZeros.copy()
                eventsGroup0Collapsed = collapseSurvivalData(eventsGroup0, GROUP_SUFFIX + str(0))
                eventsGroup1Collapsed = collapseSurvivalData(eventsGroup1, GROUP_SUFFIX + str(1))
                collapsedData = eventsGroup0Collapsed.join(eventsGroup1Collapsed, how='outer')
                collapsedData = collapsedData.fillna(0)
                eventGroups = collapsedData.filter(like=NUM_EVENT_AT_TIME_COL)
                numerLR, denomLR,newLogRankStat = calcLogRankFromEventGroups(eventGroups)
                needToShiftWindow = False
            elif curEndPosMafs.shape[0] == 1:
                #check if same value in as value out by shifting window
                if (curEndPosMafs[0] == 0 and prevStartPosMafs.shape[0]==0):
                    newLogRankStat = prevLogRankStat
                    needToShiftWindow = False
                elif prevStartPosMafs.shape[0]==0:
                    needToShiftWindow = True
                elif curEndPosMafs[0] == prevStartPosMafs[0]:
                    newLogRankStat = prevLogRankStat
                    needToShiftWindow = False
                else:
                    needToShiftWindow = True
            else:
                #so it's a monomorphic site, let's see if the start of the old window was also
                if prevStartPosMafs.shape[0]==0 or prevStartPosMafs[0] == 0:
                    newLogRankStat = prevLogRankStat
                    needToShiftWindow = False
                else:
                    if curEndPosMafs.shape[0] != 0:
                        print("Error, issue trying to set what should be zero, %s and %s" % (prevStartPos, curEndPos))
                        print(curEndPosMafs)
                        raise ValueError
            if needToShiftWindow:
                #add new MAF remove old MAF
                if curEndPosMafs.shape[0] != 0:
                    curEndPosMaf = curEndPosMafs[0]
                    if np.isnan(curEndPosMaf):
                        curEndPosMaf = 0
                else:
                    curEndPosMaf = 0
                if curEndPosMaf not in eventGroups.index:
                    if DEBUG:
                        print("adding new row for %s with MAF %s " % (curEndPos,curEndPosMaf))
                    eventGroups.loc[curEndPosMaf] = np.nan
                    eventGroups.at[curEndPosMaf, eventsGroup0Col] =0
                    eventGroups.at[curEndPosMaf, eventsGroup1Col] =1
                    if DEBUG:
                        print(eventGroups.loc[curEndPosMaf])
                else:
                    oldValue = eventGroups.loc[curEndPosMaf][eventsGroup1Col] # eventGroups.at[curEndPosMaf, eventsGroup1Col]
                    eventGroups.at[curEndPosMaf,eventsGroup1Col] = oldValue+1
                if prevStartPosMafs.shape[0]==0:
                    prevStartPosMaf = 0
                else:
                    prevStartPosMaf = prevStartPosMafs[0]
                    if np.isnan(prevStartPosMaf):
                        prevStartPosMaf = 0
                if prevStartPosMaf not in eventGroups.index:
                    print("Error, issue trying to remove previous start of window %s, MAF doesn't appear to exist in set %s" % \
                          (prevStartPos,prevStartPosMaf))
                    raise ValueError
                if DEBUG:
                    print("prevStartPosMaf = %s for eventsGroup1Col = %s" % (prevStartPosMaf,eventsGroup1Col))
                    print(eventGroups.loc[prevStartPosMaf])
                    print(eventGroups.loc[prevStartPosMaf][eventsGroup1Col])
                oldStartValue = eventGroups.loc[prevStartPosMaf][eventsGroup1Col] # eventGroups.at[prevStartPosMaf, eventsGroup1Col]
                if oldStartValue < 1:
                    print("Error, at pos = %s issue trying to remove previous start of window %s, value being set to negative " % \
                          (curPos,prevStartPos))
                    print("prevStartPosMaf = %s for eventsGroup1Col = %s" % (prevStartPosMaf, eventsGroup1Col))
                    print(eventGroups.loc[prevStartPosMaf])
                    print(eventGroups.loc[prevStartPosMaf][eventsGroup1Col])
                    print("eventGroups")
                    print(eventGroups)
                    raise ValueError
                else:
                    eventGroups.at[prevStartPosMaf, eventsGroup1Col] = oldStartValue - 1
                if eventGroups[(eventGroups[eventsGroup0Col] != 0) | (eventGroups[eventsGroup1Col] != 0)].shape[0]>0:
                    eventGroups = eventGroups[(eventGroups[eventsGroup0Col] != 0) | (eventGroups[eventsGroup1Col] != 0)]
                    eventGroups.sort_index(inplace=True)
                numerLR, denomLR, newLogRankStat = calcLogRankFromEventGroups(eventGroups)
            logRankLogPValue = calcNegLog10PvalueOneSided(logRankStat) # -1 * np.log10(stats.chi2.sf(np.square(newLogRankStat), 1))
            logRankPopTests.at[curPos, POS_COL] = curPos
            logRankPopTests.at[curPos, pop] = newLogRankStat
            logRankPopTests.at[curPos, pop + P_VALUE_SUFFIX] = logRankLogPValue
            logRankPopTests.at[curPos, pop + NUM_LR_SUFFIX] = numerLR
            logRankPopTests.at[curPos, pop + DENOM_LR_SUFFIX] = denomLR
            prevLogRankStat = newLogRankStat
    #now calculate both the stratified and weighted stratified log rank scores
    logRankPopTests[STRATIFIED_LR_COL] = logRankPopTests.apply(lambda row:
                                                               calcStratStat(row, weightedPops, useWeights=False,
                                                                             DEBUG=DEBUG),axis=1)
    logRankPopTests[STRATIFIED_LR_COL+ P_VALUE_SUFFIX] = logRankPopTests[STRATIFIED_LR_COL].apply(
                                                                    lambda x:calcNegLog10PvalueOneSided(x))
    logRankPopTests[WEIGHT_STRAT_LR_COL] =logRankPopTests.apply(lambda row:
                                                               calcStratStat(row, weightedPops, useWeights=True,
                                                                             DEBUG=DEBUG),axis=1)
    logRankPopTests[WEIGHT_STRAT_LR_COL + P_VALUE_SUFFIX] = logRankPopTests[WEIGHT_STRAT_LR_COL].apply(
                                                                    lambda x: calcNegLog10PvalueOneSided(x))
    return logRankPopTests

