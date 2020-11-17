# -------------------------------------------------------------------------------
# Name:        extractMafTableFromVcf.py
# Purpose: This is used to get the minor allele frequency (MAF) table, typically used by PCIT with the neutral
# site frequency specturm. Wrapper file for GenicIO.extractMafTableFromVcf. 
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
POP_COL = "population"

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--bedFileName", type=str,
                        help="The bed file with the chromosome, start, end, and subregions columns defining the "
                             "positions for scanning over. ")
    parser.add_argument("-o", "--outputFile", type=str,
                        default="output/neutral.tsv",
                        help="The ouput file for the MAF table being created.")
    parser.add_argument("-p", "--popFile", type=str,
                        default="input/PopFiles/allPops.txt",
                        help="File with the populations, should include both weighted and unweighted sets")
    parser.add_argument("-v", "--vcfFileName", type=str,
                        default="/data/gnomad/gnomad.genomes.r2.1.1.sites.vcf.bgz",
                        help="Input vcf file to be analyzed, must have corresponding tabix tbi file")
    parser.add_argument('--debug', dest='DEBUG', action='store_true')
    parser.set_defaults(DEBUG=False)
    start = datetime.datetime.now()
    args = parser.parse_args()
    print("Starting processing %s" % start)
    print(args)
    pops = pd.read_csv(args.popFile, sep="\t")
    GenicIO.extractMafTableFromVcf(pops[POP_COL].tolist(), args.bedFileName, args.vcfFileName, args.outputFile, args.DEBUG)
    done = datetime.datetime.now()
    elapsed = done - start
    duration=':'.join(str(elapsed).split(':')[1:])
    seconds = elapsed.total_seconds()
    print("The duration was: %s" % duration)
    print("(ie %s sec.)" % seconds)
    print("Finished processing %s" % done)


if __name__ == '__main__':
    main()


