# Population Conditional Intolerance Test (PCIT)
## [Ancestry adjustment improves genome-wide estimates of regional intolerance](https://www.biorxiv.org/content/10.1101/2020.03.05.979203v2)

Genomic regions subject to purifying selection are more likely to carry disease causing mutations. Cross species conservation is often used to identify such regions but has limited resolution to detect selection on short evolutionary timescales such as that occurring in only one species. In contrast, intolerance looks for depletion of variation relative to expectation within a species, allowing species specific features to be identified. When estimating the intolerance of noncoding sequence methods strongly leverage variant frequency distributions. As the expected distributions depend on demography, if not properly controlled for, ancestral population source may obfuscate signals of selection. We demonstrate that properly incorporating demography in intolerance estimation greatly improved variant classification (13% increase in AUC relative to comparison constraint test, CDTS; and 9% relative to conservation). We provide a genome-wide intolerance map that is conditional on demographic history that is likely to be particularly valuable for variant prioritization.


Further information to come.

-----------------------------------------------------
### Download PCIT files for test fit on gnomAD v2.1.
[Download Directory](https://upenn.box.com/v/genomeWidePCIT)  \
[Chrom 1 PCIT](https://upenn.box.com/v/pcitChr1Tests) [(tbi)](https://upenn.box.com/v/pcitChr1Tabix),
[Chrom 2 PCIT](https://upenn.box.com/v/pcitChr2Tests) [(tbi)](https://upenn.box.com/v/pcitChr2Tabix),
[Chrom 3 PCIT](https://upenn.box.com/v/pcitChr3Tests) [(tbi)](https://upenn.box.com/v/pcitChr3Tabix),  \
[Chrom 4 PCIT](https://upenn.box.com/v/pcitChr4Tests) [(tbi)](https://upenn.box.com/v/pcitChr4Tabix),
[Chrom 5 PCIT](https://upenn.box.com/v/pcitChr5Tests) [(tbi)](https://upenn.box.com/v/pcitChr5Tabix),
[Chrom 6 PCIT](https://upenn.box.com/v/pcitChr6Tests) [(tbi)](https://upenn.box.com/v/pcitChr6Tabix),  \
[Chrom 7 PCIT](https://upenn.box.com/v/pcitChr7Tests) [(tbi)](https://upenn.box.com/v/pcitChr7Tabix),
[Chrom 8 PCIT](https://upenn.box.com/v/pcitChr8Tests) [(tbi)](https://upenn.box.com/v/pcitChr8Tabix),
[Chrom 9 PCIT](https://upenn.box.com/v/pcitChr9Tests) [(tbi)](https://upenn.box.com/v/pcitChr9Tabix),  \
[Chrom 10 PCIT](https://upenn.box.com/v/pcitChr10Tests) [(tbi)](https://upenn.box.com/v/pcitChr10Tabix),
[Chrom 11 PCIT](https://upenn.box.com/v/pcitChr11Tests) [(tbi)](https://upenn.box.com/v/pcitChr11Tabix),
[Chrom 12 PCIT](https://upenn.box.com/v/pcitChr12Tests) [(tbi)](https://upenn.box.com/v/pcitChr12Tabix),  \
[Chrom 13 PCIT](https://upenn.box.com/v/pcitChr13Tests) [(tbi)](https://upenn.box.com/v/pcitChr13Tabix),
[Chrom 14 PCIT](https://upenn.box.com/v/pcitChr14Tests) [(tbi)](https://upenn.box.com/v/pcitChr14Tabix),
[Chrom 15 PCIT](https://upenn.box.com/v/pcitChr15Tests) [(tbi)](https://upenn.box.com/v/pcitChr15Tabix),  \
[Chrom 16 PCIT](https://upenn.box.com/v/pcitChr16Tests) [(tbi)](https://upenn.box.com/v/pcitChr16Tabix),
[Chrom 17 PCIT](https://upenn.box.com/v/pcitChr17Tests) [(tbi)](https://upenn.box.com/v/pcitChr17Tabix),
[Chrom 18 PCIT](https://upenn.box.com/v/pcitChr18Tests) [(tbi)](https://upenn.box.com/v/pcitChr18Tabix),  \
[Chrom 19 PCIT](https://upenn.box.com/v/pcitChr19Tests) [(tbi)](https://upenn.box.com/v/pcitChr19Tabix),
[Chrom 20 PCIT](https://upenn.box.com/v/pcitChr20Tests) [(tbi)](https://upenn.box.com/v/pcitChr20Tabix),
[Chrom 21 PCIT](https://upenn.box.com/v/pcitChr21Tests) [(tbi)](https://upenn.box.com/v/pcitChr21Tabix),  \
[Chrom 22 PCIT](https://upenn.box.com/v/pcitChr22Tests) [(tbi)](https://upenn.box.com/v/pcitChr22Tabix)


## Running PCIT

PCIT is a weighted and stratified log rank test to test for differences between the SFS within a given query window and the SFS estimated from intergenic neutral sequence across multiple ancestral populations. PCIT uses python 3.6 and standard libraries: `pysam, argparse, numpy, pandas, scipy, gzip, datetime, time.` You can find sample input files in the input directory of this repository. You can download the gnomad data [here]( https://gnomad.broadinstitute.org/downloads/).

To run PCIT, it is largely a two-step process: 1.) extracting the selecting the neutral site frequency spectrum (SFS) 2.) running the test across the region(s) of interest.

### Getting the neutral SFS

In order to get the neutral SFS first specify a bed file [like this](https://github.com/tris-10/PopCondIntolTest/blob/master/input/Neutral_r2.1.1/byChrom/chr1_CoverageFilteredNeutralPassGW.bed) over neutral regions for the data set of interest, here we use gnomAD. Additionally, specify the populations you'll be testing in the [popFile](https://github.com/tris-10/PopCondIntolTest/blob/master/input/popFiles/allPops.txt). Running the script extractMafTableFromVcf.py to get the neutral regions across each chromosome, using something like:

```
neutralDir=input/Neutral_r2.1.1/
for chrom in {1..22}; do
python src/extractMafTableFromVcf.py \
  -v data/gnomad/gnomad.genomes.r2.1.1.sites.vcf.bgz \
  -p input/popFiles/allPops.txt
  -b ${neutralDir}byChrom/chr${chrom}_CoverageFilteredNeutralPassGW.bed \
  -o ${neutralDir}byChrom/chr${chrom}_MafCoverageFilteredNeutralPassGW.txt;
done
```

Where the usage is:

```
python extractMafTableFromVcf.py --help
usage: extractMafTableFromVcf.py [-h] [-b BEDFILENAME] [-o OUTPUTFILE]
                                 [-p POPFILE] [-v VCFFILENAME] [--debug]

optional arguments:
  -h, --help            show this help message and exit
  -b BEDFILENAME, --bedFileName BEDFILENAME
                        The bed file with the chromosome, start, end, and
                        subregions columns defining the positions for scanning
                        over.
  -p POPFILE, --popFile POPFILE
                        File with the populations, should include both
                        weighted and unweighted sets
  -v VCFFILENAME, --vcfFileName VCFFILENAME
                        Input vcf file to be analyzed, must have corresponding
                        tabix tbi file
 -o OUTPUTFILE, --outputFile OUTPUTFILE
                        The ouput file for the MAF table being created.
  --debug
```

The output will be a tab delimited file with the chromosome, filter, position, reference allele, alternate allele, and the minor allele frequency (MAF) for each population specified in the popFile.  Then you can combine the files into one file across all chromosomes, an example is [here](https://github.com/tris-10/PopCondIntolTest/blob/master/input/Neutral_r2.1.1/NeutralPopMafFiltered.txt).
```
head -n1 ${neutralDir}byChrom/chr1_MafCoverageFilteredNeutralPassGW.txt >\
      ${neutralDir}NeutralPopMafFiltered.txt;
for chrom in {1..22}; do
  echo $chrom;
  grep -v CHROM ${neutralDir}byChrom/chr${chrom}_MafCoverageFilteredNeutralPassGW.txt  \
    ${neutralDir}NeutralPopMafFiltered.txt;
done
```

### Running PCIT test scan

To actually run the tests using the neutral file you just generated, you will also need to specify two population files, one that's the [weightedPopFile](https://github.com/tris-10/PopCondIntolTest/blob/master/input/popFiles/GroupsIPW.txt) used in the weighted stratified test statistic and then another one if you want to run the test separately, say on the combined populations ([like this](https://github.com/tris-10/PopCondIntolTest/blob/master/input/popFiles/Pop.txt)). The tests will be run over sliding windows centered at the position, the window is +/-(windSize/2) up and down stream of the test position. The [bedfile](https://github.com/tris-10/PopCondIntolTest/blob/master/input/genomeWideScanBeds/chrom1Step1M.bed) specifies the regions to test. For this example, every position in the first line  region (-l --lineIndex, default is -1 for all reigons).
```
chrom=1
outputDir=output/chr${chrom}/
neutralDir=input/Neutral_r2.1.1/
popFileDir=input/popFiles/
python src/PCIT.py \
  -v data/gnomad/gnomad.genomes.r2.1.1.sites.vcf.bgz \
  -b input/genomeWideScanBeds/chrom${chrom}Step1M.bed \
  -n ${neutralDir}NeutralPopMafFiltered.txt \
  --weightedPopFile ${popFileDir}/GroupsIPW.txt \
  -p ${popFileDir}/Pop.txt \
  -o ${outputDir} \
  --windowSize 100 \
  --lineIndex 1;
```

With the usage is:

```
usage: PCIT.py [-h] [-b BEDFILENAME] [-n NEUTRALFILE] [-l LINEINDEX]
               [-v VCFFILENAME] [-w WINDOWSIZE] [-p POPFILE]
               [--weightedPopFile WEIGHTEDPOPFILE] [-t TOTPOSSPANNEDNEUTRAL]
               [--dont-PASS] [-o OUTPUTDIR] [--debug]

optional arguments:
  -h, --help            show this help message and exit
  -b BEDFILENAME, --bedFileName BEDFILENAME
                        The bed file with the chromosome, start, end, and
                        subregions columns to run the scan over
  -n NEUTRALFILE, --neutralFile NEUTRALFILE
                        The neutral file name with tab seperated MAFs,
                        typically created by extractSubPopData
  -l LINEINDEX, --lineIndex LINEINDEX
                        Line index for which region to read in over bed file,
                        the idea is this allows for easy parallelization if
                        you want to split up bed regions for each job. Default
                        -1 will read in all regions.
  -v VCFFILENAME, --vcfFileName VCFFILENAME
                        Input vcf file to be analyzed (must have corresponding
                        tabix tbi file)
  -w WINDOWSIZE, --windowSize WINDOWSIZE
                        Size of window to search around position
  -p POPFILE, --popFile POPFILE
                        File with the populations, seperate from the weighted
                        set to be used in the combined analysis.
  --weightedPopFile WEIGHTEDPOPFILE
                        File with the populations and their weights
  -t TOTPOSSPANNEDNEUTRAL, --totPosSpannedNeutral TOTPOSSPANNEDNEUTRAL
                        Total number of positions spanned by neutral variants.
                        This is not the number of neutral variants but the
                        regions they represent because many sites are
                        monomorphic, so the corresponding vcf the neutral
                        sample was taken from will only have items for each
                        variant found.
  --dont-PASS           Variants must have quality PASS to be included
  -o OUTPUTDIR, --outputDir OUTPUTDIR
                        output file name or the directory for the MAF tables
                        if running over a single region
  --debug
```

This will output a file that has columns for the chromosome, position, statistics, and corresponding -log<sub>10</sub> p-values for the statistics.  Results are for each position within region(s) being scanned, there are columns for each population's log rank test statistic and the corresponding -log<sub>10</sub> p-value, along with the stratified and weighted stratified log rank test and -log<sub>10</sub> p-value. Examples for the analysis run on gnomad V2.1 can be found in the [Download Directory](https://upenn.box.com/v/genomeWidePCIT).
