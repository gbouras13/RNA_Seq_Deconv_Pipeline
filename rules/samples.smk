"""
Function for parsing the 'Bams' config and identifying samples and Bam files
"""

from itertools import chain

def samplesFromDirectory(dir):
    """Parse samples from a directory"""
    outDict = {}
    # https://stackoverflow.com/questions/11860476/how-to-unnest-a-nested-list
    samples= glob_wildcards(os.path.join(dir,'{sample}.bam'))
    samples2 = chain(*samples)
    for sample in samples2:
        outDict[sample] = {}
        bam = os.path.join(dir,f'{sample}.bam')
        if os.path.isfile(bam):
            outDict[sample]['bam'] = bam
        else:
            sys.stderr.write("\n"
                             "    FATAL: Error globbing files."
                             f"    {bam} \n"
                             "    does not exist. Ensure consistent formatting and file extensions."
                             "\n")
            sys.exit(1)
    return outDict

def parseSamplesBam(readFileDir):
    """Parse samples from a directory"""
    if os.path.isdir(readFileDir):
        sampleDict = samplesFromDirectory(readFileDir)
    else:
        sys.stderr.write("\n"
                         f"    FATAL: {readFileDir} is neither a file nor directory.\n"
                         "\n")
        sys.exit(1)
    if len(sampleDict.keys()) == 0:
        sys.stderr.write("\n"
                         "    FATAL: We could not detect any samples at all.\n"
                         "\n")
        sys.exit(1)
    return sampleDict

