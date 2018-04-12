# Simply bit of code - take the fasta files from this directory that DON'T correspond to plasmids and
# copy them to the second directory
# This is so we can pass them to the FFP program and get a tree back

from genomecrawler import pathCrawl
from Bio import SeqIO
import shutil
import os.path


def main() :
    inDir = '../polyG/campy_spp'
    outDir = '../polyG/campy_spp_fastas'

    fn = lambda x : copyIfNotPlasmid(x, outDir)

    pathCrawl(inDir, ['gb', 'gbk', 'gbf'], fn)

def copyIfNotPlasmid(target, outDir) :
    genomeRecord = next(SeqIO.parse(target, "gb"))
    if 'plasmid' in genomeRecord.features[0].qualifiers :
        return


    fapath = os.path.splitext(target)[0] + ".fasta"
    print('Copying "'+fapath+'"')
    shutil.copy2(fapath, outDir)

if __name__ == "__main__":
    main()