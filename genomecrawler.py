import os
import sys
from polygfinder import findAllPolysInStrain
from Bio import SeqIO
from natsort import natsorted
import re
import logging

import settings

def main() :
    db = crawlGenomes('c.jejuni')
    for strain in db :
        print(strain['name'], len(strain['tracts']))

def crawlGenomes(path, filter, minLengths, plugins) :
    db =[]
    fn = lambda x : addGenomeToDatabase(db, x, filter, minLengths)
    pathCrawl(path, ['gb', 'gbk', 'gbf', 'embl'], fn)
    if plugins :
        for plugin in plugins :
            if hasattr(plugin, 'annotateDb') :
                plugin.annotateDb(db)
    db = reorderDatabase(db)

    return db

def reorderDatabase(db) :
    precedence = settings.firstStrains #['NCTC 11168','81-176', '81-176 - (pTet)', '81-176 - (pVir)', 'PT14']
    speciesPrecedence = settings.firstSpecies #['C. jejuni', 'C. coli']

    def sortFunc(strain) :
        strainName = strain['name']
        species = strain['species']

        if species in speciesPrecedence :
            key = str(speciesPrecedence.index(species))
        else :
            key = str(len(speciesPrecedence))

        if strainName in precedence :
            key = key + '/'+str(precedence.index(strainName))
        else :
            key = key + '/'+str(len(precedence))
        return key+species+strainName

    return natsorted(db, key=sortFunc)

def pathCrawl(target, types, function) :
    print('Crawling', end='')
    sys.stdout.flush()
    for root, dirs, files in os.walk(target) :
        for file in files :
            ext = os.path.splitext(file)[1][1:]
            if ext not in types :
                continue

            path = os.path.join(root, file)
            function(path)
            print('.', end='')
            sys.stdout.flush()
    print('')

def createRawGenomeName(genomeRecord, path) :
    qualifiers = genomeRecord.features[0].qualifiers
    if 'strain' in  qualifiers:
        name = qualifiers['strain'][0]
    elif 'isolate' in qualifiers :
        name = qualifiers['isolate'][0]
    else :
        name = os.path.basename(path)

    if 'plasmid' in qualifiers :
        name = name + " - (" + qualifiers['plasmid'][0] + ")"

    return name

def createGenomeName(genomeRecord, path, db) :
    name = createRawGenomeName(genomeRecord, path)
    duplicates = 0
    for strain in db :
        if strain['rawName'] == name :
            duplicates += 1

    if not duplicates :
        return name

    name = name + ' (' + str(duplicates+1) + ')'

    return name

class GenomeNamer :
    def __init__(self,baseName) :
        self.baseName = baseName
        self.count = 0

    def nextName(self) :
        name = self.baseName + ' [' + str(self.count) + ']'
        self.count += 1
        return name

def addGenomeToDatabase(db, path, filter, minLengths) :
    contigs = []
    # The namer comes into play for files without clear annotation
    ext = os.path.splitext(path)[-1]
    format = 'embl' if ext.lower() == '.embl' else 'gb'
    for genomeRecord in SeqIO.parse(path, format) :
        # Skip plasmids
        if 'plasmid' in genomeRecord.features[0].qualifiers :
            continue
        if not len(genomeRecord.description) :
            genomeRecord.description = genomeRecord.id
        sortAndMapFeatures(genomeRecord)
        contigs.append({'record' : genomeRecord, 'inGenome' : path} )
    
    # If we had no non-Plasmid contigs, then we leave without doing anything
    if len(contigs) == 0 :
        logging.info('File {0} skipped, contains only plasmids'.format(path))
        return

    strain = {
        'name' : createGenomeName(genomeRecord, path, db),
        'rawName' : createRawGenomeName(genomeRecord, path),
        'path' : path,
        'contigs' : contigs }

    strain['species'] = getSpeciesName(strain)

    makeFastaVersion(strain)
    namer = GenomeNamer(strain['name'])
    polys = findAllPolysInStrain(strain, filter, minLengths, namer)
    for n, found in enumerate(polys) :
        strain['contigs'][n]['tracts'] = found

    db.append(strain)
    logging.info('File {0} added to database.'.format(path))

def makeFastaVersion(strain) :
    # Create .fasta versions of each Genbank file so we can then make blast databases from them and use bossref
    gbpath = strain['path']
    fapath = os.path.splitext(gbpath)[0] + ".fasta"
    if os.path.exists(fapath) :
        # check if the fasta file is up to date
        gbdate = os.stat(gbpath).st_mtime
        fadate = os.stat(fapath).st_mtime

        if fadate > gbdate :
            return

    SeqIO.write([contig['record'] for contig in strain['contigs']], fapath, "fasta")

def getSpeciesName(strain) :
    features = strain['contigs'][0]['record'].features
    if not features or len(features) < 1 :
        return 'Unknown sp.'
    qualifiers = features[0].qualifiers

    organism = None
    if 'organism' in qualifiers :
        organism = qualifiers['organism'][0]
    elif 'source' in qualifiers :
        organism = qualifiers['source'][0]

    species = 'Unknown sp.'
    if organism :
        match = re.match('(\S*) (\S*)', organism)
        if match :
            species = match.group(1)[0]+". "+match.group(2)

    return species

def sortAndMapFeatures(record) :
    # Sort the features to ensure they are ordered
    # But keep the source record first
    record.features = sorted(record.features, key=lambda x : (x.type != 'source', x.type, x.location.start))

    record.lookup = {}
    lastType = record.features[0].type
    first = 0
    for i, feature in enumerate(record.features) :
        if feature.type != lastType :
            record.lookup[lastType] = (first, i)
            first = i
            lastType = feature.type
    record.lookup[lastType] = (first, len(record.features))
    record.featureStarts = [f.location.start for f in record.features]

if __name__ == "__main__":
    main()
