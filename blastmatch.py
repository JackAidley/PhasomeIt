__author__ = 'Jack'

import os.path
#from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Blast import NCBIXML
#from multiprocessing import Pool
import subprocess
import featurefinder
import re
import time
#from collections import deque

import settings

# use BLAST to search for similar genes in the other strains
# and find out whether these matches are themselves phase variable
def makeBlastTemp(db, path) :
    # Create temporary query file for the blast containing every polyG tract found so we
    # can use it to blast against our genome set

    # Request temp file maybe?

    # Create a temporary fasta file containing all sequences to search
    with open(os.path.join(path,'temp/pvtracts.fa'), "w") as fastaTracts :
        for strain in db :
            # Walk the tracts
            for i, contig in enumerate(strain['contigs']) :
                for tract in contig['tracts'] :
                    if not tract['translation'] or len(tract['translation']) < 20:
                        tract['no_blast'] = True
                        continue
                    # Output the name
                    fastaTracts.write('>strain='+strain['name']+'||contig='+str(i)+'||tract='+tract['name']+"\n")
                    # Output the DNA sequence
                    fastaTracts.write(str(tract['translation'])+"\n")

def makeBlastDatabases(db) :
    # Create blast databases for everything
    for strain in db :
        # Invoke makeblastdb to create blast databases for each strain
        rootpath = os.path.splitext(strain['path'])[0]

        fapath = rootpath+'.fasta'
        dbpath = rootpath+'.nin'
        if os.path.exists(dbpath):
            # check if the db file is up to date
            dbdate = os.stat(dbpath).st_mtime
            fadate = os.stat(fapath).st_mtime

            # Check the other two dbfiles exist too
            if os.path.exists(rootpath+".nhr") and os.path.exists(rootpath+".nsq") :
                if dbdate > fadate:
                    continue

        # Invoke the makeblastdb command
        command = 'makeblastdb  -in "'+rootpath+'.fasta" -title "'+strain['name']+'" -hash_index -dbtype nucl -out "'+rootpath+'"'
        os.system(command)

def matchInPV(location, contig) :
    '''
    Find whether we have an existing PV tract that is in this location
    Need to do this first in case we have extra tracts we annotated during the crawl
    '''

    for n, tract in enumerate(contig['tracts']) :
        # See whether this match hits the location of the feature
        if location > tract['start'] and location <= tract['end'] :
            # We have a hit
            match = { 'tractNo' : n, 'locus' : tract['locus'] }
            break
    else :
        return None

    # Build tract details
    match['record'] = tract['record']
    match['function'] = tract['function']
    match['geneLength'] = tract['geneLength']

    return match

def annotateMatch(match, contig) :
    '''
    Find information about this match
    For a match that's already in the found PV tracts: locate it
    For new matches, add information about the gene
    '''

    geneRecord = match['record']
    match['geneLength'] = geneRecord.location.end - geneRecord.location.start

    # Check for it in our existing list of tracts
    for n, tract in enumerate(contig['tracts']) :
        # See whether this match hits the location of the feature
        if geneRecord == tract['record'] :
            # We have a hit
            match['tractNo'] = n
            match['locus'] = tract['locus']
            #match['tract'] = tract
            return

    # Find details about the match
    if 'product' in geneRecord.qualifiers:
        function = geneRecord.qualifiers['product'][0]
    else:
        function = 'No annotation data'

    if 'locus_tag' in geneRecord.qualifiers:
        locus = geneRecord.qualifiers['locus_tag'][0]
    else:
        locus = locus = contig['record'].name + ':' + str(geneRecord.location)

    if 'gene' in geneRecord.qualifiers:
        gene = geneRecord.qualifiers['gene'][0]
    else:
        gene = 'No Data'

    match['function'] = function
    match['locus'] = locus
    match['gene'] = gene

def makeBlastCommand(path, compareToPath, compareToName) :
    # blast all the tracts found against the database for compareTo
    dbpath = path+'/'+os.path.splitext(os.path.basename(compareToPath))[0] # + '.nin'
    querypath = path+'/temp/pvtracts.fa'
    outpath = path+'/temp/out_' + os.path.basename(compareToPath) + ".xml"

    min_evalue = 1e-6

    # Do the blast
    # Had to abandon NcbitblastnCommandline because
    # it was calling subprocess.popen and it was crashing out
    blast_cmd = 'tblastn'
    blast_cmd += ' -out "'+outpath+'" '
    blast_cmd += ' -outfmt 5' # XML format
    blast_cmd += ' -query "'+querypath+'"'
    blast_cmd += ' -db "'+dbpath+'"'
    blast_cmd += ' -evalue '+str(min_evalue)
    blast_cmd += ' -db_gencode 11' # bacterial codons
    blast_cmd += ' -matrix BLOSUM80'

    return blast_cmd

def blastAgainst(path, compareToPath, compareToName) :
    blastCmd = makeBlastCommand(path, compareToPath, compareToName)
    childProcess = subprocess.Popen(blastCmd, stdin=subprocess.PIPE,
                                     stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                     universal_newlines=True,shell=True)
    return childProcess

def getFastaName(contig) :
    description = contig['record'].description
    id = contig['record'].id
    if description.split(None, 1)[0] == id:
        return description
    else :
        return id + ' ' + description


def readBlast(db, path, compareTo) :
    outpath = os.path.join(path,'temp/out_' + os.path.basename(compareTo['path']) + ".xml")
    min_evalue = 1e-6
    min_coverage = settings.homologyCutoffTo #0.5
    min_query_coverage = settings.homologyCutoffFrom #0.4

    # Create a lookup table from contig names to the number of the contigs
    contigLookup = {}
    for i, contig in enumerate(compareTo['contigs']) :
        contigLookup[getFastaName(contig)] = i

    with open(outpath) as outHandle :
        # Crawl across all the hits
        for record in NCBIXML.parse(outHandle) :
            matches = []
            for alignment in record.alignments :

                # Find the contig of the hit
                contigNo = contigLookup[alignment.hit_def]
                contigHit = compareTo['contigs'][contigNo]

                for hsp in alignment.hsps :
                    if hsp.expect > min_evalue :
                        break

                    # find the gene for each hit
                    location = (hsp.sbjct_start + hsp.sbjct_end)//2

                    match = matchInPV(location, contigHit)
                    if match :
                        match['expect'] = hsp.expect
                        match['length'] = hsp.identities
                        match['hsp'] = hsp
                        match['contigNo'] = contigNo
                        match['contigName'] = contigHit['record'].description
                        matches.append(match)
                        continue

                    matchRecords = featurefinder.findMatchingFeatures(contigHit['record'], location, ['CDS', 'rRNA', 'gene'])
                    # store the result
                    if matchRecords :
                        match = {'record' : matchRecords[0], 'expect' : hsp.expect, 'length' : hsp.identities, 'hsp' : hsp, 'contigNo': contigNo, 'contigName' : contigHit['record'].description}
                        # Match this with an existing gene record if available
                        annotateMatch(match, contigHit)
                        matches.append(match)

            # Now, if we have matches we need to associate them with the right gene in our records
            if matches :
                # But first we want to group any matches that are to the same gene together
                # And then screen them out if the TOTAL amount matched is less than 50%
                groupedMatches = []
                lociMatched = set([m['locus'] for m in matches])
                for locus in lociMatched :
                    locusMatches = [m for m in matches if m['locus'] == locus]
                    ourMatch = dict(locusMatches[0])
                    if len(locusMatches) > 1 :
                        ourMatch['multipleHits'] = []
                        ourMatch['length'] = 0
                        for duplicate in locusMatches :
                            ourMatch['multipleHits'].append(duplicate['hsp'])
                            ourMatch['length'] = ourMatch['length'] + duplicate['length']
                        ourMatch['multipleHits'].sort(key = lambda x : x.query_start)

                    # NB: geneLength in bp, length in amino acids, hence *3
                    ourMatch['coverage'] = (ourMatch['length']*3)/ourMatch['geneLength']
                    if ourMatch['coverage'] > min_coverage :
                        groupedMatches.append(ourMatch)

                if groupedMatches :
                    for groupedMatch in groupedMatches :
                        # Get the strain and name out of the record.query
                        queryStrainName, queryContigNo, queryTractName = splitQueryName(record.query)
                        queryContig = findStrain(db, queryStrainName)['contigs'][int(queryContigNo)]
                        geneMatch = next(x for x in queryContig['tracts'] if x['name'] == queryTractName)
                        if not geneMatch :
                            print("Something has gone wrong, gene '"+record.query+"' not found")
                            continue

                        groupedMatch['queryCoverage'] = (groupedMatch['length']*3)/geneMatch['geneLength']
                        if groupedMatch['queryCoverage'] < min_query_coverage :
                            continue

                        if 'blastMatch' not in geneMatch :
                            geneMatch['blastMatch'] = {}

                        # Create a link from the query match to the subject (hit) match
                        contigLinkName = groupedMatch['contigName']
                        if contigLinkName not in geneMatch['blastMatch'] :
                            geneMatch['blastMatch'][contigLinkName] = []
                        geneMatch['blastMatch'][contigLinkName].append(groupedMatch)

                        # Do we have a link to another PV gene?
                        # If so create a bidirectional link
                        if 'tractNo' in groupedMatch :
                            otherGene = compareTo['contigs'][ groupedMatch['contigNo'] ]['tracts'][ groupedMatch['tractNo'] ]
                            geneMatch['links'].add(otherGene['uid'])
                            otherGene['links'].add(geneMatch['uid'])


def findStrain(db, strainName) :
    return next(s for s in db if s['name'] == strainName)

def findContig(db, contigName) :
    return next(c for s in db for c in s['contigs'] if c['record'].description == contigName)

def findStrainForContig(db, contigName) :
    return next(s for s in db for c in s['contigs'] if c['record'].description == contigName)

def splitQueryName(queryName) :
    parsed = re.match(r'strain=(.+)\|\|contig=(.+)\|\|tract=(.+)', queryName)
    return parsed.group(1), parsed.group(2), parsed.group(3)

def groupMatchedGenes(db) :
    '''Assign all found genes a group number'''

    print('Grouping', end='', flush=True)

    # Method
    # First pass: create bidirectional links for every gene's blast hits
    #   ---> Now done during the blast read
    # Second pass: recursively follow links to assign the lowest group number

    tractLookup = {}
    for strain in db :
        for contig in strain['contigs'] :
            for tract in contig['tracts'] :
                tractLookup[tract['uid']] = tract

    # Walk the tree, setting the group
    for strain in db :
        for contig in strain['contigs'] :
            for tract in contig['tracts'] :
                # Ungrouped genes
                if 'no_blast' in tract :
                    tract['geneGroup'] = -1
                    continue

                group = tract['geneGroup']
                if group != tract['uid'] :
                    continue
                targets = set(link for link in tract['links'])
                while targets :
                    target = tractLookup[targets.pop()]
                    if target['geneGroup'] == group :
                        continue
                    target['geneGroup'] = group
                    targets |= set(link for link in target['links'] if tractLookup[link]['uid'] != group)
                print('.', end='', flush=True)


    # Renumber the remaining groups from 0 to n
    fixupGeneGroupNumbers(db)
    print('!')

def fixupGeneGroupNumbers(db) :
    groupLookup = {-1:-1}
    nextGroup = 0
    for strain in db :
        for contig in strain['contigs'] :
            for tract in contig['tracts'] :
                oldGroup = tract['geneGroup']
                if oldGroup in groupLookup :
                    tract['geneGroup'] = groupLookup[oldGroup]
                else :
                    tract['geneGroup'] = nextGroup
                    groupLookup[oldGroup] = nextGroup
                    nextGroup += 1

def blastComparePVs(db, path) :
    
    # BLAST all the PV genes found

    print('Preparing to blast', end='')

    # First create the temporary files with the tracts in
    makeBlastTemp(db, path)
    print('.', end='')

    # And create blastdb files if needed
    makeBlastDatabases(db)
    print('.')

    # Then blast each one against all the other strains
    print('Blasting', flush=True)

    startT = time.process_time()

    repeat = 0
    targets = [(strain['path'], strain['name']) for strain in db]

    maxPoolSize = 4
    pool = []
    failed = []
    while(True) :
        # Re-fill the pool
        while len(pool) < 4 and len(targets) :
            t = targets[0]
            targets = targets[1:]
            pool.append({'name' : t[1], 'path' : t[0], 'process' : blastAgainst(path, t[0], t[1])})

        # Check the currently running processes for completion
        for process in pool :
            try :
                returncode = process['process'].wait(timeout = 0.1)
            except subprocess.TimeoutExpired :
                continue
            #returncode = process['process'].poll()
            if returncode is not None:
                # We've finished, yay!
                if returncode == 0 :
                    print('...'+process['name']+' completed')
                else :
                    print('...'+process['name']+' ('+process['path']+') failed with returncode 0x{0:02X}'.format(returncode) + ' -> '+str(process['process'].stderr.readlines()))
                    failed.append((process['path'], process['name']))

                process['process'] = None

        # Remove any removed processes from the pool
        pool = [process for process in pool if process['process']]

        # Check for completion
        if not pool and not targets :
            break

    if failed :
        print("Some strains failed to blast correctly, will have missing data")
        for t in failed :
            print('... BLAST failed for '+t[1])
        print('I can produce a .bat file you can run manually to produce the necessary files (Windows only)')
        produceBat = input('Shall I produce the .bat file? [type "yes" if you want this, or press ENTER if not]').lower()
        if produceBat == 'yes' or produceBat == 'y' :
            with open(path+'/_fixblasts.bat', 'w') as batFile :
                print('cd '+ os.getcwd(), file=batFile)
                for f in failed :
                    print(makeBlastCommand(path, f[0], f[1]), file=batFile)
                print('pause', file=batFile)
            print('Batch file produced. Find and run '+path+'/_fixblasts.bat')
            ranOk = input('Did it run okay? [type "yes" or "no"]').lower()
            if ranOk == 'y' or ranOk == 'yes' :
                failed = []



    print('Reading BLASTs', end='', flush=True)
    failedNames = [t[1] for t in failed]
    for strain in db :
        if strain['name'] not in failedNames :
            readBlast(db, path, strain)
            print('.', end='', flush=True)

    print('\n Total BLAST time: ({0:.2f}s)'.format(time.process_time()-startT))

    # Now we've got all the matching genes, we can network them into groups
    groupMatchedGenes(db)

    #from pprint import pprint
    #pprint(db)
