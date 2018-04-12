# Find all poly-G/C tracts in a genome and report them
import platform

from Bio.Seq import Seq
from Bio.Seq import reverse_complement
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import FeatureLocation

from blastmatch import getFastaName
import csv

import featurefinder
import os

# Guess what the full length product looks like
def guessFullTract(genomeRecord, geneLocation, poly) :
    """Create 0/+1/+2 tract lengths and see which gives the longest product"""

    #Start by grabbing the gene + 5000bp downstream
    start = int(geneLocation.start)
    end = int(geneLocation.end)
    if geneLocation.strand == 1 :
        seq = genomeRecord.seq[start : min(end + 10000, len(genomeRecord.seq))]
        offset = poly['at'] - start
        repeatUnit = poly['repeat']
    else :
        seq = genomeRecord.seq[max(start-10000,0) : end]
        seq = seq.reverse_complement()
        offset = end - poly['at']
        repeatUnit = reverse_complement(poly['repeat'])

    # Try varying the polyG to see what happens to the reading frame
    longest = -1
    longestExtra = 0
    longestTrans = Seq('', IUPAC.unambiguous_dna)
    for extra in range(0, 3) :
        modSeq = seq[0:offset] + Seq(repeatUnit*extra, IUPAC.unambiguous_dna) + seq[offset:]
        # Pad to multiple of 3
        modSeq = modSeq + Seq('A'*(3-(len(modSeq) % 3)), IUPAC.unambiguous_dna)
        translation = modSeq.translate(table=11, to_stop=True) # table 11 is the bacterial translation table
        if len(translation) > longest :
            longest = len(translation)
            longestExtra = extra
            longestTrans = translation

    # Now we have to deduce our true start and end positions from the longest tract length,
    # remembering to correct for the extra bases
    # Note that we only extend the annotated length, not truncate it.

    if geneLocation.strand == 1 :
        newEnd = start + longest*3 - longestExtra*len(repeatUnit)
        if newEnd > end:
            end = newEnd
    else:
        newStart = end - longest*3 + longestExtra*len(repeatUnit)
        if newStart < start:
            start = newStart

    onLength = poly['length'] + longestExtra
    if abs(onLength - poly['length']) > abs(onLength -3 - poly['length']) :
                onLength -= 3

    # Assume the longest reading frame corresponds to the ON state and store its position and size
    return {
        'onLength' : onLength,
        'longestFrame' : longest,
        'translation' : longestTrans,
        'start' : start,
        'end' : end }

# Given a gene record and a location generate an info block on it.
def buildInformationOn(genomeRecord, geneRecord, poly, namer) :

    # Work out what the current tract length is
    tractLength = poly['length']
    location = poly['at']

    if not geneRecord :
        # Return blank information
        return {
             'location' : location,
             'start' : location,
             'end'  : location + tractLength,
             'name' : 'unknown',
             'locus' : 'unknown',
             'gene' : 'unknown',
             'geneLength' : 3,
             'tract' : poly['repeat'] + str(tractLength),
             'repeatUnit' : poly['repeat'],
             'repeats' : tractLength,
             'onLength' : None,
#             'sequence' : None ,
             'translation' : None,
             'function' : 'Unknown',
             'record' : None,
             'matchScores' : {},
             'outby' : 0}

    # We cannot presume the gene record correctly identifies the full length product since frameshifts
    # So we want to begin by guessing what it should be
    fullTract = guessFullTract(genomeRecord, geneRecord.location, poly)
    start = fullTract['start']
    end = fullTract['end']

    # Find the function annotation
    if 'product' in geneRecord.qualifiers:
        function = geneRecord.qualifiers['product'][0]
    else:
        function = 'No annotation data'

    if 'locus_tag' in geneRecord.qualifiers:
        locus = geneRecord.qualifiers['locus_tag'][0]
        autoLocus = False
    else:
        autoLocus = True
        locus = namer.nextName()

    if 'gene' in geneRecord.qualifiers:
        gene = geneRecord.qualifiers['gene'][0]
    else:
        gene = 'No Data'

    # The distance from the start depends on which strand we're on
    strand = geneRecord.strand

    if not strand :
        # No strand information... bugger
        strand = 0
        repeatUnit = poly['repeat']
    elif strand == 1 :
        # On the 'sense' strand
        repeatUnit = poly['repeat']
    else :
        # On the 'antisense' strand
        # So get the reverse complement
        repeatUnit = reverse_complement(poly['repeat'])

    tract = repeatUnit+str(tractLength)

    # And we need to check for being out of the thing
    if location < start or location > end :
        # Out of the gene
        if strand == 0 :
            outby = 0
            if location < start :
                name = str(start - location) + "bp away from "+locus
            else:
                name = str(location-end) + "bp away from " + locus
        elif strand == 1 :
            if location < start:
                name = str(start - location) + "bp upstream of " + locus
                outby = -(start-location)
            else:
                name = str(location - end) + "bp downstream of " + locus
                outby = location-end
        else :
            if location < start :
                name = str(start-location) + "bp downstream of "+locus
                outby = start-location
            else:
                name = str(location - end) + "bp upstream of " + locus
                outby = -(location-end)
    else :
        name = locus
        outby = 0

    sequence = genomeRecord.seq[start:end]
    if strand == -1 :
        sequence = sequence.reverse_complement()

    return {
        'location' : location,
        'start' : start,
        'end' : end,
        'name' : name,
        'locus' : locus,
        'gene' : gene,
        'geneLength' : end-start,
        'strand' : strand,
        'tract' : tract,
        'repeatUnit' : repeatUnit,
        'repeats' : tractLength,
        'onLength' : fullTract['onLength'] if outby == 0 else None,
#        'sequence' : sequence ,
        'translation' : fullTract['translation'],
        'function' : function,
        'record' : geneRecord,
        'matchScores' : {},
        'outby' : outby,
        'autoLocus' : autoLocus}

def findAllPolysInStrain(strain, filter, minLengths, namer) :
    ''' Use bossref to find all the polys in a strain '''

    # Exclude tracts that run too close to the edge (usually artefacts of sequencing)
    edge_catch = 10

    # Work out which executable to use based on platform
    runningOn = platform.system()
    if platform.system() == 'Windows' :
        exe = 'bossref.exe'
    else :
        exe = './bossref'

    # Work out path to the fasta version
    fastaPath = os.path.splitext(strain['path'])[0] + ".fasta"
    tempPath = os.path.join(os.path.dirname(strain['path']), 'temp', os.path.basename(strain['path'])+'.tab')

    # Build the string
    arguments = ' -i "'+fastaPath+'" -o "'+tempPath+'" -c "'+' '.join(str(n) for n in minLengths)+'"'
    if len(filter) > 0 :
        arguments += ' -f '+filter

    # Run it
    result = os.system(exe+arguments)

    if (result != 0) :
        print('ERROR: Bossref failed for "'+strain['path']+'"')
        return []

    # Read in the polys for each contig
    # Note: Bossref always returns contigs in the order they appear
    contigs = strain['contigs']
    onContig = 0
    contigName = getFastaName(contigs[onContig])
    sequenceEndEdge = len(contigs[onContig]['record'].seq)-edge_catch
    polys = [[]]
    with open(tempPath) as csvFile :
        reader = csv.DictReader(csvFile, dialect=csv.excel_tab)
        for row in reader :
            if row['record'] != contigName :
                # Find the matching contig
                while True :
                    onContig += 1
                    if onContig >= len(contigs) :
                        print('ERROR: mismatching names in bossref search (internal logic failure) on file '+strain['path'])
                        return []
                    polys.append([])
                    contigName = getFastaName(contigs[onContig])
                    if  contigName == row['record'] :
                        sequenceEndEdge = len(contigs[onContig]['record'].seq)-edge_catch
                        break

            # Add to our list, excluding those that run too close to the edge.
            location = int(row['location'])
            length = int(row['length'])
            tracttype = row['tracttype']

            if location < edge_catch or location+len(tracttype)*length > sequenceEndEdge :
                continue

            polys[onContig].append({'at' : location, 'repeat' : tracttype, 'length' : length})

    for i in range(onContig, len(contigs)) :
        polys.append([])

    # Run the identification code
    found = []
    for n, contig in enumerate(contigs) :
        found.append(identifyPolys(contig['record'], polys[n], namer))

    return found

def identifyPolys(genomeRecord, polys, namer):
    ''' Identify the polys we have found '''
    found = []

    # Find the matching feature for each tract
    for poly in polys:
        location = poly['at']
        genes = featurefinder.findMatchingFeatures(genomeRecord, location, ['CDS', 'rRNA', 'gene'])

        add = None
        if len(genes) == 0:
            # Blast, no matching gene, find the closest thing we can
            # Quite a few tracts seem to be in plausible ORFs but unannotated
            # We will try and identify these and hope they group with other, annotated genes
            # we can use to discern purpose
            add = checkForORF(genomeRecord, location, poly, namer)

            # Can't identify a possible gene, assign it to the best feature
            if not add :
                _,geneRecord = featurefinder.findClosestFeature(genomeRecord, location, ['CDS', 'rRNA', 'gene'], 200)
        else:
            # Should handle this better
            geneRecord = genes[0]

        if not add :
            add = buildInformationOn(genomeRecord,geneRecord,poly, namer)

        add['geneGroup'] = add['uid'] = getUid()
        add['links'] = set()

        found.append(add)

    return found

def checkForORF(genomeRecord, location, poly, namer) :
    ''' Find a potential ORF spreading out from the given location '''

    # Need to consider all six reading-frames and potential frame shift across the tract
    # Accept the longest frame we can find

    # Conservative cutoff point for the minimum length of the identified ORF - 300aa
    cutoff = 300

    orf = None
    longest = 0
    for frame in [-3,-2,-1, 1, 2, 3] :
        newOrf = findORFinFrame(genomeRecord, location, frame, poly)
        if not newOrf :
            continue

        longestFrame = newOrf['longestFrame']
        if longestFrame > cutoff and longestFrame > longest :
            longest = longestFrame
            orf = newOrf

    if not orf :
        return None

    # Need to fill out full information for this tract?

    name = namer.nextName()
    if orf['frame'] > 0 :
        repeatUnit = poly['repeat']
    else :
        repeatUnit = reverse_complement(poly['repeat'])

    #Buildinfo
    return {
        'location' : location,
        'start' : orf['start'],
        'end' : orf['end'],
        'name' : name,
        'locus' : name,
        'gene' : 'No Data',
        'geneLength' : orf['end']-orf['start'],
        'strand' : 1 if orf['frame'] > 0 else -1,
        'tract' : repeatUnit +str(poly['length']),
        'repeatUnit' : repeatUnit,
        'repeats' : poly['length'],
        'onLength' : orf['onLength'],
#        'sequence' : sequence ,
        'translation' : orf['translation'],
        'function' : 'No data, identified by script',
        'record' : None,
        'matchScores' : {},
        'outby' : 0,
        'autoLocus' : True}

def findORFinFrame(genomeRecord, location, frame, poly) :
    # Take the sequence before, or after, the location and scan it for a potential starting site
    upstream = 900
    if frame > 0 :
        seq = genomeRecord.seq[max(location-upstream + (frame-2),0):location+(frame-2)]
    else:
        # Need reverse complement for the frame's going on the counter strand
        seq = genomeRecord.seq[location+(frame+1): location+upstream+(frame+1)]
        seq = seq.reverse_complement()

    seq = seq + Seq('A'*(3-(len(seq) % 3)), IUPAC.unambiguous_dna)
    trans = seq.translate(table=11, to_stop=True)

    # Scan back through the sequence until we find a stop codon
    # Checking for Ms along the way
    # Note that non-M start sites exist but we ignore them as they are low in occurence
    M = len(trans)
    on = len(trans)-1
    while on > 0 and trans[on] != '*' :
        if trans[on] == 'M' :
            M = on
        on -= 1

    # No M before stop codon
    if M == len(trans) :
        return None

    if frame > 0 :
        start = location - (len(trans)-M)*3 + (frame-2)
        loc = FeatureLocation(start, location, strand = 1)
    else :
        end = location + (len(trans)-M)*3 + (frame+1)
        loc = FeatureLocation(location, end, strand = -1)


    full = guessFullTract(genomeRecord,loc, poly)
    full['frame'] = frame

    return full

def getUid() :
    nextUid = getUid.counter
    getUid.counter += 1
    return nextUid
getUid.counter = 0
