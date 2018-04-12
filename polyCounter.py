# Find the length of tracts without analysis - this is so we can go lower than the usual thresholds without breaking the system

from Bio import SeqIO
from Bio import Seq
from natsort import natsorted
from genomecrawler import pathCrawl
from genomecrawler import createRawGenomeName
from Bossref import bossref
from collections import Counter

import html
import time
import re

def main() :
    targetPath = '../polyG/neisseria.both'

    cutoffs = [5, 4, 3, 2, 2, 2, 2, 2, 2]
    screen = [9, 6, 6, 5, 5, 3, 3, 3, 3]

    counts = {}
    fn = lambda x : addGenomeToCounts(x, counts, cutoffs, screen)
    pathCrawl(targetPath, ['gb', 'gbk', 'gbf'], fn)

    outputCounts(targetPath, counts)
    print("Done")

def getSpeciesName(genomeRecord) :
    features = genomeRecord.features
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

def addGenomeToCounts(path, counts, cutoffs, screen) :
    # exclude any tract that falls within this much of a contig
    edge_catch = 10

    found={}
    name = None

    # Load the genome
    for genomeRecord in SeqIO.parse(path, 'gb') :
        # Skip plasmids
        if 'plasmid' in genomeRecord.features[0].qualifiers :
            return

        if not name :
            name = createRawGenomeName(genomeRecord, path)
            species = getSpeciesName(genomeRecord)


        seq = genomeRecord.seq

        # Fire up Bossref
        finder = bossref.Bossref()
        #finder.addFilter(filter)
        finder.setMinimumLengths(cutoffs)
        if not finder.findInString(str(seq), len(seq)) :
            print("Poly G finding failed!!")
            return []

        numFound = finder.numRepeats()

        sequenceEndEdge = len(seq)-edge_catch

        for i in range(0, numFound) :
            repeat = finder.getRepeat(i)
            # Not interested in multiple-of-3 repeats
            #if repeat.unitSize % 3 == 0 :
            #    continue

            # Exclude tracts that start or end into the end of the sequence
            if repeat.at < edge_catch or repeat.at+repeat.repeats*repeat.unitSize > sequenceEndEdge :
                continue

            unit = repeat.getUnit()

            unit_rev = Seq.reverse_complement(unit)
            if unit_rev < unit :
                unit = str(unit_rev)

            length = repeat.repeats

            if unit not in found :
                found[unit] = Counter()
            found[unit][length] += 1

    # Remove any units where the longest length doesn't meet the usual screening
    removeThese = []
    for unit, lengths in found.items() :
        longest = max(lengths.keys())
        if longest < screen[len(unit)-1] :
            removeThese.append(unit)

    for unit in removeThese :
        found.pop(unit, None)

    # Add the data to the counts
    counts[name] = {'repeats' : found, 'species' : species}

def outputCounts(path, counts) :
    print("Outputting")

    styler = html.Styler()
    styler.addStyle(['table tr {vertical-align:top}', 'table tr:nth-child(odd) td{ background-color:#f4f4f4 }','table tr:nth-child(even) td{ background-color:#dddddd }'])

    with html.Doc(path+'/allTractLengths.html', title='Tract Lengths', styler=styler) as htmldoc :
        htmldoc.add("H1", "Analysis of tract lengths")
        htmldoc.addRule()
        htmldoc.add("p",  'Path: ' + path)
        htmldoc.add("p", "Created: "+time.strftime('%Y-%m-%d %H:%M:%S (%Z)'))
        htmldoc.addRule()

        def addCountTable(data) :
            # Find maximum length value
            highest = 0
            lowest = 1000
            for lengths in data.values() :
                high = max(lengths.keys())
                if high > highest :
                    highest = high
                low = min(lengths.keys())
                if low < lowest :
                    lowest = low

            units = [k for k in data.keys()]
            units.sort(key=lambda x:(len(x), x))
            with htmldoc.table() as table :
                row = ['Tract']+[str(x) for x in range(lowest, highest+1)]
                table.addHeaderRow(row)
                for unit in units :
                    row = [unit]+['.' ]*(highest-lowest+1)
                    lengths = data[unit]
                    for length, n in lengths.items() :
                        row[length-lowest+1] = n
                    table.addRow(row)

        # Start with an overall summary table
        htmldoc.add('H1', 'Totals for all genomes')
        sum = {}
        for _, units in counts.items() :
            for unit, lengths in units['repeats'].items() :
                if unit not in sum :
                    sum[unit] = Counter()
                sum[unit] += lengths

        addCountTable(sum)
        htmldoc.addRule()

        # Tract lengths for each length of repeat unit
        sum = {}

        for _, units in counts.items() :
            for unit, lengths in units['repeats'].items() :
                tag = str(len(unit))
                if tag not in sum :
                    sum[tag] = Counter()
                sum[tag] += lengths

        addCountTable(sum)
        htmldoc.addRule()

        # Species by species data
        htmldoc.add('H1', 'Totals by species')
        speciesSum = {}
        for _, units in counts.items() :
            species = units['species']
            if species not in speciesSum :
                speciesSum[species] = {}
            for unit, lengths in units['repeats'].items() :
                if unit not in speciesSum[species] :
                    speciesSum[species][unit] = Counter()
                speciesSum[species][unit] += lengths

        for species, sum in speciesSum.items() :
            htmldoc.add('H2', 'Totals for '+species)
            addCountTable(sum)
            htmldoc.addRule()

        # Then a strain breakdown for each repeat type
        htmldoc.add('H1', 'Strain-by-strain counts')
        sortedStrainNames = natsorted(counts.keys())
        for strain in sortedStrainNames :
            htmldoc.add('H2', 'Breakdown for ' + counts[strain]['species']+' '+strain)
            addCountTable(counts[strain]['repeats'])
            htmldoc.addRule()

if __name__ == '__main__' :
    main()