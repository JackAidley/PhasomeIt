'Output the data found in html format'
import os
#import io
import time
import pickle
import copy
import traceback
import logging

from collections import Counter

import html
import shutil

#from Bio import Phylo
from Bio import Seq
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import _DistanceMatrix

import settings


def outputAsHTMLFromPickle(path, plugins = None) :
    '''Read in pickled data and convert it to an HTML format'''
    logging.info('Begun HTML output')

    if not plugins :
        plugins = []

    # Create the directory if it's not already there
    if not os.path.exists(path):
        os.makedirs(path)

    if not os.path.exists(path+'/strains') :
        os.makedirs(path+'/strains')

    if not os.path.exists(path+'/groups') :
        os.makedirs(path+'/groups')

    # Copy the sorttable script to the path
    if not os.path.exists(path+'/include') :
        os.makedirs(path+'/include')

    shutil.copy('resources/sorttable.js', path+'/include')

    with HtmlCreator(path) as creator :
        creator.plugins = plugins
        try :
            creator.createStrainFiles()
        except Exception :
            logging.exception('Failed during output of strain files')

        try :
            creator.createTractSummary()
        except Exception :
            logging.exception('Failed during output of tract summary')

        try :
            creator.createGroupIndex()
        except Exception :
            logging.exception('Failed during output of group index')

        try :
            creator.createCorePhasome()
        except Exception :
            logging.exception('Failed during output of core phasome')

        try :
            for plugin in plugins :
                if hasattr(plugin, 'addToHTML') :
                    plugin.addToHTML(creator)
        except Exception :
            logging.exception('Failed during addToHTML process for plugins')

def shortenName(name, max=10, start=2, end=6) :
    if len(name) <= max :
        return name

    return name[0:start]+'..'+name[-end:]

def getRepeat(tract) :
    '''From a tract of form 'A10' or 'ATG8' get the letters'''
    for n, c in enumerate(tract) :
        if c.isdigit() :
            break
    return tract[:n]

def getRepeats(tract) :
    '''From a tract of the form 'A10', or 'ATG8' get the number'''
    for n, c in enumerate(tract) :
        if c.isdigit() :
            break
    return int(tract[n:])


class HtmlCreator :
    '''Worker class to create Html summary of the phasome data'''
    def __init__(self, path) :
        self.styler = self.createStyler()
        self.path = path
        self.creationTime = time.strftime('%Y-%m-%d %H:%M:%S (%Z)')
        self.plugins = []

    def __enter__(self) :
        # Load guide pickle
        with open(os.path.join(self.path, '../pickle', 'guide.pkl'), 'rb') as file :
            self.guide = pickle.load(file)

        # Load gene group pickle
        with open(os.path.join(self.path, '../pickle', 'groups.pkl'), 'rb') as file :
            self.geneGroups = pickle.load(file)

        # Create contig->strain and strain number lookups
        self.contigToStrain = {}
        self.strainLookup = {}
        for i, strain in enumerate(self.guide['strains']) :
            self.strainLookup[strain['name']] = i
            for contig in strain['contigs'] :
                self.contigToStrain[contig] = strain['name']

        self.index = html.Doc(self.path+"/index.html", title="Index", styler=self.styler, sortTable = 'include/sorttable.js')
        self.index.__enter__()
        self.index.add("H1", "Analysis of poly-G tracts")
        self.index.add("p", "Created: "+self.creationTime)
        self.index.addRule()
        return self

    def __exit__(self, exception_type, exception_val, trace) :
        self.index.__exit__(exception_type, exception_val, trace)

    def createStyler(self) :
        styler = html.Styler()
        styler.addStyle(['table tr {vertical-align:top}', 'table tr:nth-child(odd) td{ background-color:#f4f4f4 }','table tr:nth-child(even) td{ background-color:#dddddd }'])
        styler.addStyle(['table tr:target td{ transition:background-color 1s ease-in; -webkit-transition:background-color 1s ease-in; -moz-transition:background-color 1s ease-in; background-color: yellow; }'])
        return styler

    def loadStrainPickle(self, n) :
        with open(os.path.join(self.path, '../pickle', str(n)+'.pkl'), 'rb') as file :
            strain = pickle.load(file)
        return strain

    def createStrainFiles(self) :
        # Add to index
        print('Outputting strain files', end='', flush=True)
        allgenegroups = set()
        speciesdata = {}
        with self.index.table() as indexTable :
            indexTable.addHeaderRow(['Species', 'Strain', '#Contigs', '#PV Genes', '#Gene Groups', 'Incremental total #Gene groups'])
            # Create individual files
            for i in range(0, self.guide['numStrains']) :

                # Load the strain
                strain = self.loadStrainPickle(i)

                # Add row to index
                tracts = 0
                genegroups = set()
                for contig in strain['contigs'] :
                    tracts += len(contig['tracts'])
                    for tract in contig['tracts'] :
                        genegroups.add(tract['geneGroup'])
                allgenegroups |= genegroups

                if strain['species'] not in speciesdata :
                    speciesdata[strain['species']] = { 'genes' : 0, 'genegroups' : set(), 'numGroups' : 0, 'count' : 0 }

                speciesdata[strain['species']]['genes'] += tracts
                speciesdata[strain['species']]['genegroups'] |= genegroups
                speciesdata[strain['species']]['numGroups'] += len(genegroups)
                speciesdata[strain['species']]['count'] += 1

                indexTable.addRow([
                    strain['species'],
                    html.makeLink('strains/'+str(i)+'.html', strain['name']),
                    len(strain['contigs']),
                    tracts,
                    len(genegroups),
                    len(allgenegroups)])

                # Create the strain file for the individual strain
                logging.info('Creating strain file for '+strain['name'])
                try:
                    with html.Doc(self.path+'/strains/'+str(i)+'.html', title=strain['name'], styler=self.styler) as htmldoc :
                        htmldoc.add("H1", "Analysis of poly-G tracts in "+strain['name'])
                        htmldoc.addRule()
                        htmldoc.add("p",  'Strain name: ' + strain['name'] +'<br/>' +'Filepath: '+strain['path']+'<br />Species: '+strain['species'])
                        htmldoc.add("p", "Created: "+self.creationTime)
                        htmldoc.addRule()
                        with htmldoc.table() as table :
                            strainNames = self.guide['strainNames']
                            row = ['TractNo', 'Contig', 'Location', 'Tract', 'On length', 'Gene Group', 'Offset from gene', 'Gene', 'Length', 'Function']+strainNames
                            table.addHeaderRow(row)
                            line = 0
                            for contig in strain['contigs'] :
                                for tract in contig['tracts'] :
                                    if 'blastMatch' not in tract :
                                        matchingGenes = [""]*len(strainNames)
                                    else :
                                        if tract['blastMatch'] == 'No Data' :
                                            matchingGenes = [""]*len(strainNames)
                                        else :
                                            matchingGenes = self.listBlastMatches(tract['blastMatch'], strainNames)
                                    row = [str(line)] #[html.makeTag("a",str(line),attributes={'name':tract['name']})]
                                    row = row+[contig['record'].name, tract['location'], tract['tract'], tract['onLength']]
                                    row = row+[html.makeLink('../groups/'+str(tract['geneGroup'])+'.html', str(tract['geneGroup']))]
                                    row = row+[tract['outby'], tract['name'], str(len(tract['translation']) if 'translation' in tract and tract['translation'] else 0), tract['function']]
                                    row = row+matchingGenes
                                    table.addRow(row, attributes={'id':tract['name']})
                                    line += 1
                    print('.', end='', flush=True)
                except Exception :
                    logging.exception('Failed during creation of strain file for '+strain['name'])
                    print('!', end='', flush=True)
        print('')
        self.index.addRule()
        self.index.add('h2', 'Species data')
        with self.index.table() as speciesTable :
            speciesTable.addHeaderRow(['Species', 'Genomes', 'Total genes', 'per genome', 'Total groups', 'per genome', 'Distinct groups'])
            for name in sorted(speciesdata.keys()) :
                species = speciesdata[name]
                speciesTable.addRow(
                    [name,
                    str(species['count']),
                    str(species['genes']),
                    '{0:.1f}'.format((species['genes']/species['count']) + 0.05),
                    str(species['numGroups']),
                    '{0:.1f}'.format((species['numGroups']/species['count']) + 0.05),
                    str(len(species['genegroups']))])

        self.index.addRule()

    def createGroupIndex(self) :
        '''Create an html containing all the gene group match ups to the page'''
        print('Outputting group files', end='', flush=True)
        logging.info('Create group index')

        # Add the gene groups to the index
        self.index.add('p', html.makeLink('groups/index.html', 'Gene Groupings'))
        self.index.add('p', html.makeLink('groups/tracts.html', 'Tract length and ON/OFF by gene group'))
        self.index.addRule()

        # Matrix of gene entries to pass to tract/group lookup table
        onoffMatrix = []

        with html.Doc(
            self.path+'/groups/index.html',
            title='PV gene groups',
            styler=self.styler,
            sortTable = '../include/sorttable.js',
            includeInHeader=['addLoadEvent', 'colouredText']
        ) as htmldoc :
            htmldoc.add("H1", "PV Gene groups")
            htmldoc.addRule()
            htmldoc.add("p", "Created: "+self.creationTime)
            htmldoc.addRule()

            # Add canvas for gene group
            htmldoc.add('canvas', '', attributes={'id' : 'gg_graphic', 'width' : 800, 'height' : 800})
            htmldoc.addRule()

            # Worker function
            def geneAlreadyThere(locus, genes) :
                for gene in genes :
                    if locus == gene['locus'] :
                        return True
                return False

            # First we need to collect all the tracts together
            allTracts = []
            tractLookup = {}
            strainToSpecies = {}
            header = ['#Group', 'Name', 'Likely Function', '#PV in gene', '#Total PV', '#Total Genes']
            for i in range(0, self.guide['numStrains']) :
                strain = self.loadStrainPickle(i)
                strainToSpecies[strain['name']] = strain['species']
                header += [strain['species'] + '<br/>' + html.makeLink('../strains/'+str(i)+'.html',strain['name'])]
                for contig in strain['contigs'] :
                    for n, tract in enumerate(contig['tracts']) :
                        tract['contig'] = contig['record'].description
                        allTracts.append((strain['name'], tract))
                        tractLookup[(tract['contig'], n)] = tract

            # If including species
            #header += ['In species']

            geneCounts = []
            totalHeight = 0
            with htmldoc.table() as table :
                table.addHeaderRow(header)
                # Then go through gene group by gene group collecting all relevant genes from each strain
                addedGroups = set()
                hasNonPv = []

                # The hacky None, None tailed on the end allows us to defer the -1 ungrouped genes to the end
                for _, tract in allTracts+[(None,None)] :
                    if tract :
                        group = tract['geneGroup']
                        if group in addedGroups or group == -1 :
                            continue

                        addedGroups.add(group)
                        onoffMatrix.append([[] for i in range(0, self.guide['numStrains'])])
                        hasNonPv.append(False)
                    else :
                        group = -1


                    foundInSpecies = set()
                    genesInStrain = {}
                    tractsInGroup = []
                    for inStrain, tract in allTracts :
                        if tract['geneGroup'] != group :
                            continue

                        # Track species this is in
                        foundInSpecies.add(strainToSpecies[inStrain])

                        # For the gene group creation later
                        tract['inStrainNo'] = self.getStrainNumber(inStrain)
                        tractsInGroup.append(tract)

                        # Add this gene and any matching non-PV genes to the strain
                        if not inStrain in genesInStrain :
                            genesInStrain[inStrain] = []
                        genesInStrain[inStrain].append(
                            {
                                'locus' : tract['locus'],
                                'subgroup' : tract['subgroup'] if 'subgroup' in tract else None,
                                'linkTo' : '../strains/'+str(self.getStrainNumber(inStrain))+'.html#'+tract['locus'],
                                'attributes' : {'style' : 'background-color:' + ('#00b000' if tract['outby'] == 0 else '#ffa000')},
                                'repeat' : getRepeat(tract['tract']),
                                'tract' : tract['tract'],
                                'on' : tract['onLength'],
                                'sort' : 0 if tract['outby'] == 0 else 1,
                            })
                        if 'blastMatch' in tract :
                            for contigName, matches in tract['blastMatch'].items() :
                                strainName = self.findStrainNameForContig(contigName)
                                for match in matches :
                                    if 'tractNo' in match :
                                        # Create tract location data for the gene alignment graphics created later
                                        matched = tractLookup[(contigName, match['tractNo'])]
                                        match['repeat'] = {
                                            'location' : matched['location'],
                                            'start' : matched['start'],
                                            'end' : matched['end'],
                                            'repeats' : matched['repeats'],
                                            'repeatUnit' : matched['repeatUnit'],
                                            'strand' : matched['strand'] if 'strand' in matched else None }
                                        continue
                                    if strainName in genesInStrain and geneAlreadyThere(match['locus'], genesInStrain[strainName]) :
                                        continue
                                    if strainName not in genesInStrain :
                                        genesInStrain[strainName] = []
                                    genesInStrain[strainName].append({'locus' : match['locus'], 'attributes' : None, 'repeat' : '', 'sort' : 2})
                                    if group != -1 :
                                        hasNonPv[group] = True

                    if group != -1 :
                        groupName = str(group)
                    else :
                        groupName = 'Not in group'

                    #Create the row and groupRow lists
                    row = [html.makeLink(groupName+'.html', groupName)]
                    row += [self.geneGroups[group]['name'], self.geneGroups[group]['bestFunction']] if group != -1 else ['No data', 'No data']
                    row += ['0', '0', '0']
                    rowHead = len(row)
                    row += (['']*self.guide['numStrains'])
                    groupRow = [groupName, self.geneGroups[group]['name'] if group != -1 else 'N/A']  + [self.geneGroups[group]['bestFunction'] if group != -1 else 'No data'] + ['0', '0', '0'] + (['']*self.guide['numStrains'])

                    geneCount = [[[],[],[]] for i in range(0, self.guide['numStrains'])]
                    highestCount = 0
                    for strain, genes in genesInStrain.items() :
                        # Add the grouped genes for this strain to the row
                        x = self.getStrainNumber(strain)
                        # Sort by sort value then alphabetically
                        genes.sort(key=lambda k: (k['sort'], k['locus'].lower()))
                        for gene in genes :
                            if 'linkTo' in gene :
                                subgroup = ''
                                if gene['subgroup'] is not None :
                                    subgroup = '{'+str(gene['subgroup'])+'}'
                                row[x+rowHead] += html.makeTag('div', html.makeLink(gene['linkTo'], shortenName(gene['locus'])+subgroup), attributes=gene['attributes'])
                                groupRow[x+rowHead] += html.makeTag('div', html.makeLink('#'+gene['locus'], shortenName(gene['locus'])), attributes=gene['attributes'])
                                if group != -1 and gene['sort'] != 2:
                                    onoffMatrix[group][x].append(gene)
                            else :
                                row[x+rowHead] += shortenName(gene['locus'])+'<br />'
                                groupRow[x+rowHead] += shortenName(gene['locus'])+'<br />'
                            geneCount[x][gene['sort']].append(gene['repeat'])
                        count = sum(len(c) for c in geneCount[x])
                        if count > highestCount :
                            highestCount = count

                    # Count the types of gene found and add them to the row
                    groupRow[3] = row[3] = str( sum(len(geneCount[x][0]) for x in range(0, len(geneCount))) )
                    groupRow[4] = row[4] = str( sum(len(geneCount[x][i]) for x in range(0, len(geneCount)) for i in [0, 1]) )
                    groupRow[5] = row[5] = str( sum(len(geneCount[x][i]) for x in range(0, len(geneCount)) for i in [0, 1, 2]) )

                    totalHeight += highestCount+1

                    # Add species to the end
                    # speciesMap = {
                    #     'C. coli' : 'Ccl',
                    #     'C. concisus' : 'Ccn',
                    #     'C. curvus' : 'Ccu',
                    #     'C. fetus' : 'Cf',
                    #     'C. gracilus' : 'Cg',
                    #     'C. hominis' : 'Cho',
                    #     'C. hyointestinalis' : 'Chy',
                    #     'C. iquanorium' : 'Cig',
                    #     'C. insulaenigrae' : 'Cin',
                    #     'C. jejuni' : 'Cj',
                    #     'C. lari' : 'Cl',
                    #     'C. peloridis' : 'Cp',
                    #     'C. sp.' : 'Csp',
                    #     'C. subantarcticus' : 'Csu',
                    #     'C. ureolyticus' : 'Cu',
                    #     'C. volucris' : 'Cv' }
                    # speciesOrder = sorted(speciesMap)
                    # speciesList = []
                    # for species in speciesOrder :
                    #     if species in foundInSpecies :
                    #         speciesList.append(speciesMap[species])
                    #
                    # row += [' '.join(speciesList)]

                    # Don't add the ungrouped to the graphic
                    if group != -1 :
                        geneCounts.append(geneCount)
                        table.addRow(row)
                    else:
                        # The Not in group always stays at the bottom
                        table.addFooterRow(row)

                    self.createGroupFile(group, groupName, tractsInGroup, header, groupRow)


                    print('.', end='', flush=True)

            print('')

            htmldoc.addRule()

            # Add the javascript to update the genegroup canvas
            self.addGeneGroupGraphic(htmldoc, geneCounts, totalHeight, 'gg_graphic')

            htmldoc.add('p', '{0} of {1} gene groups have non-PV homologues ({2:.2f}%)'.format(sum(hasNonPv), len(hasNonPv), (sum(hasNonPv)/len(hasNonPv))*100))
            htmldoc.addRule()

            strainMap = {}
            strainOnMap = {}
            for strainName, tract in allTracts :
                if strainName not in strainMap :
                    strainMap[strainName] = [0]*len(addedGroups)
                    strainOnMap[strainName] = [0]*len(addedGroups)
                strainMap[strainName][tract['geneGroup']] = 1
                # Ignore out of frame tracts for the On map
                if tract['outby'] == 0 :
                    if strainOnMap[strainName][tract['geneGroup']] == '?' :
                        strainOnMap[strainName][tract['geneGroup']] = 0
                    strainOnMap[strainName][tract['geneGroup']] += (tract['repeats'] == tract['onLength'])

            # Count occurences of each group
            groupCounts = [0]*len(addedGroups)
            for groups in strainMap.values() :
                for i, present in enumerate(groups) :
                    groupCounts[i] += present

            # Categorise rarity of gene groups
            groupClass = ['']*len(addedGroups)
            classes = ['Majority (>50%)', 'Common (>15%)', 'Uncommon (>5%)', 'Rare (>2%)', 'Very rare']
            for i, count in enumerate(groupCounts) :
                prop = count/self.guide['numStrains']
                if prop >= 0.5 :
                    groupClass[i] = classes[0]
                elif prop > 0.15 :
                    groupClass[i] = classes[1]
                elif prop > 0.05 :
                    groupClass[i] = classes[2]
                elif prop > 0.02 :
                    groupClass[i] = classes[3]
                else:
                    groupClass[i] = classes[4]

            htmldoc.add('p', 'Gene group rairities by gene')
            with htmldoc.table() as table :
                table.addHeaderRow(['Species', 'Strain']+classes)
                iterTracts = iter(allTracts)
                onTract = next(iterTracts)
                totals = Counter()
                for strainName, species in zip(self.guide['strainNames'], self.guide['species']) :
                    count = Counter()
                    while onTract and onTract[0] == strainName :
                        count[groupClass[onTract[1]['geneGroup']]] += 1
                        onTract = next(iterTracts, None)
                    totals+=count
                    table.addRow([species, strainName] + [count[c] for c in classes])
                table.addFooterRow(['', 'Totals'] + [totals[c] for c in classes])

            htmldoc.add('p', 'Gene Group rarities by gene group')
            with htmldoc.table() as table :
                table.addHeaderRow(['Species', 'Strain']+classes)
                iterTracts = iter(allTracts)
                onTract = next(iterTracts)
                totals = Counter()
                for strainName, species in zip(self.guide['strainNames'], self.guide['species']) :
                    geneGroups = set()
                    while onTract and onTract[0] == strainName :
                        geneGroups.add(onTract[1]['geneGroup'])
                        onTract = next(iterTracts, None)
                    count = Counter()
                    for group in geneGroups :
                        count[groupClass[group]] += 1
                    table.addRow([species, strainName] + [count[c] for c in classes])
                count = Counter()
                for group in groupClass :
                    count[group] += 1
                table.addFooterRow(['Genegroups', 'by group'] + [count[c] for c in classes])

            htmldoc.add('p', 'Gene Group rarities by gene group and species')
            with htmldoc.table() as table :
                table.addHeaderRow(['Species']+classes)
                iterTracts = iter(allTracts)
                onTract = next(iterTracts)
                totals = Counter()
                for strainName, species in zip(self.guide['strainNames'], self.guide['species']) :
                    geneGroups = set()
                    while onTract and onTract[0] == strainName :
                        geneGroups.add(onTract[1]['geneGroup'])
                        onTract = next(iterTracts, None)
                    count = Counter()
                    for group in geneGroups :
                        count[groupClass[group]] += 1
                table.addRow([species, strainName] + [count[c] for c in classes])
                count = Counter()
                for group in groupClass :
                    count[group] += 1
                table.addFooterRow(['By group'] + [count[c] for c in classes])

            with htmldoc.tag('p') :
                htmldoc.addRaw('Majority: in >= 50% of strains<br/>')
                htmldoc.addRaw('Common: in > 15% and < 50% of strains<br/>')
                htmldoc.addRaw('Uncommon: in > 5% and <= 15% of strains<br/>')
                htmldoc.addRaw('Rare: in > 2% and <= 5% of strains<br/>')
                htmldoc.addRaw('Very rare: in <= 2% of strains<br/>')

            htmldoc.addRule()

            if self.guide['numStrains'] > 2 :
                print('...Graphing')

                def calcDistance(a,b) :
                    return sum(x!=y for x,y in zip(a,b))

                def calcDistanceAmbiguous(a,b) :
                    return sum(0 if x=='?' or y=='?' else abs(x-y) for x, y in zip(a,b))

                # Add only entries which have any tracts
                names = [n  for n in self.getStrainNames() if n in strainMap]
                distances = [[calcDistance(strainMap[a],strainMap[b]) for b in names[:n+1]] for n,a in enumerate(names)]
                numbers=[i for i in range(0, len(names))]
                htmldoc.add('p', html.makeTag('b', 'Neighbour joining tree'))
                self.runPlugins('preTree', htmldoc, 'nj')
                self.addPhylogenicTree(htmldoc, 'nj', names=names, numbers=numbers, matrix=distances)
                self.runPlugins('postTree', htmldoc, names, numbers, 'nj')

                htmldoc.addRule()

                htmldoc.add('p', html.makeTag('b', 'On/Off tree'))

                names = [n  for n in self.getStrainNames() if n in strainMap]
                distances = [[calcDistanceAmbiguous(strainOnMap[a],strainOnMap[b])*10 for b in names[:n+1]] for n,a in enumerate(names)]
                numbers=[i for i in range(0, len(names))]
                self.runPlugins('preTree', htmldoc, 'njonoff')
                self.addPhylogenicTree(htmldoc, 'nj', names=names, numbers=numbers, matrix=distances, tag='onoff')
                self.runPlugins('postTree', htmldoc, names, numbers, 'njonoff')

                htmldoc.addRule()

                subsets = []
                self.runPlugins('filteredTrees', subsets)
                if subsets :
                    htmldoc.add('p', html.makeTag('p', 'Filtered trees (nj)'))
                    for n, subset in enumerate(subsets) :
                        htmldoc.add('H3', subset['name'])
                        subnames = [names[i] for i in subset['numbers']]
                        distances = [[calcDistance(strainMap[a],strainMap[b]) for b in subnames[:n+1]] for n,a in enumerate(subnames)]
                        numbers=subset['numbers']
                        htmldoc.add('p', html.makeTag('b', 'Neighbour joining tree'))
                        self.runPlugins('preTree', htmldoc, 'nj'+str(n))
                        self.addPhylogenicTree(htmldoc, 'nj', names=subnames, numbers=numbers, matrix=distances, tag=str(n))
                        self.runPlugins('postTree', htmldoc, subnames, numbers, 'nj'+str(n))

                        distances = [[calcDistanceAmbiguous(strainOnMap[a],strainOnMap[b]) for b in subnames[:n+1]] for n,a in enumerate(subnames)]
                        numbers=subset['numbers']
                        htmldoc.add('p', html.makeTag('b', 'On/Off Neighbour joining tree'))
                        self.runPlugins('preTree', htmldoc, 'njonoff'+str(n))
                        self.addPhylogenicTree(htmldoc, 'nj', names=subnames, numbers=numbers, matrix=distances, tag='onoff'+str(n))
                        self.runPlugins('postTree', htmldoc, subnames, numbers, 'njonoff'+str(n))

                        htmldoc.addRule()

        self.createOnOffByGroup(onoffMatrix)

    def createOnOffByGroup(self, onoffMatrix) :
        ''' Create an tract length data by gene group'''
        with html.Doc(
            self.path+'/groups/tracts.html',
            title='Tract length by PV gene groups',
            styler=self.styler,
            sortTable = '../include/sorttable.js',
            includeInHeader=['addLoadEvent']
        ) as htmldoc :
            htmldoc.add("H1", "Tract lengths and ON/OFF by PV Gene groups")
            htmldoc.addRule()
            htmldoc.add("p", "Created: "+self.creationTime)
            htmldoc.addRule()

            colour_on = '#2288aa'
            colour_off = '#ee0000'
            htmldoc.add('H3', 'Predicted <span style="background-color:'+colour_on+';">On</span>/<span style="background-color:'+colour_off+';">Off</span> state')

            header = ['#Group', 'Likely Function']
            for i in range(0, self.guide['numStrains']) :
                strain = self.loadStrainPickle(i)
                header += [html.makeLink('../strains/'+str(i)+'.html',strain['name'])]

            with htmldoc.table() as table :
                table.addHeaderRow(header)
                for n, group in enumerate(onoffMatrix) :
                    row = [str(n), self.geneGroups[n]['bestFunction']] + ['']*len(group)
                    for x, genes in enumerate(group) :
                        for gene in genes :
                            if gene['on'] :
                                if getRepeats(gene['tract']) == gene['on'] :
                                    col = colour_on
                                else :
                                    col = colour_off
                            else:
                                col = '#ffffff'
                            row[x+2] += html.makeTag('div', html.makeLink(gene['linkTo'], gene['tract']), attributes={'style' : 'background-color:' + col +';'})
                    table.addRow(row)

            htmldoc.addRule()

    def addGeneGroupGraphic(self, htmldoc, geneCounts, totalHeight, targetCanvas) :
        ''' Add a graphic showing the gene groups to the htmldoc '''

        with htmldoc.tag('script', attributes={'type' : 'text/javascript'}) :
            htmldoc.addRaw('function drawGeneGroups() {')
            htmldoc.addRaw('  var canvas_element = document.getElementById("'+targetCanvas+'");')
            htmldoc.addRaw('  var canvas = canvas_element.getContext("2d");')
            htmldoc.addRaw('  var line = 0;')

            longest_name = max([len(name) for name in self.guide['displayNames']])

            htmldoc.addRaw('  var base_y = 20+10*'+str(longest_name)+';')
            htmldoc.addRaw('  var base_x = 20+10*'+str(len(str(len(geneCounts)-1)))+';')
            
            # Put the data on gene groups into the html
            htmldoc.addRaw('  var data = '+str(geneCounts)+';')

            
            htmldoc.addRaw('  var height = '+str(len(geneCounts))+';')
            htmldoc.addRaw('  var width = '+str(len(geneCounts[0]))+';')
            htmldoc.addRaw('  var col = ["#00b000", "#ffa000", "#aaaaaa"];')
            htmldoc.addRaw('  var baseCols= {"A" : "#4444dd", "C" : "#dd4444", "G" : "#44dd44", "T" : "#dddd44"};')
            htmldoc.addRaw('  var w=20;')
            htmldoc.addRaw('  var h=5;')

            # Set the canvas element to the correct size
            htmldoc.addRaw('  canvas_element.width = base_x+width*(w+1);')
            htmldoc.addRaw('  canvas_element.height = base_y+'+str(totalHeight)+'*h;')
            
            # Font for labels
            htmldoc.addRaw('  canvas.font="16px Courier New";')
            
            htmldoc.addRaw('  for(y = 0; y < height; y++) {')
            htmldoc.addRaw('    var highOff = 0;')
            # add label to the row
            htmldoc.addRaw('    if ((y%5)==0) {')
            htmldoc.addRaw('      canvas.fillStyle = "#000000";')
            htmldoc.addRaw('      canvas.textAlign="right";')
            htmldoc.addRaw('      canvas.fillText(y, base_x-10, base_y+line*h+6);')
            htmldoc.addRaw('    }')
            # draw the coloured bars
            htmldoc.addRaw('    for(x = 0; x < width; x++) {')
            htmldoc.addRaw('      var off = 0;')
            htmldoc.addRaw('      for(type = 0; type < 3; type++) {')
            htmldoc.addRaw('        if (data[y][x][type].length === 0) { continue; }')
            htmldoc.addRaw('        canvas.fillStyle = col[type];')
            htmldoc.addRaw('        canvas.fillRect(base_x+x*(w+1), base_y+(line+off)*h, w, h*data[y][x][type].length);')
            htmldoc.addRaw('        for(gene=0; gene < data[y][x][type].length; gene++) {')
            htmldoc.addRaw('          var repeat = data[y][x][type][gene];')
            htmldoc.addRaw('          for(onLetter=0; onLetter < repeat.length; onLetter++) {')
            htmldoc.addRaw('            var letter = repeat[onLetter];')
            htmldoc.addRaw('            canvas.fillStyle = baseCols[letter];')
            htmldoc.addRaw('            canvas.fillRect(base_x+x*(w+1)+onLetter*3, base_y+(line+off+gene)*h, 3, h);')
            htmldoc.addRaw('          }')
            htmldoc.addRaw('        }')
            htmldoc.addRaw('        off += data[y][x][type].length;')
            htmldoc.addRaw('      }')
            # Add a very, very faint carrier line
            htmldoc.addRaw('      if (off === 0) {')
            htmldoc.addRaw('        canvas.fillStyle = "#eeeecc";')
            htmldoc.addRaw('        canvas.fillRect(base_x+x*(w+1), base_y+(line)*h, w, h);')
            htmldoc.addRaw('      }')
            htmldoc.addRaw('      if (off > highOff) { highOff = off; }')
            htmldoc.addRaw('    }')
            htmldoc.addRaw('    line += highOff+1;')
            htmldoc.addRaw('  }')

            # Add strain name labels above the data
            htmldoc.addRaw('  var strainNames = '+str(self.guide['displayNames'])+';')
            htmldoc.addRaw('  canvas.save();')
            # Astonishingly we have to rotate the ENTIRE CANVAS to rotate the text!
            htmldoc.addRaw('  canvas.fillStyle = "#000000";')
            htmldoc.addRaw('  canvas.rotate(Math.PI/2);')
            htmldoc.addRaw('  canvas.textAlign="right";')
            htmldoc.addRaw('  for(x = 0; x < width; ++x) {')
            htmldoc.addRaw('    canvas.fillText(strainNames[x], base_y-10, -(x*(w+1)+2)-base_x);')
            htmldoc.addRaw('  }')
            htmldoc.addRaw('  canvas.restore();')
            htmldoc.addRaw('}')

            htmldoc.addRaw('window.onload = drawGeneGroups();')



    def createGroupFile(self, group, groupName, tractsInGroup, headerRow, row) :
        ''' Create an html file containing all the blast hits for a group '''
        strainNames = self.getStrainNames()

        with html.Doc(self.path+'/groups/'+groupName+'.html',
                      title='Gene group '+groupName,
                      styler=self.styler,
                      includeInHeader=['addLoadEvent', 'colouredText']) as htmldoc :
            htmldoc.add("H1", "Gene group "+groupName)
            htmldoc.add('a', '[Back to Gene Group Index]', attributes = { 'href' : 'index.html' })
            htmldoc.add("p", "Created: "+self.creationTime)
            htmldoc.addRule()

            if not tractsInGroup :
                htmldoc.add('p', 'Empty group')
                return

            # Header row with all the tracts
            with htmldoc.table() as table :
                table.addHeaderRow(headerRow)
                table.addRow(row)
            htmldoc.addRule()

            # Annotated functions
            if group != -1 :
                annotations = self.geneGroups[group]['words']
                functions = {}
                for annotation in annotations :
                    func = annotation['function'].lower()
                    if func not in functions :
                        functions[func] = {
                            'function' : annotation['function'],
                            'count' : 1,
                            'score' : annotation['score']
                        }
                    else :
                        functions[func]['count'] += 1

                funcOrdered = [(v['score'], k) for k, v in functions.items()]
                funcOrdered.sort(reverse=True)

                # Annotated functions of these genes
                htmldoc.add("H2", "Annotated functions")
                with htmldoc.table() as table :
                    table.addHeaderRow(['Function', '#Occuring', 'score'])
                    for f in funcOrdered :
                        func = functions[f[1]]
                        table.addRow([func['function'], func['count'], '{0:.2f}'.format(func['score'])])
                htmldoc.addRule()

            # (Subgroups are all or nothing so we can test the first tract)
            if 'subgroup' in tractsInGroup[0] :
                # Find subgroups and build a subgroup binary for each strain that appears at all
                subgroups = {}
                strainMap = {}
                for tract in tractsInGroup :
                    subgroup = tract['subgroup']
                    if subgroup not in subgroups :
                        subgroups[subgroup] = []
                    subgroups[subgroup].append(tract)
                    if tract['inStrainNo'] not in strainMap :
                        strainMap[tract['inStrainNo']] = set()
                    strainMap[tract['inStrainNo']].add(subgroup)

                subgroups = [subgroups[n] for n in sorted(subgroups)]

                # Binary map
                for n in strainMap :
                    present = [0]*len(subgroups)
                    for s in strainMap[n] :
                        present[s] = 1
                    strainMap[n] = present

                htmldoc.add('H2', 'Subgroups')
                subgroupHeader = ['Subgroup']+strainNames
                with htmldoc.table() as table :
                    table.addHeaderRow(subgroupHeader)
                    for n, subgroup in enumerate(subgroups) :
                        row = [str(n)]+['']*len(strainNames)
                        for tract in subgroup :
                            attributes = {'style' : 'background-color:' + ('#00b000' if tract['outby'] == 0 else '#ffa000')}
                            row[tract['inStrainNo']+1] += html.makeTag('div', html.makeLink('#'+tract['locus'], shortenName(tract['locus'])), attributes = attributes)
                        table.addRow(row)

                # On/Off scoring by subgroup
                colour_on = '#2288aa'
                colour_off = '#ee0000'
                htmldoc.add('H3', 'Predicted <span style="background-color:'+colour_on+';">On</span>/<span style="background-color:'+colour_off+';">Off</span> state')
                subgroupHeader = ['Subgroup']+strainNames
                with htmldoc.table() as table :
                    table.addHeaderRow(subgroupHeader)
                    for n, subgroup in enumerate(subgroups) :
                        row = [str(n)]+['']*len(strainNames)
                        for tract in subgroup :
                            if tract['outby'] != 0 :
                                colour = '#aaaaaa'
                            else :
                                if tract['repeats'] == tract['onLength'] :
                                    colour = colour_on
                                else :
                                    colour = colour_off
                            attributes = {'style' : 'background-color:' + colour}
                            row[tract['inStrainNo']+1] += html.makeTag('div', html.makeLink('#'+tract['locus'], shortenName(tract['locus'])), attributes = attributes)
                        table.addRow(row)
                htmldoc.addRule()

                # graph of similarity by subgroup
                def calcDistance(a,b) :
                    return sum(x!=y for x,y in zip(a,b))

                if len(subgroups) > 2 :
                    nums = sorted(strainMap)
                    distances = [[calcDistance(strainMap[a],strainMap[b]) for b in nums[:n+1]] for n,a in enumerate(nums)]
                    strainNames = self.getStrainNames()
                    names = [strainNames[n] for n in nums]
                    htmldoc.add('p', html.makeTag('b', 'Neighbour joining tree'))
                    self.runPlugins('preTree', htmldoc, 'nj')
                    self.addPhylogenicTree(htmldoc, 'nj', names, nums, distances)
                    self.runPlugins('postTree', htmldoc, names, nums, 'nj')

                    htmldoc.add('p', html.makeTag('b', 'UPGMA tree'))
                    self.runPlugins('preTree', htmldoc, 'upgma')
                    self.addPhylogenicTree(htmldoc, 'upgma', names, nums, distances)
                    self.runPlugins('postTree', htmldoc, names, nums, 'upgma')

                    htmldoc.addRule()

            onStrain = -1
            for nth, tract in enumerate(tractsInGroup) :
                if tract['inStrainNo'] != onStrain :
                    onStrain = tract['inStrainNo']
                    htmldoc.add("H2", strainNames[onStrain])
                htmldoc.add("p",
                    html.makeTag("b", strainNames[tract['inStrainNo']] + ' : ' + tract['locus']) +
                        ' ' +html.makeLink('../strains/'+str(tract['inStrainNo'])+'.html#'+tract['locus'], '[tract entry]') +
                        ' ' +html.makeLink('#top', '[top]'),
                    attributes = {'id' : tract['locus']})
                translation = tract['translation']
                tractStart, tractLength = self.getTractAminoPos(tract)
                with htmldoc.tag('p') :
                    htmldoc.add('div', 'Function: '+tract['function'], attributes={'style':'margin-left: 50px'})
                    with htmldoc.tag('div', attributes={'style':'font-family: monospace; white-space: pre; margin-left: 50px'}) :
                        step = 100
                        if translation :
                            for i in range(0, len(translation), step) :
                                if i <= tractStart < i+step :
                                    htmldoc.addRaw(' '*10+'+'+str(i).ljust(5))
                                    offset = tractStart-i
                                    htmldoc.addRaw(str(translation[i:i+offset]))
                                    htmldoc.addRaw('<span style="background-color:#0033dd;color:white;">'+str(translation[i+offset:min(i+step, i+offset+tractLength)])+'</span>')
                                    htmldoc.addRaw(str(translation[min(i+step, i+offset+tractLength):i+step])+'<br />')
                                    if offset + tractLength >= i+step :
                                        tractLength -= (i+step)-offset
                                        tractStart = i+step
                                else :
                                    htmldoc.addRaw(' '*10+'+'+str(i).ljust(5)+str(translation[i:i+step])+'<br />')
                if 'blastMatch' in tract :
                    showGraphicsByDefault = len(tractsInGroup) < 30
                    if showGraphicsByDefault :
                        self.addAlignmentGraphic(htmldoc, nth, tract)
                    htmldoc.addExpandable(groupName+'_'+str(nth)+'.html', '[show detailed alignment]')
                    with html.Doc(self.path+'/groups/'+groupName+'_'+str(nth)+'.html', styler=self.styler, includeInHeader=['addLoadEvent']) as section :
                        if not showGraphicsByDefault :
                            self.addAlignmentGraphic(section, nth, tract)
                        matchOrdering = [(self.getStrainNumber(self.findStrainNameForContig(contigName)), contigName) for contigName in tract['blastMatch']]
                        matchOrdering.sort(key=lambda x : x[0])
                        for _, contigName in matchOrdering :
                            matches = tract['blastMatch'][contigName]
                            strainName = self.findStrainNameForContig(contigName)
                            for match in matches :
                                if contigName == tract['contig'] and match['locus'] == tract['locus'] :
                                    continue
                                section.addRaw(tract['locus']+' vs: '+match['locus']+' in '+html.makeLink('../strains/'+str(self.getStrainNumber(strainName))+'.html', strainName)+ ' ('+contigName+')<br/>')
                                self.addAlignment(section, tract, match)
                else:
                    htmldoc.add('p', 'No BLAST hits')
                htmldoc.addRule()

    def addAlignmentGraphic(self, htmldoc, nth, tract) :
        ''' Add a graphic showing the alignment between this tract and its blastmatches'''
        htmldoc.add('p',html.makeTag('canvas', '', attributes={'id' : 'align_graphic_'+str(nth), 'width' : 800, 'height' : 800}))

        def h(x) :
            htmldoc.addRaw(x)

        with htmldoc.tag('script', attributes={'type' : 'text/javascript'}) :
            # Get the canvas
            h('function drawAlignment'+str(nth)+'() {')
            h('  var canvas_element = document.getElementById("align_graphic_'+str(nth)+'");')
            h('  var canvas = canvas_element.getContext("2d");')

            # Font for labels
            h('  canvas.font="10px Courier New";')

            # Calculate the canvas size needed
            longest = 0
            height = 0
            for contigName, matches in tract['blastMatch'].items() :
                for match in matches :
                    if contigName == tract['contig'] and match['locus'] == tract['locus'] :
                        continue
                    nameLength = len(self.findStrainNameForContig(contigName))+5+len(shortenName(match['locus']))
                    if nameLength > longest :
                        longest = nameLength
                    if 'multipleHits' in match :
                        height += len(match['multipleHits'])*4
                    else :
                        height += 4
                    height += 8

            textLength = longest*10
            h('  var guide_y = 20;')
            elementHeight = 20+4+height
            elementWidth = tract['geneLength'] if tract['geneLength'] < 975 else 900
            repeats = 1 if tract['geneLength'] < 975 else ((tract['geneLength']//3) + 299) // 300
            h('  canvas_element.width = 2*'+str(textLength)+'+20+'+str(elementWidth)+';')
            h('  canvas_element.height = '+str(elementHeight*repeats + 20*(repeats-1))+';')

            # Draw the guideline
            h('  var base_x = 10+'+str(textLength)+';')
            h('  canvas.fillStyle = "#000000";')
            h('  for (y=0; y<'+str(repeats)+';y++) {')
            h('    var startAt=y*300;')
            h('    var yOffset = y*'+str(elementHeight+20)+';')
            if repeats > 1 :
                h('    var endAt=Math.min((y+1)*300, '+str(tract['geneLength']//3)+');')
            else :
                h('    var endAt = '+str(tract['geneLength']//3)+';')
            h('    canvas.fillRect(base_x, guide_y+yOffset, (endAt-startAt-1)*3, 2);')
            h('    canvas.textAlign="center";')
            h('    for(i = 0; i < endAt-startAt; i+=50) {')
            h('      canvas.fillRect(base_x+i*3, yOffset+guide_y-4, 3, 4);')
            h('      canvas.fillText(startAt+i, base_x+i*3, yOffset+guide_y-6);')
            h('    }')
            h('  }')

            # Work out where our repeat should go
            if tract['outby'] == 0 :
                tractStart, tractLength = self.getTractAminoPos(tract)
                h('  canvas.fillStyle = "#0033DD";')
                if repeats > 1 :
                    h('  var y = '+str(tractStart//300)+'*'+str(elementHeight+20)+' - 6;')
                else :
                    h('  var y = -6;')
                h('  var x = base_x + '+str(tractStart)+'*3;')
                if repeats > 1 and (tractStart // 300) != ((tractStart+tractLength) // 300) :
                    h('  canvas.fillRect(x, guide_y+y, ('+str(300*((tractStart // 300)+1) - tractStart)+'-1)*3, 4);')
                    h('  canvas.fillRect(base_x, guide_y+y+'+str(elementHeight+20)+', '+str(tractLength+tractStart - 300*((tractStart // 300)+1))+'*3, 4);')
                else :
                    h('  canvas.fillRect(x, guide_y+y, '+str(tractLength)+'*3, 4);')


            # Create the alignment map
            h('  var alignmentMap = [')
            sortedNames = [(self.getStrainNumber(self.findStrainNameForContig(contigName)), contigName) for contigName in tract['blastMatch']]
            sortedNames.sort()
            for _, contigName in sortedNames :
                matches = tract['blastMatch'][contigName]
                for match in matches :
                    if contigName == tract['contig'] and match['locus'] == tract['locus'] :
                        continue
                    h('    ["'+self.findStrainNameForContig(contigName)+' ( ' + shortenName(match['locus'])+' )", "'+match['locus']+'",')
                    if 'tractNo' in match :
                        tractStart, tractLength = self.getTractAminoPos(match['repeat'])
                        if repeats > 1 and (tractStart // 300) != ((tractStart+tractLength) // 300) :
                            h('    [['+str(tractStart)+','+str(300*((tractStart // 300)+1) - tractStart)+'-1], ['+str(300*((tractStart // 300)+1))+','+str(tractLength+tractStart - 300*((tractStart // 300)+1))+']], [')
                        else :
                            h('    [['+str(tractStart)+','+str(tractLength)+']], [')
                    else :
                        h('    [] ,[')
                    if 'multipleHits' in match :
                        hsps = match['multipleHits']
                    else :
                        hsps = [match['hsp']]
                    for hsp in hsps :
                        h('      '+str(self.createAlignmentMap(hsp))+',')
                    h('    ]],')
            h('  ];')

            # Convert the map into a pretty picture

            h('  canvas.textBaseline="middle";')
            h('  var y = guide_y+6;')
            h('  var h = 4;')
            h('  var col = ["#00bb00", "#ff9900", "#cccccc"];')
            h('  for(i=0; i < alignmentMap.length; i++ ) {')
            h('    var match = alignmentMap[i];')
            h('    canvas.fillStyle = "#000000";')

            if repeats > 1 :
                h('    for (repeat = 0; repeat < '+str(repeats)+';repeat++) {')
                h('      yOffset = repeat*'+str(elementHeight+20)+';')
                h('      var startAt=repeat*300;')
                h('      var endAt=Math.min((repeat+1)*300, '+str(tract['geneLength']//3)+');')
            else :
                h('      yOffset = 0;')
                h('      startAt = 0;')
                h('      var endAt = '+str(tract['geneLength']//3)+';')
            h('      canvas.font="10px Courier New";')
            h('      canvas.textAlign="right";')
            h('      canvas.fillText(match[0], base_x-10, y+2+yOffset);')
            h('      canvas.textAlign="left";')
            h('      canvas.fillText(match[1], base_x+(endAt-startAt)*3+10, y+2+yOffset);')

            if repeats > 1 :
                h('    }')

            h('    for(j=0; j< match[3].length; j++) {')
            h('      var submatch = match[3][j];')
            h('      var x = 0;')
            h('      var nextBreak = 299;')
            h('      var yOffset = 0;')
            h('      for(k=0; k < submatch.length; k++) {')

            if repeats > 1 :
                h('        var baseEntry = submatch[k];')
                h('        if (baseEntry[0] != 3 && x+baseEntry[1] > nextBreak) {')
                # Split the entry into several
                h('          var len = baseEntry[1];')
                h('          var entries = [];')
                h('          var next = nextBreak;')
                h('          var onX = x;')
                h('          while (onX+len > next) { entries.push([baseEntry[0], next-onX]); len -= next-onX; onX=next; next += 300; }')
                h('          entries.push([baseEntry[0], len]);')
                h('        } else {')
                h('          var entries = [baseEntry];')
                h('        }')
            else :
                h('        var entries = [submatch[k]];')

            h('        for(l=0; l < entries.length; l++) {')
            h('          entry = entries[l];')
            h('          if (entry[0] < 3) {')
            h('            canvas.fillStyle = col[entry[0]];')
            h('            canvas.fillRect(base_x+x*3, yOffset+y, entry[1]*3, h);')
            h('            x+=entry[1];')
            h('          } else if (entry[0] == 3) {')
            h('            canvas.fillStyle = "#000000";')
            h('            canvas.fillRect(base_x+x*3-1, yOffset+y-2, 2, h+4);')
            h('          } else if (entry[0] == 4) {')
            h('            canvas.fillStyle = "#000000";')
            h('            canvas.fillRect(base_x+x*3, yOffset+y+1, entry[1]*3, h-2);')
            h('            x+=entry[1];')
            h('          } else {')
            h('            x+=entry[1];')
            h('          }')

            if repeats > 1 :
                h('          if (x >= nextBreak) {')
                h('            x=0;')
                h('            yOffset += '+str(elementHeight+20)+';')
                h('          }')

            h('        }')

            h('      }')
            h('      y+=h;')
            h('    }')

            # Add the position of the repeat tract for each match
            h('    canvas.fillStyle = "#0033DD";')
            h('    for (t = 0; t < match[2].length; ++t) {')
            if (repeats > 1) :
                h('       x=(match[2][t][0] % 300)*3;')
                h('       yOffset =  Math.floor(match[2][t][0]/300)*'+str(elementHeight+20)+';')
            else :
                h('       x=match[2][t][0]*3;')
                h('       yOffset = 0;')
            h('       canvas.fillRect(base_x+x, yOffset+y, match[2][t][1]*3, 2);')
            h('       canvas.fillRect(base_x+x, yOffset+y-h-2, match[2][t][1]*3, 2);')
            h('    }')

            h('    y+=2*h;')
            h('  }')

            h('}')
            # Will have multiple graphics on the page so need to call draw function for all of them
            h('addLoadEvent(drawAlignment'+str(nth)+');')


    def createAlignmentMap(self, hsp) :
        ''' Return as a coded map of the alignment '''

        # 0 = match
        # 1 = positive
        # 2 = mismatch
        # 3 = queryGap
        # 4 = subjectGap
        # 5 = nothing (e.g. before the match)

        # Convenience only
        query = hsp.query
        subject = hsp.sbjct
        align = hsp.match

        # Start by coding every letter in the alignment
        coded = [0]*len(align)
        for i,letter in enumerate(align) :
            if letter == '+' :
                coded[i] = 1
            elif letter == ' ' :
                # Need to check whether it's a mismatch or a gap
                if query[i] == '-' :
                    coded[i] = 3
                elif subject[i] == '-' :
                    coded[i] = 4
                else :
                    coded[i] = 2
            else :
                coded[i] = 0

        # Now swallow the coding into a list of similar
        map = [[5, hsp.query_start-1]]
        on = coded[0]
        run = 0
        for code in coded :
            if code == on :
                run += 1
            else :
                map.append([on, run])
                on = code
                run = 1

        # Remember to include the last match
        map.append([on, run])

        return map


    def addAlignment(self, htmldoc, tract, match) :
        '''Add an alignment to the html document, in fixed space, split in reasonable chunks'''

        geneInfo = 'Gene length: '+str(match['geneLength'])+'bp / '+str(match['geneLength']//3)+'aa'
        if 'tractNo' in match :
            geneInfo += ' PV: Yes'
        else :
            geneInfo += ' PV: No'
        htmldoc.add('span', geneInfo, attributes={'style':'margin-left: 50px'})
        if 'function' in match:
            htmldoc.addRaw('<br />')
            htmldoc.add('span', 'Function: '+match['function'], attributes={'style':'margin-left: 50px'})

        if 'multipleHits' in match :
            for hsp in match['multipleHits'] :
                self.addSingleAlignment(htmldoc, tract, match, hsp)
        else :
            self.addSingleAlignment(htmldoc, tract, match, match['hsp'])


    def addSingleAlignment(self, htmldoc, tract, match, hsp) :
        query = hsp.query
        subject = hsp.sbjct
        subject_frame = hsp.frame[1]
        align = hsp.match

        # alignment details
        htmldoc.add("p", 'Score: {0:.2f} bits: {1:.2f} e-value: {2:f}<br />length: {3:d} gaps: {4:d} id: {5:d} positives: {6:d} coverage: {7:.2f} query coverage {8:.2f}'. \
            format(hsp.score, hsp.bits, hsp.expect, hsp.align_length, hsp.gaps, hsp.identities, hsp.positives, match['coverage'], match['queryCoverage']),
            attributes={'style':'margin-left: 50px'})

        # Restore low complexity characters for display
        newQuery=['']*len(query)
        translation = tract['translation']
        at = hsp.query_start-1 # BLAST is 1-indexed
        for i, letter in enumerate(query) :
            if letter == 'X' :
                # Low complexity region replaced with Xs
                if at >= 0 and at < len(translation) :
                    newQuery[i] = translation[at].lower()
            else :
                newQuery[i] = letter

            if letter != '-' :
                at += 1

        query = ''.join(newQuery)

        # split to 100 character lines
        step = 100

        queryName = shortenName(tract['locus']).rjust(10)
        subjectName = shortenName(match['locus']).rjust(10)
        numbers = len(str(max(hsp.query_start, hsp.query_end, hsp.sbjct_start, hsp.sbjct_end)))
        prefix = ' '*(10+numbers+3)

        # *3 is because offset is in DNA bases not codons as per the alignment
        def calcOffset(i, start, end, frame) :
            if frame >= 0 :
                return i*3+start
            else :
                return end-i*3

        with htmldoc.tag('div', attributes={'style':'font-family: monospace; white-space: pre; margin-left: 50px'}) :
            for i in range(0, len(query), step) :
                right = min(step, len(query)-i)-1
                offset = calcOffset(i, hsp.sbjct_start, hsp.sbjct_end, subject_frame)
                htmldoc.addRaw(queryName+' +'+str(i+hsp.query_start).ljust(numbers+1)+query[i:i+step]+' +' + str(i+hsp.query_start+right) +'<br />')
                htmldoc.addRaw(prefix+align[i:i+step]+'<br />')
                htmldoc.addRaw(subjectName+' +'+str(offset).ljust(numbers+1)+subject[i:i+step]+' +'+str(offset+right*3)+'<br /><br />')

    def createTractSummary(self) :
        '''Create HTML file summarising tract lengths'''
        print('Outputting tract summary', end='', flush=True)

        # Helper function
        def heatMap(x, total) :
            prop = x/max(1, total)
            r=255
            g=int(255*(1-prop))
            b=int(255*(1-prop))
            data = '{0:>2d} ({1:.1f}%)'.format(x, prop*100)
            return html.makeTag(
                'div', data,
                attributes={'style' : 'background-color:#{0:02x}{1:02x}{2:02x}'.format(r,g,b)})

        def pickOrdering(seq) :
            revseq = Seq.reverse_complement(seq)
            if (seq[0] > revseq[0]) :
                return seq
            else :
                return revseq

        with html.Doc(
            self.path+'/tractLengths.html',
            title='Tract Length Summary',
            styler=self.styler,
            sortTable = 'include/sorttable.js') as htmldoc :

            htmldoc.add("H1", "Tract length Summary")
            htmldoc.addRule()
            htmldoc.add("p", "Created: "+self.creationTime)
            # Add div that we can replace with the navigation later
            htmldoc.add('a', '', attributes={'name' : 'navigation'})
            htmldoc.add('p', '<div id="navigation"></div>')
            htmldoc.addRule()

            navigation=[]

            repeats = {}
            onCrossbreak = {}
            strainNames = self.getStrainNames()
            strainCounts = []

            for strainNo in range(0, self.guide['numStrains']) :
                # Load the strain
                strain = self.loadStrainPickle(strainNo)
                strainRepeats = {}

                for contig in strain['contigs'] :
                    for tract in contig['tracts'] :

                        type = pickOrdering(tract['repeatUnit'])

                        if type not in strainRepeats :
                            strainRepeats[type] = {
                                'counts' : Counter(),
                                'onCounts' : Counter(),
                                'typeCounts' : Counter()
                            }

                        if type not in onCrossbreak :
                            onCrossbreak[type] = {'Intergenic' : Counter()}

                        if tract['outby'] == 0 :
                            # In frame, use the version in the tract direction
                            strainRepeats[type]['typeCounts'][tract['repeatUnit']] += 1
                        else :
                            strainRepeats[type]['typeCounts']['intergenic'] += 1

                        strainRepeats[type]['counts'][tract['repeats']] += 1

                        if 'onLength' in tract and tract['onLength'] != None :
                            strainRepeats[type]['onCounts'][tract['onLength']] += 1
                            if tract['onLength'] not in onCrossbreak[type] :
                                onCrossbreak[type][tract['onLength']] = Counter()
                            onCrossbreak[type][tract['onLength']][tract['repeats']] += 1
                        else :
                            onCrossbreak[type]['Intergenic'][tract['repeats']] += 1

                strainCounts.append(strainRepeats)

                for type, content in strainRepeats.items() :
                    if type not in repeats :
                        repeats[type] = copy.deepcopy(content)
                        repeats[type]['strains'] = {}
                    else :
                        for k, v in content.items() :
                            repeats[type][k] += v
                    repeats[type]['strains'][strain['name']] = copy.deepcopy(content)

                print('.', end='', flush=True)
            print('')

            # Overall repeat summary table
            # Order the repeat types
            repeatTypes = [x for x in repeats.keys()]
            repeatTypes.sort(key = lambda x : (len(x), x))

            # Find minimum and maximum repeat lengths
            shortest = 100000
            longest = 0

            totalTracts = 0

            for _, repeat in repeats.items() :
                totalTracts += sum(repeat['counts'].values())
                minRepeat = min(repeat['counts'] + repeat['onCounts'])
                maxRepeat = max(repeat['counts'] + repeat['onCounts'])
                if minRepeat < shortest :
                    shortest = minRepeat
                if maxRepeat > longest :
                    longest = maxRepeat

            # Summary table of tract types and length
            htmldoc.add('a', '', attributes={'name' : 'summary'})
            htmldoc.add("H2", html.makeLink('#navigation', '[top]')+" Summary data for Tract type and length")
            navigation.append(html.makeLink('#summary', 'Summary data for tract type and length'))
            header = ['RepeatType'] + list(range(shortest, longest+1)) + ['Genic', 'Intergenic', 'Total']
            with htmldoc.table() as table :
                table.addHeaderRow(header)
                for type in repeatTypes :
                    repeat = repeats[type]
                    row = [type]
                    total = sum(repeat['counts'].values())
                    row += [heatMap(repeat['counts'][i], total) for i in range(shortest, longest+1)]
                    row += [heatMap(total-repeat['typeCounts']['intergenic'], total), heatMap(repeat['typeCounts']['intergenic'], total), heatMap(total, totalTracts)]
                    table.addRow(row)

            htmldoc.addRule()

            # ON-length summary table
            htmldoc.add('a', '', attributes={'name' : 'onLength'})
            htmldoc.add("H2", html.makeLink('#navigation', '[top]')+" Summary data for Tract On-length (genic only)")
            navigation.append(html.makeLink('#onLength', 'Summary data for Tract On-length (genic only)'))
            header = ['RepeatType'] + list(range(shortest, longest+1)) + ['Total']
            with htmldoc.table() as table :
                table.addHeaderRow(header)
                for type in repeatTypes :
                    repeat = repeats[type]
                    total = sum(repeat['onCounts'].values())

                    if total == 0 :
                        continue

                    row = [type]
                    row += [heatMap(repeat['onCounts'][i], total) for i in range(shortest, longest+1)]
                    row += [total]
                    table.addRow(row)

            htmldoc.addRule()

            # Strain vs. type data
            htmldoc.add('a', '', attributes={'name' : 'byStrain'})
            htmldoc.add('H2', html.makeLink('#navigation', '[top]')+' Repeat types by strain')
            navigation.append(html.makeLink('#byStrain', 'Repeat types by strain'))
            with htmldoc.table() as table :
                table.addHeaderRow(['species', 'strain']+repeatTypes)
                for i, strain in enumerate(strainNames) :
                    this = strainCounts[i]
                    row = [self.guide['species'][i], strain]
                    total = sum(sum(this[t]['counts'].values()) for t in this)
                    for t in repeatTypes :
                        if t in this :
                            row += [heatMap(sum(this[t]['counts'].values()), total)]
                        else :
                            row += [heatMap(0, total)]
                    row += [total]
                    table.addRow(row)

            # Species vs. type data
            htmldoc.add('a', '', attributes={'name' : 'bySpecies'})
            htmldoc.add('H2', html.makeLink('#navigation', '[top]')+' Repeat types by species')
            navigation.append(html.makeLink('#bySpecies', 'Repeat types by species'))
            speciesCounts = {}
            for i, strainCount in enumerate(strainCounts) :
                species = self.guide['species'][i]
                if species not in speciesCounts :
                    speciesCounts[species] = {}
                for t in strainCount :
                    if t not in speciesCounts[species] :
                        speciesCounts[species][t] = Counter()
                    speciesCounts[species][t] += strainCount[t]['counts']

            speciesNames = sorted(speciesCounts)

            with htmldoc.table() as table :
                table.addHeaderRow(['species']+repeatTypes)
                for i, species in enumerate(speciesNames) :
                    this = speciesCounts[species]
                    row = [species]
                    total = sum(sum(this[t].values()) for t in this)
                    for t in repeatTypes :
                        if t in this :
                            row += [heatMap(sum(this[t].values()), total)]
                        else :
                            row += [heatMap(0, total)]
                    row += [total]
                    table.addRow(row)

            # Tract type data broken down by strain
            for nth, type in enumerate(repeatTypes) :
                repeat = repeats[type]
                htmldoc.add('a', '', attributes={'name' : 'repeats'+str(nth)})
                heading = 'Strain data for '+str(type)+'/'+str(Seq.reverse_complement(type)) + ' tracts'
                htmldoc.add('H2', html.makeLink('#navigation', '[top]')+' ' + heading)
                navigation.append(html.makeLink('#repeats'+str(nth), heading))

                # Find longest and shortest repeat
                # And also count up the totals
                shortest = 100000
                longest = 0
                totals = Counter()
                typeTotals = Counter()
                for _, strain in repeat['strains'].items() :
                    s = min(strain['counts'])
                    l = max(strain['counts'])
                    if s < shortest :
                        shortest = s
                    if l > longest :
                        longest = l
                    totals += strain['counts']
                    typeTotals += strain['typeCounts']

                with htmldoc.table() as table :
                    table.addHeaderRow(['Strain']+list(range(shortest, longest+1))+[type, Seq.reverse_complement(type), 'intergenic', 'total'])
                    for strainName in strainNames :
                        if strainName not in repeat['strains'] :
                            continue
                        strain = repeat['strains'][strainName]
                        total = sum(strain['counts'].values())
                        row = [strainName]
                        row += [heatMap(strain['counts'][i], total) for i in range(shortest, longest+1)]
                        row += [heatMap(strain['typeCounts'][i], total) for i in [type, Seq.reverse_complement(type), 'intergenic']]
                        row += [total]
                        table.addRow(row)
                    totalInStrain = sum(totals.values())
                    totalRow = ['Totals']
                    totalRow += [heatMap(totals[i], totalInStrain) for i in range(shortest, longest+1)]
                    totalRow += [heatMap(typeTotals[i], totalInStrain) for i in [type, Seq.reverse_complement(type), 'intergenic']]
                    totalRow += [totalInStrain]
                    table.addFooterRow(totalRow)

                htmldoc.addRule()

            # Crossbreak of putative ON length and actual tract length
            for nth, type in enumerate(repeatTypes) :
                crossbreak = onCrossbreak[type]
                htmldoc.add('a', '', attributes={'name' : 'putativeOn'+str(nth)})
                heading = 'On vs. Length comparison for '+str(type)+'/'+str(Seq.reverse_complement(type)) + ' tracts'
                htmldoc.add('H2', html.makeLink('#navigation', '[top]')+' ' + heading)
                navigation.append(html.makeLink('#putativeOn'+str(nth), heading))

                # Find longest and shortest repeat
                # And also count up the totals
                shortest = 100000
                longest = 0
                totals = Counter()
                for lengths in crossbreak.values() :
                    if len(lengths) == 0 :
                        continue
                    s = min(lengths)
                    l = max(lengths)
                    if s < shortest :
                        shortest = s
                    if l > longest :
                        longest = l
                    totals += lengths

                with htmldoc.table() as table :
                    table.addHeaderRow(['On Length']+list(range(shortest, longest+1))+['total'])
                    sortedLengths = sorted(crossbreak.keys(), key = lambda x : 10000 if x == 'Intergenic' else x)
                    for onLength in sortedLengths :
                        if onLength not in crossbreak :
                            continue
                        line = crossbreak[onLength]
                        total = sum(line.values())
                        row = [onLength]
                        if onLength == 'Intergenic' :
                            row += [heatMap(line[i], total) for i in range(shortest, longest+1)]
                        else :
                            row += [html.makeTag('b', heatMap(line[i], total)) if i == onLength else heatMap(line[i], total) for i in range(shortest, longest+1)]
                        row += [total]
                        table.addRow(row)
                    totalInStrain = sum(totals.values())
                    totalRow = ['Totals']
                    totalRow += [heatMap(totals[i], totalInStrain) for i in range(shortest, longest+1)]
                    totalRow += [totalInStrain]
                    table.addFooterRow(totalRow)

                htmldoc.addRule()




            # Add navigation text into existing div
            with htmldoc.tag('script', attributes={'type' : 'text/javascript'}) :
                navText = '<br/>'.join(navigation)
                navText = navText.replace('"', '\\"')
                htmldoc.addRaw('document.getElementById("navigation").innerHTML = "'+navText+'"')


        self.index.add('p', html.makeLink('tractLengths.html', 'Tract lengths'))
        with self.index.table(sortable=False) as table :
            table.addHeaderRow(['Tract Type'] + repeatTypes)
            table.addRow(['<b>Count</b>'] + [heatMap(sum(repeats[x]['counts'].values()), totalTracts) for x in repeatTypes])
        self.index.addRule()

    def getStrainNames(self) :
        return self.guide['strainNames']

    def getStrainNumber(self, strainName) :
        return self.strainLookup[strainName]

    def findStrainNameForContig(self, contigName) :
        return self.contigToStrain[contigName]

    def getTractAminoPos(self, tract) :
        tractLength = (len(tract['repeatUnit'])*tract['repeats'])//3
        if 'strand' in tract and tract['strand'] == -1 :
            # Translation started from the END of the expected region
            tractStart = (tract['end'] - tract['location'])//3 - tractLength
        else:
            tractStart = (tract['location']-tract['start'])//3
        return tractStart, tractLength

    def listBlastMatches(self, blastMatches, strains) :
        '''Convert the blastMatches into text readable entries with spaces for missing strains'''
        matchList = []

        strainMatches = set()
        for contigMatch in blastMatches.keys() :
            strainMatches.add(self.findStrainNameForContig(contigMatch))

        for n, other in enumerate(strains) :
            if other not in strainMatches :
                matchList.append('')
                continue

            strainHtml = str(n)+'.html'

            allMatches = []

            for contigName, matches in blastMatches.items() :
                if self.findStrainNameForContig(contigName) != other :
                    continue
                allMatches = allMatches + matches

            matchText = ''
            for match in allMatches :
                if len(matchText) :
                    matchText = matchText + '<br/>'
                matchText = matchText
                if 'tractNo' in match :
                    matchText = matchText+html.makeTag('span',html.makeLink(strainHtml+'#'+match['locus'], match['locus']), attributes={'title':match['locus']+' (PV, E={0:g} {1:.2f}%)'.format(match['expect'], match['coverage']*100.0)})
                else :
                    matchText = matchText+html.makeTag('span',match['locus'], attributes={'title':match['locus']+' (E={0:g} {1:.2f}%)'.format(match['expect'], match['coverage']*100.0)})
            matchList.append(matchText)

        return matchList

    def runPlugins(self, func, *args, **kwargs) :
        for plugin in self.plugins :
            call = getattr(plugin, func, None)
            if call :
                call(*args, **kwargs)

    def addPhylogenicTree(self, htmldoc, type, names, numbers, matrix, tag=None) :
        dm = _DistanceMatrix(names=['##{0:d}##'.format(i) for i in range(0, len(names))], matrix=matrix)
        constructor = DistanceTreeConstructor()
        constructorFn = getattr(constructor, type)
        try:
            tree = constructorFn(dm)
            if tag==None:
                tag=''

            def buildTreeFrom(clade) :
                if not clade.clades :
                    id = int(clade.name.replace('#', ''))
                    return 1, clade.branch_length, {'width': 1, 'branchLength' : clade.branch_length, 'id' : numbers[id]}
                width = 0
                deepest = 0
                branches = []
                for c in clade.clades :
                    w, d, b = buildTreeFrom(c)
                    if d > deepest :
                        deepest = d
                    width += w
                    branches.append(b)
                return width, deepest+clade.branch_length, {'width' : width, 'branchLength' : clade.branch_length, 'branches' : branches}

            width, depth, root = buildTreeFrom(tree.root)
            htmldoc.add('div', '', attributes={'id' : 'pretree_'+type+tag})
            self.addTreeGraphic(htmldoc, {numbers[i]:names[i] for i in range(0, len(names))}, type, depth, root, tag)
        except :
            # Problem creating tree, this is nearly always from ConstructorFn
            # Need to add the pretree to avoid errors later
            htmldoc.add('div', '', attributes={'id': 'pretree_' + type + tag})
            # Add text showing what is wrong
            htmldoc.add('p', 'Tree construction failed.')
            htmldoc.addCollapsed('<pre>'+traceback.format_exc()+'</pre>', 'Detailed error report')
            logging.exception('Tree construction failed.')

    def addTreeGraphic(self, htmldoc, names, type, depth, root, tag) :
        '''
        Add a graphical representation of this tree to the htmldoc
        New version with blobs and numbers
        '''
        def h(x) :
            htmldoc.addRaw(x)

        htmldoc.add('canvas', '', attributes={'id' : type+tag+'treeGraphic', 'width' : 1200, 'height' : 900})
        with htmldoc.tag('script', attributes={'type' : 'text/javascript'}) :
            h('function draw'+type+tag+'TreeGraphic(cols, words) {')
            h('  var canvas_element = document.getElementById("'+type+tag+'treeGraphic");')
            #h('  canvas_element.style.backgroundColor = "rgba(158, 167, 184, 0.2)";')
            h('  var canvas = canvas_element.getContext("2d");')
            h('  var tree = '+str(root)+';')
            h('  var names = '+str(names)+';')
            if settings.numberTrees :
                h('  var longestName = '+str(len(str(len(names))))+';')
            else :
                h('  var longestName = ' + str(max(len(name) for name in names.values())) + ';')
            h('  var size = 300.0 + 2.0*'+str(len(names))+';')
            h('  var originX = 500;')#size+(longestName+3)*8;')
            h('  var originY = 600;')#size+(longestName+3)*8;')
            h('  var size = originX - (longestName+3)*8;')
            h('  var longestDim = 1200;')
            h('  var shortestDim = 1000;')
            h('  canvas_element.width = shortestDim;')
            h('  canvas_element.height = longestDim;')
            h('  canvas.scale(1.0,1.0);')
            h('  var scale = size/'+str(depth)+';')
            h('  canvas.font="12px Courier New";')

            # Determine the maximum extents of the image
            h('  function extentStep(node, radius, min, max, extents) {')
            h('    var minAngle = max;')
            h('    var maxAngle = min;')
            h('    if ("id" in node) {')
            h('      var angle = (min+max)/2.0;')
            h('      minAngle = maxAngle = angle;')

            # Radius of line + blob + text
            if settings.numberTrees :
                h('      var text = (node.id + 1).toString();')
            else :
                h('      var text = names[node.id];')
                h('      if (words) {')
                h('        text += "("+words[node.id]+")";')
                h('      }')

            h('      var endRadius = radius+node["branchLength"]*scale + 24 + text.length*8;')
            h('      var cosAngle = Math.cos(angle);')
            h('      var sinAngle = Math.sin(angle);')
            h('      x = originX+cosAngle*endRadius;')
            h('      y = originY+sinAngle*endRadius;')
            # Update extents
            h('      if (x < extents[0]) extents[0] = x;')
            h('      if (x > extents[1]) extents[1] = x;')
            h('      if (y < extents[2]) extents[2] = y;')
            h('      if (y > extents[3]) extents[3] = y;')
            h('    } else {')
            h('      var branches = node["branches"];')
            h('      var step = (max - min)/node["width"];')
            h('      var on = min;')
            h('      for (var i = 0; i < branches.length; ++i) {')
            h('        var branch = branches[i];')
            h('        var max = on+(step*(branch["width"]));')
            h('        var angle;')
            h('        if ("id" in branch) {')
            h('          angle = (on+max)*0.5;')
            h('          extentStep(branch, radius, on, max, extents);')
            h('        } else {')
            h('          angle = (on+max)/2.0;')
            h('          var endRadius = radius+branch["branchLength"]*scale;')
            h('          var cosAngle = Math.cos(angle);')
            h('          var sinAngle = Math.sin(angle);')
            h('          extentStep(branch, endRadius, on, max, extents);')
            h('        }')
            h('        if (angle < minAngle) { minAngle = angle; }')
            h('        if (angle > maxAngle) { maxAngle = angle; }')
            h('        on = max;')
            h('      }')
            h('    }')
            h('    return [minAngle, maxAngle, extents];')
            h('  }')

            h('  function findOrigin(min, max, origin) {')
            h('    var over = max-origin;')
            h('    var under = origin-min;')
            h('    return (origin*2*(under/(over+under)));')
            h('  }')

            h('  function drawStep(node, radius, min, max) {')
            h('    var minAngle = max;')
            h('    var maxAngle = min;')
            h('    if ("id" in node) {')
            h('      var angle = (min+max)/2.0;')
            h('      minAngle = maxAngle = angle;')
            h('      var endRadius = radius+node["branchLength"]*scale;')
            h('      var cosAngle = Math.cos(angle);')
            h('      var sinAngle = Math.sin(angle);')
            h('      if (cols)')
            h('        canvas.strokeStyle= (node.id in cols) ? cols[node.id][0] : "#c3c3c3";')
            h('      else')
            h('        canvas.strokeStyle="#000000";')
            h('      canvas.beginPath();')
            h('      canvas.moveTo(originX+cosAngle*radius, originY+sinAngle*radius);')
            h('      canvas.lineTo(originX+cosAngle*endRadius, originY+sinAngle*endRadius);')
            h('      canvas.stroke();')
            # Draw multi-coloured lines if desired
            h('      if(cols && node.id in cols && (cols[node.id][0] != cols[node.id][1])) {')
            h('        canvas.strokeStyle= cols[node.id][1];')
            h('        canvas.setLineDash([15, 15]);')
            h('        canvas.beginPath();')
            h('        canvas.moveTo(originX+cosAngle*radius, originY+sinAngle*radius);')
            h('        canvas.lineTo(originX+cosAngle*endRadius, originY+sinAngle*endRadius);')
            h('        canvas.stroke();')
            h('        canvas.setLineDash([]);')
            h('      }')

            # Add blob
            h('      canvas.beginPath();')
            h('      if (cols)')
            h('        canvas.fillStyle= (node.id in cols) ? cols[node.id][0] : "#c3c3c3";')
            h('      else')
            h('        canvas.fillStyle="#000000";')
            h('      canvas.arc(originX+cosAngle*(endRadius+10), originY+sinAngle*(endRadius+10), 8, 0, 2*Math.PI, true);')
            h('      canvas.fill();')

            h('      if(cols && node.id in cols && (cols[node.id][0] != cols[node.id][1])) {')
            h('        canvas.beginPath();')
            h('        canvas.fillStyle= cols[node.id][1];')
            h('        var halfAngle = angle+Math.PI;')
            h('        if (angle > Math.PI/2.0 && angle < 3.0*Math.PI/2.0) halfAngle = angle;')
            h('        canvas.arc(originX+cosAngle*(endRadius+10), originY+sinAngle*(endRadius+10), 8, halfAngle, halfAngle + Math.PI, true);')
            h('        canvas.fill();')
            h('      }')

            h('      canvas.save();')
            h('      var textCols = ["#000000", "#000000"];')
            h('      if (cols)')
            h('        textCols = (node.id in cols) ? cols[node.id] : ["#c3c3c3", "#c3c3c3"];')
            h('      canvas.translate(originX+cosAngle*endRadius, originY+sinAngle*endRadius);')
            if settings.numberTrees :
                h('      var text = (node.id + 1).toString();')
            else :
                h('      var text = names[node.id];')
                h('      if (words) {')
                h('        text += "("+words[node.id]+")";')
                h('      }')
            h('      if (angle > Math.PI/2.0 && angle < 3.0*Math.PI/2.0) {')
            h('        canvas.rotate(angle+Math.PI);')
            h('        colouredText(canvas, textCols[0], textCols[1], text, -(text.length*8+24), 3);')
            h('      } else {')
            h('        canvas.rotate(angle);')
            h('        colouredText(canvas, textCols[0], textCols[1], text, 24, 3);')
            h('      }')
            h('      canvas.restore();')
            h('    } else {')
            h('      var branches = node["branches"];')
            h('      var step = (max - min)/node["width"];')
            h('      var on = min;')
            h('      for (var i = 0; i < branches.length; ++i) {')
            h('        var branch = branches[i];')
            h('        var max = on+(step*(branch["width"]));')
            h('        var angle;')
            h('        if ("id" in branch) {')
            h('          angle = (on+max)*0.5;')
            h('          drawStep(branch, radius, on, max);')
            h('        } else {')
            h('          angle = (on+max)/2.0;')
            h('          var endRadius = radius+branch["branchLength"]*scale;')
            h('          var cosAngle = Math.cos(angle);')
            h('          var sinAngle = Math.sin(angle);')
            h('          canvas.strokeStyle="#000000";')
            h('          canvas.beginPath();')
            h('          canvas.moveTo(originX+cosAngle*radius, originY+sinAngle*radius);')
            h('          canvas.lineTo(originX+cosAngle*endRadius, originY+sinAngle*endRadius);')
            h('          canvas.stroke();')
            h('          var arc = drawStep(branch, endRadius, on, max);')
            h('          canvas.strokeStyle="#000000";')
            h('          canvas.beginPath();')
            h('          canvas.arc(originX, originY, endRadius, arc[0], arc[1]);')
            h('          canvas.stroke();')
            h('        }')
            h('        if (angle < minAngle) { minAngle = angle; }')
            h('        if (angle > maxAngle) { maxAngle = angle; }')
            h('        on = max;')
            h('      }')
            h('    }')
            h('    return [minAngle, maxAngle];')
            h('  }')

            # Use the extentStep to determine extents
            # And then draw rescaled
            h('  var extents = [originX, originX, originY, originY];')
            h('  extents = extentStep(tree, 0.0, 0.0, Math.PI*2, extents)[2];')
            h('  var maximal = Math.max(extents[1]-extents[0], extents[3]-extents[2]);')
            h('  var minimal = Math.min(extents[1]-extents[0], extents[3]-extents[2]);')

            h('  if (extents[1] - extents[0] > extents[3] - extents[2]) {')
            h('    extents = [originX, originX, originY, originY];')
            h('    extents = extentStep(tree, Math.PI / 2, Math.PI / 2, Math.PI / 2 + Math.PI * 2, extents)[2];')
            h('    scale *= Math.min(longestDim / maximal, shortestDim / minimal);')
            h('    originX = findOrigin(extents[0], extents[1], originX, originX * 2);')
            h('    originY = findOrigin(extents[2], extents[3], originY, originY * 2);')
            h('    drawStep(tree, Math.PI / 2, Math.PI / 2, Math.PI / 2 + Math.PI * 2);')
            h('  } else {')
            h('    scale *= Math.min(longestDim/maximal, shortestDim/minimal);')
            h('    originX = findOrigin(extents[0], extents[1], originX);')
            h('    originY = findOrigin(extents[2], extents[3], originY);')
            h('    drawStep(tree, 0, 0, Math.PI*2);')
            h('  }')
            h('}')

            # Add the drawing event to the load events
            h('addLoadEvent(draw'+type+tag+'TreeGraphic);')

    def createCorePhasome(self) :
        '''Create a page listing the core phasome for each species'''

        # Find groups in >50% of each species for species with >2 representatives
        phasome = {}
        speciesInOrder = []
        for i in range(0, self.guide['numStrains']) :
            strain = self.loadStrainPickle(i)
            if strain['species'] not in phasome :
                phasome[strain['species']] = {'count' : Counter(), 'total' : 0}
                speciesInOrder.append(strain['species'])
            phasome[strain['species']]['total'] += 1
            done = set()
            for contig in strain['contigs'] :
                for tract in contig['tracts'] :
                    groupNo = tract['geneGroup']
                    if groupNo >= 0 and groupNo not in done:
                        phasome[strain['species']]['count'][groupNo] += 1
                        done.add(groupNo)

        # Add to index
        self.index.add('p', html.makeLink('core_phasome.html', 'Core phasome for species'))
        self.index.addRule()

        # Create HTML document
        with html.Doc(
            self.path+'/core_phasome.html',
            title='Core Phasomes',
            styler=self.styler
        ) as htmldoc :
            htmldoc.add("H1", "Core phasome for species with 2+ genomes")
            htmldoc.addRule()
            htmldoc.add("p", "Created: "+self.creationTime)
            htmldoc.addRule()

            for species in speciesInOrder :
                ours = phasome[species]
                if ours['total'] < 2 :
                    continue

                # Add header
                htmldoc.add('H2', 'Core phasome for '+species+' (n='+str(ours['total'])+')')


                # Add the phasome
                occurence = {}

                for group, count in ours['count'].items() :
                    if count > ours['total']/2 :
                        percent = count/ours['total']
                        if percent not in occurence :
                            occurence[percent] = []
                        occurence[percent].append(group)

                # Count of core phasome
                totalCore = 0
                for percent, groups in occurence.items() :
                    if percent >= 0.6 :
                        totalCore += len(groups)

                htmldoc.add('p', 'Core phasome count (groups in >= 60%) of strains: '+str(totalCore))

                if len(occurence) :
                    with htmldoc.table() as table :
                        table.addHeaderRow(['In % of strain', 'Gene group', 'Name', 'Function'])
                        for percent in sorted(occurence.keys(), reverse=True) :
                            groups = occurence[percent]
                            with table.row() as row :
                                row.addCell('{0:.1f} ({1})'.format(percent*100, len(occurence[percent])), attributes={'rowspan':len(occurence[percent])})
                                row.addCell(html.makeLink('groups/'+str(groups[0])+'.html', str(groups[0])))
                                row.addCell(self.geneGroups[groups[0]]['name'])
                                row.addCell(self.geneGroups[groups[0]]['bestFunction'])
                            for group in groups[1:] :
                                with table.row() as row :
                                    row.addCell(html.makeLink('groups/'+str(group)+'.html', str(group)))
                                    row.addCell(self.geneGroups[group]['name'])
                                    row.addCell(self.geneGroups[group]['bestFunction'])
                else :
                    htmldoc.add('p', 'No core phasome identified')

                htmldoc.addRule()