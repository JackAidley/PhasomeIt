# Find associations between the identified groups and on/off states
# And a .csv file containing metadata on the files

import csv
import os
import pickle
import html
import logging
from collections import Counter

class AssociationFinder :
    def readMetadata(self, path, csvFilename, fileColumn, columns) :
        self.meta = {}

        # Load guide pickle
        with open(os.path.join(path, 'pickle', 'guide.pkl'), 'rb') as file :
            guide = pickle.load(file)

        idToStrainNo = {}
        for n,strainName in enumerate(guide['strainNames']) :
            id = strainName.split('.')[0]
            idToStrainNo[id] = n

        try:
            idsNotInDataset = []
            with open(csvFilename, newline='') as csvFile:
                # Not using the DictReader as we want to customize slightly
                reader = csv.reader(csvFile, dialect='excel')
                header = next(reader)
                if not columns :
                    columns = range(1, len(header))
                self.tags = [header[i] for i in columns]
                self.tagValues = {t:set() for t in self.tags}
                for n,row in enumerate(reader):
                    id = row[fileColumn]
                    if id not in idToStrainNo :
                        idsNotInDataset.append(id)
                        continue
                    data = {}
                    for i in columns :
                        if i > len(row) :
                            print('Line {0} of metadata file "{1}" does not contain expected number of columns'.format(n + 1, csvFilename))
                            continue
                        value = row[i]
                        if value == '0' or value == 0 :
                            value = 'No data'
                        data[header[i]] = value
                        self.tagValues[header[i]].add(value)

                    self.meta[idToStrainNo[id]] = data

            if len(idsNotInDataset) :
                print('Warning: ids '+str(idsNotInDataset)+' appear in metadata but have not been analysed')
        except FileNotFoundError :

            print('No metadata file found, adding species as metadata')
            logging.warning('No metadata file found, adding species as metadata')
            self.addSpeciesAsMetadata(path)

        self.finaliseTags()

    def addSpeciesAsMetadata(self, path) :
        self.meta = {}

        # Load guide pickle
        with open(os.path.join(path, 'pickle', 'guide.pkl'), 'rb') as file :
            guide = pickle.load(file)

        self.tags = ['species']
        self.tagValues = {t: set() for t in self.tags}
        for i in range(0, len(guide['species'])):
            self.meta[i] = {'species': guide['species'][i]}
            self.tagValues['species'].add(guide['species'][i])

        self.finaliseTags()

    def finaliseTags(self) :
        tagValues = {}
        for t,v in self.tagValues.items() :
            v = sorted(list(v), key=lambda x : (x != 'No data', x))
            tagValues[t] = v
        self.tagValues = tagValues

        # Okay we've got the meta data, now we want to create lookups between each column
        self.lookup = {tag:{} for tag in self.tags}
        for id, isolate in self.meta.items() :
            for tag in self.tags :
                if isolate[tag] not in self.lookup[tag] :
                    self.lookup[tag][isolate[tag]] = [id]
                else:
                    self.lookup[tag][isolate[tag]] += [id]

    def loadStrainPickle(self, path, n) :
        with open(os.path.join(path, 'pickle', str(n)+'.pkl'), 'rb') as file :
            strain = pickle.load(file)
        return strain

    def findAssociations(self, path) :
        # TODO: this is not currently functional!

        if not self.meta :
            print('AssociationFinder.readMetadata must be called before findAssociations')
            return

        self.associations = {}

        # What we want to do is find gene groups that are preferentially associated with particular traits
        # Or, alternatively, genes that are preferentially on or off with particular traits
        # But first we need to get hold of our data and parse it into a suitable form
        # This repeats work done in the group file creation section which I don't like
        # Perhaps I should create this data somewhere and share it that way?

        # Load guide pickle
        with open(os.path.join(path, 'pickle', 'guide.pkl'), 'rb') as file :
            guide = pickle.load(file)

        # Load each strain and store off what gene groups it has, how many in each group and their on/off state
        geneGroups = {}
        strainsToGroups = {}
        for n, strainName in enumerate(guide['strains']) :
            strain = self.loadStrainPickle(path, n)
            for contig in strain['contigs'] :
                for tract in contig['tracts'] :
                    tractGroup = tract['geneGroup']
                    if tractGroup == -1 :
                        tractGroups = 'None'
                    if tractGroup not in geneGroups :
                        geneGroups[tractGroup] = {'exists':[], 'on':[], 'near':[], 'off':[]}
                    geneGroups[tractGroup]['exists'].append(n)
                    if tract['outby'] != 0 :
                        geneGroups[tractGroup]['near'].append(n)
                    else :
                        if tract['tract'] == tract['onLength'] :
                            geneGroups[tractGroup]['on'].append(n)
                        else :
                            geneGroups[tractGroup]['off'].append(n)

                    if n not in strainsToGroups :
                        strainsToGroups[n] = set()
                    strainsToGroups[n].add(tractGroup)

        # Produce sets for the gene groups
        for group in geneGroups.values() :
            group['set'] = set(group['exists'])

        # sub-divide by ST-complex
        # And then find most over-represented types
        # And THEN exclude ones highly represented in others
        # To give a list of over-represented
        #stHigh = {t:[] for t in self.tagValues['ST Complex']}

        #for groupNo, groupData in geneGroups.items() :


        # For each tag, find either:
        # 1. An abundance of existence or absence
        # 2. For common (> 20%) genes, an abundance of multiple copies/differential counts
        # 3. For majority (> 50%) genes, an abundance of ON genes
        #size_cutoff10 = max(len(guide['strains'])//10, 3)
        #size_cutoff5 = max(len(guide['strains'])//20, 3)
        # pvalue_cutoff = 0.05
        #
        # for groupNo, groupData in geneGroups.items() :
        #     # For each gene group check association with any of the values
        #     for tag, tagEntries in self.lookup.items() :
        #         contingency = [[], []]
        #         colNames = []
        #         for onEntry, (entryTag, entryStrains) in enumerate(tagEntries.items()) :
        #             if len(entryStrains) > size_cutoff5 :
        #                 contingency[0].append(0)
        #                 contingency[1].append(0)
        #                 colNames.append(entryTag)
        #                 for strain in entryStrains :
        #                     if strain in groupData['set'] :
        #                         contingency[0][-1] += 1
        #                     else :
        #                         contingency[1][-1] += 1
        #
        #         # Check significance of contingency table
        #         if sum(contingency[0]) < size_cutoff10 :
        #             continue
        #
        #         pvalue = stats.chi2_contingency(contingency)[1]
        #         if pvalue <= pvalue_cutoff :
        #             association = [['']+colNames]
        #             association += [['Present']+contingency[0]]
        #             association += [['Absent']+contingency[1]]
        #             self.addAssociation(tag, groupNo, association, pvalue)

    def addAssociation(self, tag, geneGroup,table, pvalue) :
        if tag not in self.associations :
            self.associations[tag] = {}
        self.associations[tag][geneGroup] = {'table' : table, 'pvalue' : pvalue }

    def addToHTML(self, creator) :
        # This page is useless, don't add it
        return
        # Add the new page to the index
        creator.index.add('p', html.makeLink('associations.html', 'Associations'))
        creator.index.addRule()

        # Create the New Page
        with html.Doc(creator.path+'/associations.html', title='Associations', styler=creator.styler) as htmldoc :
            htmldoc.add("H1", "Associations")
            htmldoc.addRule()
            htmldoc.add("p", "Created: "+creator.creationTime)
            htmldoc.addRule()

            for tag, tagData in self.associations.items() :
                htmldoc.add("H2", "Associations with: " + tag)

                for group, data in tagData.items() :
                    with htmldoc.tag('p') :
                        htmldoc.addRaw(html.makeTag('b', 'Gene group: '+html.makeLink('groups/'+str(group)+'.html', str(group)))+' with p-value '+html.sciFormat(data['pvalue'])+'<br />')
                        with htmldoc.table() as table :
                            table.addHeaderRow(data['table'][0])
                            for row in data['table'][1:] :
                                table.addRow(row)

                htmldoc.addRule()

    def preTree(self, htmldoc, prefix) :
        '''Before the tree add code for the drop-down select'''
        with htmldoc.tag('div') :
            htmldoc.addRaw('Select trait to display on graph: ')
            with htmldoc.tag('select', attributes={'onchange' : 'findAssociationsOnChange'+prefix+'(this);'}) :
                # Add all the options
                htmldoc.add('option', '----', attributes={'value' : '-1'})
                for n, tag in enumerate(self.tags) :
                    htmldoc.add('option',tag, attributes={'value' : str(n)})

    def postTree(self, htmldoc, strainNames, strainNos, prefix) :
        '''After the tree add the script to update all the values when something is chosen from the selection'''
        # Colours are the same set as the phasocols used in my thesis but with black replaced with a yellow
        colours = ['#000000', '#009280', '#ff6dce', '#4900b6', '#6db6ff', '#920000', '#aaaa29', '#6d6d6d']

        # Show un-annotated items in a faint grey
        noDataColour = '#c3c3c3'
        def nextColour(n) :
            n = n % ((len(colours)*(len(colours)+1))//2)
            if n < len(colours) :
                return [colours[n], colours[n]]
            else :
                offset = 0
                while n >= len(colours) - offset and offset < len(colours) :
                    n -= len(colours) - offset
                    offset += 1
                return [colours[n], colours[n+offset]]

        # Always use the same colours regardless of which ones are in this tree
        tagLookup = {t:{} for t in self.tags}
        for tag, states in self.tagValues.items() :
            n=0
            for state in states :
                if state != 'No data' :
                    tagLookup[tag][state] = nextColour(n)
                    n+=1
                else :
                    tagLookup[tag][state] = [noDataColour, noDataColour]

        metaLookup = {t:{} for t in self.tags}
        wordsLookup = {t:{} for t in self.tags}
        usedStates = {t:set() for t in self.tags}
        for strainNo in  strainNos :
            if strainNo in self.meta :
                for tag, state in self.meta[strainNo].items() :
                    metaLookup[tag][strainNo] = tagLookup[tag][state]
                    wordsLookup[tag][strainNo] = state
                    usedStates[tag].add(state)


        # Create key, using the tagValues ordering
        keys = {}
        longestState = 0
        for tag, states in self.tagValues.items() :
            key = [ ["#ffffff", "#ffffff", "Key for "+tag] ]
            for state in states :
                if state not in usedStates[tag] :
                    continue
                col = tagLookup[tag][state]
                key.append([col[0], col[1], state])
                if len(state) > longestState :
                    longestState = len(state)
            keys[tag] = key

        names=[s for s in strainNames]

        with htmldoc.tag('script', attributes={'type' : 'text/javascript'}) :
            htmldoc.addRaw('function findAssociationsOnChange'+prefix+'(sel) {')
            htmldoc.addRaw('  var value = sel.value;')
            htmldoc.addRaw('  var names = '+str(names)+';')
            htmldoc.addRaw('  var tags = '+str(self.tags)+';')
            htmldoc.addRaw('  var meta = '+str(metaLookup)+';')
            htmldoc.addRaw('  var words = '+str(wordsLookup)+';')
            htmldoc.addRaw('  var keys = '+str(keys)+';')
            htmldoc.addRaw('  if (value == -1) {')
            htmldoc.addRaw('    draw'+prefix+'TreeGraphic();')
            htmldoc.addRaw('    var element = document.getElementById("pretree_'+prefix+'");')
            htmldoc.addRaw('    element.innerHTML="";')
            htmldoc.addRaw('  } else {')
            htmldoc.addRaw('    var tag = tags[value];')
            htmldoc.addRaw('    draw'+prefix+'TreeGraphic(meta[tag], words[tag]);')
            htmldoc.addRaw('    var element = document.getElementById("pretree_'+prefix+'");')
            htmldoc.addRaw('    element.innerHTML=\'<canvas id="treekey_'+prefix+'" width="100" height="100"></canvas>\';')
            htmldoc.addRaw('    var canvas_element = document.getElementById("treekey_'+prefix+'");')
            htmldoc.addRaw('    var canvas = canvas_element.getContext("2d");')
            htmldoc.addRaw('    canvas_element.width = ('+str(longestState)+'*8 + 28);')
            htmldoc.addRaw('    canvas_element.height = (keys[tag].length*20 + 28);')
            #htmldoc.addRaw('    canvas.scale(3.0,3.0);')
            htmldoc.addRaw('    canvas.font="12px Courier New";')
            htmldoc.addRaw('    var y = 12;')
            htmldoc.addRaw('    function addKey(key) {')
            htmldoc.addRaw('      canvas.fillStyle = key[0];')
            htmldoc.addRaw('      canvas.beginPath();')
            htmldoc.addRaw('      canvas.arc(12, y, 8, 0, 2*Math.PI, true);')
            htmldoc.addRaw('      canvas.fill();')
            htmldoc.addRaw('      canvas.fillStyle = key[1];')
            htmldoc.addRaw('      canvas.beginPath();')
            htmldoc.addRaw('      canvas.arc(12, y, 8, Math.PI, 2*Math.PI, true);')
            htmldoc.addRaw('      canvas.fill();')
            htmldoc.addRaw('      canvas.fillStyle = "#000000";')
            htmldoc.addRaw('      canvas.fillText(key[2], 24, y+4);')
            htmldoc.addRaw('    }')
            htmldoc.addRaw('    for(var i = 0; i < keys[tag].length; ++i) { addKey(keys[tag][i]); y += 20}')
            htmldoc.addRaw('  }')
            htmldoc.addRaw('}')

    def filteredTrees(self, subsets) :
        # Create lists of strains with different ST-types
        stcounts = Counter()
        for m in self.meta.values() :
            if 'ST complex' in m :
                stcounts[m['ST complex']] +=1

        if not stcounts :
            return subsets
        stlevels = sorted(n for n in stcounts)
        for level in stlevels :
            if stcounts[level] < 3 :
                continue
            subset = {'name' : level, 'numbers' : []}
            for n,m in self.meta.items() :
                if m['ST complex'] == level :
                    subset['numbers'].append(n)
            subsets.append(subset)
        return subsets