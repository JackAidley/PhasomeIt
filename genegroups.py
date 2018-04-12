# Functions for annotating gene groups
from collections import Counter
import math
import logging

def createGroupAnnotation(db) :
    '''Add functional annotation to this gene group'''
    annotation = {}
    for strain in db :
        for contig in strain['contigs'] :
            for tract in contig['tracts'] :
                if tract['geneGroup'] not in annotation :
                    annotation[tract['geneGroup']] = {}

    for geneGroup in annotation :
        pvFunctions = []
        otherFunctions = []
        firstLocus = None
        firstGene = None
        for strain in db :
            for contig in strain['contigs'] :
                for tract in contig['tracts'] :
                    if tract['geneGroup'] != geneGroup :
                        continue

                    # Find a name to use
                    if not firstLocus and ('autoLocus' not in tract or not tract['autoLocus']) :
                        firstLocus = tract['locus']

                    if not firstGene and tract['gene'] != 'No Data' and ('autoLocus' not in tract or not tract['autoLocus']) :
                        firstGene = tract['gene']

                    pvFunctions.append(tract['function'])
                    if 'blastMatch' not in tract or tract['blastMatch'] == 'No Data' :
                        continue

                    for match in tract['blastMatch'] :
                        # Matching PV genes will be picked up elsewhere
                        if 'tractNo' in match :
                            continue
                        if 'function' in match :
                            otherFunctions.append(match['function'])
        annotation[geneGroup] = findCommonFunctions(pvFunctions, otherFunctions)
        annotation[geneGroup]['name'] = firstGene if firstGene else firstLocus if firstLocus else '#'+str(geneGroup)

    return annotation

def findCommonFunctions(pvFunctions, otherFunctions) :
    pvWords = findFunctions(pvFunctions)

    if (not pvWords) :
        return { 'bestFunction' : 'Unknown', 'words' : None }

    # If we don't have any very good hits, add in the other function keywords
    if pvWords[0]['score'] < len(pvFunctions)*0.25 :
        others = findFunctions(otherFunctions)
        if (others) :
            pvWords.append(others)
            pvWords.sort(key = lambda x : x['score'], reverse=True)

    ret = {}
    ret['bestFunction'] = pvWords[0]['function']
    ret['words'] = pvWords

    return ret

def findFunctions(inPhrases) :
    uselessWords = ['conserved', 'hypothetical', 'putative', 'protein']
    uselessPhrases = ['No data, identified by script', 'No annotation data']

    phrases = []

    for phrase in inPhrases :
        description = {
            'words' : [],
            'function' : phrase,
            'links' : [0]*len(inPhrases)
        }
        if phrase not in uselessPhrases :
            words = phrase.lower().split(' ')
            for word in words :
                if word not in uselessWords :
                    description['words'].append(word)
        phrases.append(description)

    for phrase in phrases :
        for i, other in enumerate(phrases) :
            phrase['links'][i] = countInBoth(phrase['words'], other['words'])
        phrase['score'] = sum(phrase['links'])/(5+len(phrase['words']))

    phrases.sort(key = lambda x : x['score'], reverse=True)

    return phrases if phrases else None

def countInBoth(a, b) :
    return sum((x in b) for x in a)

def addSubgroupHomologies(db) :
    ''' Further subgroup the homology groups into more similar groupings for larger groups '''
    groups = {}

    # Collect the gene groups
    for strainNo, strain in enumerate(db) :
        for contig in strain['contigs'] :
            for n, tract in enumerate(contig['tracts']) :
                geneGroup = tract['geneGroup']
                if geneGroup == -1 :
                    continue
                if geneGroup not in groups :
                    groups[geneGroup] = []
                nth = len(groups)
                groups[geneGroup].append({
                    'tract' : tract,
                    'strainNo' : strainNo,
                    'contig' : contig['record'].description,
                    'tractNo' : n})

    # Subgroup each group
    print('Subgrouping',end='',flush=True)
    for group in groups.values() :
        try:
            subgroupGroup(group)
            print('.', end='',flush=True)
        except Exception :
            # The subgrouping failed, but is not vital so carry on.
            logging.exception('Subgrouping failed for group "' + str(group) + '"')
            print('!', end='', flush=True)

    print('')

def subgroupGroup(group) :
    '''Find the subgroups within a group'''
    # If the group is small then no meaningful subgroup analysis can be done
    if len(group) < 5 :
        return

    # We DEFINE a subgroup as a group in which all member are connect by higher blast homology to each other than every link out of the group
    # i.e. For any two genes in a subgroup the greed algorithm for finding the shortest path between them does not leave the subgroup
    # And then SEEK a set of subgroups where the number of subgroups is helpful in analysing the group

    # Create a [contig][n]->tract map
    tractLookup = {}
    for nth, gene in enumerate(group) :
        if gene['contig'] not in tractLookup :
            tractLookup[gene['contig']] = {}
        tractLookup[ gene['contig'] ][ gene['tractNo'] ] = gene
        gene['distance'] = [0]*len(group)
        gene['nth'] = nth

    # Create a distance map for each
    for nth, gene in enumerate(group) :
        tract = gene['tract']
        if 'blastMatch' not in tract :
            continue
        for contigName, matches in tract['blastMatch'].items() :
            # Can have non-PV tracts in contigs we haven't indexed
            if contigName not in tractLookup :
                continue
            contigLookup = tractLookup[contigName]
            for match in matches :
                if 'tractNo' not in match :
                    continue
                distance = matchToDistance(match)
                matchedGene = contigLookup[ match['tractNo'] ]
                gene['distance'][ matchedGene['nth'] ] = max(distance,  gene['distance'][ matchedGene['nth'] ])
                matchedGene['distance'][nth] = max(distance,  matchedGene['distance'][nth])

    # Create a list of all distances in the distance maps
    allDistances = []
    for gene in group :
        for distance in gene['distance'] :
            if distance and distance not in allDistances :
                allDistances.append(distance)

    allDistances.sort()

    # Then find the largest "breakpoints" in those distances and consider them as a potential starting points for fragmenting the network
    steps = []
    for i in range(1, len(allDistances)) :
        steps.append((allDistances[i] - allDistances[i-1], i))
    steps.sort(reverse=True, key = lambda x : x[0])

    # Reduce the number of steps we look at
    if len(steps) > 100 :
        steps = steps[:100] + steps[100::2]
    if len(steps) > 200 :
        steps = steps[:200] + steps[200::2]
    if len(steps) > 300 :
        steps = steps[:300] + steps[300::4]
    if len(steps) > 500 :
        steps = steps[:500]

    # Then step down the list of largest breakpoints and see whether the related cutoffs produce a desirable breakdown of our genes
    topCutoff = allDistances[-1]
    bottomCutoff = allDistances[0]
    subgroupings = []
    for _, n in steps :
        cutoff = allDistances[n]
        # Can ignore this, will produce a fully connected graph
        if cutoff < bottomCutoff :
            continue
        if cutoff > topCutoff :
            continue
        subgroupData = subgroupsAtCutoff(group, cutoff)
        if subgroupData['numSubs'] == 1 :
            # all distances lower than this can be safely ignored
            # as they will also complete the network
            if cutoff > bottomCutoff :
                bottomCutoff = cutoff
            continue
        if subgroupData['numSubs'] >= math.ceil(len(group)/2) :
            # Too many groups, ignore this kind of cruft
            # Higher cutoffs will produce similar results
            if cutoff < topCutoff :
                topCutoff = cutoff
            continue
        # Don't add the same subgrouping multiple times
        if subgroupData['mapping'] in (s['mapping'] for s in subgroupings) :
            continue
        subgroupings.append(subgroupData)

    # No valid subgroupings found
    if not subgroupings :
        return

    # Take the ten best and refine them
    subgroupings.sort(key = lambda x : x['score'], reverse = True)
    subgroupings = subgroupings[:max(10, len(group)//10)]

    bestCutoffs = [x['cutoff'] for x in subgroupings]
    bestCutoffs.sort(reverse=True)

    for n, subgroupData in enumerate(subgroupings) :
        refined = refineSubgroup(subgroupData, [x for x in bestCutoffs if x < subgroupData['cutoff']], group)
        if refined['score'] > subgroupData['score'] :
            subgroupings[n] = refined

    subgroupings.sort(key = lambda x : x['score'], reverse = True)

    # Add the subgroup numbers to the tracts themselves
    checkValidSubgrouping(subgroupings[0]['mapping'], group)
    for n, inSub in enumerate(subgroupings[0]['mapping']) :
        group[n]['tract']['subgroup'] = inSub

def matchToDistance(match) :
    return min(match['coverage'], match['queryCoverage'])

def subgroupsAtCutoff(group, cutoff) :
    assigned = [-1]*len(group)
    onSubgroup = 0
    while True :
        # Find first unassigned
        on = 0
        while on < len(group) and assigned[on] != -1 :
            on += 1

        # All have been assigned, time to return
        if on == len(group) :
            score, subgroupScores = scoreSubgroups(assigned, group)
            return {
                'cutoff' : cutoff,
                'score' : score,
                'subgroupScores' : subgroupScores,
                'mapping' : assigned,
                'numSubs' : onSubgroup
            }


        # Spider from this one assigning the group number
        targets = set([on])
        while len(targets) != 0 :
            # Get the first one from the list
            target = targets.pop()

            if assigned[target] == -1 :
                assigned[target] = onSubgroup
                gene = group[target]
                targets = targets.union(i for i, d in enumerate( gene['distance'] ) if d >= cutoff and assigned[i] != onSubgroup)
            elif assigned[target] != onSubgroup :
                # Something went horribly, horribly wrong
                print("Subgroup homology analysis suffered catastrophic algorithmic failure. Blame (and tell!) the programmer.")
                assigned[target] = onSubgroup

        onSubgroup += 1

def scoreSubgroups(subgroups, group) :
    # We want to score this in a way that scores subgroupings that favours groupings that divide the tracts on a single strain
    strainSubgroups = {}

    for geneNo, subgroup in enumerate(subgroups) :
        gene = group[geneNo]
        strainNo = gene['strainNo']
        if strainNo not in strainSubgroups :
            strainSubgroups[strainNo] = Counter()
        strainSubgroups[strainNo][subgroup] += 1

    # Add points for
    #  1. Each strain a subgroup appears in

    # Remove points for
    #  1. Any time it appears multiple times in a strain
    #  2. Any subgroup that contains only one

    score = 0
    subgroupCount = Counter(subgroups)

    subgroupScores = [0]*len(subgroupCount)

    for counts in strainSubgroups.values() :
        for n,count in counts.items() :
            if count > 0 :
                subgroupScores[n] += 1-0.5*count

    for n, count in subgroupCount.items() :
        if count == 1 :
            subgroupScores[n] = -1

    # We prefer a moderate number of high scoring subgroups to a large number of lower scoring groups
    # So apply a modifier for the number of quality groups
    score = sum(x for x in subgroupScores)
    ratio = len([x for x in subgroupCount.values() if x > 1])/(len(group)/len(strainSubgroups))
    if ratio > 1 and score > 0 :
        score *= 1.0/math.sqrt(ratio)

    return score, subgroupScores

def refineSubgroup(subgroup, tryCutoffs, group) :
    if not tryCutoffs :
        return subgroup

    # What we're going to is "freeze" the high quality subgroups we got from the higher cutoff and try the lower cutoffs with the other subgroups to see if we can gather them up without breaking the "good" subgroups
    subgroupScores = subgroup['subgroupScores']
    assigned = list(subgroup['mapping'])

    best = subgroup
    for newCutoff in tryCutoffs :
        for n, subgroup in enumerate(assigned) :
            if subgroupScores[subgroup] <= 1.0 :
                assigned[n] = -1

        assigned = collapseMissing(assigned)

        while True :
            # Find first unassigned
            on = 0
            while on < len(assigned) and assigned[on] != -1 :
                on += 1

            # All have been assigned, time to return
            if on == len(assigned) :
                break

            # Spider from this one assigning the group number
            targets = [on]
            setToSubgroup = max(assigned)+1
            while len(targets) != 0 :
                # Get the first one from the list
                target = targets[0]
                targets = targets[1:]

                if assigned[target] != setToSubgroup :
                    if assigned[target] != -1 :
                        # We've collided with another subgroup remove it
                        assigned = [x if x != assigned[target] else -1 for x in assigned]
                    assigned[target] = setToSubgroup
                    gene = group[target]
                    for i, distance in enumerate( gene['distance'] ) :
                        if distance >= newCutoff and assigned[i] != setToSubgroup:
                            targets.append(i)

        assigned = collapseMissing(assigned)

        # It's entirely imploded
        if max(assigned) == 0 :
            return best

        score, subgroupScores = scoreSubgroups(assigned, group)

        if score > best['score'] :
            best = {
                'score' : score,
                'subgroupScores' : list(subgroupScores),
                'mapping' : list(assigned),
                'numSubs' : max(assigned)+1
            }

    return best

def collapseMissing(assigned) :
    map = {}
    on = 0
    for i in range(0, len(assigned)) :
        a = assigned[i]
        if a !=-1 :
            if a not in map :
                map[a] = on
                on += 1
            assigned[i] = map[a]

    return assigned

def checkValidSubgrouping(assigned, group) :
    # For all subgroups, there must be a path to every other element of the subgroup that
    # does not use links longer than the shortest link outside the group
    subgroups = set(assigned)

    on = 0
    invalid = False
    for subgroup in subgroups :
        # Find longest outside link
        maxOutside = 0
        for n, sub in enumerate(assigned) :
            if sub != subgroup :
                continue
            on = n
            gene = group[n]
            for i, distance in enumerate(gene['distance']) :
                if assigned[i] != subgroup and distance > maxOutside :
                    maxOutside = distance

        # Test the walk
        test = [x if x != subgroup else -1 for x in assigned]
        # Spider from this one assigning the group number
        targets = [on]
        while len(targets) != 0 :
            # Get the first one from the list
            target = targets[0]
            targets = targets[1:]

            if test[target] != subgroup :
                if test[target] != -1 :
                    # We've collided with another subgroup, our subgrouping is invalid
                    invalid = True
                test[target] = subgroup
                gene = group[target]
                for i, distance in enumerate( gene['distance'] ) :
                    if distance > maxOutside and test[i] != subgroup:
                        targets.append(i)

        if test.count(-1) :
            invalid = True

    if invalid :
        print("WARNING: Invalid subgroup generated!")

def hasOneOutlier(assigned) :
    c = Counter(assigned)
    if len(c) != 2 :
        return False
    return c[0] == 1 or c[1] == 1
