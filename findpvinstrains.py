# Find all the PV tracts in all the strains, and output them to excel
from genomecrawler import crawlGenomes
from htmlify import outputAsHTMLFromPickle
from findassociations import AssociationFinder
import blastmatch
import genegroups
import settings
from speciesannotator import SpeciesAnnotator
import logging
import loghandler

import os
import time
import pickle

from Bio.Seq import Seq

def main() :
    t=time.clock()

    #targetPath = '../polyG/c.jejuni'
    #targetPath = '../polyG/campy.new'
    #targetPath = '../PolyG/Coli_Fetus'
    #targetPath = '../polyG/campy_spp.plasmids'
    #targetPath = '../polyG/c.jejuni_test'
    #targetPath = '../polyG/c.jejuni_doylei'
    #targetPath = '../polyG/ProkkaGenomes'
    #targetPath = '../polyG/MenW'
    #targetPath = '../polyG/c.jejuni_plasmids'
    #targetPath = '../polyG/h.pylori'
    #targetPath = '../polyG/c.jejuni_six'
    #targetPath = '../polyG/c.jejuni_samtest'
    #targetPath = 'strain11168andPT14'
    #targetPath= 'strainPT14'
    #targetPath = "HappyMegan"
    #targetPath = '../polyG/c.jejuni_sam'
    #targetPath = 'c.jejuni_singlestrain'
    #targetPath = '../polyG/c.concisus'

    #targetPath = '../polyG/Haemophilus'
    #targetPath = '../polyG/Neisseria.both'
    #targetPath = '../polyG/Neisseria.test'
    #targetPath= '../polyG/MenTest'

    targetPath = '../polyG/all_st21'

    startLog(targetPath)
    logging.info('New run started from script')

    # Transfer priority species to settings
    settings.firstSpecies = ['C. jejuni']
    settings.firstStrains = ['NCTC 11168', '81-176', 'PT14', '15-537360']

    plugins = [SpeciesAnnotator(targetPath+'/_metadata.csv')]
    #findPVinFolder(targetPath, 'W10', [9, 15, 6, 5, 5, 3, 3, 3, 3], plugins) # for Neisseria
    #findPVinFolder(targetPath, 'W10', [9, 6, 10, 5, 5], plugins) # for others
    #findPVinFolder(targetPath, 'W9', [7, 6, 10, 5, 5], plugins) # For campy
    #findPVinFolder(targetPath, '', [9], plugins) # For all_st21/all_st45 test
    #reExcelPickle(targetPath, 'summary_pickle.pkl')
    #splitPickle(targetPath, 'summary_pickle.pkl')
    plugins.append(findAssociations(targetPath))
    reHTMLPickle(targetPath, plugins)

    print("Time taken: "+str(time.clock()-t)+"s")

def startLog(targetPath) :
    ''' Wrapped only to avoid duplication '''
    logging.basicConfig(
        filename=targetPath + '/_phasomeit.log',
        format = '%(asctime)s [%(levelname)s] %(name)s : %(message)s',
        filemode='w',
        level=logging.DEBUG)
    logger = logging.getLogger()
    logger.addHandler(loghandler.TrackLogSeverity())


def findPVinFolder(targetPath, filter, minLengths, plugins = None) :
    t=time.clock()

    # Create the temp directory
    if not os.path.exists(os.path.join(targetPath, 'temp')):
        os.makedirs(os.path.join(targetPath, 'temp'))

    # Find all the SSRs
    db = crawlGenomes(targetPath, filter, minLengths, plugins)

    print("Crawl time: "+str(time.clock()-t)+"s")

    t=time.clock()
    blastmatch.blastComparePVs(db, targetPath)

    # Strip everything not still in use from the data
    for strain in db :
        for contig in strain['contigs'] :
            contig['record'].seq = Seq('', contig['record'].seq.alphabet)
            contig['record'].features = None
            contig['record'].featureStarts = None
            contig['record'].lookup = None
            for tract in contig['tracts'] :
                tract.pop('record', None)
                tract.pop('sequence', None)

    if not settings.subgroupOff :
        genegroups.addSubgroupHomologies(db)
    print("Total blast analysis time: "+str(time.clock()-t)+"s")

    geneGroups = genegroups.createGroupAnnotation(db)

    print('Pickling', end='', flush=True)
    t=time.clock()

    # Store the raw data so we can re-process without having to do the heavy lifting
    # create pickle directory
    if not os.path.exists(targetPath+'/pickle') :
        os.makedirs(targetPath+'/pickle')

    # Create pickle guide
    guide = {
        'numStrains' : len(db),
        'strainNames' : [s['name'] for s in db],
        'species' : [s['species'] for s in db],
        'displayNames' : [s['species']+' '+s['name'] for s in db]}

    strains = []
    for strain in db :
        strains.append({'name' : strain['name'], 'contigs':[c['record'].description for c in strain['contigs']]})

    guide['strains']=strains

    with open(os.path.join(targetPath, 'pickle', 'guide.pkl'), 'wb') as file :
        pickle.dump(guide, file)

    # Then create pickle for each strain
    for i, strain in enumerate(db) :
        with open(os.path.join(targetPath, 'pickle', str(i)+'.pkl'), 'wb') as file :
            pickle.dump(strain, file)
        print('.', end='', flush=True)

    # And create the annotation pickle
    with open(os.path.join(targetPath, 'pickle', 'groups.pkl'), 'wb') as file :
        pickle.dump(geneGroups, file)

    print('!')

    print("Total Pickling time: "+str(time.clock()-t)+"s")

    # HTML output
    t=time.clock()
    #outputAsHTML(db, os.path.join(targetPath, 'summary_tracts'))
    #print("HTML Output time: "+str(time.clock()-t)+"s")

def splitPickle(targetPath, bigPickle) :
    print("Loading Pickle")
    t=time.clock()
    with open(os.path.join(targetPath, bigPickle), 'rb') as file :
        db = pickle.load(file)
    print("Time taken to load: "+str(time.clock()-t)+"s")

    # Strip everything not still in use from the data
    for strain in db :
        for contig in strain['contigs'] :
            contig['record'].seq = Seq('', contig['record'].seq.alphabet)
            contig['record'].features = None
            contig['record'].featureStarts = None
            contig['record'].lookup = None
            for tract in contig['tracts'] :
                tract.pop('record', None)
                tract.pop('sequence', None)

    print("Splitting Pickle", end='', flush=True)

    t=time.clock()
    # create pickle directory
    if not os.path.exists(targetPath+'/pickle') :
        os.makedirs(targetPath+'/pickle')

    # Create pickle guide
    guide = {
        'numStrains' : len(db),
        'strainNames' : [s['name'] for s in db],
        'species' : [s['species'] for s in db],
        'displayNames' : [s['species']+' '+s['name'] for s in db]}

    strains = []
    for strain in db :
        strains.append({'name' : strain['name'], 'contigs':[c['record'].description for c in strain['contigs']]})

    guide['strains']=strains

    with open(os.path.join(targetPath, 'pickle', 'guide.pkl'), 'wb') as file :
        pickle.dump(guide, file)

    # Then create pickle for each strain
    for i, strain in enumerate(db) :
        with open(os.path.join(targetPath, 'pickle', str(i)+'.pkl'), 'wb') as file :
            pickle.dump(strain, file)
        print('.', end='', flush=True)
    print('')
    print("Time taken to split: "+str(time.clock()-t)+"s")

def findAssociations(targetPath) :
    # Find associations
    #print("Finding associations")
    t = time.clock()
    finder = AssociationFinder()
    if settings.metadataCsv :
        finder.readMetadata(targetPath, os.path.join(targetPath, settings.metadataCsv), 0, settings.metadataColumns)
    else :
        finder.addSpeciesAsMetadata(targetPath)
    #finder.findAssociations(targetPath)
    #print("Association finding time: "+str(time.clock()-t)+"s")
    return finder

def reHTMLPickle(targetPath, plugins) :
    # print("Loading Pickle")
    # t=time.clock()
    # with open(os.path.join(targetPath, filename), 'rb') as file :
    #     db = pickle.load(file)
    # print("Time taken to load: "+str(time.clock()-t)+"s")


    # HTML output
    print("Outputting html")
    t=time.clock()
    outputAsHTMLFromPickle(os.path.join(targetPath, 'summary_tracts'), plugins)
    print("HTML Output time: "+str(time.clock()-t)+"s")

if __name__ == "__main__":
    main()