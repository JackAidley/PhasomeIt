# Command line parsing main entry file

import argparse
import logging
import os
import sys
import time
import loghandler

import settings
import findpvinstrains
from speciesannotator import SpeciesAnnotator

__version__ = '1.1.1'

def main() :
    '''Run phasomeit using command line arguments'''

    # Create the parser
    parser = argparse.ArgumentParser(prog='PhasomeIt', description='Identify the phasome of multiple genomes')

    # Target path
    parser.add_argument(
        '-t', '--target',
        metavar='target-folder',
        help='Folder containing all target genomes in genbank format with file extensions of .gb, .gbf, or .gbk',
        required=True)

    # Cutoffs for lengths
    parser.add_argument(
        '-c', '--cutoffs',
        metavar='repeat',
        nargs='+', # Must have at least one but can have many
        help='Cutoffs for each length, starting with 1bp repeats and going up. Enter 0 to ignore a length. Maximum of 9 allowed',
        required=True)

    parser.add_argument(
        '-f', '--filter',
        metavar='filter',
        help='Filter for types, only one accepted for now',
        required=False)

    parser.add_argument(
        '--homology',
        metavar='homology [homology]',
        nargs='+',
        help='Cutoffs for homology groups, one or two values. The first is the required match from the query to the subject, the second is the inverse. This will be calculated as 0.8*the first number if omitted. Default 0.5 0.4',
        default=[0.5, 0.4],
        type=float,
        required=False)

    parser.add_argument(
        '--firstSpecies',
        metavar='species-name',
        nargs='+',
        help='Species to place at the start of the tables',
        default=[],
        required=False)

    parser.add_argument(
        '--firstStrains',
        metavar='strain-name',
        nargs='+',
        help='Strains to place at the start of the tables, has lower precedent than any species ordering defined by firstSpecies',
        default=[],
        required=False)

    parser.add_argument(
        '--subgroupOff',
        help='Disable subgrouping. Speeds up the process and provides cleaner output.',
        action='store_false',
        default=True,
        required=False)

    parser.add_argument(
        '--numberTrees',
        help='Instead of putting names on trees, put numbers. Creates trees which are easier to read when printed.',
        action='store_true',
        default=False,
        required=False)

    parser.add_argument(
        '--metadata',
        help='The name of a csv file in the target directory containing species names in the first column and data to annotate with in the other columns. Use "--metadata-columns" to use only a subset of columns.',
        metavar='metadata-csv-filename',
        default=None,
        required=False)

    parser.add_argument(
        '--metadataColumns',
        help='One or more columns in the metadata file to use, columns are numbered from zero (default: use all columns)',
        metavar='metadata-columns',
        nargs='+',
        type=int,
        default=None,
        required=False)


    parser.add_argument(
        '--version',
        action='version',
        version = '%(prog)s '+__version__)


    args = parser.parse_args()

    # Check arguments
    if not os.path.isdir(args.target) :
        if os.path.exists(args.target) :
            print('ERROR: Specified target folder ('+args.target+') is not a folder')
        else :
            print('ERROR: Specified target folder ('+args.target+') does not exist')
        sys.exit(1)

    # Check a valid number of cutoffs
    if len(args.cutoffs) > 9 :
        print('ERROR: Only repeats up to a maximum of 9bp are supported.')
        sys.exit(1)

    # Check homology values are valid
    if len(args.homology) > 2 :
        print('ERROR: Only two homology values are accepted.')
        sys.exit(1)

    for n in args.homology :
        if n <0 or n > 1 :
            print('ERROR: homology values must be between 0 and 1')
            sys.exit(1)

    # Transfer homology values to settings
    settings.homologyCutoffTo = args.homology[0]
    settings.homologyCutoffFrom = args.homology[1] if len(args.homology) == 2 else args.homology[0]*0.8

    # Transfer priority species to settings
    settings.firstSpecies = args.firstSpecies
    settings.firstStrains = args.firstStrains

    # Transfer switches to settings
    settings.subgroupOff = args.subgroupOff
    settings.numberTrees = args.numberTrees

    # Transfer metadata settings to settings
    settings.metadataCsv = args.metadata
    settings.metadataColumns = args.metadataColumns

    t=time.clock()

    # Start the logger and add basic information
    findpvinstrains.startLog(args.target)
    logging.info('PhasomeIt '+ __version__+ ' run from command line')
    logging.info(' '.join(sys.argv))

    # Run the tool
    plugins = [SpeciesAnnotator(args.target, settings.metadataCsv)]
    findpvinstrains.findPVinFolder(args.target, args.filter, args.cutoffs, plugins)
    plugins.append(findpvinstrains.findAssociations(args.target))
    findpvinstrains.reHTMLPickle(args.target, plugins)

    print("Time taken: "+str(time.clock()-t)+"s")


    # Check what happened during the run and report if anything bad occured
    if loghandler.worstSeverity > logging.INFO :
        print('\nBEWARE: There were '+logging.getLevelName(loghandler.worstSeverity)+'(s) during the run, check the "_phasomeit.log" file for details!')
    else :
        print('\nPhasomeIt analysis completed without problems.')



if __name__ == '__main__' :
    main()
