import csv
import logging

class SpeciesAnnotator :
    def __init__(self, csvFilename, idColumn = 0) :
        self.lookup = {}
        if csvFilename == None :
            return
        try :
            with open(csvFilename, newline='') as csvFile:
                reader = csv.reader(csvFile, dialect='excel')
                header = next(reader)
                speciesCol = next((i for i,x in enumerate(header) if x.lower() == 'species'), None)
                if speciesCol == None:
                    return
                for row in reader :
                    self.lookup[row[idColumn]] = row[speciesCol]
        except FileNotFoundError :
            logging.warning('Annotation data file "'+csvFilename+'" not found. Falling back to species annotator')

    def annotateDb(self, db) :
        if not self.lookup :
            return
        for strain in db :
            strainName = strain['name']
            id = strainName.split('.')[0]
            if 'species' not in strain or strain['species'] == 'Unknown sp.' :
                if id in self.lookup :
                    strain['species'] = self.lookup[id]
