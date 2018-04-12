from Bio import Entrez
import os

Entrez.email = "jba9@leicester.ac.uk"

def main() :
    #search = searchAssembliesSpecial("Neisseria lactamica", ["meningitidis", "gonorrhoeae"])
    search=searchAssemblies('Campylobacter fetus', True)
    links = identifyLinks(search)
    downloadLinks("../PolyG/coli_fetus", links)

def searchAssemblies(organism, completeOnly = True) :
    searchText = '("' + organism + '"[Organism])'
    if completeOnly :
        searchText = searchText + 'AND ("complete genome"[Assembly Level])'
    searchHandle = Entrez.esearch(db="assembly", retmax = 250, term=searchText)
    return Entrez.read(searchHandle)

def searchAssembliesSpecial(organism, completeFor) :
#    searchText = '("'+ organism + '"[Organism] AND "complete genome"[Assembly Level])'
#    searchText += ' OR ("' + organism + '"[Organism] NOT '+' NOT '.join('"'+organism+' '+x+'"[Organism]' for x in completeFor)+')'
    searchText = '("'+ organism + '"[Organism] AND "scaffold"[Assembly Level])'
    #searchHandle = Entrez.esearch(db="assembly", retmax = 250, term=searchText)
    searchHandle = Entrez.esearch(db="assembly", retmax = 1, term=searchText)
    return Entrez.read(searchHandle)

def identifyLinks(searchRecord) :
    linkHandle = Entrez.elink(dbfrom="assembly", db="nucleotide", id=searchRecord['IdList'])
    linkRecord = Entrez.read(linkHandle)
    return linkRecord

def downloadLinks(intoDir, linkRecord) :
    # Create the directory
    if not os.path.exists(intoDir):
        os.makedirs(intoDir)

    # Create download list
    downloadIDs = []
    downloadParents = []
    for link in linkRecord :
        id = link['IdList'][0]
        linkSetDb = link['LinkSetDb']
        for linkSet in linkSetDb :
            if linkSet['LinkName'] == 'assembly_nuccore_insdc' :
                for l in linkSet['Link'] :
                    downloadIDs.append(l['Id'])
                    downloadParents.append(id)
                break
        else :
            print("Error: no INSDC assembly found for ", id)

    for i, id in enumerate(downloadIDs) :
        with open(os.path.join(intoDir, downloadParents[i] + "-" + id + ".gb"), "w") as gbOut :
            print('Downloading '+str(i+1)+'/'+str(len(downloadIDs))+' '+downloadParents[i] + '-' + id)
            fetchHandle = Entrez.efetch(db="nuccore", id=id, rettype="gb", retmode="text")
            gbOut.write(fetchHandle.read())
            fetchHandle.close()

# Remove auto-run
if __name__ == "__main__":
    main()