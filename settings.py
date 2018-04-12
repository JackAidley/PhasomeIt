# Global settings and set functions

# Turn off the normal output, and only show errors
silent = False

# Cutoff for homology
homologyCutoffTo = 0.5
homologyCutoffFrom = 0.4

# Strains to put first
firstSpecies = []
firstStrains = []

# Switches
subgroupOff = False
numberTrees = False

# metadata
metadataCsv = None
metadataColumns = None

# Function that prints or not based on the state of silent
def output(*args, **kwargs) :
    if not silent :
        print(*args, **kwargs)
