"""Functions to identify features at a location in a genome"""
from bisect import bisect_left

def findMatchingFeatures(record, location, types) :
    for type in types :
        features = findMatchingFeatureOfType(record, location, type)
        if features :
            return features

    return []

def findMatchingFeatureOfType(record, location, type) :
    """Find and return all features that span the given location"""
    if type not in record.lookup :
        return []

    lookup = record.lookup[type]
    features = record.features

    # bisect will find first record it can possible be inside, scan right to find any others that contain it
    lo = lookup[0]
    hi = lookup[1]
    # Note the '-1' at the end is required because bisect_left returns the first element larger than ours.
    n = max(0, bisect_left(record.featureStarts, location, lo=lo, hi=hi)-1)

    ret = []
    for feature in features[n:lookup[1]] :
        if feature.location.start > location :
            break
        if location in feature :
            ret.append(feature)

    return ret


def findClosestFeature(record, location, types, directional_bias = 0) :
    # Find the closest for each type
    found = [(findClosestFeatureOfType(record, location, t, directional_bias), n) for n,t in enumerate(types)]
    # Sort the by how close they are, but maintain preference order in case of a tie
    found.sort(key=lambda x : (abs(x[0][0]), x[1]))
    return found[0][0]

def findClosestFeatureOfType(record, location, type, directional_bias = 0) :
    """ Find the closest* feature to a point (assumes not in a feature)
        * - after applying a bias towards features facing away from the point
    """
    if type not in record.lookup :
           return (1e9, None)

    lookup = record.lookup[type]
    features = record.features

    # bisect will find first record it can possible be inside, scan right to find any others that contain it
    lo = lookup[0]
    hi = lookup[1]
    # Note the '-1' at the end is required because bisect_left returns the first element larger than ours.
    n = max(0, bisect_left(record.featureStarts, location, lo=lo, hi=hi)-1)

    closest = None
    distance = len(record.seq)

    for i in range(max(n-2, lo), min(n+3, hi)) :
        f = features[i]

        start = f.location.start
        end = f.location.end
        if (location < start):
            fdist = start - location
            if f.strand == -1 :
                fdist += directional_bias
        else:
            fdist = location - end
            if f.strand == 1 :
                fdist += directional_bias

        if (fdist < distance):
            closest = f
            distance = fdist

    return (distance, closest)