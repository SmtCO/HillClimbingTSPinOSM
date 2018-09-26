'''
AI_OSM
'''
import fiona
fiona.supported_drivers
import ogr
import mplleaflet
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.path import Path
from shapely.geometry import *
from datetime import datetime
from math import sqrt
import random


tstart = datetime.now()
print tstart

def hillclimb(init_function, move_operator, objective_function, max_evaluations):
    # Hillclimb until either max_evaluations is reached or we are at a local optima

    best = init_function()
    best_score = objective_function(best)

    num_evaluations = 1

    while num_evaluations < max_evaluations:
        # Examine moves around our current position
        move_made = False
        for next in move_operator(best):
            if num_evaluations >= max_evaluations:
                break

            # See if this move is better than the current
            next_score = objective_function(next)
            num_evaluations += 1
            if next_score > best_score:
                best = next
                best_score = next_score
                move_made = True
                break  # depth first search

        if not move_made:
            break  # couldn't find a better move (must be at a local maximum)

    return (num_evaluations, best_score, best)


def hillclimb_and_restart(init_function, move_operator, objective_function, max_evaluations):
    # Repeatedly hillclimb until max_evaluations is reached

    best = None
    best_score = 0

    num_evaluations = 0

    while num_evaluations < max_evaluations:
        remaining_evaluations = max_evaluations - num_evaluations

        evaluated, score, found = hillclimb(init_function, move_operator, objective_function, remaining_evaluations)

        num_evaluations += evaluated
        if score > best_score or best is None:
            best_score = score
            best = found

    return (num_evaluations, best_score, best)


def reversed_sections(tour):
    # Generator to return all possible variations where the section between two cities are swapped
    for i, j in all_pairs(len(tour)):
        if i != j:
            copy = tour[:]
            if i < j:
                copy[i:j + 1] = reversed(tour[i:j + 1])
            else:
                copy[i + 1:] = reversed(tour[:j])
                copy[:j] = reversed(tour[i + 1:])
            if copy != tour:  # no point returning the same tour
                if copy[0] == 0:  # tour must start from location 0
                    yield copy
                elif copy[-1] == 0:
                    yield copy[::-1]


def swapped_locations(tour):
    # Generator to create all possible variations where two cities have been swapped
    for i, j in all_pairs(len(tour)):
        if i < j and i != 0:
            copy = tour[:]
            copy[i], copy[j] = tour[j], tour[i]
            yield copy


def all_pairs(size, shuffle=random.shuffle):
    # Generates all i,j pairs for i,j from 0-size uses shuffle to randomise (if provided)
    r1 = range(0, size)
    r2 = range(0, size)
    if shuffle:
        shuffle(r1)
        shuffle(r2)
    for i in r1:
        for j in r2:
            if j != 0:
                yield (i, j)

def cartesian_matrix(coords):
    # Create a distance matrix for the city coords that uses straight line distance
    matrix = {}
    for i, (x1, y1) in enumerate(coords):
        for j, (x2, y2) in enumerate(coords):
            dx, dy = x1 - x2, y1 - y2
            dist = sqrt(dx * dx + dy * dy)
            matrix[i, j] = dist
    return matrix


def tour_length(matrix, tour):
    # Total up the total length of the tour based on the distance matrix
    total = 0
    num_locations = len(tour)
    for i in range(num_locations):
        j = (i + 1) % num_locations
        loc_i = tour[i]
        loc_j = tour[j]
        if j != 0:  # travel route is acyclic
            total += matrix[loc_i, loc_j]
    return total


def init_random_tour(tour_length):
    tour = range(1, tour_length) # tour=range(tour_length)
    random.shuffle(tour)
    return [0] + tour

# Read OSM
driver = ogr.GetDriverByName('OSM')
data = driver.Open('mapsunumroma.osm')
# or for Ankara
# data = driver.Open('mapsunumankara.osm')
layer = data.GetLayer('points')

features = [x for x in layer]
print len(features)

data_list = []
coordsmaps = []
coordsTSP = []
coordsfile = open('coords.txt', 'w')

for feature in features:
    data = feature.ExportToJson(as_object=True)
    coords = data['geometry']['coordinates']
    shapely_geo = Point(coords[0], coords[1])
    name = data['properties']['name']
    highway = data['properties']['highway']
    other_tags = data['properties']['other_tags']

    if other_tags and ('amenity') in other_tags:
        feat = [x for x in other_tags.split(',') if 'amenity' in x][0]
        amenity = feat[feat.rfind('>') + 2:feat.rfind('"')]
        if (amenity == 'place_of_worship'):
            coordsmaps.append([name, coords[0], coords[1]])
            coordsTSP.append([float(coords[0]), float(coords[1])])
            coordsfile.write("%s%s %f%s %f\n" % (name, ',', coords[0], ',', coords[1]))
    else:
        amenity = None
    data_list.append([name, highway, amenity, shapely_geo])

    if other_tags and ('tourism') in other_tags:
        feat = [x for x in other_tags.split(',') if 'tourism' in x][0]
        tourism = feat[feat.rfind('>') + 2:feat.rfind('"')]

        if (tourism == 'museum' or tourism == 'artwork'):
            coordsmaps.append([name, coords[0], coords[1]])
            coordsTSP.append([float(coords[0]), float(coords[1])])
            coordsfile.write("%s%s %f%s %f\n" % (name, ',', coords[0], ',', coords[1]))
    else:
        tourism = None
    data_list.append([name, highway, tourism, shapely_geo])

    if other_tags and ('leisure') in other_tags:
        feat = [x for x in other_tags.split(',') if 'leisure' in x][0]
        leisure = feat[feat.rfind('>') + 2:feat.rfind('"')]

        if (leisure == 'theme_park'):
            coordsmaps.append([name, coords[0], coords[1]])
            coordsTSP.append([float(coords[0]), float(coords[1])])
            coordsfile.write("%s%s %f%s %f\n" % (name, ',', coords[0], ',', coords[1]))
    else:
        leisure = None
    data_list.append([name, highway, leisure, shapely_geo])

    if other_tags and ('historic') in other_tags:
        feat = [x for x in other_tags.split(',') if 'historic' in x][0]
        historic = feat[feat.rfind('>') + 2:feat.rfind('"')]

        if (historic == 'monument'):
            coordsmaps.append([name, coords[0], coords[1]])
            coordsTSP.append([float(coords[0]), float(coords[1])])
            coordsfile.write("%s%s %f%s %f\n" % (name, ',', coords[0], ',', coords[1]))
    else:
        historic = None
    data_list.append([name, highway, historic, shapely_geo])

coordsfile.close()

gdf = gpd.GeoDataFrame(data_list, columns=['Name', 'Highway', 'AmenityTourismLeisureLanduse', 'geometry'],
                       crs={'init': 'epsg:4326'}).to_crs(epsg=3310)
gdf.tail()

places = gdf[gdf.AmenityTourismLeisureLanduse.isin(['place_of_worship', 'museum', 'artwork', 'theme_park', 'monument'])]

print places

move_operator = reversed_sections # move_operator = swapped_locations
max_iterations = 40000

# Input file and store coordinates
# coords = [(37.7768016, -122.4169151),(37.7860105, -122.4025377),(37.7821494, -122.4058960),(37.7689269, -122.4029053),(37.7768800, -122.3911496),(37.7706628, -122.4040139),(37.7870361, -122.4039444),(37.7507903, -122.3877184),(37.7914417, -122.3927229),(37.8672841, -122.5010216)]
coords = coordsTSP
#[(32.824690, 39.972573), (32.824972, 39.972576), (32.820081, 39.969206), (32.819826, 39.969249), (32.818039, 39.972891), (32.824037, 39.970472)]

# Initialize random tour
init_function = lambda: init_random_tour(len(coords))
# tour = init_random_tour(len(coords))
# Construct distance matrix
matrix = cartesian_matrix(coords)
# matrix=distance_matrix(coords)

# Define objective_function
objective_function = lambda tour: -tour_length(matrix, tour)

# Execute hill_climb
iterations, score, best = hillclimb_and_restart(init_function, move_operator, objective_function, max_iterations)
shortest = []
print 'objective_function'
print objective_function
for loc in best:
    shortest.append(loc + 1)

# Write txt and graph for TSP input
fig, ax = plt.subplots(figsize=(10, 10))
createfile = open('test.txt', 'w')
for i, row in places.iterrows():
    x = row['geometry'].x
    y = row['geometry'].y
    name = row['Name']
    print("%s%s %f%s %f" % (name, ',', x, ',', y))
    createfile.write("%s%s %f%s %f\n" % (name, ',', x, ',', y))
    plt.annotate(row['Name'], xy=(x, y), size=13, xytext=(0, 5), textcoords='offset points')
    plt.plot(x, y, 'o', color='#f16824')
    ax.set(aspect=1)
createfile.close()
plt.title('Input', size=16)
print plt.show()

# Nodes and instanceSize are passed into main() using another program
# Just gave them default values for this example
# The Node lookup table.

openfile = open('test.txt', 'r')
Nodes = {}
Allplaces = []
Allx = []
Ally = []
number = 1
for i in openfile:
    text = i.split(', ')  # split the input text based on space & store in the list 'text'
    Nodes[str(number)] = float(text[1]), float(text[1])  # assign the 1st item to key and 2nd item to value of the dictionary
    LineNumber = str(number)
    lineName = text[0]
    lineGx = float(text[1])
    lineGy = float(text[2])
    Allplaces.append([LineNumber, lineName, lineGx, lineGy])
    Allx.append(lineGx)
    Ally.append(lineGy)
    number += 1

print Allplaces
openfile.close()


ShortLine = []

i = 0
for harf in shortest:
    ShortLine.append(harf)
    i += 1
x = []
y = []
j = 0

print ShortLine
coordsmapTSP = []

while j < len(ShortLine):
    d = ShortLine[j]-1
    x.append([Allx[d]])
    y.append([Ally[d]])
    coordsmapTSP.append(coordsmaps[d])
    print coordsmapTSP
    j += 1

fig = plt.figure()

ax1 = fig.add_subplot(111)

ax1.set_title("Output")
ax1.plot(x, y, c='r', label='the data')

ax1.legend()

plt.show()


coords = []
coordsx = []
coordsy = []
verts = []


i = 0
print coordsmapTSP
for line in coordsmapTSP:
    aux = line
    plt.plot(aux[1], aux[2], 'go-', markersize=10)
    coords.append(aux)
    coordsx.append(aux[1])
    coordsy.append(aux[2])
    verts.append((float(coordsx[i]), float(coordsy[i])))
    path = Path(verts)
    print("Cluster: " + str(i) + "\n" + str(verts))
    i += 1
    patches.PathPatch(path)
    plt.plot(coordsx, coordsy, '-', lw=2, color='blue', ms=10)
mplleaflet.show()

tend = datetime.now()
print tend
