#-*- coding: utf-8 -*-
# R. Gautier 2014
# Basic functions for packing defects determination



# Euclidian distance without square root
def dist(coord1, coord2):
    return ((coord1[0]-coord2[0])**2 +
           (coord1[1]-coord2[1])**2 +
           (coord1[2]-coord2[2])**2)

# distance for one axis without square root
def dist_oneAxis(axe1,axe2):
    return (axe1-axe2)**2

