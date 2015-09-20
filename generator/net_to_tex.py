#! /usr/bin/python
import sys
import time
import random
import os
import math


f = open(sys.argv[1])

lines = f.readlines()

nodes = []
links = []
statics = []
bases = []
relays = []


def distance((x1,y1),(x2,y2)):
    return ((x1 - x2)**2 + (y1 - y2)**2)**0.5


def printSet(Id, l):
    if len(l) == 0:
        print "\def\%s{}"%(Id)
        return
    print "\def\%s{"%(Id)
    for i in range(0,len(l)-1):
        print "%.2f/%.2f/%s,"%(l[i][1][0],l[i][1][1],l[i][0])
    i = len(l)-1
    if i >= 0:
        print "%.2f/%.2f/%s}"%(l[i][1][0],l[i][1][1],l[i][0])
    print "\n"

maxx = (0,0)
minx = (1000,1000)
for line in lines:
    s = line.split()
    if s[0] == "node":
        nid = s[1]
        nx = float(s[2])
        ny = float(s[3])
        nodes.append((nid,(nx,ny)))
        maxx = max((nx,ny),maxx)
        minx = min((nx,ny),minx)
        if nid[0] == "s":
            statics.append((nid,(nx,ny)))
        if nid[0] == "b":
            bases.append((nid,(nx,ny)))
        if nid[0] == "r":
            relays.append((nid,(nx,ny)))

    if s[0] == "route":
        i = int(s[1])
        j = int(s[2])
        links.append((nodes[i][0],nodes[j][0]))



printSet("S",statics)
printSet("R",relays)
printSet("B",bases)

if len(sys.argv) > 2:
    if sys.argv[2] == "-r": #Find the links
        tx_range = float(sys.argv[3])
    for i in range(len(nodes)):
        n1 = nodes[i]
        for j in range(i+1,len(nodes)):
            n2 = nodes[j]
            if distance(n1[1],n2[1]) < tx_range:
                links.append((n1[0],n2[0]))


if len(links)==0:
    print "\def\links{}"
else:
    print "\def\links{"
    for l in links[:-1]:
        (i,j) = l
        print "%s/%s,"%(i,j)
    (i,j) = links[-1]
    print "%s/%s}"%(i,j)

area = max(maxx[0],maxx[1])
print "\def\A{%d}"%(math.ceil(area+5))

print "% Number of links = ", len(links)
print "% Number of nodes = ", str(len(statics) + len(bases) + len(relays))
print "% Number of relays = ", len(relays)
