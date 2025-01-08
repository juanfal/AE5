#!/usr/bin/env python3
# encoding: utf-8
# ae5.py
# Carlos Villagrasa, Javier Falgueras
# juanfc 2019-02-16
# ae5 2024-01-22

__version__ = 0.51 # 2024-02-08
# Predation management added 2023-11-14



# TODO:
# - Quitar el DistVal (ver antes)
#  Añadir a json:
#       --numGen int
#       --varia
#       --setRandomSeed int
#       --verbose
#       --saveExcel
#       --outDir 'str'
#       --outFName 'str'
#       --NumberOfCells int
#       --NumberOfRsrcsInEachCell int
#       --Distribution 'str'
#       --DistType 'str'
#       --DistVal int
#       --species 'str'
#       -p
#       --redirectStdout
#       --swallow
#       --saveWorld
#       --noZero 'str'
#       initFile







# ############################################################################# #
#                                 IMPORTS                                       #
# ############################################################################# #

# os
# sys
# json
# textwrap
# pprint
# random
# pathlib
# argparse
# datetime
# pip3 install xlsxwriter
# pip3 install matplotlib
# functools
# subprocess
# urllib

# Added to help installation but repeated later

try:
    import matplotlib
except ImportError:
    print("You need to install the library matplotlib. Copy next line, paste and execute it:")
    print("pip3 install matplotlib")
try:
    import xlsxwriter
except ImportError:
    print("You need to install the library xlsxwriter. Copy next line, paste and execute it:")
    print("pip3 install xlsxwriter")

# end



import os
import sys
import json
import textwrap
from pprint import pprint
import random

from pathlib import Path
import argparse
from datetime import datetime, timedelta

import numpy as np
import numpy.random as npr
import matplotlib.pyplot as plt

# format for printing numpy arrays
np.set_printoptions(formatter={'int': '{: 7d}'.format})

import xlsxwriter
# sys.path.append('/usr/local/lib/python3.7/site-packages')
# import xlrd

NUMPY_EXPERIMENTAL_ARRAY_FUNCTION=0


# # for dynplot
# import matplotlib.pyplot as plt
# from functools import partialmethod


# ############################################################################# #
#                                 CONSTANTS                                     #
# ############################################################################# #
INDIVIDUAL      = 0 # they consume df
ACTOR           = 1 # they consume df+if of recipient
RECIPIENT       = 2 # they consume df (it doesn't consume directly)
RECIPR_INTRA    = 3 # they consume (intra 2x df+if) (inter df1+if1+df2+if2)
RECIPR_INTER    = 4 # they consume (intra 2x df+if) (inter df1+if1+df2+if2)
NUMBER_OF_FORMS = 5

AS_INDIVIDUAL = -1
AS_RECIPIENT  = -2

# Types of distribution
NEIGHBOURS_DISTRIBUTION = 0 # randomly distribute among neighbours around
RANDOM_GLOBAL_AVG       = 1 # random change around global average
RANDOM_GLOBAL_BY_CELLS  = 2 # RANDOM_GLOBAL_AVG leaving empty cells


COMMANDS = "comandos.txt"
COMMANDSLEFT = "comandosRestan.txt"
TO_ABORT_EXECUTION_FILENAME = "ACORTAR"
TIMERECORDFILE = "_TIEMPOS.txt"


# ############################################################################# #
#                                 CLASSES                                       #
# ############################################################################# #
from classes import CQueue, dynplot



# ############################################################################# #
#                                MAIN SUBPROGRAMS                               #
# ############################################################################# #

def newWorld():
    """Build our world.
     First index:  cell index
     Second index: species index
    So we have:
        gWorld[iCell, iSpecies]
    RETURNS
      gWorld, gStatsAnt, gStatsPost
    """
    return (np.zeros((gNumberOfCells,  gNumberOfSpecies), dtype=int),
            np.zeros((gNumberOfSpecies, NUMBER_OF_FORMS), dtype=int),
            np.zeros((gNumberOfSpecies, NUMBER_OF_FORMS), dtype=int)
            )

def randomDist(nitems, nCells, distVal):
    """Returns a list of nitems numbers randomly distributed in nCells"""

    delta = distVal/100
    average = nitems / nCells
    # print(nCells)
    r = np.sort(npr.randint(0, nitems+1, nCells-1))
    r = np.concatenate((r,[nitems])) - np.concatenate(([0],r))
    # by Javi :)
    r = np.around(r * delta + (1-delta) * average).astype(int)
    # compensate rounding simply adding the difference
    # to the first cells
    dif = nitems - r.sum()
    # TODO Javi: hacer una lista de i y permutarla.
    rcell = npr.permutation(range(nCells))
    i = 0
    while dif != 0:
        if r[rcell[i]] + dif >= 0:
            r[rcell[i]] += dif
            dif = 0
        else:
            dif += r[rcell[i]]
            r[rcell[i]] = 0
        i += 1
    npr.shuffle(r)
    return r


def randomDistByCells(nItems, nCells, distVal):
    destCells = int(nCells * (1-distVal/100)+0.5)
    if destCells <= 1:  # when distVal is close to 100 -> all items into 1 cell
        r = np.zeros(nCells, dtype=int)
        rndDesti = npr.randint(0, nCells)
        r[rndDesti] = nItems # all to a only random cell
    else:
        # print(f"nItems, nCells, distVal {nItems}, {nCells}, {distVal}.  DA: {destCells} {nCells} {distVal}")
        r = np.concatenate((randomDist(nItems, destCells, distVal),
                            np.zeros(nCells-destCells, dtype=int)))
        npr.shuffle(r)

    return r



def doInitialDistributionOf1Species(iSpecies):
    n = gConf["species"][iSpecies]["NumberOfItems"]
    meanPerCell = n // gNumberOfCells
    remainder   = n %  gNumberOfCells
    gWorld[ :, iSpecies] = [meanPerCell for _ in range(gNumberOfCells)]
    # the excess of items are equally given out to the first cells
    # it will give this remainder number (0..gNumberOfCells)
    if remainder:
        for _ in range(remainder):
            gWorld[npr.randint(0, gNumberOfCells), iSpecies] += 1


def doInitialDistribution():
    """NOT anymore: In fact this is unnecessary as first step in the generation process
    will do it (again) distribute all equally
    doDistribute() left for the end of generations loop"""

    # If a world file is there, take it!, but only when the number of cells fits

    if gNumberOfCells == checkAndCountPrevWorldNumberOfCells(gWorldCompFileName):
        # adding the possibility of changing the number of cells
        # First: compute the previous number of them
        #
        gWorld[:,:] = [list(map(int, line.split())) for line in open(gWorldCompFileName, encoding="utf-8")]
        # 2020-05-26 fixed the situation in which the read world comes
        # a number of items of one species but the user wants a global
        # new one, number which comes from either the init json or command args or both
        for iSpecies in range(gNumberOfSpecies):
            n = gConf["species"][iSpecies]["NumberOfItems"]
            if gWorld[ :, iSpecies].sum() != n:
                doInitialDistributionOf1Species(iSpecies)
    else:
        # Simplest, averaged distribution
        # Should we call doDistribute() at the end of this?
        for iSpecies in range(gNumberOfSpecies):
            doInitialDistributionOf1Species(iSpecies)
        # TODO: Shall we keep this?
        doDistribute()


def doDistributeOf1Species(iSpecies):
    distType, distVal = getDist(iSpecies)
    if distType == NEIGHBOURS_DISTRIBUTION:     # NEIGHBOURS_DISTRIBUTION
        tempDist = np.zeros((gNumberOfCells), dtype=int)
        distLen = int(gNumberOfCells*distVal/100/2)
        for iCell in range(gNumberOfCells):
            nitems = gWorld[iCell, iSpecies]
            # TODO: if the number of items nitems
            # is large, we could consider directly writing the
            # average on each cell
            for _ in range(nitems):
                dist = random.randint(-distLen, distLen)
                # dist = npr.randint(-distLen, distLen+1)
                p = iCell+dist
                # if p < 0 python automatically wraps
                if p >= gNumberOfCells:
                    p %= gNumberOfCells
                tempDist[p] += 1
        gWorld[ :, iSpecies] = tempDist


    elif distType == RANDOM_GLOBAL_AVG:         # RANDOM_GLOBAL_AVG: x_a \delta + xm (1-\delta)
        nitems = gWorld[ :, iSpecies].sum()
        wildDist = randomDist(nitems, gNumberOfCells, distVal)
        gWorld[ :, iSpecies] = wildDist
    else:                                       # RANDOM_GLOBAL_BY_CELLS
        nitems = gWorld[ :, iSpecies].sum()
        wildDist = randomDistByCells(nitems, gNumberOfCells, distVal)
        gWorld[ :, iSpecies] = wildDist
    # print("----------> ", nitems, delta, average)
    # print("wildDist despu: ", wildDist.sum(), [wildDist[i] for i in range(gNumberOfCells) if wildDist[i] > 0])
    # print("-----> ", nitems - wildDist.sum())


def doDistribute():
    for iSpecies in range(gNumberOfSpecies):
        doDistributeOf1Species(iSpecies)

def doGrouping():
    """Form the groups following the group partner"""
    for iCell in range(gNumberOfCells):
        n = 0
        # printv("Cell grouping %d" % n)
        n += 1
        # stillPossibleGroupInCell = False
        permutedListOrigGroups = npr.permutation(gWithPartnerList)
        # print("Grouping permutation list:", permutedListOrigGroups)
        for iSpecies in permutedListOrigGroups:
            phenFlex = gConf["species"][iSpecies]["PhenotypicFlexibility"]
            i_partnerList = iGetPartnerList(iSpecies)      # list of i partners for the groups
            i_groupedList = iGetGroupStartingInList(iSpecies)  # i already formed group that iSpecies starts
            # number of iSpecies items in that iCell
            ni = gWorld[iCell, iSpecies]

            # item to group with
            ni_partnerList = np.array([ gWorld[iCell, i_partner] for i_partner in i_partnerList ])

            # Phenotypic complexity tells the % (0-1) from the total
            # existing particular species, that should be grouped, so we
            # compare the current amount of items (ungrouped)
            # and check if there is room for more grouping
            ni_toGroup = int(ni * phenFlex)
            ni_total_partners = ni_partnerList.sum()
            if ni_toGroup < ni_total_partners:
                ni_toGroupList = np.trunc(ni_partnerList/ni_total_partners * ni_toGroup).astype(int)
            else:
                ni_toGroupList = ni_partnerList

            ni_feasible = ni_toGroupList.sum()

            # print("ni_toGroup: %d with %d, ni_feasible: %d" % (ni_toGroupList, i_groupedList, ni_feasible))
            if ni_feasible > 0:
                for i, ni_partnerListi in enumerate(i_partnerList):
                    n = ni_toGroupList[i]
                    if ni_partnerListi  == iSpecies:
                        n //= 2
                    gWorld[ iCell,iSpecies] -= n
                    gWorld[iCell, ni_partnerListi] -= n

                    gWorld[iCell, i_groupedList[i]] += n

#
def doAssociation(iCell, queue, predatoryOffspring):
    """ performs the association and immediately moves the new associated couple
    to the queue ready to be shuffled and then to go eating
    Participants in the association leave the cell both, after eating in queue
    they will be born as individuals """
    #
    listOfAssocOrigs = npr.permutation(gListOfAssociationOrigs)
    for iOrig in listOfAssocOrigs:
        ni_Orig = gWorld[iCell, iOrig]

        iAssocList = npr.permutation(iAssociatedTargetGetList(iOrig))
        # iAssocList = iAssociatedTargetGetList(iOrig)

        # build list of populations
        ni_assocList = [ gWorld[iCell, i] for i in iAssocList ]
        ni_assocListSum = sum(ni_assocList)

        if ni_Orig >= ni_assocListSum:
            ni_toAssocList = ni_assocList
        else:
            ni_toAssocList = [ (n * ni_Orig)//ni_assocListSum for n in ni_assocList ]
            # print("Hola", ni_assocList, end='      ') Avisado
            # print(ni_toAssocList)

        # Kind of association: IRAR: Individual, Recipient, Actor, Reciprocal
        for iAssoc, ni_Assoc in zip(iAssocList, ni_toAssocList):
            # assert queue.free() >= ni_Assoc, f"NO QUEPO (1) origs !!! {ni_Assoc} -> {queue.free()}"
            # TODO Javi: dirFitOri  = gConf["species"][iOrig]["DirectOffspring"]
            #            indirFitOri = gConf["species"][iOrig]["IndirectOffspring"]
            #            dirFitAss   = gConf["species"][iAssoc]["DirectOffspring"]
            #            indirFitAss = gConf["species"][iAssoc]["IndirectOffspring"]
            if iOrig == iAssoc:
                # Reciprocal intra A><A
                ni_Assoc -= ni_Assoc % 2
                queue.push(ni_Assoc * [[iOrig, iAssoc]])
                gWorld[iCell, iOrig] -= ni_Assoc
                gStatsAnt[iOrig, RECIPR_INTRA] += ni_Assoc
                #TODO Javi: gPlotGREED[iOrig, RECIPR_INTRA] == DirFitOri

#CAMBIAR
            elif iOrig in iAssociatedTargetGetList(iAssoc):
                # Reciprocal inter-specific A><B
                # assert queue.free() >= ni_Assoc, f"NO QUEPO (2) intra/inter-specific!!! {ni_Assoc} -> {queue.free()}"

                if gConf["species"][iOrig]["IndirectOffspring"] >= 0 and \
                   gConf["species"][iAssoc]["IndirectOffspring"] >= 0:
                    queue.push(ni_Assoc * [[iOrig, iAssoc]])
                    gWorld[iCell, iOrig] -= ni_Assoc
                    queue.push(ni_Assoc * [[iAssoc, iOrig]])
                    gWorld[iCell, iAssoc] -= ni_Assoc
                else:
                    queue.push(ni_Assoc * [[iOrig, AS_INDIVIDUAL]])
                    gWorld[iCell, iOrig] -= ni_Assoc
                    queue.push(ni_Assoc * [[iAssoc, AS_INDIVIDUAL]])
                    gWorld[iCell, iAssoc] -= ni_Assoc
                    predatoryOffspring[iOrig] += ni_Assoc*gConf["species"][iAssoc]["IndirectOffspring"]
                    predatoryOffspring[iAssoc] += ni_Assoc*gConf["species"][iOrig]["IndirectOffspring"]


                gStatsAnt[iOrig, RECIPR_INTER] += ni_Assoc
                gStatsAnt[iAssoc, RECIPR_INTER] += ni_Assoc
                # TODO Javi: gPlotGREED[iOrig, RIN] += ni_Assoc * (DirFitOri+indFitAssoc-abs(indFitOri))
                #            gPlotGREED[iAssoc, RIN] += ni_Assoc * (DirFitAssoc+indFitOrig-abs(indFitAssoc))
            else:
                # either A>B or A>B>C
                ##assert queue.free() >= ni_Assoc, "NO QUEPO (3) as recipient!!! "+str(ni_Assoc)+" -> "+str(queue.free()




                if gConf["species"][iOrig]["IndirectOffspring"] >= 0:
                    queue.push(ni_Assoc * [[iOrig, iAssoc]])
                    gWorld[iCell, iOrig] -= ni_Assoc
                    queue.push(ni_Assoc * [[iAssoc, AS_RECIPIENT]])
                    gWorld[iCell, iAssoc] -= ni_Assoc
                else:
                    queue.push(ni_Assoc * [[iOrig, AS_INDIVIDUAL]])
                    gWorld[iCell, iOrig] -= ni_Assoc
                    queue.push(ni_Assoc * [[iAssoc, AS_RECIPIENT]])
                    gWorld[iCell, iAssoc] -= ni_Assoc
                    predatoryOffspring[iAssoc] += ni_Assoc*gConf["species"][iOrig]["IndirectOffspring"]

#
                gStatsAnt[iOrig, ACTOR] += ni_Assoc
                gStatsAnt[iAssoc, RECIPIENT] += ni_Assoc
                # TODO Javi: gPlotGREED[iOrig, ACTOR] == DirFitOri-abs(indFitOri)
                #            gPlotGREED[iAssoc, RCT] += ni_Assoc * (DirFitAssoc+indFitOrig-abs(indFitAssoc))

            # TODO Javi: gPlotGREED[iSpecies, RECIPR_INTER] == gPlotGREED[iSpecies, RIN]// gStatsAnt[iSpecies, RECIPR_INTER]
            #            gPlotGREED[iSpecies, RECIPIENT]== gPlotGREED[iSpecies, RIN]// gStatsAnt[iSpecies, RECIPIENT]
            #           gPlotSTD = np.around(np.std(gStatAnt[:,:], axis=0), decimals=1)




def doAssociationQueueAndConsume():
    """ performs the association and immediately moves the new associated
    to the queue ready to shuffle and then eat
    Participants in the association leave the cell both. Then individuals """
    #

    if gArgs["varia"]:
        newDirFit  = np.zeros((gNumberOfSpecies), dtype=[("sum", "int"), ("N", "int")])

    predatoryOffspring = np.zeros(gNumberOfSpecies, dtype=int)
    for iCell in range(gNumberOfCells):
        # ENQUEUEING
        # queue the index of item and its receptor
        # the second number next to the index is:
        #     -2 for recipient AS_RECIPIENT
        #     -1 for individual AS_INDIVIDUAL
        # if second number >= 0
        #   If second number == first => receptor == orig => reciprocal intraspecific
        #   If second number != first => receptor != orig => reciprocal interspecific
        nTotalInCell = gWorld[iCell,:].sum()
        if nTotalInCell == 0:
            continue
        # print(f"Cell[{iCell}]: {nTotalInCell}")
        #@ Habría que reservar sólo los de direct > 0 o indirect > 0,
        queue = CQueue(nTotalInCell)   # size of the queue

        # ADD Pairs of associated to the queue
        predatoryOffspring[:] = 0
        doAssociation(iCell, queue, predatoryOffspring)

        # ADD SPARE INDIVIDUALS TO THE QUEUE
        for iSpecies in range(gNumberOfSpecies):
            n = gWorld[iCell, iSpecies]
            if n > 0:
                ##assert queue.free() >= n, "NO QUEPO (4) resto !!! "+n+" -> "+str(queue.free())
                queue.push( n * [[iSpecies, AS_INDIVIDUAL]] )
                gWorld[iCell, iSpecies] = 0
                gStatsAnt[iSpecies, INDIVIDUAL] += n
                # TODO Javi: gPlotGREED[iSpecies, INDIVIDUAL] == gConf[iSpecies]["DirectOffspring"]
            gWorld[iCell, iSpecies] = predatoryOffspring[iSpecies]




        # CONSUME!
        # printv("queue size before eating: ", len(queue))
        ##assert nTotalInCell == len(queue), "Houston we've got a problem"
        queue.shuffle()
        rsrc = gConf["NumberOfRsrcsInEachCell"]
        q = iter(queue)
        leftover = 0
        while q.toIter() > 0 and ((not gArgs["swallow"] and rsrc >= gMinEating) or leftover == 0):
            couple = next(q)
            iOrig = couple[0]
            iDest = couple[1]

            dirFit   = gConf["species"][iOrig]["DirectOffspring"]
            indirFit = gConf["species"][iOrig]["IndirectOffspring"]

            rsrcToEat = 0
            if iDest < 0:  # iOrig is INDIVIDUAL or RECIPIENT

                ##TODO Javi?
                rsrcToEat = dirFit
                leftover = rsrcToEat
                #rsrcToEat = dirFit + abs(indirFit) # taking advantage of the indirFit for himself

                if rsrc >= rsrcToEat:
                    rsrc -= rsrcToEat
                    leftover = 0
                    gWorld[iCell, iOrig] += rsrcToEat
                    if iDest == AS_INDIVIDUAL:
                        gStatsPost[iOrig, INDIVIDUAL] += 1
                    else:
                        gStatsPost[iOrig, RECIPIENT] += 1
            else:
                # giving indirect offspring
                # indirFit = gConf["species"][iOrig]["IndirectOffspring"]
                if gArgs["varia"]:
                    ##TODO only for those with list of assoc and after
                    #      seeing the form
                    fitVarLimit   = gConf["species"][iOrig]["FitnessVariationLimit"]
                    dirFit, indirFit = fitnessVariations(dirFit, indirFit, fitVarLimit)
                    newDirFit[iOrig]["sum"] += dirFit * dirFit
                    newDirFit[iOrig]["N"] += dirFit

                # The others (ACTOR, and RECIPR_INTRA/RECIPR_INTER) have to give IndirectOffspring
                if indirFit > 0:
                    rsrcToEat = dirFit + indirFit
                else:
                    rsrcToEat = dirFit
                leftover = rsrcToEat
                if rsrc >= rsrcToEat:
                    rsrc -= rsrcToEat
                    leftover = 0
                    gWorld[iCell, iOrig] += dirFit
                    #gWorld[iCell, iDest] = noNeg( gWorld[iCell, iDest] + indirFit )
                    gWorld[iCell, iDest] += indirFit #JAVI
                    if iOrig == iDest:
                        gStatsPost[iOrig, RECIPR_INTRA] += 1
                    elif iOrig in iAssociatedTargetGetList(iDest): # is iOrig among the assoc of assoc
                        gStatsPost[iOrig, RECIPR_INTER] += 1
                    else:
                        gStatsPost[iOrig, ACTOR] += 1





    if gArgs["varia"]:
        for iSpecies in range(gNumberOfSpecies):
            # CHANGE THE GLOBAL DirectOffspring gConf parameter
            if newDirFit[iSpecies]["N"]:
                incFit = gConf["species"][iSpecies]["DirectOffspring"] + gConf["species"][iSpecies]["IndirectOffspring"]
                dirFit = int(round(newDirFit[iSpecies]["sum"] / newDirFit[iSpecies]["N"]))
                gConf["species"][iSpecies]["DirectOffspring"] = dirFit
                gConf["species"][iSpecies]["IndirectOffspring"] = incFit - dirFit

def doUngroup():
    for iOrig in gWithPartnerList:
        iGroupList = iGetGroupStartingInList(iOrig)
        for iGroup in iGroupList:
            phenFlex = gConf["species"][iGroup]["PhenotypicFlexibility"]
            for iCell in range(gNumberOfCells):

                # number of iGroup items in that iCell
                ni = gWorld[iCell,iGroup] # form index is 0 always

                ni_unGroup = int(round(ni * (1.0 - phenFlex)))

                if ni_unGroup > 0:
                    iPartner = iGetPartnerFromOrigAndGroup(iOrig, iGroup)
                    gWorld[iCell,   iOrig] += ni_unGroup
                    gWorld[iCell,iPartner] += ni_unGroup
                    gWorld[iCell,  iGroup] -= ni_unGroup


# ############################################################################# #
#                                     TOOLS                                     #
# ############################################################################# #

def isInt(anyNumberOrString):
    try:
        int(anyNumberOrString)
        return True
    except ValueError:
        return False

def iFrom_id(theId):
    """Returns index of species from its id, or -1"""
    i = gNumberOfSpecies - 1
    while i >= 0 and gConf["species"][i]["id"] != theId:
        i -= 1
    return i

def iListFrom_idList(idList):
    """Returns index of species from its id, or -1"""
    return [iFrom_id(toFindId) for toFindId in idList]

def noNeg(n):
    return max(n, 0)

def checkNegs():
    for iCell in range(gNumberOfCells):
        for iSpecies in range(gNumberOfSpecies):
            if gWorld[iCell, iSpecies] < 0:
                print(f"************************  gWorld[{iCell}, {iSpecies}] NEGATIVE {gWorld[iCell, iSpecies]}")
                # print("************************  gWorld[%d, %d] NEGATIVE %d" % (iCell, iSpecies, gWorld[iCell, iSpecies]))


def add_cont(s):
    if s.endswith('.json'):
        s = s[:s.rfind('.')]
    if s.startswith('-') or s.endswith('.py') or s.endswith('_cont'):
        return s

    return s + '_cont'


# ------------ Groups

def iGetPartnerList(iSpecies):
    "returns a list of indexes of the partners for grouping"
    toFindIdList = gConf["species"][iSpecies]["GroupPartners"]
    return [iFrom_id(toFindId) for toFindId in toFindIdList]

def iGetGroupStartingInList(iSpecies):
    idOrig = gConf["species"][iSpecies]["id"]
    toFindIList = [idOrig + '|' + idPartner for idPartner in gConf["species"][iSpecies]["GroupPartners"]]
    return [iFrom_id(toFindId) for toFindId in toFindIList]

def iGetPartnerFromOrigAndGroup(iOrig, iGroup):
    idOrig = gConf["species"][iOrig]["id"]
    idGroup = gConf["species"][iGroup]["id"]
    return iFrom_id(idGroup[len(idOrig)+1:])

def getListOfOrigGroups():
    theList = [i for i in range(gNumberOfSpecies) if len(gConf["species"][i]["GroupPartners"]) ]

    # Puts complex groups first
    theList.sort(key = lambda i:
            gConf["species"][i]["id"].count('|') + \
                str(gConf["species"][i]["GroupPartners"]).count('|'), reverse=True)
    return theList

# ------------ Association

def collectAssociationActors():
    """Return a list of indexes of origs that are actually associated with
    some other Now when no associates, is [] """

    return [i for i in range(gNumberOfSpecies) if len(gConf["species"][i]["AssociatedSpecies"]) ]

def iAssociatedTargetGetList(iSpecies):
    '''list of ids of its associates '''
    return iListFrom_idList(gConf["species"][iSpecies]["AssociatedSpecies"])



# ------------ Others

def getDist(iSpecies):
    """from 45n or 50n returns the integer and the type
    n: NEIGHBOURS_DISTRIBUTION
    r: RANDOM_GLOBAL_AVG
    """
    if "Distribution" in gConf["species"][iSpecies] or "Distribution" in gConf:
        if "Distribution" in gConf["species"][iSpecies]:
            dist = gConf["species"][iSpecies]["Distribution"]
        else:
            dist = gConf["Distribution"]

        if dist.endswith("n"): # average neighbours distribution
            distType = NEIGHBOURS_DISTRIBUTION
        elif dist.endswith("r"):
            distType = RANDOM_GLOBAL_AVG
        else:
            distType = RANDOM_GLOBAL_BY_CELLS

        distVal = int(dist[:-1])
        if distVal > 100 or distVal < 0:
            print("Error in Distribution val (should be 0<=x<=100. Is:", distVal)
            sys.exit(1)



    return distType, distVal

def fitnessVariations(direct, indirect, fitVarLimit):
    """Returns a pair of new rand int direc-indirect fitnesses inside the
    fitVarLimit"""

    tot = direct + abs(indirect)

    # the pairs ordered in the smoothest possible way
    directRange = list(range(tot, -1, -1)) + list(range(0, tot))
    indireRange = list(range(0, tot + 1)) + list(range(-tot, 0, 1))
    pairsOfFitness = list(zip(directRange, indireRange))

    currentI = pairsOfFitness.index((direct, indirect))
    lpairsOfFitness = len(pairsOfFitness)

    # 4 cases: fitVarLimit > Maxfitness, upper-overflow, lower-underflow, inside
    if 2*fitVarLimit >= lpairsOfFitness:
        bottom = 0
        topp = lpairsOfFitness
    elif currentI + fitVarLimit >= lpairsOfFitness:
        extraOver = currentI + fitVarLimit - (lpairsOfFitness - 1)
        bottom = currentI - fitVarLimit - extraOver
        topp = lpairsOfFitness
    elif currentI - fitVarLimit < 0:
        extraBelow = fitVarLimit - currentI
        bottom = 0
        topp = currentI + fitVarLimit + extraBelow
    else:
        bottom = currentI - fitVarLimit
        topp = currentI + fitVarLimit + 1

    newI = int(npr.randint(bottom, topp))
    return pairsOfFitness[newI]

# def ranking(l, value):
#     rank = unique(l.flatten())
#     return 1+where(rank==value)[0][0]

def checkAndCountPrevWorldNumberOfCells(fName):
    if not os.path.isfile(fName):
        return 0
    with open(fName, encoding="utf-8") as f:
        return sum(1 for _ in f)

def neededFoodFor(iSpecies):
    return gConf["species"][iSpecies]["DirectOffspring"]
    # dirFit   = gConf["species"][iSpecies]["DirectOffspring"]
    # indirFit = gConf["species"][iSpecies]["IndirectOffspring"]
    # return dirFit + abs(indirFit)

def minNeededFood():
    return min(neededFoodFor(i) for i in range(gNumberOfSpecies))



def initExcel():
    global gExcelCellHeader, gExcelCellID, gExcelWorkbook, gExcelWorksheet

    # # Check if there already is a collective Excel
    # if gGlobalExcel:
    #     gExcelWorkbook = xlsxwriter.Workbook(gGlobalExcel)
    #     if os.path.isfile(gGlobalExcel):
    #         prexcel = xlrd.open_workbook(gGlobalExcel, encoding="utf-8")
    #         sheets = prexcel.sheets()
    #         # run through the sheets and store sheets in workbook
    #         # this still doesn't write to the file yet
    #         for sheet in sheets: # write data from old file
    #             newSheet = gExcelWorkbook.add_worksheet(sheet.name)
    #             for row in range(sheet.nrows):
    #                 for col in range(sheet.ncols):
    #                     newSheet.write(row, col, sheet.cell(row, col).value)

    #     gExcelWorksheet = gExcelWorkbook.add_worksheet(gOutFNameBase)
    # else:
    #     gExcelOut = os.path.join(gOutDir, gOutFNameBase + ".xlsx")
    #     gExcelWorkbook = xlsxwriter.Workbook(gExcelOut)
    #     gExcelWorksheet = gExcelWorkbook.add_worksheet(gOutFNameBase)
    gExcelWorkbook = xlsxwriter.Workbook(gExcelOut)
    gExcelWorksheet = gExcelWorkbook.add_worksheet()

    gExcelCellHeader = gExcelWorkbook.add_format({
                                           'align':'center', 'bg_color': '#CCCCFF', 'bold': True})
    gExcelCellID = gExcelWorkbook.add_format({'align':'center', 'bg_color': '#33FF99', 'bold': True})

    gExcelWorkbook.set_properties({
        'title':    'Evolutionary Automata',
        'subject':  gInitConfFile + '_' + gThedatetimeStamp,
        'author':   'Javier Falgueras, Juan Falgueras, Andrés Moya',
        'manager':  'Javier Falgueras, Juan Falgueras, Andrés Moya',
        'company':  'Univ. Valencia + Univ. Málaga',
        'category': 'Research Thesis',
        'keywords': 'Evolution Automata',
        'comments': " ".join(sys.argv)})

def saveExcel(numGen):
    """Saving in Excel"""
    # id, NumberOfItems, DirectOffspring, IndirectOffspring,
    # AssociatedSpecies, FitnessVariationLimit, INDIVIDUAL, ACTOR, RECIPIENT,
    # RECIPR_INTRA/RECIPR_INTER
    if 'gExcelCellHeader' not in globals():
        initExcel()
    txtOutName = os.path.join(gOutDir, gOutFNameBase + ".txt")
    txtOut = open(txtOutName, "a", encoding="utf-8")

    globalsHeader =  ["NCel", "RsCel", "Dst"]
    iSpecHeader =  ["ID", "dst", "D", "I", ">", "Fv", "Gr", "Ph", "Drσ", "IND", "2", "ACT", "2", "RNT", "2", "RIN", "2", "REX", "2"]
    globalsHeaderLen = len(globalsHeader)
    iSpecHeaderLen = len(iSpecHeader)

    if numGen == 1:
        gExcelWorksheet.write_row(0,0, globalsHeader, gExcelCellHeader)
    gExcelWorksheet.write_row(numGen,0, [gConf["NumberOfCells"], gConf["NumberOfRsrcsInEachCell"], gConf["Distribution"]])
    print(gConf["NumberOfCells"], "\t", gConf["NumberOfRsrcsInEachCell"],"\t", gConf["Distribution"], "\t", file=txtOut, end="", sep="")


    standardDistribution = np.around(np.std(gWorld[:,:], axis=0), decimals=1)

    if numGen == 1:
        for iSpecies in range(gNumberOfSpecies):
            gExcelWorksheet.write_row(0,globalsHeaderLen+iSpecies*iSpecHeaderLen, iSpecHeader, gExcelCellHeader)


    for iSpecies in range(gNumberOfSpecies):
        gExcelWorksheet.write(numGen, globalsHeaderLen+iSpecies*iSpecHeaderLen,
                 gConf["species"][iSpecies]["id"], gExcelCellID)
        print(gConf["species"][iSpecies]["id"] + "\t", file=txtOut, end="", sep="")

        if "Distribution" in gConf["species"][iSpecies]:
            speciesDist = gConf["species"][iSpecies]["Distribution"]
        else:
            speciesDist = ""

        toWrite =[
             speciesDist,
             gConf["species"][iSpecies]["DirectOffspring"],
             gConf["species"][iSpecies]["IndirectOffspring"],
             ",".join(gConf["species"][iSpecies]["AssociatedSpecies"]),
             gConf["species"][iSpecies]["FitnessVariationLimit"],
             ",".join(gConf["species"][iSpecies]["GroupPartners"]),
             gConf["species"][iSpecies]["PhenotypicFlexibility"],
             standardDistribution[iSpecies],
             gStatsAnt[ iSpecies,INDIVIDUAL],
             gStatsPost[iSpecies,INDIVIDUAL],
             gStatsAnt[ iSpecies,ACTOR],
             gStatsPost[iSpecies,ACTOR],
             gStatsAnt[ iSpecies,RECIPIENT],
             gStatsPost[iSpecies,RECIPIENT],
             gStatsAnt[ iSpecies,RECIPR_INTRA],
             gStatsPost[iSpecies,RECIPR_INTRA],
             gStatsAnt[ iSpecies,RECIPR_INTER],
             gStatsPost[iSpecies,RECIPR_INTER]
        ]
        gExcelWorksheet.write_row(numGen, globalsHeaderLen+iSpecies*iSpecHeaderLen+1, toWrite)
        print("\t".join(map(str,toWrite)), file=txtOut, end="", sep="")
        if iSpecies < gNumberOfSpecies-1:
            print("\t", file=txtOut, end="")

        if numGen == 1:
            ori = iSpecies*iSpecHeaderLen + 1 + globalsHeaderLen
            end = ori + iSpecHeaderLen - 2
            gExcelWorksheet.set_column(ori, end, None, None, {'level': 1, 'hidden': True})

    print(file=txtOut)

def saveWorld(genNumber):
    ext = ".world.txt"
    fname = gNewWorldCompFileName[:-10]
    gendirname = fname + "_Generations"
    if not os.path.isdir(gendirname):
        os.makedirs(gendirname)
    with open(f"{gendirname}/{genNumber:04d}{ext}" , 'w', encoding="utf-8") as outfile:
        for iCell in range(gNumberOfCells):
            print("\t".join(map(str, gWorld[iCell, :])), file=outfile)


def saveConf():
    """Save conf in a new _cont.json file
    if not there, or rewrite previous _cont.json file if was the initial conf loaded
    Added: save corresponding final gWorld state"""

    # gThedatetimeStamp = datetime.now().strftime("%Y%m%d-%H%M%S.%f")


    # save the conf, but ala!, the number of items is totally different
    for iSpecies in range(gNumberOfSpecies):
        gConf["species"][iSpecies]["NumberOfItems"] = int(gWorld[:, iSpecies].sum())
    with open(gNewConfCompFileName, 'w', encoding="utf-8") as outfile:
        json.dump(gConf, outfile, sort_keys = True, indent = 4,
                   ensure_ascii = False)
    # save the matrix of cells state for possible re-reading
    with open(gNewWorldCompFileName, 'w', encoding="utf-8") as outfile:
        for iCell in range(gNumberOfCells):
            print("\t".join(map(str, gWorld[iCell, :])), file=outfile)

def checkConf(conf):
    """Verifies conf returning inconsistencies or an empty string"""
    # TODO
    # verify gConf so
    #   - groups agree
    #
    problems = ""
    mustBeDefined = []

    items = {"NumberOfCells", "NumberOfRsrcsInEachCell", "species"}

    itemsspecies = {"id", "Distribution", "NumberOfItems", "DirectOffspring",
        "GroupPartners", "PhenotypicFlexibility", "AssociatedSpecies",
        "IndirectOffspring", "FitnessVariationLimit"}

    if not items.issubset(conf.keys()):
        problems += f"Some essential item(s) is not defined\n\t Check the next items are all there:\n\t{list(items)}\n"

    # Check values for Distribution, if there
    if "Distribution" in conf and conf["Distribution"][-1] not in "rhn":
        problems += "Distribution global value must end either in r or n\n"


    previousIds = set([])
    longListSpecies = len(conf["species"])
    if not isinstance(conf["species"], list) or longListSpecies == 0:
        problems += "Species must be a list with species inside\n"
    for i in range(longListSpecies):
        ##TODO check the IndirectOffspring and variability are 0 when
        #      the list of associated is empty
        confi = conf["species"][i]
        if "Distribution" in confi and confi["Distribution"][-1] not in "rhn":
            problems += f"Distribution for species {i} must end either in r or n\n"


        theId = confi["id"]
        # partnerId = confi["GroupPartners"]
        if theId in previousIds:
            problems += f"Id '{theId}' REPEATED\n"
        else:
            previousIds.add(theId)
        # if partnerId != "":
        #     if iGetPartnerList(i) == -1:
        #         problems += "Species %d, id: '%s' has not the partner species '%s' in the conf\n" % \
        #             (i, theId, partnerId)
        #     if iGetGroupStartingInList(i) == -1:
        #         problems += "Group id: '%s' is not configured\n" % \
        #             (theId + '|' + partnerId)
        for aId in theId.split('|'):
            if iFrom_id(aId) == -1:
                problems += f"Component id: '{aId}' from group '{theId}' is not configured\n"

        mustBeDefined += [ theId+'|'+ i    for i in confi["GroupPartners"]]


        # print("***********", confi.keys())
        diffKeys = set(confi.keys()) - itemsspecies
        if len(diffKeys):
            problems += f"""\
Some essential item(s) of species '{theId}' (number {i}) is not defined

    Check syntax: {diffKeys}"""

    # pprint(mustBeDefined)
    for theID in set(mustBeDefined):
        if iFrom_id(theID) == -1:
            problems += f"The species '{theID}' should be defined as it is a potential group\n"

    if problems:
        print(f"\nERROR(s)\nIn the configuration file '{gInitConfFile}' or in the command line arguments added there are inconsistencies:\n\n{problems}\n\n")
        sys.exit(1)

def completeDataFileNames(fromfilename, inDirectory):
    conf = os.path.join(inDirectory, fromfilename + ".json")
    world = os.path.join(inDirectory, fromfilename + ".world.txt")
    return conf, world


def anyZeros(l):
    i = 0
    # print(f"zero[{l[i]}]? {gWorld[:,l[i]].sum(axis=0)}")
    while i < len(l) and gWorld[:,l[i]].sum(axis=0) > 0:
        i += 1
    return i < len(l)


def printv(*args):
    if gArgs["verbose"]:
        if isinstance(args[0], str):
            print(*args)
        else:
            pprint(args)

# ######################## #
# ARGS AND CONF PROCESSING #
# ######################## #

# (to abort execution create a file "{TO_ABORT_EXECUTION_FILENAME}" in this very directory)
# References:

def defineAndGetCommandLineArgs():
    """Parses and returns command line arguments"""

    theArgParser = argparse.ArgumentParser(description="""
                     * Evolutionary Automata *

            Automata to evolve life through generations

To follow the examples described in the Article this program is made for:

Once the program is installed and running, to try each of the examples you can
use the name of the file in 'data' (no need to use the name of the `data`
directory) like: 'Paper2/2Plank1.json'

Or you can add specific parameters like the ones in the Article this program is
the base for.  In this case the parameters parameters for each case can be
copy–pasted into the Terminal from the table “Terminal Command or Rank:” found
in the Supplementary Material.

For example, for the 'Paradox of the plankton' one of the examples is:

  python3 ae5.py Paper1/1Plank4 --numGen=200 --saveExcel --setRandomSeed=1 --verbose








""",
                    formatter_class=argparse.RawTextHelpFormatter)


    # initial configuration to load from a file
    theArgParser.add_argument(
        "initFile",
        nargs="?",
        default="defaultInit",
        help=textwrap.dedent("""\
        When provided it sets the file with the initial configuration.
        All the configuration files are inside the './data' directory.
        Do not add the .extension to the file name.
        If not provided, it defaults to 'defaultInit'
        HINT:
            You can use
                egoism/test1
            as init file and then, _cont files will be generated inside that data/egoism/
            folder""")
    )

    # Number of generations to run
    theArgParser.add_argument(
        "--numGen", type=int,
        default=10,
        metavar="int",
        help="Sets the number of generations to run, default: 10")

    # We want phenotypic variability
    theArgParser.add_argument(
        "--varia", help=textwrap.dedent("""\
        It takes the FitnessVariationLimit for each species
        and changes in each generation the Direct/Indirect fitnesses
        inside that limits: currentValue+- FitnessVariationLimit
        following a list of pairs.
        False if not provided"""),
        action='store_true')

    # Setting random seed to 0 to repeat pseudorandom values
    theArgParser.add_argument(
        "--setRandomSeed", type=int,
        default=-1,
        metavar="int",
        help=textwrap.dedent("""\
        Sets the random seed to a fixed initial value so
          to repeat same random sequences. If not, each
          running will start with different seeds"""),
        )

    # We want phenotypic variability
    theArgParser.add_argument(
        "--verbose", help=textwrap.dedent("""\
        Gives as much detailed information as it can"""),
        default=False,
        action='store_true')

    # # We want egoism multilevel selection
    # theArgParser.add_argument(
    #     "--egoism", help=textwrap.dedent("""\
    #     Considers egoism of each item in multilevel selection"""),
    #     action='store_true')

    # We want to save final status
    theArgParser.add_argument(
        "--saveExcel", help=textwrap.dedent("""\
        Save stats in 'Excel' file and in a txt file.
        See Excel iSpecHeader for meaning of txt columns
        It is set anyway to True
            if any --outFName or outDir is set"""),
        action='store_true')

    # Directory for the results.
    # If no /start relative to ./results
    theArgParser.add_argument(
        "--outDir", type=str,
        metavar="'str'",
        help=textwrap.dedent("""\
        Specifies another than the default output directory
        where to save the .txt and .xlsx files with global outputs.
        If the path of the output dir starts in /
            is considered an ABSOLUTE path.
        else
            is considered a path relative to ./results/
        It sets --saveExcel to True
        Example:
           --outDir=assocTests
           --outDir=/absoluteDir
           --outDir=more/than/one/level""")
    )

    # Filename for the results.
    # If none provided, the default original
    # conf+date is taken
    theArgParser.add_argument(
        "--outFName", type=str,
        metavar="'str'",
        help=textwrap.dedent("""\
        Specifies another than the default output filename
        (initFilename+date)
        where to save the .txt and .xlsx files with global outputs.
        It sets --saveExcel to True
        It appends a new sheet after to the Excel file. The name
        of the Excel file is the name of the output directory
        Example:
            --outFNname=assocTest20""")
    )

    theArgParser.add_argument(
        "--NumberOfCells", type=int,
        default=argparse.SUPPRESS, metavar="int",
        help="Sets the total number of cells in our world")

    theArgParser.add_argument(
        "--NumberOfRsrcsInEachCell", type=int,
        default=argparse.SUPPRESS, metavar="int",
        help="Sets the items of resources in each cell for each generation (int)")


    theArgParser.add_argument(
        "--Distribution", type=str,
        metavar="'str'",
        default=argparse.SUPPRESS,
        help=textwrap.dedent("""\
        It allows specify the desired type of distribution between:

          a) RANDOM_GLOBAL_AVG          from 0r to 100r (with suffix r)
                           to distribute all the elements among all the cells
                           0r    means every cell has the average, the same value
                           100r: means a total random number in each cell

                           All the items are taken from the original cells,
                           their original positions are forgotten and then
                           distributed randomly

                           The randomization is done splitting the queue of
                           elements to assign in random lengths and those
                           lengths are the number of items in each cell

          b) RANDOM_GLOBAL_BY_CELLS     from 0h to 100h (with suffix h)
                           to distribute all the elements among _some_ cells, it
                           excludes a proportion of cells, which are left
                           without items

                           With higher h we have more empty cells and it
                           quickly fills fewer cells with larger amount
                           of items in each one

                           0h    all the cells are meant to fill and they will
                                    get the same average number of items
                           30h   it leaves 30%% of cells empty
                           100h  this extreme value is used to put all the
                                    items in one cell

                           The way the items are sent to the cells is the
                           RANDOM_GLOBAL_AVG (r) way.  Here we use the very same
                           value of the percentage of cells to leave empty that
                           goes with h as the level of randomization on those
                           destination cells.  This strengthen the sense of
                           randomization as h (=r) grows.  So
                                30h means:
                                    30%% of cells left empty
                                    70%% of cells filled in a 30r way


          c) NEIGHBOURS_DISTRIBUTION     from 0n to 100n (with suffix n)
                           to distribute each item in each cell among the cells
                           around.  The number n used is a percentage of the
                           number of total cells, so [-n/2, +n/2] is the
                           destination segment around the original cell, in a
                           circular way at the limits.  The distribution of the
                           original cell items into this segment of cells is
                           done randomly choosing the destination cell with
                           randint for each item. The list of cells around is
                           considered circular. 0 means no cells around are
                           considered, so cells are isolated and are not
                           redistributed

        In general:
          0r   means no random, it is equivalent to 100n
                distributed in case (c)
          100r means a _totally_ random number into all the cells

          0h   means no random, it is equivalent to 100n
                distributed in case (c)
          30h  means 30%% of cells will be left empty, the other 70%%
                will be filled with a (r, we use here 30r) distribution
          100h means all the items into 1 cell

            """)
    )


    theArgParser.add_argument(
        "--species", type=str,
        default=argparse.SUPPRESS, metavar="'str'",
        help=textwrap.dedent("""\
    Sets one or more species parameters
    Surround everything with " "
    No spaces

    Parameters you can change/add in .json conf:

        "id": "B",
        "NumberOfItems": 100000,
        "DirectOffspring": 0,           # int >= 0
        "GroupPartners": ["B", "C"],
        "PhenotypicFlexibility": 1.0,   # 0-1

        "AssociatedSpecies":   [],
        "IndirectOffspring": 0,         # int
        "FitnessVariationLimit": 0      # int >= 0

    Parameters you can change/add as command line parameter:

        NumberOfItems=100000;
        DirectOffspring=0;
        PhenotypicFlexibility=1.0;
        IndirectOffspring=0;
        FitnessVariationLimit=0;


    Put first the 0-index of the species:
        --species="0;DirectOffspring=8;1;IndirectOffspring=1"
    Changes the DirectOffspring of the first species and the
    IndirectOffspring of the second.

    A series of indexes separated by semicolons ; is valid:
        --species="0;3;1;DirectOffspring=8"
    changes DirectOffspring of species 0, 3 and 1

    A Python standard range (no 'steps' here) is valid:
        --species="0;8:12;DirectOffspring=8"
    as with Python 8:12 means 8;9;10;11 (12 is not included)

    When an index is negative, is suppose from the end:
        --species="0;3;1;8:-2;DirectOffspring=8"
    where 8:-2  is  8,9,10,11,12,13  if there were 15 species.

    An empty index in ranges means a extreme:
        --species=":;DirectOffspring=8"
    means change all species
    """))


    theArgParser.add_argument('-p', action='store_true',
    help=textwrap.dedent("""\
    Adding -p to the parameters makes a matplotlib
    graph to appear and dynamically get updated
    as generations go on
    """))

    # We want to redirect stdout to
    theArgParser.add_argument("--redirectStdout",
    action='store_true', default=False,
    help=textwrap.dedent("""\
        Redirect stdout to outFName_stdout.txt""")
    )


    theArgParser.add_argument('--swallow',
    action='store_true', default=False,
    help=textwrap.dedent("""\
    the item in spite of not having enough food left takes
    that insufficient food and prevents the rest of queue from eating it""")
    )

    theArgParser.add_argument('--saveWorld',
    action='store_true', default=False,
    help=textwrap.dedent("""\
    save the state of the World after each generation in ./data/...""")
    )


    theArgParser.add_argument(
        "--noZero", type=str,
        default="", metavar="'str'",
        help=textwrap.dedent("""\
        Name or give the index of the species that when disappear, the execution has to finish.
        If -, then any zero in any species will stop the execution.

        When it stops by this criteria no enter is required and
        the graphic output will disappear at the end

        For example:
            --noZero="gamA0,gamB0"
            --noZero=1,0,2
            --noZero=-""")
    )

    return vars(theArgParser.parse_args())

def readInitConfFile(fileName):
    """Load config from init file"""
    with open(fileName, encoding="utf-8") as f:
        conf = json.load(f)

    return conf

def replaceAndAddDefaultsArgsInConf(conf, args):
    """Replace the conf args with the args provided in the command line"""
    def parseSpeciesArgs(s):

        lenSpecies = len(conf["species"])
        l = s.split(';')
        indxsRead = True
        for arg in l:
            if isInt(arg) or ':' in arg:   # int (+-) or range
                if indxsRead:
                    ns = []
                    indxsRead = False
                if isInt(arg):
                    ns.append(int(arg))
                else:
                    if arg == ':':
                        first = 0
                        last = lenSpecies
                    elif arg.startswith(':'):
                        first = 0
                        last = int(arg[1:])
                    elif arg.endswith(':'):
                        first = int(arg[1:])
                        last = lenSpecies
                    else:
                        first, last = map(int,arg.split(':'))
                    if first < 0: first += lenSpecies
                    if last  < 0: last  += lenSpecies
                    ns += list(range(first, last))
                    printv(ns)
            else:
                indxsRead = True
                k, val=arg.split('=')
                for n in ns:
                    if k in conf["species"][n]:
                        # if it was an existing key, adapt
                        t = type(conf["species"][n][k])
                        conf["species"][n][k] = t(val)
                    else:
                        # if not, create
                        conf["species"][n][k] = val

    # substitute init file conf specific parameters
    # provided through the args in the command line
    for param in args:
        if param in conf:
            if param == "species":
                parseSpeciesArgs(args["species"])
            else:
                if conf[param] is not None:
                    t = type(conf[param])
                    conf[param] = t(args[param])
        else:
            conf[param] = args[param]
    mandatory = ["NumberOfItems",
                 "DirectOffspring"]
    defaults = {"IndirectOffspring": 0,
                "AssociatedSpecies": [],
                "GroupPartners": []
                }

    for species in conf["species"]:
        problems = set(mandatory) - set(species.keys())
        if not problems:
            for tod in set(defaults) - set(species.keys()):
                species[tod] = defaults[tod]
        else:
            print(f"""
ERROR(s) in the init file '{gInitConfFile}.json'
         or in command line arguments. Missing values:

For species {species["id"]} no value(s) for: {", ".join(problems)}
""")
            sys.exit(2)
    return conf

def initPlot():
    fig, ax = plt.subplots()
    ax.set_autoscaley_on(True)

    seqPobPrev = np.array([gWorld[:,:].sum(axis=0)])
    for i in range(gNumberOfSpecies):
        plt.plot([0],  seqPobPrev[:,i])
    _ = ax.set_title('Population for each species')
    ax.yaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter('{x:,.0f}'))

    ax.grid(visible=True, which='both', color='0.65', linestyle='-')
    plt.ion()
    plt.show()

    return fig, ax, seqPobPrev

def addPlot(fig, ax, seqPobPrev, genNumber):
    seqPobPrev = np.append(seqPobPrev, [gWorld[:,:].sum(axis=0)], axis=0)
    for i in range(gNumberOfSpecies):
        if genNumber > 1:
            plt.plot(range(genNumber+1),  seqPobPrev[:,i], color = f'C{i}')
        else:
            plt.plot(range(genNumber+1),  seqPobPrev[:,i], label=gConf["species"][i]["id"], color = f'C{i}')

    ax.legend(loc='upper left')
    _ = ax.set_xlabel(f'Generation: {genNumber}')
    ax.set_title('Population of each species')

    plt.draw()
    plt.pause(0.0005)


    return fig, ax, seqPobPrev



def timeRecording(tNow, nGens, totGens, td, fromArgs):
    if not os.path.isfile(TIMERECORDFILE):
        with open(TIMERECORDFILE, "w", encoding="utf-8") as f:
            print("Fecha    -Hora:    h:mm:segs      [ngen/tot],  dirFile", file=f)
            print("================== ============== ============ ========", file=f)

    # days, hours, minutes, seconds = td.days, td.seconds // 3600, td.seconds % 3600 / 60.0, td.seconds % 60.0
    # hours += days*24
    with open(TIMERECORDFILE, "a", encoding="utf-8") as f:
        print(f"{{{tNow.strftime('%Y%m%d-%H%M%S')}}}: {td} [{nGens:4}/{totGens}], {fromArgs}", file=f)
        # print(f"{{{tNow}}}: {hours:2}h {minutes:2}m {seconds:2}s [{nGens:4}/{totGens}], {fromArgs}", file=f)

def build_defaultjson(inDir):
    with open(os.path.join(inDir, 'defaultInit.json'), 'w', encoding="utf-8") as f:
        f.write("""{
    "NumberOfCells": 15,
    "NumberOfRsrcsInEachCell": 10000,
    "Distribution": "100r",
    "species": [
        {"id": "A",
            "NumberOfItems": 100,
            "DirectOffspring": 5,

            "GroupPartners":   ["A"],
            "PhenotypicFlexibility": 0.0,

            "AssociatedSpecies":   [],
            "IndirectOffspring": 2,
            "FitnessVariationLimit": 0
        },
        {"id": "A|A",
            "NumberOfItems": 0,
            "DirectOffspring": 0,

            "GroupPartners":   [],
            "PhenotypicFlexibility": 0.0,

            "AssociatedSpecies":   [],
            "IndirectOffspring": 2,
            "FitnessVariationLimit": 0
        }
    ]
}
""")



######################################################################
#                                                                    #
#                                            88                      #
#                                            ""                      #
#                                                                    #
#            88,dPYba,,adPYba,   ,adPPYYba,  88  8b,dPPYba,          #
#            88P'   "88"    "8a  ""     `Y8  88  88P'   `"8a         #
#            88      88      88  ,adPPPPP88  88  88       88         #
#            88      88      88  88,    ,88  88  88       88         #
#            88      88      88  `"8bbdP"Y8  88  88       88         #
#                                                                    #
#                                                                    #
######################################################################

# ####### #
# GLOBALS #
# ####### #

gOnLinux = False


startTime = datetime.now()
# on Linux to check time limits
REDZONETIME = None
if gOnLinux:
    REDZONETIME = startTime + timedelta(days=6, hours=18)

# INPUT
# GET CONF FROM .json INIT FILE AND THEN FROM ARGS
# init conf files (json and world)
gArgs = defineAndGetCommandLineArgs()
gInitConfFile = gArgs["initFile"]  # base file name

gInDir  = "../data"
gOutDir = "../results"

# A -> A.json, A.world.txt
gInitConfCompName, gWorldCompFileName = completeDataFileNames(gInitConfFile, gInDir)
# 1. A.json -> A_cont.json -> A_cont.json
if gInitConfFile.endswith('_cont'):
    newInitFileNameBase = gInitConfFile
else:
    newInitFileNameBase = gInitConfFile + '_cont'

# 2. …and newWorld and newConf cont out
gNewConfCompFileName, gNewWorldCompFileName = completeDataFileNames(newInitFileNameBase, gInDir)

gConf = replaceAndAddDefaultsArgsInConf(readInitConfFile(gInitConfCompName), gArgs)

gNumberOfSpecies = len(gConf["species"])


# SET UP CONVENIENT GLOBALS
#
gNumberOfCells   = gConf["NumberOfCells"]

gToSaveExcel     = gArgs["saveExcel"]
gGlobalExcel = None
gListOfAssociationOrigs = collectAssociationActors() # list of species starting association


gWithPartnerList = getListOfOrigGroups()




# SPECIES FOR WHICH ZERO MEANS END
#
noZeroListi = []
if gConf["noZero"]:
    if gConf["noZero"] == '-':
        noZeroListi = list(range(gNumberOfSpecies))
    else:
        noZeroListi = gConf["noZero"].split(',')
        try:
            noZeroListi = list(map(int,noZeroListi))
        except Exception:
            noZeroListi = iListFrom_idList(noZeroListi)

# print(f"noZeroListi: {noZeroListi}")


# OUTPUT DIRS AND FILENAMES
#
if gArgs["outDir"]:
    gToSaveExcel = True
    if gArgs["outDir"].startswith("/"):
        gOutDir = gArgs["outDir"]
    else:
        gOutDir = os.path.join(gOutDir, gArgs["outDir"])
    gGlobalExcel = os.path.join(gOutDir, "global.xlsx")

if not os.path.isdir(gInDir):
    Path(gInDir).mkdir()
    build_defaultjson(gInDir)

Path(gOutDir).mkdir(parents=True, exist_ok=True)
if not os.path.isdir(gOutDir):
    Path(gOutDir).mkdir(parents=True)

gInitConfFileNoSlashes = gInitConfFile.replace('/', '-')
gThedatetimeStamp = startTime.strftime("%Y%m%d-%H%M%S")
gOutFNameBase = gInitConfFileNoSlashes + '_' + gThedatetimeStamp


if gArgs["outFName"]:
    gToSaveExcel = True
    gGlobalExcel = os.path.join(gOutDir, gOutDir + ".xlsx")
    gOutFNameBase = gArgs["outFName"]

if gToSaveExcel:
    gExcelOut = os.path.join(gOutDir, gOutFNameBase + ".xlsx")
    if os.path.isfile(gExcelOut):
        baseName = gExcelOut[:gExcelOut.rfind('.')]
        i = 2
        while os.path.isfile(baseName+'.xlsx'):
            baseName = f"{baseName}_{i}"
            i += 1
        gExcelOut = baseName + '.xlsx'

if gConf['redirectStdout']:
    sys.stdout=open(os.path.join(gOutDir, gOutFNameBase + "_stdout.txt"), "a", encoding="utf-8")

if gConf['p']:
    gGlobalPDF = os.path.join(gOutDir, gInitConfFileNoSlashes + '_' + gThedatetimeStamp + ".pdf")
    if gArgs["verbose"]:
        print(f"PDF to save: {gGlobalPDF}")


if gOnLinux:
    print("Executing\n", " ".join(sys.argv))
checkConf(gConf)

printv("gWithPartnerList:", gWithPartnerList)
printv(gConf)
printv(gArgs)

myNumberOfProcess = ""
myDirectoryOfProcess = ""
for a in sys.argv:
    if a.startswith("--outFName"):
        myNumberOfProcess = a[a.index('=')+1:]
        if isInt(myNumberOfProcess):
            myNumberOfProcess = str(int(a[a.index('=')+1:]))
    if a.startswith("--outDir"):
        myDirectoryOfProcess = a[a.index('=')+1:]
myIDOfProcess = myDirectoryOfProcess + myNumberOfProcess



# exit(0)

# ########## #
# LETS DO IT #
# ########## #

# To avoid total random, we admit setting a specific seed
if gArgs["setRandomSeed"] != -1:
    npr.seed(int(gArgs["setRandomSeed"]))










# if "egoism" in gArgs:
#     calcEgoism()

if not gArgs["swallow"]:
    gMinEating = minNeededFood()

noNegger = np.vectorize(noNeg) #JAVI



# RESTART POINT
thereWasZeros = True
while thereWasZeros: # Sorry, we need this now to repeat the whole experiment

    gWorld, gStatsAnt, gStatsPost = newWorld()
    doInitialDistribution()

    # print("     %s" % " ".join([gConf["species"][i]["id"] for i in range(gNumberOfSpecies)]))

    print(f"\n\tUsing configuration:\n\t{gInDir}/{gInitConfFile}.json\n")
    print(" "*5, end='')
    for i in range(gNumberOfSpecies):
        print( f'{gConf["species"][i]["id"][:7]:8s}', end='')
    print()

    if gConf['p']:
        fig, ax, seqPobPrev = initPlot()

    lastTime = None
    if gOnLinux:
        lastTime = startTime






    # ########### #
    #  MAIN LOOP  #
    # ########### #

    for genNumber in range(1, gArgs["numGen"]+1):
        gStatsAnt.fill(0)
        gStatsPost.fill(0)

        doGrouping()

        print(f"{genNumber:3d}: {gWorld[:,:].sum(axis=0)} Tot: {gWorld[:,:].sum():8d}")

        if gConf['p']:
            fig, ax, seqPobPrev = addPlot(fig, ax, seqPobPrev, genNumber)


        if anyZeros(noZeroListi):
            thereWasZeros = True
            break

        doAssociationQueueAndConsume()
        gWorld = noNegger(gWorld) #JAVI
        if gArgs["verbose"]:
            checkNegs()


        doUngroup()
        doDistribute()

        if gArgs["saveWorld"]:
            saveWorld(genNumber)

        if gToSaveExcel:
            saveExcel(genNumber)

        if os.path.isfile(TO_ABORT_EXECUTION_FILENAME) or myIDOfProcess and \
           os.path.isfile(TO_ABORT_EXECUTION_FILENAME  +  myIDOfProcess):
            print(f"""
    BROKE BECAUSE THERE IS A FILE "{TO_ABORT_EXECUTION_FILENAME}"
    Remove that file to avoid this to be repeated:
            rm {TO_ABORT_EXECUTION_FILENAME}..
    """)
            break

        if gOnLinux:
            timeNow = datetime.now()
            deltaTimeGen = timeNow - lastTime
            lastTime = timeNow
            timeRecording(timeNow, genNumber, gArgs["numGen"], deltaTimeGen, myIDOfProcess)
            # [x for x in sys.argv if x.startswith("--outDir=")][-1].split('=')[-1])
            if timeNow > REDZONETIME:
                break

        sys.stdout.flush()

    thereWasZeros = len(noZeroListi) > 0 and anyZeros(noZeroListi)
    # if thereWasZeros:
    #     if 'plt' in locals() and gConf['p']:
    #         print(type(plt))
    #         print(plt)
    #         print("Deleting previous plot")
    #         plt.close()

# END OF WHILE TRUE


# END RUNNING ####################################


saveConf()

if gToSaveExcel:
    gExcelWorkbook.close()

if gConf['redirectStdout']:
    sys.stdout.close()


if gConf['p']:
    plt.savefig(gGlobalPDF)
    if genNumber >= gArgs["numGen"]:
        # it got broken by a zero
        input("Press Enter key  ⮐  to exit")

newArgv = " ".join(sys.argv)

if gOnLinux and genNumber < gArgs["numGen"]:
    genLeft = gArgs["numGen"] - genNumber
    theArgv = map(add_cont, sys.argv)

    newArgv = " ".join(theArgv) + ' --numGen=' + str(genLeft)
    with open(COMMANDSLEFT, "a", encoding="utf-8") as f:
        print(newArgv, file=f)
