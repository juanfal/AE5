#!/usr/bin/env python3
# encoding: utf-8
# rae5.py
# juanfc 2019-07-08
#        2020-02-17     Simple copy from ae3
#
# Builds the ranges for each parameter and make the corresponding
# ae5.py calls

# 2020-07-15
#   Fallos desde el principio.  Se construía un comandos.txt y lanzaba sbatch,
#   lo cual funciona si no tocas el comandos.txt.  Al renombrar el
#   comandos.txt, cuando finalmente entraba el sbatch,  no veía el
#   comandos.txt.

#   Arreglado del tirón creando un comandos.personalizado para cada rae5 que
#   se lanza, y así no se toca más

# TODO: --species can have more than one range… so it has to be expanded by itself
#       naming the output folder and each file

# SOLUTION:
#       --species has been treated especially.  As an special parameter, sorry
#       its ranges do not admit , yet
#       rae5.py --algo=[1:3] --uno=[a] lacagaste --none=nada --species=1,a=[2:4],b=[2:6]

# FOR LISTS INSIDE RANGES
#
#                               USE ;
#
#       rae5.py --algo=[1:3] --uno=[a] lacagaste --none=nada --species="1,2,a=[2:4;7],b=[2:6]"

# WORKING!:
#       rae5.py --outDir=kkkkkkkk 0619Amen2 --species="-1,IndirectOffspring=[-8:-4]"
# --numGen int
# --verbose
# --saveExcel

# i - initFile
# r - --setRandomSeed int

# v - --varia
# n - --NumberOfCells int
# s - --NumberOfRsrcsInEachCell int
# d - --Distribution 'str'
#     --species 'str'
#     "NumberOfItems": 100,
#     "DirectOffspring": 5,

#     "GroupPartners":   "A",
#     "PhenotypicFlexibility": 0.0,

#     "AssociatedSpecies":   "",
#     "IndirectOffspring": 2,
#     "FitnessVariationLimit": 0





# --outDir 'str'
# --outFName 'str'


# print(sys.argv)
# status = os.system('./ae5.py 0619Amen2 --species="-1,IndirectOffspring=-8"')
# print(status)

# print(buildList("1,2,4:8,100,1000:1010"))
# print(expandArg("--algo=[1,2,4:8,100,1000:1010]"))
# print(expandArg("--algo=[1,2,3]"))
# print(expandArg("--algo=[a]"))
# print(expandArg("--algo=a=[1:6],b=[-1,8:10]"))
# print(expandArg("--algo=nada"))

# print(ranges)
# ranges = [expandArg("--algo=[1:3]"), expandArg("--uno=[a]"), expandArg("lacagaste"), expandArg("--none=nada"), expandArg("--vars=[1,2,4:8]")]
# --algo=[1:3] --uno=[a] lacagaste --none=nada --vars=[1,2,4:8]


import sys
import os
import subprocess
from pathlib import Path
import shutil
import argparse
import re
from datetime import datetime
import time


def frange(x, y, jump):
  while x < y:
    yield x
    x += jump

def buildList(pyExpr):
    r = []
    tmp = pyExpr.split(",")
    for item in tmp:
        if item.count(':') == 1:
            first, last = item.split(':')
            if '.' in first or '.' in last:
                first = float(first)
                last = float(last)
            else:
                first = int(first)
                last = int(last)
            r +=  list(map(str,list(frange(first, last))))
        elif item.count(':') == 2:
            first, last, step = item.split(':')
            if '.' in first or '.' in last or '.' in step:
                first = float(first)
                last = float(last)
                step = float(step)
            else:
                first = int(first)
                last = int(last)
                step = int(step)
            r +=  list(map(str,list(frange(first, last, step))))
        else:
            r.append(item)

    return r

def expandArg(anArg):
    r = []
    if "[" not in anArg or anArg.startswith("GroupPartners"):
        r = ["", [anArg]]
    else:
        r = anArg.strip()[:-1].split("[")
        # print(r)
        r = [r[0], buildList(r[1])]

    return r

def build(base, ranges, sep):
    tsep = ""
    commandList = [base]
    for theName, theRange in ranges:
        temp = []
        for i in theRange:
            for prev in commandList:
                temp += [prev + tsep + theName + str(i)]
        commandList = temp
        tsep = sep
    return commandList


def quoteSpecies(commandList):
    r = []
    for command in commandList:
        r.append( re.sub(r'--species=([^ ]+)', r'--species="\1"', command))
    return r


def fileNumbering(n):
    return "%04d" % n

outDir = datetime.now().strftime('%Y%m%d-%H%M%S')  # easying no range execs

# def commandsBack(cmds):
#     os.rename(cmds, "comandos_"+str(outDir)+"_"+datetime.now().strftime('%Y%m%d-%H%M%S')+".txt")


if len(sys.argv) > 1 and "-h" == sys.argv[1]:
    print("""Examples:
          rae5.py --outDir=kkkkkkkk 0619Amen2 --species="-1;IndirectOffspring=[-8:-4,2];0;DirectOffspring=[1,2]"
          rae5.py --outDir=kkkkkkkk 0619Amen2 --NumberOfCells=[100:10000:100] --species="-1;IndirectOffspring=[-8:-4,2];0;DirectOffspring=[1,2]"
          rae5.py --outDir=kkkkkkkk -t 0619Amen2 --NumberOfCells=[100:103] --species="-1;IndirectOffspring=[-8:-4,2];0;DirectOffspring=[1,2]"
          rae5.py --outDir=kkkkkkkk -t 0619Amen2 --NumberOfCells="[100:110:3,8,10:14:2]" --species="-1;IndirectOffspring=[-8:-4:2,2];0;DirectOffspring=[1,2]"


          You must express an output directory in the first parameter:
                    --outDir=nameOfTheDirectory
          if it it starts in / it'd be an absolute path
          in other case, it will be created inside 'results' directory
          UPDATE: if none provided it defaults to a default out dir, to ease no-range rae

          You can add the argument:

                --time=7-0            # for allotting the maximum 7 days
                --time=0:30:00        # for allotting 30 minutes

          IMPORTANT:
              Observe in the last examples:
              To express a list of different values or ranges, separate them with ;
              but to prevent the Terminal shell from an undesired interpretation, surround it in ';' or ";"
    """)
    exit(0)


onLinux = sys.platform == "linux" or sys.platform == "linux2"

COMMANDSLEFT = "comandosRestan.txt"
TIMESLURM="7-0"

if onLinux:
    if os.path.isfile(COMMANDSLEFT):
        print("****** RESUMING PREVIOUS UNFINISHED EXECUTION *****")
        print()
        with open(COMMANDSLEFT) as f:
            n = sum(1 for _ in f)
        x = subprocess.call(["sbatch",
                        "--array=1-"+str(n),
                        "--time="+TIMESLURM,
                        "arbatch.sh",
                        COMMANDSLEFT
                        ])
        # if not x:
        #     shutil.move(COMMANDSLEFT, COMMANDSLEFT+'_'+datetime.now().strftime('%Y%m%d-%H%M%S')+'.txt')
        sys.exit(x)




jsonFile = ''
testMode = False
ranges = []
for arg in sys.argv[1:]:

    if "-t" == arg:
        testMode = True
        continue

    if not arg.startswith("-"): # json ??
        jsonFile = arg
        if jsonFile.endswith('.json'):
            jsonFile = jsonFile[:jsonFile.rfind('.')]
        continue

    if arg.startswith("--time="):
        _, TIMESLURM = arg.split("=")
        continue

    if arg.startswith("--outDir="):
        _, outDir = arg.split("=")
        continue

    if not arg.startswith("--species"):
        ranges.append(expandArg(arg))
    else:
        #print("*****************", arg[10:])
        ranges.append(["", build("--species=", map(expandArg, arg[10:].split(";")), ";" )])

commandList = build("./ae5.py ", ranges, " ")
commandList = quoteSpecies(commandList)


if not outDir:
    print("You must give an output directory in parameters:")
    print("          --outDir=nameOfTheDirectory")
    print("if the dir starts in / it'd be an absolute path")
    print("in other case, it will be created inside 'results' directory")
    exit(1)

outDir = Path(outDir)
if not outDir.root:
    finalOutDir = Path("results") / outDir
if not testMode:
    if not os.path.isdir(str(finalOutDir)):
        finalOutDir.mkdir(parents=True)
    # finalOutDir.mkdir(parents=True, exist_ok=True)

outListFName = finalOutDir / "_list.txt"

n = 1

if not testMode:
    listing = open(str(outListFName), "w")

ae5CommandsFile = "comandos_"+str(outDir)+"_"+datetime.now().strftime('%Y%m%d-%H%M%S')+".txt"

comandos = open(ae5CommandsFile, "w")  #  adding lines for next sbatch arbatch

for command in commandList:
    fname = fileNumbering(n)
    command += " --outDir=" + str(outDir)
    command += " --outFName=" + fname
    command += " --setRandomSeed=1"
    command += " --redirectStdout"
    # command += " --sendTelegram"
    if jsonFile:
        newjsonFile = jsonFile+"_"+fname
        command += " "+newjsonFile
        if jsonFile.startswith('/'):
            realjsonFile = jsonFile
            realnewjsonFile = newjsonFile
        else:
            realjsonFile = "data/"+jsonFile
            realnewjsonFile = "data/"+newjsonFile

        if testMode:
            print("to copy: "+ realjsonFile+'.json'+" -> "+realnewjsonFile+'.json')
        else:
            try:
                shutil.copy2(realjsonFile+'.json', realnewjsonFile+'.json')
            except Exception as e:
                print("ERROR: File do not exist: "+jsonFile+'.json')


    print(command)
    if not testMode:
        print(command, file=listing)
        print(command, file=comandos)
        # os.system(command)

    n += 1

if not testMode:
    listing.close()
comandos.close()

if not testMode:
    if onLinux:
        # print(["sbatch", "--array=1-"+str(len(commandList)), "arbatch.sh"])
        subprocess.call(["sbatch",
                        "--array=1-"+str(len(commandList)),
                        "--time="+TIMESLURM,
                        "arbatch.sh",
                        ae5CommandsFile
                        ])
        # print("Dando 3 segundos…")
        # time.sleep(3)


