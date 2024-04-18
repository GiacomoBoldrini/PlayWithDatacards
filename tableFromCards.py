#!/usr/bin/env python

import re
from sys import argv
import os.path
from optparse import OptionParser
from math import sqrt,fabs
parser = OptionParser()
parser.add_option("-s", "--stat",   dest="stat",          default=False, action="store_true")  # ignore systematic uncertainties to consider statistical uncertainties only
parser.add_option("-S", "--force-shape", dest="shape",    default=False, action="store_true")  # ignore systematic uncertainties to consider statistical uncertainties only
parser.add_option("-a", "--asimov", dest="asimov",  default=False, action="store_true")
parser.add_option("-m", "--mass", dest="mass",  default=125, type="float")


(options, args) = parser.parse_args()
options.bin = True # fake that is a binary output, so that we parse shape lines
options.noJMax = False
options.nuisancesToExclude = ''
options.allowNoSignal = True

from DatacardParser import *

DC = parseCard(file(args[0]), options)
nuisToConsider = [ y for y in DC.systs ]

errors = {}
for nuis in nuisToConsider:
    if nuis[2] == 'gmN': gmN = nuis[3][0]
    else               : gmN = 0
    for channel in nuis[4]:
        #print channel
        if channel not in errors.keys(): errors[channel] = {}
        for process in nuis[4][channel]:
            newError = 0.
            if nuis[2] == 'gmN': gmN = nuis[3][0]
            if nuis[4][channel][process] == 0: continue
            #print nuis[2],gmN
            if gmN != 0:
                newError = nuis[4][channel][process] * sqrt(gmN) / DC.exp[channel][process]
            else:
                #print nuis[4][channel][process]
                if not isinstance ( nuis[4][channel][process], float ) :
                    # [0.95, 1.23]
                #if len(nuis[4][channel][process]) == 2 :
                    newError = fabs((nuis[4][channel][process][1]-nuis[4][channel][process][0])/2.)   # symmetrized
                else : 
                    newError = fabs(1-nuis[4][channel][process])
            if process in errors[channel].keys():
                errors[channel][process] += newError*newError
            else:
                errors[channel][process] = newError*newError

for channel in errors:
    for process in errors[channel]:
        errors[channel][process] = sqrt(errors[channel][process])

for x in DC.exp:
    for y in DC.exp[x]:
        print "%10s %10s %10.2f +/- %10.2f (%10.2f \\%%)" % (x,y,DC.exp[x][y],DC.exp[x][y]*errors[x][y],errors[x][y]*100)




size = "footnotesize"

print "\n"
print "========================="
print "\n latex style \n"

print "\\begin{table}[h!]\\begin{center}"
print ("\\%s{\\begin{tabular}{" % size)

print ("c|"),
for channel in DC.exp:
    print ("c |"),
print "} \\hline"

print (" "),
for channel in DC.exp:
    print ("& %13s " % channel.replace('_', '-')),
print ("\\\\ \\hline")

signals     = DC.list_of_signals()
backgrounds = DC.list_of_backgrounds()

totsig    = {}
errtotsig = {}
totbkg    = {}
errtotbkg = {}

for s in signals :
    print (" %13s " % s.replace('_', '-')),
    for channel in DC.exp:
        if s in DC.exp[channel].keys(): # possible that some backgrounds appear only in some channels
            print (" & %10.2f +/- %10.2f (%10.0f \\%%) " % (DC.exp[channel][s],DC.exp[channel][s]*errors[channel][s],errors[channel][s]*100)),
            if channel not in    totsig.keys():    totsig[channel] = 0.0
            if channel not in errtotsig.keys(): errtotsig[channel] = 0.0
            totsig[channel]    = totsig[channel]    + DC.exp[channel][s]
            errtotsig[channel] = errtotsig[channel] + (DC.exp[channel][s]*errors[channel][s] * DC.exp[channel][s]*errors[channel][s])
            #print " <<<  sqrt errtotsig[",channel,"] = ", sqrt(errtotsig[channel]) , ">>> "
        else :
            print (" & - "),

    print ("\\\\")


print ("\\hline")
print (" %13s " % "signal"),
for channel in DC.exp:
    errtotsig[channel] = sqrt(errtotsig[channel])
    print (" & %10.2f +/- %10.2f (%10.0f \\%%) " % (totsig[channel],errtotsig[channel],errtotsig[channel]/totsig[channel]*100)),
print ("\\\\")
print ("\\hline")


for b in backgrounds :
    print (" %13s " % b.replace('_', '-')),
    for channel in DC.exp:
        if b in DC.exp[channel].keys(): # possible that some backgrounds appear only in some channels
            print (" & %10.2f +/- %10.2f (%10.0f \\%%) " % (DC.exp[channel][b],DC.exp[channel][b]*errors[channel][b],errors[channel][b]*100)),
            if channel not in    totbkg.keys():    totbkg[channel] = 0.0
            if channel not in errtotbkg.keys(): errtotbkg[channel] = 0.0
            totbkg[channel]    = totbkg[channel]    + DC.exp[channel][b]
            errtotbkg[channel] = errtotbkg[channel] + (DC.exp[channel][b]*errors[channel][b] * DC.exp[channel][b]*errors[channel][b])
        else :
            print (" & - "),

    print ("\\\\")


print ("\\hline")
print (" %13s " % "background"),
for channel in DC.exp:
    errtotbkg[channel] = sqrt(errtotbkg[channel])
    print (" & %10.2f +/- %10.2f (%10.0f\\%%) " % (totbkg[channel],errtotbkg[channel],errtotbkg[channel]/totbkg[channel]*100)),
print ("\\\\")
print ("\\hline")



print "\\end{tabular}"
print "}"
print "\\end{center}"


print "\n\n\n"


print "\\end{table}"


print "========================="
print "\n\n\n"


