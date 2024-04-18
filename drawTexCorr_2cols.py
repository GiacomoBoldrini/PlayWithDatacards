#!/usr/bin/env python

#
# original code starting from https://github.com/amassiro/HiggsAnalysis-CombinedLimit
#

import re
import os.path
from math import *
from optparse import OptionParser


os.environ['TERM'] = 'linux'


parser = OptionParser()
parser.add_option("-f", "--format",  dest="format",   default="html",          type="string",  help="Format for output number")
parser.add_option("-m", "--mass",    dest="mass",     default=0,               type="float",  help="Higgs mass to use. Will also be written in the Workspace as RooRealVar 'MH'.")
parser.add_option("-D", "--dataset", dest="dataname", default="data_obs",      type="string",  help="Name of the observed dataset")
parser.add_option("-a", "--all",     dest="all",      default=False, action='store_true',  help="Report all nuisances (default is only lnN)")
parser.add_option("", "--noshape",   dest="noshape",  default=False, action='store_true',  help="Counting experiment only (alternatively, build a shape analysis from combineCards.py -S card.txt > newcard.txt )")
parser.add_option("-o", "--output",  dest="output",   default="minitable.tex", type="string",  help="Summary table")
parser.add_option("", "--blind",     dest="blind",    default=False, action='store_true',  help="blind (default is false)")
parser.add_option("", "--legend",    dest="legend",   default="",              type="string",  help="Addition information in caption")

                      

(options, args) = parser.parse_args()
options.stat = False
options.bin = True # fake that is a binary output, so that we parse shape lines
options.out = "tmp.root"
options.fileName = args[0]
options.cexpr = False
options.fixpars = False
options.libs = []
options.verbose = 0
options.poisson = 0
options.nuisancesToExclude = []
options.noJMax = True
options.allowNoSignal = True
options.modelparams = []
options.optimizeTemplateBins = False
options.format = "html"

# import ROOT with a fix to get batch mode (http://root.cern.ch/phpBB3/viewtopic.php?t=3198)
import sys
sys.argv = [ '-b-']
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")

from HiggsAnalysis.CombinedLimit.DatacardParser import *
from HiggsAnalysis.CombinedLimit.ShapeTools     import *
if options.fileName.endswith(".gz"):
    import gzip
    file = gzip.open(options.fileName, "rb")
    options.fileName = options.fileName[:-3]
else:
    file = open(options.fileName, "r")

DC = parseCard(file, options)
if not DC.hasShapes: DC.hasShapes = True
MB = ShapeBuilder(DC, options)
if not options.noshape: MB.prepareAllShapes()

def commonStems(list, sep="_"):
    hits = {}
    for item in list:
        base = ""
        for token in item.split(sep):
            base = base + sep + token if base else token
            if base not in hits: hits[base] = 0
            hits[base] += 1
    veto = {}
    for k,v in hits.iteritems():
        pieces = k.split(sep)
        for i in xrange(1, len(pieces)):
            k2 = "_".join(pieces[:-i])
            if hits[k2] == v: 
                veto[k2] = True
            else:
                veto[k] = True
    ret = []
    for k,v in hits.iteritems():
       if k not in veto: ret.append((k,v))
    ret.sort()
    return ret 
    
report = {}; errlines = {}
for (lsyst,nofloat,pdf,pdfargs,errline) in DC.systs:
    print(lsyst, pdf, len(errline))
    #if not options.all and pdf != "lnN": continue
    if not len(errline) : 
       print("OH NO")
       continue
    types = []
    minEffect, maxEffect = 999.0, 1.0
    processes = {}
    channels  = []
    errlines[lsyst] = errline
    for b in DC.bins:
        types.append(pdf)
        channels.append(b)
        for p in DC.exp[b].iterkeys():
            if errline[b][p] == 0: continue
            processes[p] = True
            if "shape" in pdf :
                vals = [] 
                try:
                   objU,objD,objC = MB.getShape(b,p,lsyst+"Up"), MB.getShape(b,p,lsyst+"Down"), MB.getShape(b,p)
                   if objC.InheritsFrom("TH1"): valU,valD,valC =  objU.Integral(), objD.Integral(), objC.Integral()
                   elif objC.InheritsFrom("RooDataHist"): valU,valD,valC =  objU.sumEntries(), objD.sumEntries(), objC.sumEntries()

                except:
                   objC = None
                   valC = 0 

                if valC!=0: 
                        errlines[lsyst][b][p] = "%.3f/%.3f"%(valU/valC,valD/valC)
                        vals.append(valU/valC)
                        vals.append(valD/valC)
                else: 
                        errlines[lsyst][b][p] = "NAN/NAN"
                        vals.append(1.)
                        vals.append(1.)
            else: vals = errline[b][p] if type(errline[b][p]) == list else [ errline[b][p] ]
            for val in vals:
                if val < 1 and val !=0 : val = 1.0/val
                minEffect = min(minEffect, val)
                maxEffect = max(maxEffect, val)
    channelsShort = commonStems(channels)
    types = ",".join(set(types))
    report[lsyst] = { 'channels':channelsShort, 'bins' : channels, 'processes': sorted(processes.keys()), 'effect':("%5.3f"%minEffect,"%5.3f"%maxEffect), 'types':types }

# Get list
names = report.keys() 
if "brief" in options.format:
    names = [ k for (k,v) in report.iteritems() if len(v["bins"]) > 1 ]
# alphabetic sort
names.sort()

print(names)
# now re-sort by category (preserving alphabetic sort inside)
namesCommon = [ n for n in names if re.match(r"(pdf_|QCD|lumi|UE|BR).*", n) ]
namesCMS1   = [ n for n in names if re.match(r"CMS_(eff|scale|fake|res|trig).*", n) ]
namesCMS2   = [ n for n in names if re.match(r"CMS_.*", n) and n not in namesCMS1 ]
namesRest   = [ n for n in names if n not in namesCommon and n not in namesCMS1 and n not in namesCMS2 ]
names = namesCommon + namesCMS1 + namesCMS2 + namesRest


# channels = ["VBSOS", "VBSSS_WZ", "VBS4l", "VBSWV", "VBSZV", "VBSZZ2l2nu", "VBSSS_tau"]
channels = ["VBSOS", "VBSSS_WZ", "VBS4l", "VBSWV", "VBSZV"]

n_tables = 20

import numpy as np

names_dict = np.array_split(names, n_tables)

cols = " & ".join("\\textbf{" + i.replace("_", "\_") + "}" for i in channels) + " & ".join("\\textbf{" + i.replace("_", "\_") + "}" for i in channels)
print("\\begin{center}")
print("   \\renewcommand{\\arraystretch}{0.8}")
print("   \\begin{tabular}{ |c|c|" + " ".join("c|" for i in channels) + " }")
print("   \\hline")
print("   \\textbf{Nuisance} & \\textbf{type} & " + " & ".join("\\textbf{" + i.replace("_", "\_") + "}" for i in channels) + " \\textbf{Nuisance} & \\textbf{type} & " + " & ".join("\\textbf{" + i.replace("_", "\_") + "}" for i in channels) +  " \\\ ")
print("   \\hline")
for idx,names_list in enumerate(names_dict):

   
   idx_ = 0
   for nuis in names_list:
           val = report[nuis]
           type_ = val["types"]
           if "?" in type_: type_ = "shape"
           if idx_ == 0:
              l = "   " + nuis.replace("_", "\_") + " & " + type_ + " & "
           else:
              l += " & " +  "   " + nuis.replace("_", "\_") + " & " + type_ + " & "

           dings = [0]*len(channels)
           for x in sorted(val["bins"]):
              if "VBSOS" in x: chan = "VBSOS"
              if "VBSjjlnu" in x: chan = "VBSWV"
              if "VBSjjll" in x: chan = "VBSZV"
              #if "ZZ" in x and "ZZ2e2nu" not in x: chan = "ZZ"
              #if "ZZ2e2nu" in x: chan = "ZZ2e2nu"
              if "VBS4l" in x: chan = "VBS4l"
              if "VBSSS_WZ" in x: chan = "VBSSS_WZ"
              if "VBSSS_tau" in x: chan = "VBSSS_tau"
              if "VBSZZ2l2nu" in x: chan = "VBSZZ2l2nu"
             

              if len(", ".join(["%s(%s)"%(k,v) for (k,v) in errlines[nuis][x].iteritems() if v != 0])) != 0:
                 dings[channels.index(chan)] = 1
        
           for j in dings:
              if j == 1: l += " \\checkmark & "
              else: l += " - & "
   
           l = l[:-2]
           # l += " \\\ "

           idx_ += 1
           if idx_ == 2:
              l += " \\\ "
              print(l)
              idx_ = 0
  
   
   
   
print("   \\hline")
print("   \\end{tabular}")
print("\\end{center}")
