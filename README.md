PlayWithDatacards
=================

Tools for datacards

Where:

    /afs/cern.ch/user/a/amassiro/Limit/PlayWithDatacards


# Modification tools

look here:

    https://github.com/amassiro/ModificationDatacards



# Transform datacard into a table


    python tableFromCards.py  /afs/cern.ch/user/a/amassiro/public/xLatinos/ww/WWDFcut0jet/hww-19.36fb.mH125.of_0j_shape.txt

e.g.

    for ww
    bash scripts/doDatacard2TableWW.sh



# Combine cards

    cd /afs/cern.ch/user/a/amassiro/scratch0/VBF/Limit/CMSSW_6_1_0/src
    export SCRAM_ARCH=slc5_amd64_gcc462
    cmsenv
    combineCards.py ww_sf_0j=/afs/cern.ch/user/a/amassiro/public/xLatinos/ww/WWSFcut0jet/hww-19.36fb.mH125.sf_0j_shape.txt   \
                    ww_of_0j=/afs/cern.ch/user/a/amassiro/public/xLatinos/ww/WWDFcut0jet/hww-19.36fb.mH125.of_0j_shape.txt   \
                    ww_sf_1j=/afs/cern.ch/user/a/amassiro/public/xLatinos/ww/WWSFcut1jet/hww-19.36fb.mH125.sf_1j_shape.txt   \
                    ww_of_1j=/afs/cern.ch/user/a/amassiro/public/xLatinos/ww/WWDFcut1jet/hww-19.36fb.mH125.of_1j_shape.txt > \
                    /afs/cern.ch/user/a/amassiro/public/xLatinos/ww/hww-19.36fb.mH125.txt

