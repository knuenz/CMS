#MERGE TTree / RooDataSets

import os
import array
import sys
import re

import ROOT as R

###############################################
## define output file name
def mergeChain(outfile, files, pathOnFs):
    ## choose Name of TCHAIN, should be the same as input TTree name
    chain = R.TChain("data","data")
    print "---- Start processing TTrees ----"
    ## merge trees into chain
    for i in range( len(files) ):
        print "-- Processing " + files[i]
        chain.Add( pathOnFs + files[i] )
     
        print "-- -- Added file: ", files[i]
            
    print "----------"
    print "----- ( I am merging your ", len(files)," TTrees now...PATIENCE...) -----"
    chain.Merge(outfile)
    print "----- merging...SUCCEEDED... ;) -----"

print "---- (Start Script) ----"    

## choose input files
print "---- load Files list ----"

files_1 = [
#your other files here

'TTree_Upsilon3S_Onia2MuMu-v6_10_1_Z0i.root',
'TTree_Upsilon3S_Onia2MuMu-v6_11_1_bG0.root',
'TTree_Upsilon3S_Onia2MuMu-v6_12_1_7oA.root',
'TTree_Upsilon3S_Onia2MuMu-v6_13_1_9Zz.root',
'TTree_Upsilon3S_Onia2MuMu-v6_14_1_yUP.root',
'TTree_Upsilon3S_Onia2MuMu-v6_15_1_KDH.root',
'TTree_Upsilon3S_Onia2MuMu-v6_16_1_elz.root',
'TTree_Upsilon3S_Onia2MuMu-v6_17_1_ugb.root',
'TTree_Upsilon3S_Onia2MuMu-v6_18_1_OzZ.root',
'TTree_Upsilon3S_Onia2MuMu-v6_19_1_flW.root',
'TTree_Upsilon3S_Onia2MuMu-v6_1_1_OWC.root',
'TTree_Upsilon3S_Onia2MuMu-v6_20_1_1FG.root',
'TTree_Upsilon3S_Onia2MuMu-v6_21_1_lDk.root',
'TTree_Upsilon3S_Onia2MuMu-v6_2_1_oVU.root',
'TTree_Upsilon3S_Onia2MuMu-v6_3_1_bT0.root',
'TTree_Upsilon3S_Onia2MuMu-v6_4_1_Njq.root',
'TTree_Upsilon3S_Onia2MuMu-v6_5_1_G4d.root',
'TTree_Upsilon3S_Onia2MuMu-v6_6_1_1Am.root',
'TTree_Upsilon3S_Onia2MuMu-v6_7_1_erI.root',
'TTree_Upsilon3S_Onia2MuMu-v6_8_1_BFB.root',
'TTree_Upsilon3S_Onia2MuMu-v6_9_1_MFT.root'


]

#files_2 = [
#your other files here
#]

#files_3 = [
#your other files here
#]

## folder name of Input files
#pathOnFs_1 = "rfio:///"
#pathOnFs_2 = "rfio:///" 
#pathOnFs_3 = "rfio:///"

###############################################
#
# mergeChain(outfile, files, pathOnFs)
#
###############################################


#mergeChain("TTree_hltReport_1.root", files_1, pathOnFs_1)
#mergeChain("TTree_hltReport_2.root", files_2, pathOnFs_2)
#mergeChain("TTree_hltReport_3.root", files_3, pathOnFs_3)

###############################################
#files = ["TTree_hltReport_1.root",
#         "TTree_hltReport_2.root",
#         "TTree_hltReport_3.root"]


###############################################
#mergeChain(NameOfOutFile, ListOfInputFiles, pathOnFileSystem)
pathOnFs_1 = "file://Upsilon3S/"
mergeChain("TTree_MCprompt_Upsilon3S_Fall10.root", files_1, pathOnFs_1)

print "----- (Tree created+saved) -----"

###############################################
## delete TChain objects
#print "----- (delete TChain objects) -----"
#os.system("rm *_tree.root")

## done   
print "----- (THE END) -----"


