#
# Read chic SelEvents in python 
# and dump RooDataSet, make histos
#
# $Id: chibanalysis_new.py,v 1.1 2011/12/10 01:34:22 giordano Exp $

from ROOT import TTree, TFile,TH1F,TH2F, gROOT, TMath
from ROOT import RooDataSet,RooRealVar,RooArgSet
#import ChiSquaredProbability
from optparse import OptionParser

import sys
import ROOT

from DataFormats.FWLite import Events, Handle

#from ChibInputFileList import inputfilesA, inputfilesB 

from math import sqrt, pi as M_PI

def deltaPhi( phi1,  phi2):
    result = phi1 - phi2;
    while (result > M_PI):
        result -= 2*M_PI;
    while (result <= -M_PI):
        result += 2*M_PI;
    return result;

def deltaR2(eta1, phi1,  eta2,  phi2): 
    deta = eta1 - eta2;
    dphi = deltaPhi(phi1, phi2);
    return deta*deta + dphi*dphi;

def deltaR( eta1,  phi1,  eta2,  phi2):
    return sqrt(deltaR2 (eta1, phi1, eta2, phi2));

def main(options,args):

    inputfiles = [ 'default', 'default' ]
#    inputfiles = [ 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_2011A.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_2011B.root' ]
#    inputfiles = [ 'file:/scratch/knuenz/Polarization/RootInput/Chib/chib_2011A_gtV21A_2082pb.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chib_2011B_gtV21A_2613pb.root' ]
#    inputfiles = [ 'file:/scratch/knuenz/Polarization/RootInput/Chib/chicstep5_chibV1_A.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chicstep5_chibV1_B.root' ]
#    inputfiles = [ 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_v13d_A_17Jan.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_v13d_B_17Jan.root' ]
#    inputfiles = [ 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_chibV1_A_20Jan.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_chibV1_B_20Jan.root' ]
#    inputfiles = [ 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_pi0Rej_false_A_27jan.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_pi0Rej_false_B_27jan.root' ]
    inputfiles = [ 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_pi0Rej_false_A_27jan.root' ]

#    gROOT.ProcessLine('.L /afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/CommonTools/Statistics/src/ChiSquaredProbability.cc+')
#    gROOT.ProcessLine('.L /afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/CommonTools/Statistics/src/IncompleteGammaComplement.cc+')
#    from ROOT import ChiSquaredProbability
#    from ROOT import IncompleteGammaComplement
    Y1Smass0=9.46
    Y2Smass0=10.023
    Y3Smass0=10.355
    Ymass_a=0.058
    Ymass_b=0.047
    Ymass_c=0.22
            
    events= Events(inputfiles)
    
    
    # candidate Chi
    chicCand_h  = Handle ("vector<reco::CompositeCandidate>")
    chicCand_l = ( "chib","chicCompCand","chibs5" )
    
    #candidate J/Psi
    jpsiCand_h  = Handle ("vector<pat::CompositeCandidate>")
    jpsiCand_l = ( "chib","jpsiCompCand" ,"chibs5")
    
    jpsiCandPAT_h  = Handle ("vector<pat::CompositeCandidate>")
    jpsiCandPAT_l = ( "onia2MuMuPatTrkTrk","" )
    
    
    #candidate gamma
    gammaCand_h  = Handle ("vector<reco::CompositeCandidate>")
    gammaCand_l = ( "chib","gammaCompCand" ,"chibs5")
    
    convCand_h = Handle ("vector<reco::Conversion>")
    convCand_l = ( "allConversions","" )
    
    
    # Create histograms, etc.
    gROOT.SetBatch()        # don't pop up canvases
    gROOT.SetStyle('Plain') # white background
    
    #fill a RooDataSet
    invm1S     = RooRealVar("invm1S",   "invm1S",9,50)
    invm2S     = RooRealVar("invm2S",   "invm2S",9,50)
    invm3S     = RooRealVar("invm3S",   "invm3S",9,50)
    jpsipt   = RooRealVar("jpsipt", "jpsipt",0,100)
    jpsimass   = RooRealVar("jpsimass", "jpsimass",8,12)
    gammapt   = RooRealVar("gammapt", "gammapt",0,100)
    deltaRChiJpsi  = RooRealVar("deltaRChiJpsi","deltaRChiJpsi",0,1)
    deltaRJpsig    = RooRealVar("deltaRJpsig","deltaRJpsig",0,1)
    jpsieta  = RooRealVar("jpsieta","jpsieta",-5,5)
    ctpv     = RooRealVar("ctpv",   "ctpv",-5,5)
    ctpverr  = RooRealVar("ctpverr","ctpverr",0,10)
    ctbs     = RooRealVar("ctbs",   "ctbs",-5,5)
    ctbserr  = RooRealVar("ctbserr","ctbserr",0,10)
    ctpvsig     = RooRealVar("ctpvsig",   "ctpvsig",0,1000)
    ctbssig     = RooRealVar("ctbssig",   "ctbssig",0,1000)
    weight   = RooRealVar("weight", "weight",0,1)
    Rconv    = RooRealVar("Rconv","Rconv",0,100)
    jpsiVprob    = RooRealVar("jpsiVprob","jpsiVprob",0,1.5)
    Y1Smass_nSigma    = RooRealVar("Y1Smass_nSigma","Y1Smass_nSigma",-100,100)
    Y2Smass_nSigma    = RooRealVar("Y2Smass_nSigma","Y2Smass_nSigma",-100,100)
    Y3Smass_nSigma    = RooRealVar("Y3Smass_nSigma","Y3Smass_nSigma",-100,100)
    jpsipx   = RooRealVar("jpsipx", "jpsipx",-100,100)
    jpsipy   = RooRealVar("jpsipy", "jpsipy",-100,100)
    jpsipz   = RooRealVar("jpsipz", "jpsipz",-100,100)
    gammapx   = RooRealVar("gammapx", "gammapx",-100,100)
    gammapy   = RooRealVar("gammapy", "gammapy",-100,100)
    gammapz   = RooRealVar("gammapz", "gammapz",-100,100)
    vertexChi2ProbGamma   = RooRealVar("vertexChi2ProbGamma", "vertexChi2ProbGamma",0,1)
    
    argSet = RooArgSet()
    argSet.add(invm1S)
    argSet.add(invm2S)
    argSet.add(invm3S)
    argSet.add(jpsipt)
    argSet.add(jpsieta)
    argSet.add(jpsimass)
    argSet.add(gammapt)
    argSet.add(deltaRJpsig)
    argSet.add(deltaRChiJpsi)
    argSet.add(ctpv)
    argSet.add(ctpverr)
    argSet.add(ctpvsig)
    argSet.add(ctbs)
    argSet.add(ctbserr)
    argSet.add(ctbssig)
    argSet.add(Rconv)
    argSet.add(jpsiVprob)
    argSet.add(Y1Smass_nSigma)
    argSet.add(Y2Smass_nSigma)
    argSet.add(Y3Smass_nSigma)

    argSet.add(jpsipx)
    argSet.add(jpsipy)
    argSet.add(jpsipz)
    argSet.add(gammapx)
    argSet.add(gammapy)
    argSet.add(gammapz)
    argSet.add(vertexChi2ProbGamma)


    rds    = RooDataSet("d","d",argSet,"weight")
    
    
    hmass    = TH1F("hmass","hmass",400,9.0,25)
    Ymass    = TH1F("Ymass","Ymass",400,9.0,13)
    hmass_DeltaRChiJpsi    = TH2F("hmass_DeltaRChiJpsi","hmass_DeltaRChiJpsi",10,0,0.5,400,9.0,13)
    hmass_DeltaRJpsig    = TH2F("hmass_DeltaRJpsig","hmass_DeltaRJpsig",80,0,4,400,9.0,13)
    hmass_DeltaPhiJpsig    = TH2F("hmass_DeltaPhiJpsig","hmass_DeltaPhiJpsig",80,0,4,400,9.0,13)
    hmass_DeltaEtaJpsig    = TH2F("hmass_DeltaEtaJpsig","hmass_DeltaEtaJpsig",80,0,4,400,9.0,13)
    hmass_gammapt    = TH2F("hmass_gammaPt","hmass_gammaPt",80,0,4,400,9.0,13)
    hgammapt_DeltaRJpsig    = TH2F("hgammaPt_DeltaRJpsig","hgammaPt_DeltaRJpsig",80,0,4,100,0,10)
    hgammapt = TH1F("gpt","gpt",100,0,10)
    hct      = TH1F("ct","ct",100,-0.5,0.5)
    hmass_Rconv    = TH2F("hmass_Rconv","hmass_Rconv",100,0,50,400,9.0,13)
    
    hGammaP    = TH1F("hGammaP","hGammaP",100,-15,15)
    hUpsP    = TH1F("hUpsP","hUpsP",100,-50,50)
    hCosAlphaP    = TH1F("hCosAlphaP","hCosAlphaP",100,-1,1)
    
    # loop over events
    ncands=0
    nevents=0
    ncoll=0
    nSel=0
    for event in events:
        nevents+=1
        print 'nevents', nevents
        #print dir(event._event.getTFile().GetName())
        #print dir(event._event.time())
        
        #print "run %d \t lumi %d \t orbit %d \t file %s" % (event._event.getRun().run(),event._event.luminosityBlock(), event._event.orbitNumber(),event._event.getTFile().GetName())
        #break
    
        
        #print event._event.getEvent()
        
        #print event
        #if not event.getByLabel (chicCand_l,chicCand_h) : continue
    
    
        event.getByLabel (chicCand_l,chicCand_h)
        event.getByLabel (jpsiCand_l,jpsiCand_h)
        event.getByLabel (jpsiCandPAT_l,jpsiCandPAT_h)
        event.getByLabel (gammaCand_l,gammaCand_h)
        event.getByLabel (convCand_l,convCand_h)
        
        if not chicCand_h.isValid() : continue
        
        chicColl    =  chicCand_h.product()
        jpsiColl    =  jpsiCand_h.product()
        jpsiCollPAT =  jpsiCandPAT_h.product()
        gammaColl   =  gammaCand_h.product()
        convColl    =  convCand_h.product()
        
        #loop over candidates
        for i in range(chicColl.size()):
            ncands+=1
            print 'ncands', ncands
    
            chicCand    = chicColl[i]
            jpsiCand    = jpsiColl[i]
            jpsiCandPAT = jpsiCollPAT[0]
            gammaCand   = gammaColl[i]
            convCand    = convColl[0]
    
            #match gamma candidate and conversion candidate crudely
            dptmin = 9999
            for c in convColl:
                ncoll+=1
                print 'ncoll', ncoll
                diff = sqrt(convCand.pairMomentum().perp2()) - gammaCand.pt()
                if diff < dptmin:
                    dptmin = diff
                    convCand = c
             #print 'dptmin',dptmin    
            
             #print chicCand.mass() , jpsiCand.mass(), gammaCand.pt(), 
            #if abs(chicCand.eta()) > 1.1 : continue
#            if abs(jpsiCand.eta()) > YrapCut : continue
            #reject conversions outside pixles
    
#            if jpsiCandPAT.userFloat('ppdlPV') > YLTCut : continue
            
#            if convCand.conversionVertex().position().rho() > 12 : continue
#            if convCand.conversionVertex().position().rho() < 2 : continue
    
            #if gammaCand.pt() < 1.0 :continue
#            if gammaCand.pt() < 1.0 :continue

    
#            if jpsiCand.mass() > YmassMin or jpsiCand.mass() < YmassMax : continue
            
            deltaR_Jpsig = deltaR(jpsiCand.eta(),jpsiCand.phi(),gammaCand.eta(),gammaCand.phi())
            deltaPhi_Jpsig = deltaPhi(jpsiCand.phi(),gammaCand.phi())
            deltaEta_Jpsig = abs(jpsiCand.eta()-gammaCand.eta())
            deltaR_ChiJpsi = deltaR(jpsiCand.eta(),jpsiCand.phi(),chicCand.eta(),chicCand.phi())
    
            #if deltaR_Jpsig > 1.3 : continue
    
    
            Qval = chicCand.mass()-jpsiCand.mass()
    
            #if mass < 10.43 and mass > 10.42:
                #print "\n--------------\nmass : %f\tgammapt : %f\tjpsiPt :%f \t run %d \t lumi %d \t orbit %d \t file %s"  %(mass, gammaCand.pt(), jpsiCand.pt(),event._event.getRun().run(),event._event.luminosityBlock(), event._event.orbitNumber(),event._event.getTFile().GetName())
    
            #print "\n--------------\nrun %d \t event %015d \t lumi %06d \t orbit %010d \t file %s"  %(event._event.getRun().run(),event._event.id().event(),event._event.luminosityBlock(), event._event.orbitNumber(),event._event.getTFile().GetName())
    
            #print "\n--------------\nrun %d \t event %015d \t mass : %f\tgammapt : %f\tjpsiPt :%f \t file %s"  %(event._event.getRun().run(),event._event.id().event(),mass, gammaCand.pt(), jpsiCand.pt(),event._event.getTFile().GetName())
            nSel+=1
                    
            Ymass.Fill(jpsiCand.mass())
            hmass.Fill(Qval)
            hgammapt.Fill(gammaCand.pt())
            hgammapt_DeltaRJpsig.Fill(deltaR_Jpsig,gammaCand.pt())
            hct.Fill(jpsiCandPAT.userFloat('ppdlPV'))
            hmass_DeltaRJpsig.Fill(deltaR_Jpsig,Qval)
            hmass_DeltaRChiJpsi.Fill(deltaR_ChiJpsi,Qval)
            hmass_DeltaPhiJpsig.Fill(deltaPhi_Jpsig,Qval)
            hmass_DeltaEtaJpsig.Fill(deltaEta_Jpsig,Qval)
            hmass_gammapt.Fill(gammaCand.pt(),Qval)
            hmass_Rconv.Fill(convCand.conversionVertex().position().rho(),Qval)
            
            gammaP=sqrt(gammaCand.px()*gammaCand.px()+gammaCand.py()*gammaCand.py()+gammaCand.pz()*gammaCand.pz())
            UpsP=sqrt(jpsiCand.px()*jpsiCand.px()+jpsiCand.py()*jpsiCand.py()+jpsiCand.pz()*jpsiCand.pz())
            hGammaP.Fill(gammaP)
            hUpsP.Fill(UpsP)
            hCosAlphaP.Fill((gammaCand.px()*jpsiCand.px()+gammaCand.py()*jpsiCand.py()+gammaCand.pz()*jpsiCand.pz())/(gammaP*UpsP))
            
            #fill RooDataSet
            #if mass > 3.2 and mass < 4.0:
    
            sigma=Ymass_a+Ymass_b*(abs(jpsiCand.y())-Ymass_c)
            if abs(jpsiCand.y())<Ymass_c:
                sigma=Ymass_a
            
            
            deltaRChiJpsi.setVal(deltaR_ChiJpsi)
            deltaRJpsig.setVal(deltaR_Jpsig)  
            invm1S.setVal( Qval + Y1Smass0 )
            invm2S.setVal( Qval + Y2Smass0 )
            invm3S.setVal( Qval + Y3Smass0 )
            jpsipt.setVal(jpsiCand.pt())
            jpsimass.setVal(jpsiCand.mass())
            jpsieta.setVal(jpsiCand.y())
            ctpv.setVal(jpsiCandPAT.userFloat('ppdlPV'))
            ctpverr.setVal(jpsiCandPAT.userFloat('ppdlErrPV'))
            ctbs.setVal(jpsiCandPAT.userFloat('ppdlBS'))
            ctbserr.setVal(jpsiCandPAT.userFloat('ppdlErrBS'))
            gammapt.setVal(gammaCand.pt())
            Rconv.setVal(convCand.conversionVertex().position().rho())
            ctpvsig.setVal(abs(jpsiCandPAT.userFloat('ppdlPV'))/jpsiCandPAT.userFloat('ppdlErrPV'))
            ctbssig.setVal(abs(jpsiCandPAT.userFloat('ppdlBS'))/jpsiCandPAT.userFloat('ppdlErrBS'))
            jpsiVprob.setVal(jpsiCandPAT.userFloat('vProb'))
            Y1Smass_nSigma.setVal((jpsiCand.mass()-Y1Smass0)/sigma)
            Y2Smass_nSigma.setVal((jpsiCand.mass()-Y2Smass0)/sigma)
            Y3Smass_nSigma.setVal((jpsiCand.mass()-Y3Smass0)/sigma)
            vertexChi2ProbGamma.setVal(TMath.Prob(convCand.conversionVertex().chi2(),int(convCand.conversionVertex().ndof())))
#            vertexChi2ProbGamma.setVal(1)
            jpsipx.setVal(jpsiCand.px())
            jpsipy.setVal(jpsiCand.py())
            jpsipz.setVal(jpsiCand.pz())
            gammapx.setVal(gammaCand.px())
            gammapy.setVal(gammaCand.py())
            gammapz.setVal(gammaCand.pz())
 
            weight.setVal(1)
 
            rds.add(argSet)
            
            print TMath.Prob(convCand.conversionVertex().chi2(),int(convCand.conversionVertex().ndof())), "from ", int(convCand.conversionVertex().ndof()), " and " ,convCand.conversionVertex().ndof(), " and ", convCand.conversionVertex().chi2()
    
        if ncands%1000==0:
            print ncands
            #break
    
    
    print "Number of chib candidates:          ", ncands 
    print "Number of selected chib candidates: ", nSel
    
    outdataset= TFile("rooDS_"+str(options.cutName)+".root","recreate")
    
    outdataset.cd()
    rds.Write()
    tree = rds.tree()
    tree.SetName("tree")
    tree.Write()
    
    hlist = [hmass,Ymass,hgammapt,hct, hgammapt_DeltaRJpsig , hmass_DeltaRJpsig ,hmass_DeltaPhiJpsig ,hmass_DeltaEtaJpsig , hmass_gammapt, hmass_DeltaRChiJpsi, hmass_Rconv,hCosAlphaP,hGammaP,hUpsP]
    
    for h in hlist : h.Write()

if __name__ == "__main__":
    parser = OptionParser(description='%prog : Chib Invariant Mass Fitter.',usage='chic-unbinnedlh.py --options')
    parser.add_option('--cutName',dest='cutName',default='Def',help='The name of your cut')
    (options,args) = parser.parse_args()

    main(options,args)
