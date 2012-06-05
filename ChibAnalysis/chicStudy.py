from ROOT import TTree, TFile,TH1F,TH2F, gROOT , TRandom2, TMath
from ROOT import RooDataSet,RooRealVar,RooArgSet
import sys
from DataFormats.FWLite import Events, Handle
from math import sqrt, pi as M_PI
from ROOT import TLorentzVector

from ROOT import TTree, TFile,TH1F,TH2F, gROOT, TMath
from ROOT import RooDataSet,RooRealVar,RooArgSet
from optparse import OptionParser

global ncands
ncands=0

invm1S     	= RooRealVar("invm1S",   "invm1S",2.8,65)
invm2S     	= RooRealVar("invm2S",   "invm2S",2.8,65)
invm3S     	= RooRealVar("invm3S",   "invm3S",2.8,65)

Y1Smass_nSigma  = RooRealVar("Y1Smass_nSigma","Y1Smass_nSigma",-100,100)
Y2Smass_nSigma  = RooRealVar("Y2Smass_nSigma","Y2Smass_nSigma",-100,100)
Y3Smass_nSigma  = RooRealVar("Y3Smass_nSigma","Y3Smass_nSigma",-100,100)

jpsimass   	= RooRealVar("jpsimass", "jpsimass",0,12)
jpsiVprob   	= RooRealVar("jpsiVprob","jpsiVprob",0,1.5)

jpsipt     	= RooRealVar("jpsipt", "jpsipt",0,200)
gammapt    	= RooRealVar("gammapt", "gammapt",0,100)

deltaRChiJpsi   = RooRealVar("deltaRChiJpsi","deltaRChiJpsi",0,10)
deltaRJpsig     = RooRealVar("deltaRJpsig","deltaRJpsig",0,10)

jpsieta     	= RooRealVar("jpsieta","jpsieta",-5,5)
gammaeta    	= RooRealVar("gammaeta", "gammaeta",-5,5)

ctpv        	= RooRealVar("ctpv",   "ctpv",-5,5)
ctpverr     	= RooRealVar("ctpverr","ctpverr",0,10)
ctbs        	= RooRealVar("ctbs",   "ctbs",-5,5)
ctbserr     	= RooRealVar("ctbserr","ctbserr",0,10)
ctpvsig     	= RooRealVar("ctpvsig",   "ctpvsig",0,1000)
ctbssig     	= RooRealVar("ctbssig",   "ctbssig",0,1000)
weight      	= RooRealVar("weight", "weight",0,1)

Rconv       	= RooRealVar("Rconv","Rconv",0,100)

jpsipx   	= RooRealVar("jpsipx", "jpsipx",-100,100)
jpsipy   	= RooRealVar("jpsipy", "jpsipy",-100,100)
jpsipz   	= RooRealVar("jpsipz", "jpsipz",-100,100)
gammapx   	= RooRealVar("gammapx", "gammapx",-100,100)
gammapy   	= RooRealVar("gammapy", "gammapy",-100,100)
gammapz   	= RooRealVar("gammapz", "gammapz",-100,100)

Q    		= RooRealVar("Q", "Q",0,100)
ndofConv   	= RooRealVar("ndofConv", "ndofConv",0,100)
chi2Conv   	= RooRealVar("chi2Conv", "chi2Conv",0,100)
vertexChi2ProbGamma   = RooRealVar("vertexChi2ProbGamma", "vertexChi2ProbGamma",0,1)

UpsIndex   	= RooRealVar("UpsIndex", "UpsIndex",0,10000000)
GammaIndex   	= RooRealVar("GammaIndex", "GammaIndex",0,500)


Pi0Mass     	= RooRealVar("Pi0Mass", "Pi0Mass",0,5)
alphaPi0     	= RooRealVar("alphaPi0", "alphaPi0",0,1)

vtxdz        	= RooRealVar("vtxdz"      , "vtxdz"      ,-15,15)
vtxerrdz        = RooRealVar("vtxerrdz"   , "vtxerrdz"   ,  0, 15)
vtxNsigmadz     = RooRealVar("vtxNsigmadz", "vtxNsigmadz",-15,15)

RunNb     = RooRealVar("RunNb", "RunNb",0,1e9)
EventNb     = RooRealVar("EventNb", "EventNb",0,1e15)

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

def sqr(a):
    return a*a

def checkDistance(p, maxDist, vtx):
    dist=sqr(p.vertex().z()-vtx.position().z())+sqr(p.vertex().y()-vtx.position().y())+sqr(p.vertex().x()-vtx.position().x())
    if dist<maxDist['d']:
        maxDist['d']=dist
        return True
    return False

def AssignVtx(aJpsiCand, aPhotonCand, privtxs):
    if privtxs.size()==1 :
        #print "only 1 privtx"
        return (True, privtxs[0], privtxs[0])

    #print "privtx size" , privtxs.size()
    
    # const reco::Vertex *JpsiVtx=0, *PhotonVtx=0;
    dist_JpsiVtx={'d':100000.0}
    dist_PhotonVtx={'d':100000.0}
    # Boost_foreach(const reco::Vertex& vtx, *privtxs.product()){
    idx=0
    for vtx in privtxs:
        if checkDistance(aJpsiCand,dist_JpsiVtx,vtx) : iJpsiVtx=idx
        if checkDistance(aPhotonCand,dist_PhotonVtx,vtx) : iPhotonVtx=idx
        idx+=1

        
    #print "distances " , dist_JpsiVtx['d'] , dist_PhotonVtx['d']    
    JpsiVtx=privtxs[iJpsiVtx]
    PhotonVtx=privtxs[iPhotonVtx]
    #print "the Jpsi has vertex rho, phi, z " , JpsiVtx.position().rho(), JpsiVtx.position().phi(), JpsiVtx.position().z()
    #print "the Phot has vertex rho, phi, z " , PhotonVtx.position().rho(), PhotonVtx.position().phi(), PhotonVtx.position().z()

    if iJpsiVtx==iPhotonVtx :
        #print "vertex indexes " , iPhotonVtx    
        return (True, JpsiVtx, PhotonVtx)
    else:
        return (False, JpsiVtx, PhotonVtx)
    
def checkVtxCompatibility(aJpsiCand, aPhotonCand, privtxs):

    (result, JpsiVtx, PhotonVtx) = AssignVtx(aJpsiCand, aPhotonCand, privtxs)
    
    if result==True:
        return (result, JpsiVtx, PhotonVtx) 
    
    #print "covariance " , JpsiVtx.covariance(0,0) , JpsiVtx.covariance(2,2)  
    
    dist2=(JpsiVtx.position()-PhotonVtx.position()).Mag2();
    sigma2=0;
    DX=[(JpsiVtx.position().x()-PhotonVtx.position().x()),
        (JpsiVtx.position().y()-PhotonVtx.position().y()),
        (JpsiVtx.position().z()-PhotonVtx.position().z())
        ]
    
    i=0
    while i<3:
        sigma2+=DX[i]*DX[i]*(JpsiVtx.covariance(i,i)+PhotonVtx.covariance(i,i))
        i+=1

       
        #   print "sigma2 " << sigma2 << " " << DX[i] << " " <<  (JpsiVtx->covariance(i,i)+PhotonVtx->covariance(i,i)) << std::endl;
    
    sigma2/= dist2;
    
    print  "distance among vertexes %f sigma %f" %( dist2 , sigma2 )
    
    if dist2/sigma2>25 :
        print "rejected pair "
        return (False, JpsiVtx, PhotonVtx)
    return (True, JpsiVtx, PhotonVtx)


def f_dz(conv,  refVtx) :
    vtx = conv.conversionVertex();
    
    mom = conv.refittedPairMomentum();

    dz = (vtx.z()-refVtx.z()) - ((vtx.x()-refVtx.x())*mom.x()+(vtx.y()-refVtx.y())*mom.y())/mom.rho() * mom.z()/mom.rho();
    return dz;


def f_dzError(conv, refVtx):
    
    vtx = conv.conversionVertex();
    
    mom = conv.refittedPairMomentum();
    
    sigmadx2 = vtx.covariance(0,0)+refVtx.covariance(0,0);
    sigmady2 = vtx.covariance(1,1)+refVtx.covariance(1,1);
    sigmadz2 = vtx.covariance(2,2)+refVtx.covariance(2,2);
    
    pt  = mom.rho();
    px  = mom.x();
    py  = mom.y();
    pz  = mom.z();

    sigmadz2_B = sqr(px/pt*pz/pt)*sigmadx2 + sqr(py/pt*pz/pt)*sigmady2;

    return sqrt(sigmadz2 + sigmadz2_B);


def checkNearestVertex(conv, vtx):

    _dz= f_dz(conv,vtx)
    _dzError=f_dzError(conv,vtx)

    # this is just a check
    _sigmasTkVtxComp=3
    #if(abs(_dz/_dzError) > _sigmasTkVtxComp):
        #print '---->'
        #print "vtx z" , vtx.z() ," \t dz " ,_dz ," \t " ,_dzError ," \t dz/dzErr " ,_dz/_dzError

    return(_dz, _dzError, _dz/_dzError);


def findNearestVertex(conv, vtxColl):
    '''Find the vertex nearest to the conversion'''


    minSigmas=100000
    _sigmasTkVtxComp=100
   
    for vtx in vtxColl:
        _dz= f_dz(conv,vtx)
        _dzError=f_dzError(conv,vtx)
        if abs(_dz/_dzError) > _sigmasTkVtxComp and abs(_dz/_dzError)< minSigmas:
            minSigmas = abs(_dz/_dzError)
            theVtx = vtx

    if (minSigmas < _sigmasTkVtxComp):
        return (True, theVtx)
    else:
        return (False, 0)    
    
def bookHistos():


    global argSet
    argSet = RooArgSet()
    argSet.add(invm1S)
    argSet.add(invm2S)
    argSet.add(invm3S)
    argSet.add(jpsipt)
    argSet.add(jpsieta)
    argSet.add(jpsimass)
    argSet.add(gammapt)
    argSet.add(gammaeta)
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
    argSet.add(Q)
    argSet.add(ndofConv)
    argSet.add(chi2Conv)

    argSet.add(Pi0Mass)
    argSet.add(alphaPi0)

    argSet.add(vtxdz)
    argSet.add(vtxerrdz)
    argSet.add(vtxNsigmadz)

    argSet.add(UpsIndex)
    argSet.add(GammaIndex)

    argSet.add(RunNb)
    argSet.add(EventNb)

    global rds
    rds    = RooDataSet("d","d",argSet,"weight")



    global hlist
    hlist = []

def saveHistos():
    outdataset= TFile("rooDS_"+str(options.fileName)+".root","recreate")
    outdataset.cd()
    rds.Write()
    tree = rds.tree()
    tree.SetName('atree')
    tree.Write()
    for h in hlist : h.Write()

def fillHistos(chicCand, jpsiCand, jpsiCandPAT, gammaCand, convCand, z_tupla, nearestMass, i, nevents, event, alpha):

    deltaR_Jpsig = deltaR(jpsiCand.eta(),jpsiCand.phi(),gammaCand.eta(),gammaCand.phi())
    deltaPhi_Jpsig = deltaPhi(jpsiCand.phi(),gammaCand.phi())
    deltaEta_Jpsig = abs(jpsiCand.eta()-gammaCand.eta())
    deltaR_ChiJpsi = deltaR(jpsiCand.eta(),jpsiCand.phi(),chicCand.eta(),chicCand.phi())
    
    
    
    Qval = chicCand.mass()-jpsiCand.mass()
    
    Y1Smass0=3.096916
    Y2Smass0= 3.68609 
    Y3Smass0=10.3552
    Ymass_a=0.058
    Ymass_b=0.047
    Ymass_c=0.22
    
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
    gammaeta.setVal(gammaCand.y())
    Rconv.setVal(convCand.conversionVertex().position().rho())
    ctpvsig.setVal(abs(jpsiCandPAT.userFloat('ppdlPV'))/jpsiCandPAT.userFloat('ppdlErrPV'))
    ctbssig.setVal(abs(jpsiCandPAT.userFloat('ppdlBS'))/jpsiCandPAT.userFloat('ppdlErrBS'))
    jpsiVprob.setVal(jpsiCandPAT.userFloat('vProb'))
    Y1Smass_nSigma.setVal((jpsiCand.mass()-Y1Smass0)/sigma)
    Y2Smass_nSigma.setVal((jpsiCand.mass()-Y2Smass0)/sigma)
    Y3Smass_nSigma.setVal((jpsiCand.mass()-Y3Smass0)/sigma)
    vertexChi2ProbGamma.setVal(TMath.Prob(convCand.conversionVertex().chi2(),int(convCand.conversionVertex().ndof())))
    jpsipx.setVal(jpsiCand.px())
    jpsipy.setVal(jpsiCand.py())
    jpsipz.setVal(jpsiCand.pz())
    gammapx.setVal(gammaCand.px())
    gammapy.setVal(gammaCand.py())
    gammapz.setVal(gammaCand.pz())
    Q.setVal(Qval)
    ndofConv.setVal(convCand.conversionVertex().ndof())
    chi2Conv.setVal(convCand.conversionVertex().chi2())
    
    weight.setVal(1)

    Pi0Mass.setVal(nearestMass)

    vtxdz.setVal(z_tupla[0])
    vtxerrdz.setVal(z_tupla[1])
    vtxNsigmadz.setVal(z_tupla[2])

    GammaIndex.setVal(i)
    UpsIndex.setVal(nevents)
    
    RunNb.setVal(event._event.id().run())
    EventNb.setVal(event._event.id().event())
    
    alphaPi0.setVal(alpha)
    
    rds.add(argSet)


def evalDiPhoton(conv, pf):

    pi0Cand = conv.p4() + pf.p4()
    return pi0Cand

def rejectPF(pf):
    ''' apply same cuts as in analyzer '''

    etaMinGammaForPi0Rejection = 0
    etaMaxGammaForPi0Rejection = 3
    
    if pf.particleId() == 4 and abs(pf.eta()) >= etaMinGammaForPi0Rejection and abs(pf.eta()) < etaMaxGammaForPi0Rejection:
        return False
    else:
        return True


def vertexCompatibility(jpsiCand,gammaCand,convCand,vtxColl):
    # Assign to the Jpsi and to the Photon (compCandidate) the vertex from the vtxCollection
    # this is needed to access the covariance matrix that is stored in the vertex object
    # but not in the vertex method of the compCandidate
    (result, JpsiVtx, PhotonVtx) = AssignVtx(jpsiCand,gammaCand,vtxColl)

    # Find the nearest vertex among the conversion object and the vtxCollection
    # This is a simple check that the nearest vertex is the one already found in the analyzer 
    # (resultFind, convVtx) =  findNearestVertex(convCand,vtxColl)
    # if resultFind==True and convVtx.z()!=PhotonVtx.z() :
    #     print "problem " , convVtx.z() , PhotonVtx.z() 
    #     raise Exception

    # Find the nearest vertex among the conversion object and the Jpsi
    # Finding this again allows to cut stronger than before (that was at 5 sigmas)
    # NB: the approach is different respect to the one in the analyzer. There we
    #     cut on the distance among the conversion and the nearest primary vertex
    #     then we check also that this prvtx is compatible with the Jpsi
    #     Here we do the two things at the same time
    z_tupla=checkNearestVertex(convCand,JpsiVtx)    
    return z_tupla       

def analyzeEvent(args):

    loopOnChi(args)

def loopOnChi(args):
    
    for i in range(args['chicColl'].size()):
        global ncands
        ncands+=1
                
        chicCand    = args['chicColl'][i]
        jpsiCand    = args['jpsiColl'][i]
        jpsiCandPAT = args['jpsiCollPAT'][0]
        gammaCand   = args['gammaColl'][i]

        #match gamma candidate and conversion candidate crudely
        dptmin = 9999
        for c in args['convColl']:
            diff = abs(c.refittedPairMomentum().rho() - gammaCand.pt())
            if diff < dptmin:
                dptmin = diff
                convCand = c


        z_tupla = vertexCompatibility(jpsiCand,gammaCand,convCand,args['vtxColl'])
        
        pi0_tupla = pi0Cut(gammaCand,args['pfColl'])
        nearestMass = pi0_tupla[0]
        alpha = pi0_tupla[1]

        fillHistos(chicCand, jpsiCand, jpsiCandPAT, gammaCand, convCand, z_tupla, nearestMass, i, args['nevents'], args['event'] , alpha)
        
def pi0Cut(conv,pfColl):
    
    NearestMass = 0
    deltaMass = 10000
    Pi0Mass = 0.134

    for pf in pfColl:
        if rejectPF(pf):
            continue

        pi0Cand = evalDiPhoton(conv, pf)
        if abs(pi0Cand.M()-Pi0Mass) < deltaMass:
            NearestMass = pi0Cand.M()
            deltaMass = abs(pi0Cand.M()-Pi0Mass)  
	        
            E1 = conv.energy()
            E2 = pf.energy() 
            alpha = abs( E1-E2 ) / ( E1+E2 )         

    return (NearestMass,alpha);


def loopOnEvents():

    events= Events(inputfiles)

    pfCand_h = Handle ("vector<reco::PFCandidate>")
    pfCand_l = ( "particleFlow","" )

    # candidate Chi
    chicCand_h  = Handle ("vector<reco::CompositeCandidate>")
    chicCand_l = ( "chic","chicCompCand","chics5" )
    
    #candidate J/Psi
    jpsiCand_h  = Handle ("vector<pat::CompositeCandidate>") #sometimes reco:: oder pat::
    jpsiCand_l = ( "chic","jpsiCompCand" ,"chics5")
    
    jpsiCandPAT_h  = Handle ("vector<pat::CompositeCandidate>")
    jpsiCandPAT_l = ( "onia2MuMuPatTrkTrk","" )
    
    #candidate gamma
    gammaCand_h  = Handle ("vector<reco::CompositeCandidate>")
    gammaCand_l = ( "chic","gammaCompCand" ,"chics5")
    
    convCand_h = Handle ("vector<reco::Conversion>")
    convCand_l = ( "allConversions","" )
    
    vtxCand_h = Handle('vector<reco::Vertex>')
    vtxCand_l = ('offlinePrimaryVertices')

    events= Events(inputfiles)

    # loop over events
    nevents = 0
    for event in events:
        nevents += 1

        if nevents % 100 == 0 :
            print 'Events Processed ' , nevents

        #if nevents > 1000:
            #break
    
        #print "run %d \t lumi %d \t orbit %d \t file %s" % (event._event.getRun().run(),event._event.luminosityBlock(), event._event.orbitNumber(),event._event.getTFile().GetName())
    

        event.getByLabel (pfCand_l,pfCand_h)
        event.getByLabel (chicCand_l,chicCand_h)
        event.getByLabel (jpsiCand_l,jpsiCand_h)
        event.getByLabel (jpsiCandPAT_l,jpsiCandPAT_h)
        event.getByLabel (gammaCand_l,gammaCand_h)
        event.getByLabel (convCand_l,convCand_h)
        event.getByLabel (vtxCand_l,vtxCand_h)

        if not chicCand_h.isValid() : continue
        if vtxCand_h.isValid() != True :
            print "no vertex in the event "
            continue

        args={}   
        args['chicColl']    =  chicCand_h.product()
        args['jpsiColl']    =  jpsiCand_h.product()
        args['jpsiCollPAT'] =  jpsiCandPAT_h.product()
        args['gammaColl']   =  gammaCand_h.product()
        args['convColl']    =  convCand_h.product()
        args['vtxColl']     =  vtxCand_h.product()
        args['pfColl']      =  pfCand_h.product()
        args['nevents']      =  nevents
        args['event']      =  event

        analyzeEvent(args)

    print "Number of events:          ", nevents 
    print "Number of chic candidates:          ", ncands 


if __name__  == '__main__':
    
    parser = OptionParser(description='%prog : Chic Invariant Mass Fitter.',usage='chic-unbinnedlh.py --options')
    parser.add_option('--fileName',dest='fileName',default='Def',help='The name of your output file')
    (options,args) = parser.parse_args()

    # Create histograms, etc.
    gROOT.SetBatch()        # don't pop up canvases
    gROOT.SetStyle('Plain')

#    inputfiles = [ 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_05Apri112_A001.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_05Apri112_A002.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_05Apri112_A003.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_05Apri112_A004.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_05Apri112_A005.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_05Apri112_A006.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_05Apri112_A007.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_05Apri112_A008.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_05Apri112_A009.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_05Apri112_A010.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_05Apri112_A011.root' ]
#    inputfiles = [ 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_05Apri112_B.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_05Apri112_B001.root',  'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_05Apri112_B002.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_05Apri112_B003.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_05Apri112_B004.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_05Apri112_B005.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_05Apri112_B006.root']

#    inputfiles = [ 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_2011A001.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_2011A002.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_2011A003.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_2011A004.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_2011A005.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_2011A006.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_2011A007.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_2011A008.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_2011A009.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_2011A010.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_2011A011.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_2011A012.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_2011A013.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_2011A014.root' ]
#    inputfiles = [ 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_2011B001.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_2011B002.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_2011B003.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_2011B004.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_2011B005.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_2011B006.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_2011B007.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_2011B008.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_2011B009.root' ]
#    inputfiles = [ '/tmp/knuenz/runA/chic_SelEvts_100_1_9mq.root']

#    inputfiles = [ 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_100_1_9mq.root' ]

#    inputfiles = [ 'file:/tmp/mahmad/Psi2S_2011B/chicStep5_Psi2S_02May_Part1000.root',
#    'file:/tmp/mahmad/Psi2S_2011B/chicStep5_Psi2S_02May_Part1001.root',
#    'file:/tmp/mahmad/Psi2S_2011B/chicStep5_Psi2S_02May_Part1002.root',
#    'file:/tmp/mahmad/Psi2S_2011B/chicStep5_Psi2S_02May_Part1003.root',
#    'file:/tmp/mahmad/Psi2S_2011B/chicStep5_Psi2S_02May_Part1004.root',
#    'file:/tmp/mahmad/Psi2S_2011B/chicStep5_Psi2S_02May_Part1005.root',
#    'file:/tmp/mahmad/Psi2S_2011B/chicStep5_Psi2S_02May_Part1006.root',
#    'file:/tmp/mahmad/Psi2S_2011B/chicStep5_Psi2S_02May_Part1007.root',
#    'file:/tmp/mahmad/Psi2S_2011B/chicStep5_Psi2S_02May_Part1008.root',
#    'file:/tmp/mahmad/Psi2S_2011B/chicStep5_Psi2S_02May_Part1009.root',
#    'file:/tmp/mahmad/Psi2S_2011B/chicStep5_Psi2S_02May_Part1010.root',
#    'file:/tmp/mahmad/Psi2S_2011B/chicStep5_Psi2S_02May_Part2000.root',
#    'file:/tmp/mahmad/Psi2S_2011B/chicStep5_Psi2S_02May_Part2001.root',
#    'file:/tmp/mahmad/Psi2S_2011B/chicStep5_Psi2S_02May_Part2002.root',
#    'file:/tmp/mahmad/Psi2S_2011B/chicStep5_Psi2S_02May_Part2003.root',
#    'file:/tmp/mahmad/Psi2S_2011B/chicStep5_Psi2S_02May_Part2004.root',
#    'file:/tmp/mahmad/Psi2S_2011B/chicStep5_Psi2S_02May_Part2005.root',
#    'file:/tmp/mahmad/Psi2S_2011B/chicStep5_Psi2S_02May_Part2006.root',
#    'file:/tmp/mahmad/Psi2S_2011B/chicStep5_Psi2S_02May_Part2007.root',
#    'file:/tmp/mahmad/Psi2S_2011B/chicStep5_Psi2S_02May_Part2008.root',
#    'file:/tmp/mahmad/Psi2S_2011B/chicStep5_Psi2S_02May_Part2009.root',
#    'file:/tmp/mahmad/Psi2S_2011B/chicStep5_Psi2S_02May_Part2010.root',
#    'file:/tmp/mahmad/Psi2S_2011B/chicStep5_Psi2S_02May_Part2011.root']
 
    inputfiles = [ 'file:/tmp/mahmad/Psi2S_2011A/chicStep5_Psi2S_2011A_Part_1000.root',
    'file:/tmp/mahmad/Psi2S_2011A/chicStep5_Psi2S_2011A_Part_1001.root',
    'file:/tmp/mahmad/Psi2S_2011A/chicStep5_Psi2S_2011A_Part_1002.root',
    'file:/tmp/mahmad/Psi2S_2011A/chicStep5_Psi2S_2011A_Part_1003.root',
    'file:/tmp/mahmad/Psi2S_2011A/chicStep5_Psi2S_2011A_Part_1004.root',
    'file:/tmp/mahmad/Psi2S_2011A/chicStep5_Psi2S_2011A_Part_1005.root',
    'file:/tmp/mahmad/Psi2S_2011A/chicStep5_Psi2S_2011A_Part_1006.root',
    'file:/tmp/mahmad/Psi2S_2011A/chicStep5_Psi2S_2011A_Part_1007.root',
    'file:/tmp/mahmad/Psi2S_2011A/chicStep5_Psi2S_2011A_Part_1008.root',
    'file:/tmp/mahmad/Psi2S_2011A/chicStep5_Psi2S_2011A_Part_2000.root',
    'file:/tmp/mahmad/Psi2S_2011A/chicStep5_Psi2S_2011A_Part_2001.root',
    'file:/tmp/mahmad/Psi2S_2011A/chicStep5_Psi2S_2011A_Part_2002.root',
    'file:/tmp/mahmad/Psi2S_2011A/chicStep5_Psi2S_2011A_Part_2003.root',
    'file:/tmp/mahmad/Psi2S_2011A/chicStep5_Psi2S_2011A_Part_2004.root',
    'file:/tmp/mahmad/Psi2S_2011A/chicStep5_Psi2S_2011A_Part_2005.root',
    'file:/tmp/mahmad/Psi2S_2011A/chicStep5_Psi2S_2011A_Part_2006.root',
    'file:/tmp/mahmad/Psi2S_2011A/chicStep5_Psi2S_2011A_Part_2007.root',
    'file:/tmp/mahmad/Psi2S_2011A/chicStep5_Psi2S_2011A_Part_2008.root',
    'file:/tmp/mahmad/Psi2S_2011A/chicStep5_Psi2S_2011A_Part_2009.root',
    'file:/tmp/mahmad/Psi2S_2011A/chicStep5_Psi2S_2011A_Part_3000.root',
    'file:/tmp/mahmad/Psi2S_2011A/chicStep5_Psi2S_2011A_Part_3001.root',
    'file:/tmp/mahmad/Psi2S_2011A/chicStep5_Psi2S_2011A_Part_3002.root',
    'file:/tmp/mahmad/Psi2S_2011A/chicStep5_Psi2S_2011A_Part_3003.root',
    'file:/tmp/mahmad/Psi2S_2011A/chicStep5_Psi2S_2011A_Part_3004.root',
    'file:/tmp/mahmad/Psi2S_2011A/chicStep5_Psi2S_2011A_Part_3005.root',
    'file:/tmp/mahmad/Psi2S_2011A/chicStep5_Psi2S_2011A_Part_3006.root',
    'file:/tmp/mahmad/Psi2S_2011A/chicStep5_Psi2S_2011A_v4_2000.root',
    'file:/tmp/mahmad/Psi2S_2011A/chicStep5_Psi2S_2011A_v4_2001.root',
    'file:/tmp/mahmad/Psi2S_2011A/chicStep5_Psi2S_2011A_v4_2002.root',
    'file:/tmp/mahmad/Psi2S_2011A/chicStep5_Psi2S_2011A_v4_2003.root',
    'file:/tmp/mahmad/Psi2S_2011A/chicStep5_Psi2S_2011A_v4_2004.root']
 
 
    #inputfiles = [ 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_100_1_9mq.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_101_1_fcB.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_102_1_xBf.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_103_1_R6A.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_103_2_PpM.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_104_1_81b.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_105_1_S14.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_106_1_hqi.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_107_1_nyx.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_108_1_oVs.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_109_1_dW3.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_10_1_0nW.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_110_1_ByH.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_111_1_bkY.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_112_1_dti.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_113_1_Cfu.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_113_2_F6e.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_114_1_guC.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_115_1_z4n.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_116_1_jm4.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_117_1_cQv.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_118_1_w5Q.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_119_1_C3M.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_11_1_Knm.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_120_1_I0N.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_121_2_O9k.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_122_1_qjO.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_123_1_FqB.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_124_1_AZt.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_125_1_V6F.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_126_1_PQr.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_127_1_5k9.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_128_2_a8P.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_129_1_kMN.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_129_2_qKs.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_12_1_E5N.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_130_1_k77.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_131_1_Ngc.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_132_1_Uqt.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_133_5_Cqh.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_134_1_6XG.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_135_1_f0x.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_136_1_0wq.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_136_2_NFt.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_136_2_mv0.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_137_1_Cgc.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_138_1_wik.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_139_1_Ya6.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_13_1_BgA.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_140_1_Uy3.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_141_1_5qH.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_142_1_cPA.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_143_1_t8l.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_144_1_qfW.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_144_3_eic.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_145_1_1R0.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_146_1_omu.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_147_1_r8s.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_148_1_uxk.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_148_3_vcO.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_149_1_M5y.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_14_1_kJe.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_150_1_rmf.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_151_2_7yi.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_152_2_QN2.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_153_1_KCu.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_154_1_195.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_155_3_hF8.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_156_1_D15.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_157_1_4vI.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_158_1_dWK.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_159_1_SWp.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_15_1_vEE.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_160_1_jeo.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_161_1_Qa5.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_161_1_Vib.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_162_1_W94.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_163_4_FYp.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_164_1_ncG.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_164_2_Ilg.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_165_2_9cF.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_166_1_57j.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_167_2_f6d.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_168_1_zZ3.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_169_1_0h4.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_16_1_zKE.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_170_1_Kn4.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_170_2_zgm.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_171_1_Bfi.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_172_1_IS5.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_173_1_dcR.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_174_2_jlm.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_175_3_KI7.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_176_1_ul7.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_177_1_ZX9.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_178_2_0tU.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_179_2_lMD.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_17_1_IBZ.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_180_1_S0A.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_181_1_yIx.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_182_1_dag.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_183_1_qFA.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_184_1_uJy.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_185_1_uWs.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_186_1_sHk.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_187_1_uLA.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_188_1_1jO.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_189_1_6aM.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_18_1_SGl.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_190_1_l1B.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_191_1_KC6.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_192_2_0Kf.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_193_2_Cke.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_194_1_yKJ.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_195_1_nC0.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_196_1_Oey.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_197_1_OQs.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_198_1_pp8.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_199_1_g1i.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_19_1_7m2.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_1_1_6e6.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_200_1_Rxq.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_201_1_R1z.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_202_1_UJd.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_203_3_dTn.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_204_1_hL3.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_205_1_1IB.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_206_1_Yww.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_207_1_3vR.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_208_1_xqk.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_209_1_Ijn.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_20_1_S0N.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_210_1_tAg.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_211_1_PJi.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_212_1_Fkk.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_213_1_EXx.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_214_1_PcO.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_215_1_qTH.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_216_2_EGX.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_217_1_Ff8.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_218_1_KZn.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_219_1_DrA.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_219_2_3VN.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_21_1_NYg.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_220_3_KA3.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_221_1_hel.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_222_3_c11.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_223_1_Uqc.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_224_1_Pqp.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_225_3_fXC.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_226_1_3fY.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_227_1_3hO.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_228_1_3mh.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_229_2_thx.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_22_1_4z1.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_230_1_vLu.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_231_2_3kT.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_232_1_3x0.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_233_1_o4t.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_234_1_iCS.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_235_1_T9O.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_236_2_R5x.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_237_1_7Dm.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_238_1_Rn3.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_239_1_hf1.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_23_1_6YY.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_240_1_Opo.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_241_2_pvk.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_242_1_yKm.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_243_1_hhl.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_244_1_paq.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_245_1_BVN.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_246_1_YuD.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_247_1_AiB.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_248_1_Ieh.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_249_2_RWu.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_24_1_jCR.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_250_1_RR8.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_251_1_Rsa.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_252_1_9xc.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_253_1_Qe4.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_254_1_kjr.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_255_1_9ya.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_256_1_7st.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_257_1_AWh.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_25_1_lmS.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_26_2_ptA.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_27_1_nwl.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_28_1_fDN.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_29_1_2Fw.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_2_1_TFp.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_30_1_8D8.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_31_1_j6D.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_31_2_Q4U.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_32_1_wZL.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_33_1_qD3.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_34_1_lNg.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_35_1_FZb.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_36_1_hDV.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_37_1_iDY.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_38_1_8ZO.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_39_1_0Q2.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_3_1_jNQ.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_40_1_V8c.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_41_1_06e.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_42_1_fjt.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_43_2_jC2.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_44_1_SUw.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_44_2_2TN.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_45_1_ZYV.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_46_1_U3H.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_47_1_aOT.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_48_1_xT3.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_49_1_L60.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_4_1_j8D.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_50_1_ULO.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_51_1_Rxc.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_51_2_PXG.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_52_2_cUL.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_53_1_GPj.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_53_2_YgB.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_54_1_dLe.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_55_1_Huh.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_56_3_bJl.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_57_1_9yn.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_58_1_sy8.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_59_1_Crn.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_59_2_HhH.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_5_1_ZGI.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_60_1_7BB.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_61_1_KWl.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_62_1_tiE.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_63_1_zmz.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_64_1_ilo.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_65_2_2GZ.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_66_1_mMm.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_66_3_04K.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_67_2_BPv.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_68_1_JSw.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_68_2_KPm.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_69_1_8Dd.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_6_1_VG7.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_70_1_me6.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_71_2_uNr.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_72_1_RpX.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_73_1_Ood.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_74_1_AYh.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_75_1_BRK.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_76_1_Pbt.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_77_1_b5H.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_78_1_pRC.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_79_1_1MZ.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_7_1_Qqi.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_80_1_EL8.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_81_1_JjJ.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_82_1_rrJ.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_83_1_89M.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_84_2_Juu.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_84_3_iyi.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_85_1_Nek.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_86_1_3OA.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_87_1_XJ7.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_88_1_qbZ.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_89_1_LLn.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_8_1_fNL.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_90_2_LAM.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_91_1_SMi.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_91_2_hY3.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_92_1_r9p.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_93_1_wxe.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_94_1_ng5.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_95_1_Nby.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_96_1_aQc.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_97_2_iLp.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_97_2_roz.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_98_1_VNM.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_99_2_Iug.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011A/chic_SelEvts_9_1_VxR.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_100_1_qqX.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_101_1_F2B.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_102_1_qL6.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_103_2_O0O.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_104_2_PiW.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_105_1_n5d.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_106_1_2jv.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_107_1_mEd.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_108_1_4Vq.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_109_2_V7L.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_10_1_mZT.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_110_1_Hlh.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_111_1_RU1.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_112_2_efP.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_113_1_c9z.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_114_2_Re1.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_115_1_SQX.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_116_1_toQ.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_117_1_BCB.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_118_1_Ncl.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_119_1_LRB.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_11_1_5gu.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_120_1_rUd.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_121_1_DyM.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_122_1_gOP.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_123_1_mnd.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_124_1_hgE.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_125_1_h2s.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_126_1_f4m.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_127_1_uLd.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_128_1_L0K.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_129_2_H7H.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_12_1_jJP.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_130_1_wHk.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_131_2_cL4.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_132_1_Nb1.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_133_1_mWg.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_134_1_oeZ.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_135_1_ZPB.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_136_1_Hoj.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_137_1_euU.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_138_1_dAx.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_139_1_s7G.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_13_1_CEb.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_140_1_LEw.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_141_1_zgd.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_142_1_dgX.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_143_1_6Mb.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_144_1_TWd.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_145_1_V89.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_146_1_g0r.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_147_1_70Z.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_148_2_7LY.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_149_1_0g0.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_149_2_3V5.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_14_1_KJO.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_150_1_ZmY.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_151_1_Qv8.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_152_1_M5W.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_153_1_KE4.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_154_1_yY1.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_155_1_rpT.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_156_1_4Zt.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_157_1_6nM.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_158_1_Ean.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_159_2_5IB.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_15_1_kQE.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_160_1_sP8.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_161_1_IVp.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_162_1_qtp.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_164_1_stw.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_165_1_P27.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_166_1_pt3.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_167_1_fCo.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_168_1_sjG.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_169_1_HhR.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_16_1_19c.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_170_1_0pM.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_171_1_NTj.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_172_1_QkH.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_173_1_J9f.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_174_1_yX7.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_175_2_wmt.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_176_1_gTj.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_177_1_qtH.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_178_1_UVX.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_178_2_Hks.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_179_1_c6Z.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_17_1_927.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_180_1_Hnb.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_181_1_7wn.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_182_1_M9F.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_183_1_Akq.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_184_1_fp9.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_185_1_exM.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_186_1_1tS.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_187_1_rCP.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_188_1_qeF.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_189_1_fP8.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_18_1_k9G.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_190_2_Io8.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_191_3_3qL.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_192_1_gM9.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_193_1_tep.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_194_2_ffN.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_195_1_JUl.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_196_1_LHS.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_197_1_cZt.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_198_1_bgF.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_199_1_TBh.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_19_1_Gcy.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_1_1_yUQ.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_200_1_3Ty.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_201_1_Lni.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_202_1_sQu.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_203_1_JtR.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_20_1_9Dy.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_21_1_37Q.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_22_1_6YR.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_23_1_kvK.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_24_1_ESZ.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_25_1_e2F.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_26_1_Spq.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_27_1_wFE.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_28_1_1yC.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_29_1_NN4.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_2_1_CIh.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_30_1_8Co.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_31_2_R7j.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_32_1_pmJ.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_33_1_PPn.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_34_2_uVP.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_35_1_hVI.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_36_2_EtU.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_37_1_pks.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_38_1_lGF.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_39_1_4Qb.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_3_1_6sh.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_40_1_Nun.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_41_1_0DC.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_42_1_ycD.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_43_1_b8B.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_44_1_5bq.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_45_1_TGu.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_46_1_C75.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_47_1_Nwz.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_48_1_Qgz.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_49_1_cOZ.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_4_1_0D5.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_50_1_ihz.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_51_1_A5v.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_52_2_mlA.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_53_1_GaG.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_53_2_cxS.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_54_1_scm.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_55_1_t6R.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_56_1_H95.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_57_1_uEt.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_57_2_Z6Y.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_58_1_5Kg.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_59_3_IRA.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_5_1_p6j.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_60_1_hnR.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_61_2_Xp9.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_61_3_HoC.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_62_1_7OM.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_63_1_18X.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_64_1_smJ.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_65_1_EFB.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_66_1_MLU.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_67_1_4K1.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_68_1_95o.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_69_1_l50.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_6_1_7GL.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_70_1_8k4.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_71_1_T3I.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_72_1_9zm.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_73_1_cZo.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_74_1_1U7.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_75_1_X5F.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_76_1_TXS.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_77_1_BTO.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_78_1_Ujx.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_79_1_a1b.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_7_1_nG4.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_80_1_4qf.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_81_1_C7M.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_82_1_HAE.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_83_2_SlL.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_84_2_XtP.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_85_1_wq0.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_86_1_SZ2.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_87_2_iBV.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_88_1_2N0.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_89_1_IgN.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_8_1_gUI.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_90_1_B6z.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_91_1_LU7.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_92_1_yMt.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_93_1_oS9.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_94_1_JmQ.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_95_1_VqX.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_96_1_x2Q.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_97_1_7Tq.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_98_1_Csn.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_99_1_DzY.root', 'rfio:/castor/cern.ch/user/m/mahmad/Chib/v13dp1/Conv1/01May12/2011B/chic_SelEvts_9_1_s2K.root' ]
    bookHistos()
    loopOnEvents()
    saveHistos()
