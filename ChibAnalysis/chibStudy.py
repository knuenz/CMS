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

invm1S     	= RooRealVar("invm1S",   "invm1S",9,65)
invm2S     	= RooRealVar("invm2S",   "invm2S",9,65)
invm3S     	= RooRealVar("invm3S",   "invm3S",9,65)

Y1Smass_nSigma  = RooRealVar("Y1Smass_nSigma","Y1Smass_nSigma",-100,100)
Y2Smass_nSigma  = RooRealVar("Y2Smass_nSigma","Y2Smass_nSigma",-100,100)
Y3Smass_nSigma  = RooRealVar("Y3Smass_nSigma","Y3Smass_nSigma",-100,100)

jpsimass   	= RooRealVar("jpsimass", "jpsimass",8,12)
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
    if(abs(_dz/_dzError) > _sigmasTkVtxComp):
        print '---->'
        print "vtx z" , vtx.z() ," \t dz " ,_dz ," \t " ,_dzError ," \t dz/dzErr " ,_dz/_dzError

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
    
    Y1Smass0=9.46030
    Y2Smass0=10.02326
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
    chicCand_l = ( "chib","chicCompCand","chibs5" )
    
    #candidate J/Psi
    jpsiCand_h  = Handle ("vector<pat::CompositeCandidate>") #sometimes reco:: oder pat::
    jpsiCand_l = ( "chib","jpsiCompCand" ,"chibs5")
    
    jpsiCandPAT_h  = Handle ("vector<pat::CompositeCandidate>")
    jpsiCandPAT_l = ( "onia2MuMuPatTrkTrk","" )
    
    #candidate gamma
    gammaCand_h  = Handle ("vector<reco::CompositeCandidate>")
    gammaCand_l = ( "chib","gammaCompCand" ,"chibs5")
    
    convCand_h = Handle ("vector<reco::Conversion>")
    convCand_l = ( "allConversions","" )
    
    vtxCand_h = Handle('vector<reco::Vertex>')
    vtxCand_l = ('offlinePrimaryVertices')

    events= Events(inputfiles)

    # loop over events
    nevents = 0
    for event in events:

        if nevents % 100 == 0 :
            print 'Events Processed ' , nevents

        #if nevents > 1000:
            #break
    
        #print "run %d \t lumi %d \t orbit %d \t file %s" % (event._event.getRun().run(),event._event.luminosityBlock(), event._event.orbitNumber(),event._event.getTFile().GetName())
    

        nevents += 1
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
    print "Number of chib candidates:          ", ncands 


if __name__  == '__main__':
    
    parser = OptionParser(description='%prog : Chib Invariant Mass Fitter.',usage='chic-unbinnedlh.py --options')
    parser.add_option('--fileName',dest='fileName',default='Def',help='The name of your output file')
    (options,args) = parser.parse_args()

    # Create histograms, etc.
    gROOT.SetBatch()        # don't pop up canvases
    gROOT.SetStyle('Plain')

#    inputfiles = [ 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_05Apri112_A001.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_05Apri112_A002.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_05Apri112_A003.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_05Apri112_A004.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_05Apri112_A005.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_05Apri112_A006.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_05Apri112_A007.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_05Apri112_A008.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_05Apri112_A009.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_05Apri112_A010.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_05Apri112_A011.root' ]
#    inputfiles = [ 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_05Apri112_B.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_05Apri112_B001.root',  'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_05Apri112_B002.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_05Apri112_B003.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_05Apri112_B004.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_05Apri112_B005.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_05Apri112_B006.root']

#    inputfiles = [ 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_2011A001.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_2011A002.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_2011A003.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_2011A004.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_2011A005.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_2011A006.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_2011A007.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_2011A008.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_2011A009.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_2011A010.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_2011A011.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_2011A012.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_2011A013.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_2011A014.root' ]
#    inputfiles = [ 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_2011B001.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_2011B002.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_2011B003.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_2011B004.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_2011B005.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_2011B006.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_2011B007.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_2011B008.root', 'file:/scratch/knuenz/Polarization/RootInput/Chib/chibstep5_2011B009.root' ]
    inputfiles = [ 'file:/afs/hephy.at/scratch/k/knuenz/ErnestArea/CMSSW_4_2_3/src/chib_SelEvts.root']
     
    bookHistos()
    loopOnEvents()
    saveHistos()