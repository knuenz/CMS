Float_t PtMin = 0;
Float_t PtMax = 40;
Float_t rapMin = -2.3;
Float_t rapMax = 2.3;
Float_t JpsictMin = -1;
Float_t JpsictMax = 2.5;
double JpsiMassMin=2.7;
double JpsiMassMax=3.5;

int Angular_Binning_costh = nBinsCosT/2;
int Angular_Binning_phi = nBinsPhiPol/2;

const int numberofpTbins=1;//kNbPTBins;
const int numberofRapbins=1;//kNbRapForPTBins;
int NpTBinStart=4;
int NRapBinStart=1;

Int_t AnalyticIndex = 0; //0...without analytical component fits
Int_t SBComparisonIndex = 0; //0...no SB Comparison


RooRealVar JpsiMass("JpsiMass","M [GeV]",JpsiMassMin,JpsiMassMax);
RooRealVar JpsiRap("JpsiRap","#nu",rapMin,rapMax);
RooRealVar Jpsict("Jpsict","l_{J/#psi} [mm]",JpsictMin,JpsictMax);
RooRealVar JpsiPt("JpsiPt","pT [GeV]",PtMin,PtMax);
RooRealVar costh_CS("costh_CS","cos #theta_{CS}",-1,1);
RooRealVar phi_CS("phi_CS","#phi_{CS} [deg]",0,360);
RooRealVar costh_HX("costh_HX","cos#theta_{HX}",-1,1);
RooRealVar phi_HX("phi_HX","#phi_{HX} [deg]",0,360);
RooRealVar MCType_idx("MCType_idx","MCType_idx",0,2);//0=PR,1=NP,2=BK
RooRealVar JpsiType_idx("JpsiType_idx","JpsiType_idx",0,1.5);//0=GG,1=GT,2=TT
RooRealVar MCweight("MCweight","MCweight",0,1.1);


double pi = TMath::Pi();

RooPlot* JpsiMassframeMCTruthtestCB;
RooPlot* JpsictframePromptAnal;
RooPlot* JpsictframePromptAnalChi2;
RooPlot* JpsictframeNonPromptAnal;
RooPlot* JpsictframeNonPromptAnalChi2;
RooPlot* JpsictframeBackgroundAnalRealDataSBR;
RooPlot* JpsictframeBackgroundAnalRealDataSBRchi2;
RooPlot* JpsictframeBackgroundAnalRealDataSBL;
RooPlot* JpsictframeBackgroundAnalRealDataSBLchi2;
RooPlot* phi_CSframe;
RooPlot* costh_CSframe;
RooPlot* phi_HXframe;
RooPlot* costh_HXframe;
RooPlot* JpsiMassframe;
RooPlot* JpsiMassframeChi2;
RooPlot* JpsictframePrompt;
RooPlot* JpsictframePromptAnal;
RooPlot* JpsictframePromptAnalChi2;
RooPlot* JpsictframeBackground;
RooPlot* JpsictframeBackgroundAnal;
RooPlot* JpsictframeBackgroundAnalChi2;
RooPlot* JpsictframeNonPrompt;
RooPlot* JpsictframeNonPromptAnal;
RooPlot* JpsictframeNonPromptAnalChi2;
RooPlot* JpsictframeDataChi2;
RooPlot* JpsictframeSB;
RooPlot* JpsictframeChi2;
RooPlot* Jpsictframe;

double pTMin;
double pTMax;
double RapMin;
double RapMax;

double pTmean[numberofpTbins],errpTmean[numberofpTbins];
double yrap1[numberofpTbins],eyrap1[numberofpTbins],yrap2[numberofpTbins],eyrap2[numberofpTbins],yrap3[numberofpTbins],eyrap3[numberofpTbins];
double yrap1MCTRUTH[numberofpTbins],eyrap1MCTRUTH[numberofpTbins],yrap2MCTRUTH[numberofpTbins],eyrap2MCTRUTH[numberofpTbins],yrap3MCTRUTH[numberofpTbins],eyrap3MCTRUTH[numberofpTbins];
double ypol_CS_1[numberofpTbins],eypol_CS_1[numberofpTbins],ypol_HX_1[numberofpTbins],eypol_HX_1[numberofpTbins],ypol_CS_2[numberofpTbins],eypol_CS_2[numberofpTbins],ypol_HX_2[numberofpTbins],eypol_HX_2[numberofpTbins],ypol_CS_3[numberofpTbins],eypol_CS_3[numberofpTbins],ypol_HX_3[numberofpTbins],eypol_HX_3[numberofpTbins];
double ypol_CS_1_rap2[numberofpTbins],eypol_CS_1_rap2[numberofpTbins];
double ypol_HX_1_rap2[numberofpTbins],eypol_HX_1_rap2[numberofpTbins],ypol_CS_2_rap2[numberofpTbins],eypol_CS_2_rap2[numberofpTbins],ypol_HX_2_rap2[numberofpTbins],eypol_HX_2_rap2[numberofpTbins],ypol_CS_3_rap2[numberofpTbins],eypol_CS_3_rap2[numberofpTbins],ypol_HX_3_rap2[numberofpTbins],eypol_HX_3_rap2[numberofpTbins];
double ypol_CS_1_rap3[numberofpTbins],eypol_CS_1_rap3[numberofpTbins], ypol_HX_1_rap3[numberofpTbins],eypol_HX_1_rap3[numberofpTbins],ypol_CS_2_rap3[numberofpTbins],eypol_CS_2_rap3[numberofpTbins],ypol_HX_2_rap3[numberofpTbins],eypol_HX_2_rap3[numberofpTbins],ypol_CS_3_rap3[numberofpTbins],eypol_CS_3_rap3[numberofpTbins],ypol_HX_3_rap3[numberofpTbins],eypol_HX_3_rap3[numberofpTbins];

double BFRACTION;
double ERRBFRACTION;

RooDataSet *reddata, *reddataMCtest;

TLegend *BFractionLegend;
TLegend *SBLegend;
TLegend *GlobalLegend;
TLegend *BackgroundLegendAnalDataSBL;
TLegend *BackgroundLegendAnalDataSBR;
TLegend *BackgroundLegendAnal;
TLegend *BackgroundLegend;
TLegend *NonPromptLegendAnal;
TLegend *NonPromptLegend;
TLegend *PromptLegendAnal;
TLegend *PromptLegend;
TLegend *MassLegend;

TH1* legendMassdata;
TH1* legendMassDashed;
TH1* legendMassDashedred;
TH1* legendMassDashedgreen;
TH1* legendMass;
TH1* legendMassdatared;
TH1* legendMassdatablue;
TH1* legendMassdataStyle;
TH1* legendMassdataredStyle;
TH1* legendMassdatablueStyle;
TH1* legendMassdatagreenStyle;
TH1* legendMassdatayellowStyle;

RooDataSet* realdata000;

TH1* polfunchisto_CS;
TH1* datahistoangular_CS;
TH1* polfunchisto_HX;
TH1* datahistoangular_HX;

TFile *fIn2;
TFile *fIn3;
TFile *fIn4;
TTree* treeData;
TTree* treeDataTest;
RooDataSet *realdata00;
RooDataSet *realdata;

char reducestr[200];
char reducestr2[200];
char reducestr3[200];
char reducestr4[200];
char reducestr5[200];
char reducestr6[200];
char splitMC[200];
char splitMCtest[200];
char reduceMassData[200];
char reduceSignalwindowData[200];

RooDataSet *MCModelDataSet;
RooDataSet *MCTestDataSet;

RooRealVar lambda;
RooRealVar NbkgMass;
RooExponential expo;
RooExtendPdf eexpo;
RooRealVar CBm;
RooRealVar CBs;
RooRealVar CBa;
RooRealVar CBn;
RooCBShape CBMass;
RooGaussian gauss;
RooRealVar fracgauss;
RooAddPdf gaussCB;
RooRealVar NsigMassCB;
RooExtendPdf eCBMass;
RooAddPdf sumJpsiMass2;

RooRealVar ctmeanSig1;
RooRealVar ctsigmaSig1;
RooRealVar ctcoefSig1;
RooGaussModel ctgaussSig1;
RooExtendPdf ctegaussSig1;
RooRealVar ctsigmaSig2;
RooRealVar ctcoefSig2;
RooGaussModel ctgaussSig2;
RooExtendPdf ctegaussSig2;
RooRealVar ctmeanSig3;
RooRealVar ctsigmaSig3;
RooRealVar ctcoefSig3;
RooGaussModel ctgaussSig3;
RooExtendPdf ctegaussSig3;
RooAddModel sumJpsictPrompt;
RooRealVar ctcoefPrompt;
RooExtendPdf esumJpsictPrompt;

RooDataHist *realhData;
RooDataSet *hltIndex0PR, *hltIndex0NP, *hltIndex0BK, *hltIndex0PRMCtest, *hltIndex0NPMCtest, *hltIndex0BKMCtest, *SignalDataSet;
RooDataHist *hData;
RooDataHist *hDataMCtest;

RooDataHist *hltIndex0NPhist;
RooDataHist *hltIndex0PRhist;
RooDataHist *hltIndex0BKhist;

RooFitResult* fitresMass;

double NsigTot;double errNsigTot;double NbkgTot;double errNbkgTot;

double errfracgauss;
double fracgaussdouble;
double resol;
double errresol;
double meanSig1double;
double CBmean;
double errCBmean;
double JpsiMass3SigmaMin;
double JpsiMass3SigmaMax;
double JpsiMass6SigmaMin;
double JpsiMass6SigmaMax;

RooFitResult* fitresctPrompt;
double ctmeanSig1double;
double ctsigmaSig1double;
double ctsigmaSig2double;
double ctsigmaSig3double;
double ctcoefSig1double;
double ctcoefSig2double;
double ctcoefSig3double;
double ctcoefSigdoubleall;
double ctcoefSig1doublenorm;
double ctcoefSig2doublenorm;
double ctcoefSig3doublenorm;

RooDataSet *hltIndex0PRredMCtest;
RooDataSet *hltIndex0NPredMCtest;
RooDataSet *hltIndex0BKredMCtest;
RooDataHist *hDatared;
RooDataHist *hDataredMCtest;
RooDataSet *redreddata;
RooDataSet *hltIndex0PRred,*hltIndex0NPred,*hltIndex0BKred;
RooDataSet *redreddataMCtest;
RooDataSet *hltIndex0PRredMCtest,*hltIndex0NPredMCtest,*hltIndex0BKredMCtest;
RooDataSet *realdataBkg;
RooDataHist *realdataBkghist;
RooDataSet *realdataSBR;
RooDataHist *realdataSBRhistfit;
double NSBR;
RooDataSet *realdataSBL;
RooDataHist *realdataSBLhistfit;
double NSBL;
double SBscale;
RooDataHist *realdataSBLscale;
TH1* realdataSBLscale0;
RooDataHist *realdataSBLhist;
RooHistPdf histrealSBL;
RooRealVar coefctrealSBL;
RooDataHist *realdataSBRhist;
RooHistPdf histrealSBR;
RooRealVar coefctrealSBR;
RooAddPdf histrealSB;
RooDataSet *MCdataBkg;
RooDataHist *MCdataBkghist;
RooDataSet *MCdataSBL;
RooDataSet *MCdataSBR;
RooDataHist *MCdataBkghistctau;
RooHistPdf histMCSB;
RooDataSet *redrealdata;
RooDataHist *redrealhData;
RooDataSet *hltIndex0PRredcut;
RooDataSet *hltIndex0NPredcut;
RooDataSet *hltIndex0PRredcutMCtest;
RooDataSet *hltIndex0NPredcutMCtest;
RooDataHist *hltIndex0PRredhistMCtest;
RooDataHist *hltIndex0NPredhistMCtest;
RooDataHist *hltIndex0BKredhistMCtest;

RooDataHist *hltIndex0PRhistbin;
RooHistPdf histPrompt;
RooRealVar coefctPrompt;
RooExtendPdf ehistPrompt;
RooDataHist *hltIndex0BKhistbin;
RooHistPdf histBackground;
RooRealVar coefctBackground;
RooExtendPdf ehistBackground;
RooDataHist *hltIndex0NPhistbin;
RooHistPdf histNonPrompt;
RooRealVar coefctNonPrompt;
RooExtendPdf ehistNonPrompt;
RooExtendPdf ehistMCSB;
RooExtendPdf ehistrealSB;
RooDataHist *hltIndex0PRredhist;
RooDataHist *hltIndex0BKredhist;

RooRealVar ctmeanSig1fixnotfix;
RooRealVar ctmeanSig1fix;
RooRealVar ctsigmaSig1fix;
RooRealVar ctcoefSig1fix;
RooGaussModel ctgaussSig1fixnotfix;
RooGaussModel ctgaussSig1fix;
RooRealVar ctsigmaSig2fix;
RooRealVar ctcoefSig2fix;
RooGaussModel ctgaussSig2fix;
RooGaussModel ctgaussSig2fixnotfix;
RooRealVar ctmeanSig3fix;
RooRealVar ctsigmaSig3fix;
RooRealVar ctcoefSig3fix;
RooGaussModel ctgaussSig3fix;
RooGaussModel ctgaussSig3fixnotfix;
RooAddModel sumJpsictPromptfix;
RooAddModel sumJpsictPromptfixnotfix;
RooRealVar ctcoefPromptfix;
RooExtendPdf esumJpsictPromptfix;
RooExtendPdf esumJpsictPromptfixnotfix;
RooRealVar tau2NonPrompt;
RooDecay NonPrompt2;
RooRealVar coefctNonPrompt2;
RooExtendPdf eNonPrompt2;
RooRealVar tau3NonPrompt;
RooDecay NonPrompt3;
RooRealVar coefctNonPrompt3;
RooExtendPdf eNonPrompt3;
RooRealVar tau4NonPrompt;
RooDecay NonPrompt4;
RooRealVar coefctNonPrompt4;
RooExtendPdf eNonPrompt4;
RooAddPdf NonPromptctTOT;
RooRealVar tau2;
RooDecay bkg2;
RooRealVar coefctbkg2;
RooExtendPdf ebkg2;
RooRealVar tau3;
RooDecay bkg3;
RooRealVar coefctbkg3;
RooExtendPdf ebkg3;
RooRealVar tau4;
RooDecay bkg4;
RooRealVar coefctbkg4;
RooExtendPdf ebkg4;
RooAddPdf bkgctTOT;
RooExtendPdf ebkgctTOT;
RooRealVar coefctNonPromptanal;
RooExtendPdf eNonPromptctTOT;


RooAbsReal* integralSignal;
double integralNsigRange;
RooAbsReal* integralbkg;
double integralNbkgRange;
double NsigRange;
double errNsigTotscaled;
double NbkgRange;
double errNbkgTotscaled;
double sigOvbkg;
double errsigOvbkg;


RooFitResult* fitresctPromptTemplate;
RooFitResult* fitresctNonPromptTemplate;
RooFitResult* fitresctBackgroundTemplate;
RooFitResult* fitresctNonPrompt;
RooFitResult* fitresctBackground;
RooFitResult* fitresGlobal;
RooFitResult* fitresctBackgroundRealDataSBR;
RooFitResult* fitresctBackgroundRealDataSBL;

RooAddPdf sumJpsict;

double chi2ctGlobal;

double lambda_theta_CS_double;
double lambda_phi_CS_double;
double lambda_thetaphi_CS_double;
double errlambda_theta_CS_double;
double errlambda_phi_CS_double;
double errlambda_thetaphi_CS_double;

double lambda_theta_HX_double;
double lambda_phi_HX_double;
double lambda_thetaphi_HX_double;
double errlambda_theta_HX_double;
double errlambda_phi_HX_double;
double errlambda_thetaphi_HX_double;

int nArray;
double xploterror;
double xplot;

int nbin;
int numberofbins;
int npT;
double numberofprocessedbins;
int nRap=1;

double NPromptForPol;double NNonPromptForPol;
RooRealVar NPromptForPolVar;
RooRealVar NNonPromptForPolVar;

RooDataHist *AngularPRhistMCtest_CS;
RooDataHist *AngularNPhistMCtest_CS;
RooDataHist *AngularBKhistMCtest_CS;
RooDataHist *AngularPRhistMCtest_HX;
RooDataHist *AngularNPhistMCtest_HX;
RooDataHist *AngularBKhistMCtest_HX;

double ctresol;
double DecayLengthCut;
RooDataSet *AngularPromptMCtest;
RooDataHist *AngularPromptMCtest_hist_CS;RooDataHist *AngularPromptMCtest_hist_HX;

RooDataSet *AngularPromptRealData;
RooDataHist *AngularPromptRealData_hist_CS;RooDataHist *AngularPromptRealData_hist_HX;

RooFitResult* fitresthetaphi_CS;
RooFitResult* fitresthetaphi_HX;

RooRealVar lambda_theta_CS;
RooRealVar lambda_phi_CS;
RooRealVar lambda_thetaphi_CS;
RooGenericPdf polfunc_CS;
RooExtendPdf epolfunc_CS;
RooProdPdf CSPDF;
//RooAddPdf CSPDF;
RooHistPdf AngularHistPdf_CS_Bkg;

RooRealVar lambda_theta_HX;
RooRealVar lambda_phi_HX;
RooRealVar lambda_thetaphi_HX;
RooGenericPdf polfunc_HX;
RooExtendPdf epolfunc_HXs;
RooProdPdf HXPDF;
//RooAddPdf HXPDF;
RooHistPdf AngularHistPdf_HX_Bkg;

double bfractionMCTRUTH;
double bfractionMCTRUTHrange;

FILE *outputFile2;

double ERRBFRACTIONRANGE;
double BFRACTIONRANGE;
double NPromptoverNNonPrompt;
double errNNonPromptTotscaled;
double errNPromptTotscaled;

double NPROMPT;
double ERRNPROMPT;
double NNONPROMPT;
double ERRNNONPROMPT;
char outputfilename[200];
double chi2Mass;
double chi2ctGlobal;
double chi2ctPrompt;
double chi2ctNonPrompt;
double chi2ctBackground;
double chi2ctPromptAnal;
double chi2ctPromptAnal;
double chi2ctGlobalData;

double NbkgTot;
double errNbkgTot;
double NPromptSig1;
double NPromptSig2;
double NPromptSig3;
double errNPromptSig1;
double errNPromptSig2;
double errNPromptSig3;
double NPromptTot;
double errNPromptTot;
double NNonPromptTot;
double errNNonPromptTot;
double NPromptTottemp ;
double errNPromptTottemp;
double NPromptTotfix;
double errNPromptTotfix;
double NPromptSig1fix;
double NPromptSig2fix;
double NPromptSig3fix;
double errNPromptSig1fix;
double errNPromptSig2fix;
double errNPromptSig3fix;
double NPromptTotfixadd;
double errNPromptTotfixadd;
double NNonPrompt2;
double NNonPrompt3;
double NNonPromptges;
double NNonPromptTotanal;
double errNNonPromptTotanal;
double LifeTimeCutMin;
double LifeTimeCutMax;
RooAbsReal* integralPromptFit;
RooAbsReal* integralNonPromptFit;
double NPromptRange;
double NNonPromptRange;

double reddatasum;
double realdata000sum;
double realdata00sum;
double realdatasum;
double reddataMCtestsum;
double hltIndex0PRsum;
double hltIndex0NPsum;
double hltIndex0BKsum;
double hltIndex0PRMCtestsum;
double hltIndex0NPMCtestsum;
double hltIndex0BKMCtestsum;
double redreddatasum;
double hltIndex0PRredsum;
double hltIndex0NPredsum;
double hltIndex0BKredsum;
double hltIndex0PRredcutsum;
double hltIndex0NPredcutsum;
double hltIndex0PRredMCtestsum;
double hltIndex0NPredMCtestsum;
double hltIndex0BKredMCtestsum;
double hltIndex0PRredcutMCtestsum;
double hltIndex0NPredcutMCtestsum ;
double realdataBkgsum;
double realdataSBLsum;
double realdataSBRsum;
double MCdataBkgsum;
double MCdataSBLsum;
double MCdataSBRsum;
double redrealdatasum;
double sigovbkgMCTRUTH;
FILE *outputFile;

TFile *fAcc;
TH2F *hAcc2D_pol_pT_rap[kNbFrames][kNbPTBins+4][kNbRapForPTBins+1];

double NBkgRange;
double SigOvBkgRange;

double invariantlambda_CS;
double invariantlambda_HX;
double invariantF_CS;
double invariantF_HX;
double invariantlambda_CS_rap1[numberofpTbins];
double invariantF_CS_rap1[numberofpTbins];
double invariantlambda_HX_rap1[numberofpTbins];
double invariantF_HX_rap1[numberofpTbins];
double invariantlambda_CS_rap2[numberofpTbins];
double invariantF_CS_rap2[numberofpTbins];
double invariantlambda_HX_rap2[numberofpTbins];
double invariantF_HX_rap2[numberofpTbins];
double invariantlambda_CS_rap3[numberofpTbins];
double invariantF_CS_rap3[numberofpTbins];
double invariantlambda_HX_rap3[numberofpTbins];
double invariantF_HX_rap3[numberofpTbins];


double errinvariantlambda_CS;
double errinvariantlambda_HX;
double errinvariantF_CS;
double errinvariantF_HX;
double errinvariantlambda_CS_rap1[numberofpTbins];
double errinvariantF_CS_rap1[numberofpTbins];
double errinvariantlambda_HX_rap1[numberofpTbins];
double errinvariantF_HX_rap1[numberofpTbins];
double errinvariantlambda_CS_rap2[numberofpTbins];
double errinvariantF_CS_rap2[numberofpTbins];
double errinvariantlambda_HX_rap2[numberofpTbins];
double errinvariantF_HX_rap2[numberofpTbins];
double errinvariantlambda_CS_rap3[numberofpTbins];
double errinvariantF_CS_rap3[numberofpTbins];
double errinvariantlambda_HX_rap3[numberofpTbins];
double errinvariantF_HX_rap3[numberofpTbins];

double chi2_phi_HX;
double chi2_costh_HX;
double chi2_phi_CS;
double chi2_costh_CS;

RooDataSet *MCPROMPT;
RooDataSet *MCNONPROMPT;
RooDataSet *MCBACKGROUND;

RooDataSet *redMCPROMPT;
RooDataSet *redMCNONPROMPT;
RooDataSet *redMCBACKGROUND;

RooDataHist *redMCPROMPThist;
RooDataHist *redMCNONPROMPThist;
RooDataHist *redMCBACKGROUNDhist;

RooDataSet *redredMCPROMPT;
RooDataSet *redredMCNONPROMPT;
RooDataSet *redredMCBACKGROUND;
