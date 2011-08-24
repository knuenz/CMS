namespace ToyMC{

double ScenarioSig [3][9]={{-1,-0.5,0,0.5,1,0,0,0.5,1},{0,0,0,0,0,0.5,-0.5,-0.75,-1},{0,0,0,0,0,0,0,0,0}};//lamth_Signal,lamph_Signal,lamthph_Sigmal

double ScenarioBkg [3][9]={{-1,-0.5,0,0.5,1,0,0,0.5,1},{0,0,0,0,0,0.5,-0.5,-0.75,-1},{0,0,0,0,0,0,0,0,0}};////lamth_Bkg,lamph_Bkg,lamthph_Bkg


int MarkerStyle[3] = {0, 20, 21};
int MarkerStyle2[3] = {0, 24, 25};
int MarkerStyle3[3] = {0, 26, 32};

int MarkerColor[3] = {0, 601, 632};
int MarkerColor2[3] = {0, 632, 632};
int MarkerColor3[3] = {0, 418, 632};


const int nPtBins=8;
const int nRapBins=2;

double ptRange[nPtBins+1]={0,5,10,15,20,25,30,50,100};

double ptCentre[nPtBins]={2.5,7.5,12.5,17.5,22.5,27.5,40,75};

double ptCentreErr[nPtBins]={2.5,2.5,2.5,2.5,2.5,2.5,10,25};

int numEvents[nRapBins][nPtBins]={{5722, 46093,27989,10212,3459,1405,1036,89},{15483, 63466,30923,10423,3606,1500,1076,103}};

double fracBackground[nRapBins][nPtBins]={{0.179,0.167,0.120,0.085,0.087,0.072,0.136,0.184},{0.272,0.218,0.156,0.116,0.105,0.103,0.129,0.214}};

const int nEffs=3;
const int FidCuts=3;

}
