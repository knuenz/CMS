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

int numEvents[nRapBins][nPtBins]={{87520,621900,361600,131800,46650,17700,13240,1165},{160500,637300,305600,103500,35230,13220,10050,886}};

double fracBackground[nRapBins][nPtBins]={{0.17,0.17,0.12,0.08,0.08,0.09,0.12,0.17},{0.27,0.22,0.15,0.11,0.11,0.11,0.13,0.19}};

double EffCorrFrac[nRapBins][nPtBins]={{0.24,0.25,0.43,0.57,0.66,0.72,0.78,0.87},{0.24,0.25,0.43,0.57,0.65,0.71,0.77,0.88}};

}
