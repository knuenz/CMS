// root -l -b -q addWeights.C
#include <string>
//#include <pair>
#include <vector>
void addWeights()
{
  std::vector<std::string> files;
  std::vector<float> weights;
  std::string base_path("/scratch/knuenz/Polarization/RootInput/");
  std::string treename("data");
  
  //Scale all samples to the sample with least amount of int lum (Min Bias)

  double MCfraction =  0.898735;
  double lumi = 2*156.63;

  files.push_back(base_path+std::string("TTree_red_LambdaBerr.root"));
  weights.push_back(float(10./(lumi*MCfraction)));

  files.push_back(base_path+std::string("TTree_red_LambdaBerr_pseudo.root"));
  weights.push_back(float(10./(lumi*(1-MCfraction))));
/*
  files.push_back(base_path+std::string("Spring10/nonPromptJpsi/TTree_pol_Mu0Track0Jpsi_BpToJPsiMuMu.root"));
  weights.push_back(float(1.0525/36.699));

  files.push_back(base_path+std::string("Spring10/nonPromptJpsi/TTree_pol_Mu0Track0Jpsi_BsToJPsiMuMu.root"));
  weights.push_back(float(1.0525/33.498));

  files.push_back(base_path+std::string("Spring10/nonPromptJpsi/TTree_pol_Mu0Track0Jpsi_LambdaBToJPsiMuMu.root"));
  weights.push_back(float(1.0525/156.638));

  files.push_back(base_path+std::string("Spring10/promptJpsi/TTree_pol_Mu0Track0Jpsi_MCprompt.root"));
  weights.push_back(float(1.0525/13.708));

  files.push_back(base_path+std::string("Spring10/MinBias/TTree_pol_Mu0Track0Jpsi_MCMinBias.root"));
  weights.push_back(float(1.0));
*/
  for(Int_t i = 0; i < files.size(); ++i)
    {
      std::cout<<"Opening: "<< files[i];
      std::cout << " and loading TTree, " << treename << std::endl;
      TFile* currentFile = new TFile(files[i].c_str(),"update");

      TTree *theTree = (TTree*)currentFile->Get(treename.c_str());

      Int_t numEntries = (Int_t)theTree->GetEntries();

      Float_t weight;

      TBranch *newBranch = theTree->Branch("MCweight",&weight,"MCweight/F");

      std::cout << "Weight is: " << weights[i] << std::endl;
      std::cout << "Editing " << numEntries << " entries." << std::endl;
      for(Int_t j=0; j<numEntries; ++j)
	{
	  weight=weights[i];
	  theTree->GetEntry(j);

	  newBranch->Fill();
	}

      theTree->Write("",TObject::kOverwrite);
      currentFile->Close();
      delete currentFile;
    }


}
