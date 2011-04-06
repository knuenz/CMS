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
  
//  0.926,0.645,0.785,0.158,4.344 for  B0,B+,Bs,Lambdab,prompt
  double lumiB0 = 0.926;
  double lumiBs = 0.785;
  double lumiBp = 0.645;
  double lumiLambdaB = 0.158;

//  double luminorm = 100;

  files.push_back(base_path+std::string("TTree_prep_B0_.root"));
  weights.push_back(float(lumiB0));

  files.push_back(base_path+std::string("TTree_prep_Bs_.root"));
  weights.push_back(float(lumiBs));

  files.push_back(base_path+std::string("TTree_prep_Bp_.root"));
  weights.push_back(float(lumiBp));

  files.push_back(base_path+std::string("TTree_prep_LambdaB_.root"));
  weights.push_back(float(lumiLambdaB));








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
