#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

#include "NanoAODTree.h"

using namespace TMVA;


float jet_NHF;
float jet_CHF;
float jet_NEMF;
float jet_NHM;
float jet_CHM;



TMVA::Reader* BookJetMVAReader(std::string basePath, std::string weightFileName)
{
   TMVA::Reader* reader = new TMVA::Reader("!Color:!Silent");

   reader->AddVariable("jet_NHF", &jet_NHF);
   reader->AddVariable("jet_CHF", &jet_CHF);
   reader->AddVariable("jet_NEMF", &jet_NEMF);
   reader->AddVariable("jet_NumNeutralParticles", &jet_NHM);
   reader->AddVariable("jet_CHM", &jet_CHM);

   reader->BookMVA("BDT method", basePath+"/"+weightFileName);

   return reader;
}





int main(int argc, char** argv) {

  string status_sample = *(argv + 1);
  bool isMC = false;
  bool isData = false;
  if (status_sample.compare("mc") == 0) isMC = true;
  if (status_sample.compare("data") == 0) isData = true;


  string output;
  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--output") {
      if (i + 1 < argc) {
	output = argv[i+1];
	break;
      } else {
	std::cerr << "--output option requires one argument." << std::endl;
	return 1;
      }      
    }  
  }
  if(output==""){
    std::cerr << "--output argument required" << std::endl;
    return 1;
  }
    
  

  string input;
  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--input") {
      if (i + 1 < argc) {
	input = argv[i+1];
	break;
      } else {
	std::cerr << "--intput option requires one argument." << std::endl;
	return 1;
      }      
    }  
  }
  if(input==""){
    std::cerr << "--input argument required" << std::endl;
    return 1;
  }


  bool overwrite;
  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--overwrite") {
      overwrite = true;
      break;
    }
  }
 

  TFile* f_new = TFile::Open(output.c_str());
  if(f_new!=0 && !overwrite){
    cout<<output<<" already exists, please delete it before converting again"<<endl;
    return 0;
  }
  f_new = TFile::Open(output.c_str(),"RECREATE");


  TChain* oldtree = new TChain("Events");
  oldtree->Add(input.c_str());
  NanoAODTree* tree = new NanoAODTree(oldtree);

  TTree* tree_new=new TTree("BDT_tree","BDT_tree");

  //New branches

  int _nJet;
  float _Jet_BDT_sig1000_light[kJetMax];
  float _Jet_BDT_sig1000_b[kJetMax];
  float _Jet_BDT_sig1_light[kJetMax];
  float _Jet_BDT_sig1_b[kJetMax];

  tree_new->Branch("nJet",             &_nJet,             "nJet/I");
  tree_new->Branch("Jet_BDT_sig1000_light", &_Jet_BDT_sig1000_light, "Jet_BDT_sig1000_light[nJet]/F");
  tree_new->Branch("Jet_BDT_sig1000_b",     &_Jet_BDT_sig1000_b,     "Jet_BDT_sig1000_b[nJet]/F");
  tree_new->Branch("Jet_BDT_sig1_light",    &_Jet_BDT_sig1_light,    "Jet_BDT_sig1_light[nJet]/F");
  tree_new->Branch("Jet_BDT_sig1_b",        &_Jet_BDT_sig1_b,        "Jet_BDT_sig1_b[nJet]/F");

  
  TMVA::Reader* MVA_sig1000_light_reader = BookJetMVAReader("BDT_weights","/BDT_sig1000_light_noExtVar_BDT.weights.xml");
  TMVA::Reader* MVA_sig1000_b_reader = BookJetMVAReader("BDT_weights","/BDT_sig1000_b_noExtVar_BDT.weights.xml");
  TMVA::Reader* MVA_sig1_light_reader = BookJetMVAReader("BDT_weights","/BDT_sig1_light_noExtVar_BDT.weights.xml");
  TMVA::Reader* MVA_sig1_b_reader = BookJetMVAReader("BDT_weights","/BDT_sig1_b_noExtVar_BDT.weights.xml");

  int nentries = tree->GetEntries();
  cout<<"Nentries="<<nentries<<endl;

  for (int iEntry = 0; iEntry < nentries ; iEntry++){

    tree->GetEntry(iEntry);

   if(iEntry%10000==0) cout<<"Entry #"<<iEntry<<" "<< int(100*float(iEntry)/nentries)<<"%"<<endl;
  
    _nJet = tree->nJet;

    for(unsigned int i_jet=0; i_jet<_nJet;i_jet++){

      //FIXME Remove the rawFactor when nanoAOD is fixed
      jet_NHF = tree->Jet_neHEF[i_jet]/(1-tree->Jet_rawFactor[i_jet]);
      jet_CHF = tree->Jet_chHEF[i_jet]/(1-tree->Jet_rawFactor[i_jet]);
      jet_NEMF = tree->Jet_neEmEF[i_jet]/(1-tree->Jet_rawFactor[i_jet]);
      jet_NHM = tree->Jet_nNeutralConst[i_jet];
      jet_CHM = tree->Jet_nChargedConst[i_jet];      

      _Jet_BDT_sig1000_light[i_jet] = MVA_sig1000_light_reader->EvaluateMVA("BDT method");
      _Jet_BDT_sig1000_b[i_jet] = MVA_sig1000_b_reader->EvaluateMVA("BDT method");
      _Jet_BDT_sig1_light[i_jet] = MVA_sig1_light_reader->EvaluateMVA("BDT method");
      _Jet_BDT_sig1_b[i_jet] = MVA_sig1_b_reader->EvaluateMVA("BDT method");
        
    }

    tree_new->Fill();

  }

  
  f_new->cd();
  tree_new->AddFriend("Events",input.c_str());

  tree_new->Write();
  f_new->Close();

  return 0;
}
