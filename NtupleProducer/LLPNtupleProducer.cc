#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"

#include "NanoAODTree.h"



int main(int argc, char** argv) {

  string status_sample = *(argv + 1);
  bool isMC = false;
  bool isData = false;
  if (status_sample.compare("mc") == 0) isMC = true;
  if (status_sample.compare("data") == 0) isData = true;
  string output = *(argv + 2);
  
  string input = *(argv + 3);

  bool overwrite = false;
  if(argv[3])
    overwrite = *(argv + 4);

  TFile* f_new = TFile::Open(output.c_str());
  if(f_new!=0 && !overwrite){
    cout<<output<<" already exists, please delete it before converting again"<<endl;
    return 0;
  }
  f_new = TFile::Open(output.c_str(),"RECREATE");


  TChain* oldtree = new TChain("Events");
  oldtree->Add(input.c_str());
  
  NanoAODTree* tree = new NanoAODTree(oldtree);

  TTree* tree_new=new TTree("LLPtree","LLPtree");
  tree_new->AddFriend(oldtree);

  //New branches
  int _nMuon;
  bool _Muon_isPresel[kMuonMax];
  int _nMuon_presel;
  int _Muon_presel_index[kMuonMax];

  int _nElectron;
  bool _Electron_isPresel[kElectronMax];
  int _nElectron_presel;
  int _Electron_presel_index[kElectronMax];

  int _nJet;
  bool _Jet_isPresel[kJetMax];
  int _nJet_presel;
  int _Jet_presel_index[kJetMax];  

  tree_new->Branch("nMuon",             &_nMuon,             "nMuon/I");
  tree_new->Branch("Muon_isPresel",     &_Muon_isPresel,     "Muon_isPresel[nMuon]/O");
  tree_new->Branch("nMuon_presel",      &_nMuon_presel,      "nMuon_presel/I");
  tree_new->Branch("Muon_presel_index", &_Muon_presel_index, "_Muon_presel_index[nMuon_presel]/I");

  tree_new->Branch("nElectron",             &_nElectron,             "nElectron/I");
  tree_new->Branch("Electron_isPresel",     &_Electron_isPresel,     "Electron_isPresel[nElectron]/O");
  tree_new->Branch("nElectron_presel",      &_nElectron_presel,      "nElectron_presel/I");
  tree_new->Branch("Electron_presel_index", &_Electron_presel_index, "Electron_presel_index[nElectron_presel]/I");

  tree_new->Branch("nJet",             &_nJet,             "nJet/I");
  tree_new->Branch("Jet_isPresel",     &_Jet_isPresel,     "Jet_isPresel[nElectron]/O");
  tree_new->Branch("nJet_presel",      &_nJet_presel,      "nJet_presel/I");
  tree_new->Branch("Jet_presel_index", &_Jet_presel_index, "Jet_presel_index[nElectron_presel]/I");



  for (int iEntry = 0; iEntry < tree->GetEntries() ; iEntry++){


    tree->GetEntry(iEntry);

    if(iEntry%1000==0) cout<<"Entry #"<<iEntry<<endl;

    _nMuon = tree->nMuon;
    _nMuon_presel = 0;

    for(unsigned int i_mu=0; i_mu<_nMuon;i_mu++){

      bool ispresel = tree->Muon_pt[i_mu]>5;
      ispresel &= abs(tree->Muon_eta[i_mu])<2.4;
      ispresel &= abs(tree->Muon_dxy[i_mu])<0.05;
      ispresel &= abs(tree->Muon_dz[i_mu])<0.1;
      ispresel &= tree->Muon_pfRelIso04_all[i_mu]<0.4;
      ispresel &= tree->Muon_mediumId[i_mu];   
      
      _Muon_isPresel[i_mu] = ispresel;

      if(ispresel){
	_Muon_presel_index[_nMuon_presel] = i_mu;
	_nMuon_presel++;
      }

    }


    _nElectron = tree->nElectron;
    _nElectron_presel = 0;

    for(unsigned int i_ele=0; i_ele<_nElectron;i_ele++){

      bool ispresel = tree->Electron_pt[i_ele]>7;
      ispresel &= abs(tree->Electron_eta[i_ele])<2.5;
      ispresel &= abs(tree->Electron_dxy[i_ele])<0.05;
      ispresel &= abs(tree->Electron_dz[i_ele])<0.1;
      ispresel &= tree->Electron_pfRelIso03_all[i_ele]<0.4;
      ispresel &= tree->Electron_mvaSpring16GP_WP80[i_ele];
      ispresel &= tree->Electron_convVeto[i_ele];
      ispresel &= tree->Electron_lostHits[i_ele]<=1;

      _Electron_isPresel[i_ele] = ispresel;

      if(ispresel){
	_Electron_presel_index[_nElectron_presel] = i_ele;
	_nElectron_presel++;
      }

    }

    
    _nJet = tree->nJet;
    _nJet_presel = 0;

    for(unsigned int i_jet=0; i_jet<_nJet;i_jet++){

      bool ispresel = tree->Jet_pt[i_jet]>25;
      ispresel &= abs(tree->Jet_eta[i_jet])<2.4;
      ispresel &= tree->Jet_jetId[i_jet]>0;

      //Lepton veto
      int muId1 = tree->Jet_muonIdx1[i_jet];
      int muId2 = tree->Jet_muonIdx2[i_jet];
      int eleId1 = tree->Jet_electronIdx1[i_jet];
      int eleId2 = tree->Jet_electronIdx2[i_jet];
      
      if(muId1>0 && _Muon_isPresel[muId1]) ispresel = false;
      else if(muId2>0 && _Muon_isPresel[muId2]) ispresel = false;
      else if(eleId1>0 && _Electron_isPresel[eleId1]) ispresel = false;
      else if(eleId2>0 && _Electron_isPresel[eleId2]) ispresel = false;

      _Jet_isPresel[i_jet] = ispresel;

      if(ispresel){
	_Jet_presel_index[_nJet_presel] = i_jet;
	_nJet_presel++;
      }

    }
            

    tree_new->Fill();

  }

  
  f_new->cd();

  tree_new->GetListOfFriends()->Remove(tree_new->GetListOfFriends()->At(0));
  tree_new->AddFriend("Events",input.c_str());

  tree_new->Write();
  f_new->Close();

  return 0;
}