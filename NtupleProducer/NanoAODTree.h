// Author: T. Strebler (IC)
// Date:   19 January 2018
//
// Wrapper for NanoAOD tree
// Can be either included in a interpreted macro or compiled in c++
// (use `root-config --glibs --cflags`)
//
// Create the nanoAODTree object from the pointer to the tree, then access the stored objects from it
// Common TTree functions GetEntry (entry), GetEntries() are implemented



#ifndef NANOAODTREE_H
#define NANOAODTREE_H

#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>

const int kMuonMax = 100;
const int kElectronMax = 100;
const int kJetMax = 100;

using namespace std;

class NanoAODTree {
public :
   TChain          *_tree; 

   //Old NanoAOD branches used
   int run;
   int luminosityBlock;
   int event;

   int nMuon;
   float Muon_pt[kMuonMax];
   float Muon_eta[kMuonMax];
   float Muon_phi[kMuonMax];
   float Muon_mass[kMuonMax];
   float Muon_dxy[kMuonMax];
   float Muon_dz[kMuonMax];
   float Muon_pfRelIso04_all[kMuonMax];
   bool Muon_mediumId[kMuonMax];

   int nElectron;
   float Electron_pt[kElectronMax];
   float Electron_eta[kElectronMax];
   float Electron_phi[kElectronMax];
   float Electron_mass[kElectronMax];
   float Electron_dxy[kElectronMax];
   float Electron_dz[kElectronMax];
   float Electron_pfRelIso03_all[kElectronMax];
   bool Electron_mvaSpring16GP_WP80[kElectronMax];
   bool Electron_convVeto[kElectronMax];
   int Electron_lostHits[kElectronMax];

   int nJet;
   float Jet_pt[kJetMax];
   float Jet_eta[kJetMax];
   int Jet_jetId[kJetMax];
   int Jet_muonIdx1[kJetMax];
   int Jet_muonIdx2[kJetMax];
   int Jet_electronIdx1[kJetMax];
   int Jet_electronIdx2[kJetMax];
   float Jet_chHEF[kJetMax];
   float Jet_neHEF[kJetMax];
   float Jet_neEmEF[kJetMax];   
   int Jet_nChargedConst[kJetMax];
   int Jet_nNeutralConst[kJetMax];
   float Jet_rawFactor[kJetMax];
   float Jet_btagCMVA[kJetMax];
   int Jet_puId[kJetMax];
   float Jet_qgl[kJetMax];
   int Jet_partonFlavour[kJetMax];


   // methods
   NanoAODTree (TChain* tree);
   ~NanoAODTree();
   void Init(TChain* tree);
   Int_t GetEntry(int entry);
   Long64_t GetEntries();
   TChain* GetTree();

};



NanoAODTree::NanoAODTree (TChain* tree)
{
    Init(tree);
}


NanoAODTree::~NanoAODTree() {}



void NanoAODTree::Init(TChain* tree)
{

  // Set branch addresses and branch pointers
  if (!tree) return;
  _tree = tree;  
  _tree->SetMakeClass(1); // needed especially when compiling
  
  _tree->SetBranchAddress("run",&run);
  _tree->SetBranchAddress("luminosityBlock",&luminosityBlock);
  _tree->SetBranchAddress("event",&event);

  _tree->SetBranchAddress("nMuon",&nMuon);  
  _tree->SetBranchAddress("Muon_pt",&Muon_pt);  
  _tree->SetBranchAddress("Muon_eta",&Muon_eta);
  _tree->SetBranchAddress("Muon_phi",&Muon_phi);
  _tree->SetBranchAddress("Muon_mass",&Muon_mass);
  _tree->SetBranchAddress("Muon_dxy",&Muon_dxy);
  _tree->SetBranchAddress("Muon_dz",&Muon_dz);
  _tree->SetBranchAddress("Muon_pfRelIso04_all",&Muon_pfRelIso04_all);
  _tree->SetBranchAddress("Muon_mediumId",&Muon_mediumId);

  _tree->SetBranchAddress("nElectron",&nElectron);  
  _tree->SetBranchAddress("Electron_pt",&Electron_pt);  
  _tree->SetBranchAddress("Electron_eta",&Electron_eta);
  _tree->SetBranchAddress("Electron_phi",&Electron_phi);
  _tree->SetBranchAddress("Electron_mass",&Electron_mass);
  _tree->SetBranchAddress("Electron_dxy",&Electron_dxy);
  _tree->SetBranchAddress("Electron_dz",&Electron_dz);
  _tree->SetBranchAddress("Electron_pfRelIso03_all",&Electron_pfRelIso03_all);
  _tree->SetBranchAddress("Electron_mvaSpring16GP_WP80",&Electron_mvaSpring16GP_WP80);
  _tree->SetBranchAddress("Electron_convVeto",&Electron_convVeto);
  _tree->SetBranchAddress("Electron_lostHits",&Electron_lostHits);

  _tree->SetBranchAddress("nJet",&nJet);  
  _tree->SetBranchAddress("Jet_pt",&Jet_pt);  
  _tree->SetBranchAddress("Jet_eta",&Jet_eta);
  _tree->SetBranchAddress("Jet_jetId",&Jet_jetId);  
  _tree->SetBranchAddress("Jet_muonIdx1",&Jet_muonIdx1);
  _tree->SetBranchAddress("Jet_muonIdx2",&Jet_muonIdx2);
  _tree->SetBranchAddress("Jet_electronIdx1",&Jet_electronIdx1);
  _tree->SetBranchAddress("Jet_electronIdx2",&Jet_electronIdx2);
  _tree->SetBranchAddress("Jet_chHEF",&Jet_chHEF);
  _tree->SetBranchAddress("Jet_neHEF",&Jet_neHEF);
  _tree->SetBranchAddress("Jet_neEmEF",&Jet_neEmEF);
  _tree->SetBranchAddress("Jet_nChargedConst",&Jet_nChargedConst);
  _tree->SetBranchAddress("Jet_nNeutralConst",&Jet_nNeutralConst);
  _tree->SetBranchAddress("Jet_rawFactor",&Jet_rawFactor);
  _tree->SetBranchAddress("Jet_btagCMVA",&Jet_btagCMVA);
  _tree->SetBranchAddress("Jet_puId",&Jet_puId);
  _tree->SetBranchAddress("Jet_qgl",&Jet_qgl);
  _tree->SetBranchAddress("Jet_partonFlavour",&Jet_partonFlavour);

}


Int_t NanoAODTree::GetEntry(int entry)
{
    return _tree->GetEntry(entry);
} 

Long64_t NanoAODTree::GetEntries()
{
    return _tree->GetEntries();
}

TChain* NanoAODTree::GetTree()
{
    return _tree;
}

#endif
