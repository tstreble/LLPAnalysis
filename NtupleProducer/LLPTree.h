// Author: T. Strebler (IC)
// Date:   19 January 2018
//
// Wrapper for LLP tree
// Can be either included in a interpreted macro or compiled in c++
// (use `root-config --glibs --cflags`)
//
// Create the LLPTree object from the pointer to the tree, then access the stored objects from it
// Common TTree functions GetEntry (entry), GetEntries() are implemented


#ifndef LLPTREE_H
#define LLPTREE_H

#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>


using namespace std;

class LLPTree {
public :
   TChain          *_tree; 

   //Old LLP branches used  
   int nJet_sel;
   int Jet_sel_index[kJetMax];

   // methods
   LLPTree (TChain* tree);
   ~LLPTree();
   void Init(TChain* tree);
   Int_t GetEntry(int entry);
   Long64_t GetEntries();
   TChain* GetTree();

};



LLPTree::LLPTree (TChain* tree)
{
    Init(tree);
}


LLPTree::~LLPTree() {}



void LLPTree::Init(TChain* tree)
{

  // Set branch addresses and branch pointers
  if (!tree) return;
  _tree = tree;  
  _tree->SetMakeClass(1); // needed especially when compiling  

  _tree->SetBranchAddress("nJet_sel",&nJet_sel);  
  _tree->SetBranchAddress("Jet_sel_index",&Jet_sel_index);

}


Int_t LLPTree::GetEntry(int entry)
{
    return _tree->GetEntry(entry);
} 

Long64_t LLPTree::GetEntries()
{
    return _tree->GetEntries();
}

TChain* LLPTree::GetTree()
{
    return _tree;
}

#endif
