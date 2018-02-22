#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TGraphAsymmErrors.h"
#include "TH2F.h"

#include "NanoAODTree.h"
#include "LLPTree.h"


map<pair<int,int>, vector<vector<int> > > all_combinations; //key=(nJet,kJet), element=all possible combinations of kJet among nJet
map<pair<int,int>, vector<vector<int> > > all_combinations_compl; //key=(nJet,kJet), element=all complementary of possible combinations of kJet among nJet

vector<vector<int> > getAll_combinations(int nJet, int kJet){

  vector<vector<int> > combinations;
  if(nJet==0 || kJet==0 || kJet>nJet)
    return combinations;

  pair<int,int> key = make_pair(nJet,kJet);

  if(all_combinations[key].size()>0)
    return all_combinations[key];
  
  else if(kJet==1){
    cout<<"nJet="<<nJet<<" kJet="<<kJet<<endl;
    for(int i=0; i<nJet; i++) combinations.push_back({i});
    all_combinations[key] = combinations;
  }
  
  else{
    cout<<"nJet="<<nJet<<" kJet="<<kJet<<endl;
    vector<vector<int> > combinations_km1 = getAll_combinations(nJet,kJet-1);
    for(unsigned int i=0; i<combinations_km1.size(); i++){      
      for(int j=combinations_km1[i].back()+1; j<nJet; j++){
	vector<int> new_comb = combinations_km1[i];
	new_comb.push_back(j);
	combinations.push_back(new_comb);
      }
    }
    all_combinations[key] = combinations;
  }

  return combinations;

}


vector<vector<int> > getAll_combinations_compl(int nJet, int kJet){

  vector<vector<int> > combinations_compl;
  if(nJet==0 || kJet==0 || kJet>nJet)
    return combinations_compl;

  pair<int,int> key = make_pair(nJet,kJet);

  if(all_combinations_compl[key].size()>0)
    return all_combinations_compl[key];

  vector<vector<int> > combinations = getAll_combinations(nJet,kJet);
  
  for(unsigned int i_comb; i_comb<combinations.size(); i_comb++){
    vector<int> comb_compl;
    for(int i=0; i<nJet; i++){
      bool in_comb=false;
      for(auto j : combinations[i_comb]){
	if(i==j){
	  in_comb=true;
	  break;
	}
      }
      if(in_comb) continue;
      comb_compl.push_back(i);
    }
    combinations_compl.push_back(comb_compl);
  }

  all_combinations_compl[key] = combinations_compl;
  return combinations_compl;

}



float fakerate_from_TGraph(TGraphAsymmErrors* graph, float pt){

  int n = graph->GetN();

  double x, y;
  graph->GetPoint(0,x,y);
  if(pt<x)
    return y;

  graph->GetPoint(n-1,x,y);
  if(pt>x)
    return y;
  

  for(int i=0; i<n;i++){

    graph->GetPoint(i,x,y);
    double err_low_x = graph->GetErrorXlow(i);
    double err_high_x = graph->GetErrorXhigh(i);

    if( (x-err_low_x) <= pt && pt < (x+err_high_x) )
      return y;

  }

  return -1.;

}




float fakerate_from_TH2(TH2F* graph, float pt, float eta){

  float absEta = fabs(eta);
  int nbins_x = graph->GetXaxis()->GetNbins();
  float pt_sat = min(pt,float(0.9999999*graph->GetXaxis()->GetBinLowEdge(nbins_x+1)));
  
  int bin = graph->FindBin(pt_sat,absEta);
  float fake_rate = graph->GetBinContent(bin);

  if(fake_rate>0)
    return fake_rate;

  return -1.;

}





int main(int argc, char** argv) {

  bool isMC = false;

  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--mc") {
      isMC = true;
      break;
    }  
  }

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



  string input_LLP;
  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--input_LLP") {
      if (i + 1 < argc) {
	input_LLP = argv[i+1];
	break;
      } else {
	std::cerr << "--intput_LLP option requires one argument." << std::endl;
	return 1;
      }      
    }  
  }
  if(input_LLP==""){
    std::cerr << "--input_LLP argument required" << std::endl;
    return 1;
  }



  bool overwrite;
  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--overwrite") {
      overwrite = true;
      break;
    }
  }
 


  string fakerate_file_b;
  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--fakerate_file_b") {
      if (i + 1 < argc) {
	fakerate_file_b = argv[i+1];
	cout<<"Using "<<fakerate_file_b<<" for b-jets fake rate"<<endl;
	break;
      } else {
	std::cerr << "--fakerate_file_b option requires one argument." << std::endl;
	return 1;
      }      
    }  
  }
  if(fakerate_file_b==""){
    std::cerr << "--fakerate_file_b argument required" << std::endl;
    return 1;
  }


  string fakerate_file_light;
  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--fakerate_file_light") {
      if (i + 1 < argc) {
	fakerate_file_light = argv[i+1];
	cout<<"Using "<<fakerate_file_light<<" for light jets fake rate"<<endl;
	break;
      } else {
	std::cerr << "--fakerate_file_light option requires one argument." << std::endl;
	return 1;
      }      
    }  
  }
  if(fakerate_file_light==""){
    std::cerr << "--fakerate_file_light argument required" << std::endl;
    return 1;
  }


  TFile* f_new = TFile::Open(output.c_str());
  if(f_new!=0 && !overwrite){
    cout<<output<<" already exists, please delete it before converting again"<<endl;
    return 0;
  }
  f_new = TFile::Open(output.c_str(),"RECREATE");


  TFile* f_fakerate_b = TFile::Open(fakerate_file_b.c_str());
  TH2F* gr_fakerate_high_CMVA = (TH2F*)f_fakerate_b->Get("gr_CMVA_high");
  gr_fakerate_high_CMVA->SetDirectory(0);

  //TGraphAsymmErrors* gr_fakerate_high_CMVA_b = (TGraphAsymmErrors*)f_fakerate_b->Get("gr_CMVA_high_b");
  //TGraphAsymmErrors* gr_fakerate_high_CMVA_ec = (TGraphAsymmErrors*)f_fakerate_b->Get("gr_CMVA_high_ec"); 

  f_fakerate_b->Close();
  
  TFile* f_fakerate_light = TFile::Open(fakerate_file_light.c_str());
  TH2F* gr_fakerate_low_CMVA = (TH2F*)f_fakerate_light->Get("gr_CMVA_low");
  gr_fakerate_low_CMVA->SetDirectory(0);

  //TGraphAsymmErrors* gr_fakerate_low_CMVA_b = (TGraphAsymmErrors*)f_fakerate_light->Get("gr_CMVA_low_b");
  //TGraphAsymmErrors* gr_fakerate_low_CMVA_ec = (TGraphAsymmErrors*)f_fakerate_light->Get("gr_CMVA_low_ec");

  f_fakerate_light->Close();


  TChain* oldLLPtree = new TChain("LLPtree");
  oldLLPtree->Add(input_LLP.c_str());
  LLPTree* LLPtree = new LLPTree(oldLLPtree);
  LLPtree->GetEntries();

  TChain* oldtree = new TChain("Events");
  oldtree->Add(input.c_str());
  NanoAODTree* tree = new NanoAODTree(oldtree);  


  TTree* tree_new=new TTree("FakeRate_tree","FakeRate_tree");
 
  //New branches
  int _nJet;
  int _nJet_sel;
  float _Jet_FR[kJetMax];

  float _EventWeight_FR_JetBin[kJetMax];
  float _EventWeight_FR;

  tree_new->Branch("nJet",             &_nJet,             "nJet/I");
  tree_new->Branch("nJet_sel",         &_nJet_sel,         "nJet_sel/I");
  tree_new->Branch("Jet_FR",           &_Jet_FR,           "Jet_FR[nJet]/F");

  tree_new->Branch("EventWeight_FR_JetBin",           &_EventWeight_FR_JetBin,           "EventWeight_FR_JetBin[nJet_sel]/F");
  tree_new->Branch("EventWeight_FR",           &_EventWeight_FR,           "EventWeight_FR/F");



  for (int iEntry = 0; iEntry < tree->GetEntries() ; iEntry++){

    tree->GetEntry(iEntry);
    LLPtree->GetEntry(iEntry);


    if(iEntry%1000==0) 
      cout<<"Entry #"<<iEntry<<endl;

    _nJet = tree->nJet;

    for(unsigned int i_jet=0; i_jet<_nJet;i_jet++){

      float pt = tree->Jet_pt[i_jet];
      float eta = tree->Jet_eta[i_jet];
      float CMVA = tree->Jet_btagCMVA[i_jet];

      _Jet_FR[i_jet] = -1;

      if(CMVA>0.4432){
	/*if(abs(eta)<1.479)
	  _Jet_FR[i_jet] = fakerate_from_TGraph(gr_fakerate_high_CMVA_b,pt);
	else if(abs(eta)>1.479 && abs(eta)<2.4)
	_Jet_FR[i_jet] = fakerate_from_TGraph(gr_fakerate_high_CMVA_ec,pt);*/
	_Jet_FR[i_jet] = fakerate_from_TH2(gr_fakerate_high_CMVA,pt,eta);
      }

      else{
	/*if(abs(eta)<1.479)
	  _Jet_FR[i_jet] = fakerate_from_TGraph(gr_fakerate_low_CMVA_b,pt);
	else if(abs(eta)>1.479 && abs(eta)<2.4)
	_Jet_FR[i_jet] = fakerate_from_TGraph(gr_fakerate_low_CMVA_ec,pt);*/
	_Jet_FR[i_jet] = fakerate_from_TH2(gr_fakerate_low_CMVA,pt,eta);
      }

    }


    _nJet_sel = LLPtree->nJet_sel;
    _EventWeight_FR = 0;
    
    for(int kJet = 1; kJet<=_nJet_sel; kJet++){
      
      _EventWeight_FR_JetBin[kJet-1] = 0;      

      vector<vector<int> > combinations = getAll_combinations(_nJet_sel,kJet);
      vector<vector<int> > combinations_compl = getAll_combinations_compl(_nJet_sel,kJet);

      for(unsigned int i_comb = 0; i_comb<combinations.size(); i_comb++){

	float EventWeight_FR_comb = 1;
	for(auto i_jet : combinations[i_comb]){
	  int i_jet_sel = LLPtree->Jet_sel_index[i_jet];	  
	  EventWeight_FR_comb *= _Jet_FR[i_jet_sel];
	}
	for(auto i_jet : combinations_compl[i_comb]){
	  int i_jet_sel = LLPtree->Jet_sel_index[i_jet];	  
	  EventWeight_FR_comb *= (1-_Jet_FR[i_jet_sel]);
	}
	
	_EventWeight_FR_JetBin[kJet-1] += EventWeight_FR_comb;

      }

      _EventWeight_FR += _EventWeight_FR_JetBin[kJet-1];


    }

    tree_new->Fill();

  }

  
  f_new->cd();
  tree_new->AddFriend("Events",input.c_str());
  tree_new->Write();
  f_new->Close();

  return 0;


}
