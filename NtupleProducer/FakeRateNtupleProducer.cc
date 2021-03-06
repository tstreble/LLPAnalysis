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






map<int,float> fakerate_from_TH2(TH2F* graph, float pt, float eta){

  float absEta = fabs(eta);
  int nbins_x = graph->GetXaxis()->GetNbins();
  float pt_sat = min(pt,float(0.9999999*graph->GetXaxis()->GetBinLowEdge(nbins_x+1)));
  
  int bin = graph->FindBin(pt_sat,absEta);
  float fake_rate = graph->GetBinContent(bin);
  float fake_rate_up = fake_rate + graph->GetBinError(bin);
  float fake_rate_down = max(fake_rate - graph->GetBinError(bin),0.);

  if(fake_rate_down<=0. && graph->GetXaxis()->FindBin(pt)>0 && graph->GetYaxis()->FindBin(absEta)<=graph->GetYaxis()->GetNbins()){
    cout<<"Negative FR down, please investigate"<<endl;
    cout<<"pt="<<pt<<" eta="<<eta<<endl;
  }

  map<int,float> FR;
  FR[0] = fake_rate;
  FR[+1] = fake_rate_up;
  FR[-1] = fake_rate_down;
  
  if(FR[0]>0)
    return FR;

  FR[0] = -1;
  FR[+1] = -1;
  FR[-1] = -1;

  return FR;

}





map<TString,float> fakerate_from_FRreader(pair<TH2F*,map<TString,TGraphAsymmErrors*> > reader, float pt, float eta){

  TH2F* histo = reader.first;
  float absEta = fabs(eta);
  int nbins_x = histo->GetXaxis()->GetNbins();
  int nbins_y = histo->GetYaxis()->GetNbins();
  float pt_sat = min(pt,float(0.9999999*histo->GetXaxis()->GetBinLowEdge(nbins_x+1)));
  

  map<TString,TGraphAsymmErrors*> graphs = reader.second;

  int bin = histo->FindBin(pt_sat,absEta);
  map<TString,float> FR;
  for(auto& gr : graphs) FR[gr.first] = gr.second->Eval(bin);

  if(FR["nomin"]>0 && pt>histo->GetXaxis()->GetBinLowEdge(1) && absEta<histo->GetYaxis()->GetBinLowEdge(nbins_y+1))
    return FR;

  for(auto& gr : graphs) FR[gr.first] = -1;

  return FR;

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


  string fakerate_file_light_quark;
  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--fakerate_file_light_quark") {
      if (i + 1 < argc) {
	fakerate_file_light_quark = argv[i+1];
	cout<<"Using "<<fakerate_file_light_quark<<" for light jets quark fake rate"<<endl;
	break;
      } else {
	std::cerr << "--fakerate_file_light_quark option requires one argument." << std::endl;
	return 1;
      }      
    }  
  }
  if(fakerate_file_light_quark==""){
    std::cerr << "--fakerate_file_light_quark argument required" << std::endl;
    return 1;
  }



  string fakerate_file_light_gluon;
  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--fakerate_file_light_gluon") {
      if (i + 1 < argc) {
	fakerate_file_light_gluon = argv[i+1];
	cout<<"Using "<<fakerate_file_light_gluon<<" for light jets gluon fake rate"<<endl;
	break;
      } else {
	std::cerr << "--fakerate_file_light_gluon option requires one argument." << std::endl;
	return 1;
      }      
    }  
  }
  if(fakerate_file_light_gluon==""){
    std::cerr << "--fakerate_file_light_gluon argument required" << std::endl;
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
  map<TString,TGraphAsymmErrors*> gr_fakerate_high_CMVA_unfold;
  gr_fakerate_high_CMVA_unfold["nomin"] = (TGraphAsymmErrors*)f_fakerate_b->Get("gr_CMVA_high_unfold");
  gr_fakerate_high_CMVA_unfold["up"] = (TGraphAsymmErrors*)f_fakerate_b->Get("gr_CMVA_high_unfold_Up");
  gr_fakerate_high_CMVA_unfold["down"] = (TGraphAsymmErrors*)f_fakerate_b->Get("gr_CMVA_high_unfold_Down");
  gr_fakerate_high_CMVA_unfold["ptUp"] = (TGraphAsymmErrors*)f_fakerate_b->Get("gr_CMVA_high_unfold_ptUp");
  gr_fakerate_high_CMVA_unfold["ptDown"] = (TGraphAsymmErrors*)f_fakerate_b->Get("gr_CMVA_high_unfold_ptDown");
  f_fakerate_b->Close();
  pair<TH2F*,map<TString,TGraphAsymmErrors*> > FR_reader_high_CMVA = make_pair(gr_fakerate_high_CMVA,gr_fakerate_high_CMVA_unfold);
  
  TFile* f_fakerate_light_quark = TFile::Open(fakerate_file_light_quark.c_str());
  TH2F* gr_fakerate_low_CMVA_quark = (TH2F*)f_fakerate_light_quark->Get("gr_CMVA_low_QGL_high");
  gr_fakerate_low_CMVA_quark->SetDirectory(0);
  map<TString,TGraphAsymmErrors*> gr_fakerate_low_CMVA_quark_unfold;
  gr_fakerate_low_CMVA_quark_unfold["nomin"] = (TGraphAsymmErrors*)f_fakerate_light_quark->Get("gr_CMVA_low_QGL_high_unfold");
  gr_fakerate_low_CMVA_quark_unfold["up"] = (TGraphAsymmErrors*)f_fakerate_light_quark->Get("gr_CMVA_low_QGL_high_unfold_Up");
  gr_fakerate_low_CMVA_quark_unfold["down"] = (TGraphAsymmErrors*)f_fakerate_light_quark->Get("gr_CMVA_low_QGL_high_unfold_Down");
  gr_fakerate_low_CMVA_quark_unfold["ptUp"] = (TGraphAsymmErrors*)f_fakerate_light_quark->Get("gr_CMVA_low_QGL_high_unfold_ptUp");
  gr_fakerate_low_CMVA_quark_unfold["ptDown"] = (TGraphAsymmErrors*)f_fakerate_light_quark->Get("gr_CMVA_low_QGL_high_unfold_ptDown");
  f_fakerate_light_quark->Close();
  pair<TH2F*,map<TString,TGraphAsymmErrors*> > FR_reader_low_CMVA_quark = make_pair(gr_fakerate_low_CMVA_quark,gr_fakerate_low_CMVA_quark_unfold);


  TFile* f_fakerate_light_gluon = TFile::Open(fakerate_file_light_gluon.c_str());
  TH2F* gr_fakerate_low_CMVA_gluon = (TH2F*)f_fakerate_light_gluon->Get("gr_CMVA_low_QGL_low");
  gr_fakerate_low_CMVA_gluon->SetDirectory(0);
  map<TString,TGraphAsymmErrors*> gr_fakerate_low_CMVA_gluon_unfold;
  gr_fakerate_low_CMVA_gluon_unfold["nomin"] = (TGraphAsymmErrors*)f_fakerate_light_gluon->Get("gr_CMVA_low_QGL_low_unfold");
  gr_fakerate_low_CMVA_gluon_unfold["up"] = (TGraphAsymmErrors*)f_fakerate_light_gluon->Get("gr_CMVA_low_QGL_low_unfold_Up");
  gr_fakerate_low_CMVA_gluon_unfold["down"] = (TGraphAsymmErrors*)f_fakerate_light_gluon->Get("gr_CMVA_low_QGL_low_unfold_Down");
  gr_fakerate_low_CMVA_gluon_unfold["ptUp"] = (TGraphAsymmErrors*)f_fakerate_light_gluon->Get("gr_CMVA_low_QGL_low_unfold_ptUp");
  gr_fakerate_low_CMVA_gluon_unfold["ptDown"] = (TGraphAsymmErrors*)f_fakerate_light_gluon->Get("gr_CMVA_low_QGL_low_unfold_ptDown");
  f_fakerate_light_gluon->Close();
  pair<TH2F*,map<TString,TGraphAsymmErrors*> > FR_reader_low_CMVA_gluon = make_pair(gr_fakerate_low_CMVA_gluon,gr_fakerate_low_CMVA_gluon_unfold);


  TChain* oldLLPtree = new TChain("LLP_tree");
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
  float _Jet_FR_b_up[kJetMax];
  float _Jet_FR_b_down[kJetMax];
  float _Jet_FR_q_up[kJetMax];
  float _Jet_FR_q_down[kJetMax];
  float _Jet_FR_g_up[kJetMax];
  float _Jet_FR_g_down[kJetMax];
  float _Jet_FR_b_pt_up[kJetMax];
  float _Jet_FR_b_pt_down[kJetMax];
  float _Jet_FR_q_pt_up[kJetMax];
  float _Jet_FR_q_pt_down[kJetMax];
  float _Jet_FR_g_pt_up[kJetMax];
  float _Jet_FR_g_pt_down[kJetMax];

  float _EventWeight_FR_JetBin[kJetMax];
  float _EventWeight_FR_JetBin_b_up[kJetMax];
  float _EventWeight_FR_JetBin_b_down[kJetMax];
  float _EventWeight_FR_JetBin_q_up[kJetMax];
  float _EventWeight_FR_JetBin_q_down[kJetMax];
  float _EventWeight_FR_JetBin_g_up[kJetMax];
  float _EventWeight_FR_JetBin_g_down[kJetMax];
  float _EventWeight_FR_JetBin_b_pt_up[kJetMax];
  float _EventWeight_FR_JetBin_b_pt_down[kJetMax];
  float _EventWeight_FR_JetBin_q_pt_up[kJetMax];
  float _EventWeight_FR_JetBin_q_pt_down[kJetMax];
  float _EventWeight_FR_JetBin_g_pt_up[kJetMax];
  float _EventWeight_FR_JetBin_g_pt_down[kJetMax];

  float _EventWeight_FR;
  float _EventWeight_FR_b_up;
  float _EventWeight_FR_b_down;
  float _EventWeight_FR_q_up;
  float _EventWeight_FR_q_down;
  float _EventWeight_FR_g_up;
  float _EventWeight_FR_g_down;
  float _EventWeight_FR_b_pt_up;
  float _EventWeight_FR_b_pt_down;
  float _EventWeight_FR_q_pt_up;
  float _EventWeight_FR_q_pt_down;
  float _EventWeight_FR_g_pt_up;
  float _EventWeight_FR_g_pt_down;

  tree_new->Branch("nJet",             &_nJet,             "nJet/I");
  tree_new->Branch("nJet_sel",         &_nJet_sel,         "nJet_sel/I");
  tree_new->Branch("Jet_FR",           &_Jet_FR,           "Jet_FR[nJet]/F");
  tree_new->Branch("Jet_FR_b_up",      &_Jet_FR_b_up,        "Jet_FR_b_up[nJet]/F");
  tree_new->Branch("Jet_FR_b_down",    &_Jet_FR_b_down,      "Jet_FR_b_down[nJet]/F");
  tree_new->Branch("Jet_FR_q_up",      &_Jet_FR_q_up,        "Jet_FR_q_up[nJet]/F");
  tree_new->Branch("Jet_FR_q_down",    &_Jet_FR_q_down,      "Jet_FR_q_down[nJet]/F");
  tree_new->Branch("Jet_FR_g_up",      &_Jet_FR_g_up,        "Jet_FR_g_up[nJet]/F");
  tree_new->Branch("Jet_FR_g_down",    &_Jet_FR_g_down,      "Jet_FR_g_down[nJet]/F");
  tree_new->Branch("Jet_FR_b_pt_up",   &_Jet_FR_b_pt_up,     "Jet_FR_b_pt_up[nJet]/F");
  tree_new->Branch("Jet_FR_b_pt_down", &_Jet_FR_b_pt_down,   "Jet_FR_b_pt_down[nJet]/F");
  tree_new->Branch("Jet_FR_q_pt_up",   &_Jet_FR_q_pt_up,     "Jet_FR_q_pt_up[nJet]/F");
  tree_new->Branch("Jet_FR_q_pt_down", &_Jet_FR_q_pt_down,   "Jet_FR_q_pt_down[nJet]/F");
  tree_new->Branch("Jet_FR_g_pt_up",   &_Jet_FR_g_pt_up,     "Jet_FR_g_pt_up[nJet]/F");
  tree_new->Branch("Jet_FR_g_pt_down", &_Jet_FR_g_pt_down,   "Jet_FR_g_pt_down[nJet]/F");

  tree_new->Branch("EventWeight_FR_JetBin",           &_EventWeight_FR_JetBin,           "EventWeight_FR_JetBin[nJet_sel]/F");
  tree_new->Branch("EventWeight_FR_JetBin_b_up",      &_EventWeight_FR_JetBin_b_up,      "EventWeight_FR_JetBin_b_up[nJet_sel]/F");
  tree_new->Branch("EventWeight_FR_JetBin_b_down",    &_EventWeight_FR_JetBin_b_down,    "EventWeight_FR_JetBin_b_down[nJet_sel]/F");
  tree_new->Branch("EventWeight_FR_JetBin_q_up",      &_EventWeight_FR_JetBin_q_up,      "EventWeight_FR_JetBin_q_up[nJet_sel]/F");
  tree_new->Branch("EventWeight_FR_JetBin_q_down",    &_EventWeight_FR_JetBin_q_down,    "EventWeight_FR_JetBin_q_down[nJet_sel]/F");
  tree_new->Branch("EventWeight_FR_JetBin_g_up",      &_EventWeight_FR_JetBin_g_up,      "EventWeight_FR_JetBin_g_up[nJet_sel]/F");
  tree_new->Branch("EventWeight_FR_JetBin_q_down",    &_EventWeight_FR_JetBin_g_down,    "EventWeight_FR_JetBin_g_down[nJet_sel]/F");
  tree_new->Branch("EventWeight_FR_JetBin_b_pt_up",   &_EventWeight_FR_JetBin_b_pt_up,      "EventWeight_FR_JetBin_b_pt_up[nJet_sel]/F");
  tree_new->Branch("EventWeight_FR_JetBin_b_pt_down", &_EventWeight_FR_JetBin_b_pt_down,    "EventWeight_FR_JetBin_b_pt_down[nJet_sel]/F");
  tree_new->Branch("EventWeight_FR_JetBin_q_pt_up",   &_EventWeight_FR_JetBin_q_pt_up,      "EventWeight_FR_JetBin_q_pt_up[nJet_sel]/F");
  tree_new->Branch("EventWeight_FR_JetBin_q_pt_down", &_EventWeight_FR_JetBin_q_pt_down,    "EventWeight_FR_JetBin_q_pt_down[nJet_sel]/F");
  tree_new->Branch("EventWeight_FR_JetBin_g_pt_up",   &_EventWeight_FR_JetBin_g_pt_up,      "EventWeight_FR_JetBin_g_pt_up[nJet_sel]/F");
  tree_new->Branch("EventWeight_FR_JetBin_q_pt_down", &_EventWeight_FR_JetBin_g_pt_down,    "EventWeight_FR_JetBin_g_pt_down[nJet_sel]/F");


  tree_new->Branch("EventWeight_FR",           &_EventWeight_FR,           "EventWeight_FR/F");
  tree_new->Branch("EventWeight_FR_b_up",      &_EventWeight_FR_b_up,        "EventWeight_FR_b_up/F");
  tree_new->Branch("EventWeight_FR_b_down",    &_EventWeight_FR_b_down,      "EventWeight_FR_b_down/F");
  tree_new->Branch("EventWeight_FR_q_up",      &_EventWeight_FR_q_up,        "EventWeight_FR_q_up/F");
  tree_new->Branch("EventWeight_FR_q_down",    &_EventWeight_FR_q_down,      "EventWeight_FR_q_down/F");
  tree_new->Branch("EventWeight_FR_g_up",      &_EventWeight_FR_g_up,        "EventWeight_FR_g_up/F");
  tree_new->Branch("EventWeight_FR_g_down",    &_EventWeight_FR_g_down,      "EventWeight_FR_g_down/F");
  tree_new->Branch("EventWeight_FR_b_pt_up",   &_EventWeight_FR_b_pt_up,     "EventWeight_FR_b_pt_up/F");
  tree_new->Branch("EventWeight_FR_b_pt_down", &_EventWeight_FR_b_pt_down,   "EventWeight_FR_b_pt_down/F");
  tree_new->Branch("EventWeight_FR_q_pt_up",   &_EventWeight_FR_q_pt_up,     "EventWeight_FR_q_pt_up/F");
  tree_new->Branch("EventWeight_FR_q_pt_down", &_EventWeight_FR_q_pt_down,   "EventWeight_FR_q_pt_down/F");
  tree_new->Branch("EventWeight_FR_g_pt_up",   &_EventWeight_FR_g_pt_up,     "EventWeight_FR_g_pt_up/F");
  tree_new->Branch("EventWeight_FR_g_pt_down", &_EventWeight_FR_g_pt_down,   "EventWeight_FR_g_pt_down/F");


  int nentries = tree->GetEntries();
  cout<<"Nentries="<<nentries<<endl;

  for (int iEntry = 0; iEntry < tree->GetEntries() ; iEntry++){

    tree->GetEntry(iEntry);
    LLPtree->GetEntry(iEntry);

    if(iEntry%10000==0) cout<<"Entry #"<<iEntry<<" "<< int(100*float(iEntry)/nentries)<<"%"<<endl;


    _nJet = tree->nJet;

    for(unsigned int i_jet=0; i_jet<_nJet;i_jet++){

      float pt = tree->Jet_pt[i_jet];
      float eta = tree->Jet_eta[i_jet];
      float CMVA = tree->Jet_btagCMVA[i_jet];
      float qgl = tree->Jet_qgl[i_jet];

      map<TString,float> FR;

      if(CMVA>0.4432){

	FR = fakerate_from_FRreader(FR_reader_high_CMVA,pt,eta);	

	_Jet_FR[i_jet] = FR["nomin"];
	_Jet_FR_b_up[i_jet] = FR["up"];
	_Jet_FR_b_down[i_jet] = FR["down"];
	_Jet_FR_q_up[i_jet] = FR["nomin"];
	_Jet_FR_q_down[i_jet] = FR["nomin"];
	_Jet_FR_g_up[i_jet] = FR["nomin"];
	_Jet_FR_g_down[i_jet] = FR["nomin"];
	_Jet_FR_b_pt_up[i_jet] = FR["ptUp"];
	_Jet_FR_b_pt_down[i_jet] = FR["ptDown"];
	_Jet_FR_q_pt_up[i_jet] = FR["nomin"];
	_Jet_FR_q_pt_down[i_jet] = FR["nomin"];
	_Jet_FR_g_pt_up[i_jet] = FR["nomin"];
	_Jet_FR_g_pt_down[i_jet] = FR["nomin"];

      }

      else{
	if(qgl>0.5){

	  FR = fakerate_from_FRreader(FR_reader_low_CMVA_quark,pt,eta);	

	  _Jet_FR[i_jet] = FR["nomin"];
	  _Jet_FR_b_up[i_jet] = FR["nomin"];
	  _Jet_FR_b_down[i_jet] = FR["nomin"];
	  _Jet_FR_q_up[i_jet] = FR["up"];
	  _Jet_FR_q_down[i_jet] = FR["down"];
	  _Jet_FR_g_up[i_jet] = FR["nomin"];
	  _Jet_FR_g_down[i_jet] = FR["nomin"];
	  _Jet_FR_b_pt_up[i_jet] = FR["nomin"];
	  _Jet_FR_b_pt_down[i_jet] = FR["nomin"];
	  _Jet_FR_q_pt_up[i_jet] = FR["ptUp"];
	  _Jet_FR_q_pt_down[i_jet] = FR["ptDown"];
	  _Jet_FR_g_pt_up[i_jet] = FR["nomin"];
	  _Jet_FR_g_pt_down[i_jet] = FR["nomin"];
  
	}
	else{

	  FR = fakerate_from_FRreader(FR_reader_low_CMVA_gluon,pt,eta);	

	  _Jet_FR[i_jet] = FR["nomin"];
	  _Jet_FR_b_up[i_jet] = FR["nomin"];
	  _Jet_FR_b_down[i_jet] = FR["nomin"];
	  _Jet_FR_q_up[i_jet] = FR["nomin"];
	  _Jet_FR_q_down[i_jet] = FR["nomin"];
	  _Jet_FR_g_up[i_jet] = FR["up"];
	  _Jet_FR_g_down[i_jet] = FR["down"];
	  _Jet_FR_b_pt_up[i_jet] = FR["nomin"];
	  _Jet_FR_b_pt_down[i_jet] = FR["nomin"];
	  _Jet_FR_q_pt_up[i_jet] = FR["nomin"];
	  _Jet_FR_q_pt_down[i_jet] = FR["nomin"];
	  _Jet_FR_g_pt_up[i_jet] = FR["ptUp"];
	  _Jet_FR_g_pt_down[i_jet] = FR["ptDown"];

	}
      }


      
    }


    _nJet_sel = LLPtree->nJet_sel;
    _EventWeight_FR = 0;
    _EventWeight_FR_b_up = 0;
    _EventWeight_FR_b_down = 0;
    _EventWeight_FR_q_up = 0;
    _EventWeight_FR_q_down = 0;
    _EventWeight_FR_g_up = 0;
    _EventWeight_FR_g_down = 0;
    _EventWeight_FR_b_pt_up = 0;
    _EventWeight_FR_b_pt_down = 0;
    _EventWeight_FR_q_pt_up = 0;
    _EventWeight_FR_q_pt_down = 0;
    _EventWeight_FR_g_pt_up = 0;
    _EventWeight_FR_g_pt_down = 0;
    
    for(int kJet = 1; kJet<=_nJet_sel; kJet++){
      
      _EventWeight_FR_JetBin[kJet-1] = 0;      
      _EventWeight_FR_JetBin_b_up[kJet-1] = 0;      
      _EventWeight_FR_JetBin_b_down[kJet-1] = 0;      
      _EventWeight_FR_JetBin_q_up[kJet-1] = 0;      
      _EventWeight_FR_JetBin_q_down[kJet-1] = 0;    
      _EventWeight_FR_JetBin_g_up[kJet-1] = 0;      
      _EventWeight_FR_JetBin_g_down[kJet-1] = 0;   
      _EventWeight_FR_JetBin_b_pt_up[kJet-1] = 0;      
      _EventWeight_FR_JetBin_b_pt_down[kJet-1] = 0;      
      _EventWeight_FR_JetBin_q_pt_up[kJet-1] = 0;      
      _EventWeight_FR_JetBin_q_pt_down[kJet-1] = 0;    
      _EventWeight_FR_JetBin_g_pt_up[kJet-1] = 0;      
      _EventWeight_FR_JetBin_g_pt_down[kJet-1] = 0;        

      vector<vector<int> > combinations = getAll_combinations(_nJet_sel,kJet);
      vector<vector<int> > combinations_compl = getAll_combinations_compl(_nJet_sel,kJet);

      for(unsigned int i_comb = 0; i_comb<combinations.size(); i_comb++){

	float EventWeight_FR_comb = 1;
	float EventWeight_FR_comb_b_up = 1;
	float EventWeight_FR_comb_b_down = 1;
	float EventWeight_FR_comb_q_up = 1;
	float EventWeight_FR_comb_q_down = 1;
	float EventWeight_FR_comb_g_up = 1;
	float EventWeight_FR_comb_g_down = 1;
	float EventWeight_FR_comb_b_pt_up = 1;
	float EventWeight_FR_comb_b_pt_down = 1;
	float EventWeight_FR_comb_q_pt_up = 1;
	float EventWeight_FR_comb_q_pt_down = 1;
	float EventWeight_FR_comb_g_pt_up = 1;
	float EventWeight_FR_comb_g_pt_down = 1;


	for(auto i_jet : combinations[i_comb]){
	  int i_jet_sel = LLPtree->Jet_sel_index[i_jet];	  
	  EventWeight_FR_comb *= _Jet_FR[i_jet_sel];
	  EventWeight_FR_comb_b_up *= _Jet_FR_b_up[i_jet_sel];
	  EventWeight_FR_comb_b_down *= _Jet_FR_b_down[i_jet_sel];
	  EventWeight_FR_comb_q_up *= _Jet_FR_q_up[i_jet_sel];
	  EventWeight_FR_comb_q_down *= _Jet_FR_q_down[i_jet_sel];
	  EventWeight_FR_comb_g_up *= _Jet_FR_g_up[i_jet_sel];
	  EventWeight_FR_comb_g_down *= _Jet_FR_g_down[i_jet_sel];
	  EventWeight_FR_comb_b_pt_up *= _Jet_FR_b_pt_up[i_jet_sel];
	  EventWeight_FR_comb_b_pt_down *= _Jet_FR_b_pt_down[i_jet_sel];
	  EventWeight_FR_comb_q_pt_up *= _Jet_FR_q_pt_up[i_jet_sel];
	  EventWeight_FR_comb_q_pt_down *= _Jet_FR_q_pt_down[i_jet_sel];
	  EventWeight_FR_comb_g_pt_up *= _Jet_FR_g_pt_up[i_jet_sel];
	  EventWeight_FR_comb_g_pt_down *= _Jet_FR_g_pt_down[i_jet_sel];
	}
	for(auto i_jet : combinations_compl[i_comb]){
	  int i_jet_sel = LLPtree->Jet_sel_index[i_jet];	  
	  EventWeight_FR_comb *= (1-_Jet_FR[i_jet_sel]);
	  EventWeight_FR_comb_b_up *= (1-_Jet_FR_b_up[i_jet_sel]);
	  EventWeight_FR_comb_b_down *= (1-_Jet_FR_b_down[i_jet_sel]);
	  EventWeight_FR_comb_q_up *= (1-_Jet_FR_q_up[i_jet_sel]);
	  EventWeight_FR_comb_q_down *= (1-_Jet_FR_q_down[i_jet_sel]);
	  EventWeight_FR_comb_g_up *= (1-_Jet_FR_q_up[i_jet_sel]);
	  EventWeight_FR_comb_g_down *= (1-_Jet_FR_q_down[i_jet_sel]);
	  EventWeight_FR_comb_b_pt_up *= (1-_Jet_FR_b_pt_up[i_jet_sel]);
	  EventWeight_FR_comb_b_pt_down *= (1-_Jet_FR_b_pt_down[i_jet_sel]);
	  EventWeight_FR_comb_q_pt_up *= (1-_Jet_FR_q_pt_up[i_jet_sel]);
	  EventWeight_FR_comb_q_pt_down *= (1-_Jet_FR_q_pt_down[i_jet_sel]);
	  EventWeight_FR_comb_g_pt_up *= (1-_Jet_FR_q_pt_up[i_jet_sel]);
	  EventWeight_FR_comb_g_pt_down *= (1-_Jet_FR_q_pt_down[i_jet_sel]);
	}
	
	_EventWeight_FR_JetBin[kJet-1] += EventWeight_FR_comb;
	_EventWeight_FR_JetBin_b_up[kJet-1] += EventWeight_FR_comb_b_up;
	_EventWeight_FR_JetBin_b_down[kJet-1] += EventWeight_FR_comb_b_down;
	_EventWeight_FR_JetBin_q_up[kJet-1] += EventWeight_FR_comb_q_up;
	_EventWeight_FR_JetBin_q_down[kJet-1] += EventWeight_FR_comb_q_down;
	_EventWeight_FR_JetBin_g_up[kJet-1] += EventWeight_FR_comb_g_up;
	_EventWeight_FR_JetBin_g_down[kJet-1] += EventWeight_FR_comb_g_down;
	_EventWeight_FR_JetBin_b_pt_up[kJet-1] += EventWeight_FR_comb_b_pt_up;
	_EventWeight_FR_JetBin_b_pt_down[kJet-1] += EventWeight_FR_comb_b_pt_down;
	_EventWeight_FR_JetBin_q_pt_up[kJet-1] += EventWeight_FR_comb_q_pt_up;
	_EventWeight_FR_JetBin_q_pt_down[kJet-1] += EventWeight_FR_comb_q_pt_down;
	_EventWeight_FR_JetBin_g_pt_up[kJet-1] += EventWeight_FR_comb_g_pt_up;
	_EventWeight_FR_JetBin_g_pt_down[kJet-1] += EventWeight_FR_comb_g_pt_down;

      }

      _EventWeight_FR += _EventWeight_FR_JetBin[kJet-1];
      _EventWeight_FR_b_up += _EventWeight_FR_JetBin_b_up[kJet-1];
      _EventWeight_FR_b_down += _EventWeight_FR_JetBin_b_down[kJet-1];
      _EventWeight_FR_q_up += _EventWeight_FR_JetBin_q_up[kJet-1];
      _EventWeight_FR_q_down += _EventWeight_FR_JetBin_q_down[kJet-1];
      _EventWeight_FR_g_up += _EventWeight_FR_JetBin_g_up[kJet-1];
      _EventWeight_FR_g_down += _EventWeight_FR_JetBin_g_down[kJet-1];
      _EventWeight_FR_b_pt_up += _EventWeight_FR_JetBin_b_pt_up[kJet-1];
      _EventWeight_FR_b_pt_down += _EventWeight_FR_JetBin_b_pt_down[kJet-1];
      _EventWeight_FR_q_pt_up += _EventWeight_FR_JetBin_q_pt_up[kJet-1];
      _EventWeight_FR_q_pt_down += _EventWeight_FR_JetBin_q_pt_down[kJet-1];
      _EventWeight_FR_g_pt_up += _EventWeight_FR_JetBin_g_pt_up[kJet-1];
      _EventWeight_FR_g_pt_down += _EventWeight_FR_JetBin_g_pt_down[kJet-1];

    }

    tree_new->Fill();

  }

  
  f_new->cd();
  tree_new->AddFriend("Events",input.c_str());
  tree_new->Write();
  f_new->Close();

  return 0;


}
