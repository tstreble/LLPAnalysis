#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1F.h>
#include <TF1.h>

#include<algorithm>

#include "Helpers.C"

using namespace std;

float luminosity = 35867; //35.867 fb-1

map<int,map<int,map<int,float> > > XS_SUSY; //First index ctau, 2nd mgno, 3rd mchi

float XS_ttbar_DL = 87.3;
float XS_ttbar_SLfromT = 182;
float XS_ttbar_SLfromTbar = 182;
float XS_ttbar = 831.76;

float XS_ZNuNu_BJets = 48.71;


vector<float> XS;

vector<TString> tree = {"LLP_tree","BDT_tree","Events"};
vector<TString> tree_FR = {"LLP_tree","BDT_tree","Events","FakeRate_tree"};

//MC

map<int,vector<TString> > filename_T1qqqq_ctau;
map<int,TString> filename_norm_T1qqqq_ctau;


vector<TString> filename_ttbar = {"/vols/cms/tstreble/LL_ntuples/TTJets_incl/LLPntuples/TTJets_incl_LLP_*",
				  "/vols/cms/tstreble/LL_ntuples/TTJets_incl/BDTntuples/TTJets_incl_BDT_*",
				  "/vols/cms/tstreble/LL_ntuples/TTJets_incl/ntuples_postproc/TTJets_incl_NANO_*",
                                  "/vols/cms/tstreble/LL_ntuples/TTJets_incl/FRntuples/TTJets_incl_FR_TTbar_EMu_sig1000_b_cut0.10_*.root"};
vector<TString> filename_ZNuNu_BJets = {"/vols/cms/tstreble/LL_ntuples/DYBJets_ZNuNu/LLPntuples/DYBJets_ZNuNu_LLP.root",
					"/vols/cms/tstreble/LL_ntuples/DYBJets_ZNuNu/BDTntuples/DYBJets_ZNuNu_BDT.root",
					"/vols/cms/tstreble/LL_ntuples/DYBJets_ZNuNu/ntuples_postproc/DYBJets_ZNuNu_NANO_Skim.root"};

TString filename_norm_ttbar = "/vols/cms/tstreble/LL_ntuples/TTJets_incl/PUntuples/TTJets_incl_NANO_*";
TString filename_norm_ZNuNu_BJets = "/vols/cms/tstreble/LL_ntuples/DYBJets_ZNuNu/PUntuples/DYBJets_ZNuNu_NANO_Skim.root";


vector<vector<TString> > filename_MC;
vector<TString> filename_norm_MC;
vector<TString> MC_sample_name;
vector<TString> MC_sample_cut;
vector<TString> MC_sample_cut_norm;


vector<vector<TString> > filename_data;



// Systematics
vector<TString> fake_syst_names;



void initialize(){

  XS.clear();
  filename_MC.clear();
  filename_norm_MC.clear();
  MC_sample_name.clear();
  MC_sample_cut.clear();
  MC_sample_cut_norm.clear();

  filename_data.clear();

  for(unsigned int mgno = 600; mgno<1200; mgno+=200){
    for(unsigned int mchi = 0; mchi<500; mchi+=200){    
      XS_SUSY[1000][mgno][mchi] = 1/luminosity;
    }
  }
  
  filename_T1qqqq_ctau[1000] = {"/vols/cms/tstreble/LL_ntuples/T1qqqq/LLPntuples/T1qqqq_ctau-1000_LLP.root",
					    "/vols/cms/tstreble/LL_ntuples/T1qqqq/BDTntuples/T1qqqq_ctau-1000_BDT.root",
					    "/vols/cms/tstreble/LL_ntuples/T1qqqq/PUntuples/T1qqqq_ctau-1000_Skim.root"};

  filename_norm_T1qqqq_ctau[1000] = "/vols/cms/tstreble/LL_ntuples/T1qqqq/PUntuples/T1qqqq_ctau-1000_Skim.root";
  


  for(auto& it1 : XS_SUSY){
    int ctau = it1.first;
    for(auto& it2 : XS_SUSY[ctau]){
      int mgno = it2.first;
      for(auto& it3 : XS_SUSY[ctau][mgno]){
	int mchi = it3.first;

	XS.push_back(XS_SUSY[ctau][mgno][mchi]);
	filename_MC.push_back(filename_T1qqqq_ctau[ctau]);
	filename_norm_MC.push_back(filename_norm_T1qqqq_ctau[ctau]);
	MC_sample_name.push_back(Form("T1qqqq_ctau%i_mgno%i_mchi%i",ctau,mgno,mchi));
	MC_sample_cut.push_back(Form("Sum$(GenPart_pdgId==1000021 && GenPart_mass>%i && GenPart_mass<%i)>0 && Sum$(GenPart_pdgId==1000022 && GenPart_mass>%i && GenPart_mass<%i)>0",mgno-1,mgno+1,mchi-1,mchi+1));
	MC_sample_cut_norm.push_back(Form("Sum$(GenPart_pdgId==1000021 && GenPart_mass>%i && GenPart_mass<%i)>0 && Sum$(GenPart_pdgId==1000022 && GenPart_mass>%i && GenPart_mass<%i)>0",mgno-1,mgno+1,mchi-1,mchi+1));

      }
    }
  }



  XS.push_back(XS_ttbar);
  filename_MC.push_back(filename_ttbar);
  filename_norm_MC.push_back(filename_norm_ttbar);
  MC_sample_name.push_back("ttbar");
  MC_sample_cut.push_back("1");
  MC_sample_cut_norm.push_back("1");


  /*XS.push_back(XS_ZNuNu_BJets);
  filename_MC.push_back(filename_ZNuNu_BJets);
  filename_norm_MC.push_back(filename_norm_ZNuNu_BJets);
  MC_sample_name.push_back("ZNuNu_BJets");
  MC_sample_cut.push_back("1");
  MC_sample_cut_norm.push_back("1");*/


  fake_syst_names.push_back("b_up");
  fake_syst_names.push_back("b_down");
  fake_syst_names.push_back("b_pt_up");
  fake_syst_names.push_back("b_pt_down");
  fake_syst_names.push_back("q_up");
  fake_syst_names.push_back("q_down");
  fake_syst_names.push_back("q_pt_up");
  fake_syst_names.push_back("q_pt_down");
  fake_syst_names.push_back("g_up");
  fake_syst_names.push_back("g_down");
  fake_syst_names.push_back("g_pt_up");
  fake_syst_names.push_back("g_pt_down");

}






void datacard_maker(TString var1, int nbin, float xmin, float xmax,
		    TString cut_cat, TString jet_sel, TString file, 
		    bool syst=false, 
		    TString var1_JecUp = "", TString cut_cat_JecUp = "", TString jet_sel_JecUp="",
		    TString var1_JecDown = "", TString cut_cat_JecDown = "", TString jet_sel_JecDown=""){

  TString var = var1 + "*(" + var1 + Form("<=%f) + %f*(",xmax,0.999*xmax) + var1 + Form(">%f)",xmax);
  cout<<var<<endl;

  initialize();

  
  TFile* f_new = TFile::Open(file,"RECREATE");

  f_new->cd();


  vector<float> norm_MC;

  for(unsigned i_MC=0;i_MC<MC_sample_name.size();i_MC++){

    cout<<MC_sample_name[i_MC]<<endl;

    
    TH1F* h_MC_norm = single_plot(filename_norm_MC[i_MC],
				  "Events",
				  MC_sample_cut_norm[i_MC],
				  "genWeight*puWeight",
				  3,0,3);
    
    norm_MC.push_back(h_MC_norm->Integral());

  }

  cout<<"OK norm"<<endl;


  //////////////////////////////////////////
  //////////////////////////////////////////
  /////////          NOMINAL        ////////
  //////////////////////////////////////////
  //////////////////////////////////////////
  
  
  for(unsigned i_MC=0;i_MC<MC_sample_name.size();i_MC++){

    cout<<MC_sample_name[i_MC]<<endl;

    TH1F* h_MC = single_plot(filename_MC[i_MC],
			     tree,
			     var,
			     "genWeight*puWeight*("+MC_sample_cut[i_MC] + " && " + cut_cat + " && Sum$(" + jet_sel + ")>0)",
			     nbin,xmin,xmax);

    h_MC->Scale(luminosity*XS[i_MC]/norm_MC[i_MC]);
    h_MC->SetNameTitle("x_"+MC_sample_name[i_MC],"x_"+MC_sample_name[i_MC]);
    makeBinContentsPositive(h_MC,true);
    h_MC->Write();

    if(MC_sample_name[i_MC].Contains("T1qqqq")) continue;


    TH1F* h_MC_FR = single_plot(filename_MC[i_MC],
				tree_FR,
				var,
				"EventWeight_FR/(1-EventWeight_FR)*genWeight*puWeight*("+MC_sample_cut[i_MC] + " && " + cut_cat + "&& Sum$(" + jet_sel + ")==0)",
				nbin,xmin,xmax);

    h_MC_FR->Scale(luminosity*XS[i_MC]/norm_MC[i_MC]);
    h_MC_FR->SetNameTitle("x_"+MC_sample_name[i_MC]+"_fake","x_"+MC_sample_name[i_MC]+"_fake");
    makeBinContentsPositive(h_MC,true);
    h_MC_FR->Write();

  }


  //Fakes
    
  TH1F* h_fakes = (TH1F*)f_new->Get("x_ttbar_fake")->Clone();

  h_fakes->SetNameTitle("x_data_fakes","x_data_data_fakes");
  makeBinContentsPositive(h_fakes,true);
  h_fakes->Write();

  //Data

  TH1F* h_data_obs = (TH1F*)f_new->Get("x_ttbar")->Clone();
  
  h_data_obs->SetNameTitle("x_data_obs","x_data_obs");
  h_data_obs->Write();

  if(!syst){
    f_new->Close();
    return;
  }



  cout<<"OK nominal"<<endl;
  
  
  //////////////////////////////////////////
  //////////////////////////////////////////
  /////////       SYSTEMATICS       ////////
  //////////////////////////////////////////
  //////////////////////////////////////////

  TString var_JecUp = var1_JecUp + "*(" + var1_JecUp + Form("<=%f) + %f*(",xmax,0.999*xmax) + var1_JecUp + Form(">%f)",xmax);
  TString var_JecDown = var1_JecDown + "*(" + var1_JecDown + Form("<=%f) + %f*(",xmax,0.999*xmax) + var1_JecDown + Form(">%f)",xmax);


  for(unsigned i_MC=0;i_MC<MC_sample_name.size();i_MC++){


    TH1F* h_MC_up = single_plot(filename_MC[i_MC],
				tree,
				var_JecUp,
				"genWeight*puWeight*("+MC_sample_cut[i_MC] + " && " + cut_cat_JecUp + " && Sum$(" + jet_sel_JecUp + ")>0)",
				nbin,xmin,xmax);

    h_MC_up->Scale(luminosity*XS[i_MC]/norm_MC[i_MC]);
    h_MC_up->SetNameTitle("x_"+MC_sample_name[i_MC]+"_jUp","x_"+MC_sample_name[i_MC]+"_jUp");
    makeBinContentsPositive(h_MC_up,true);
    h_MC_up->Write();
  
    TH1F* h_MC_down = single_plot(filename_MC[i_MC],
				  tree,
				  var_JecDown,
				  "genWeight*puWeight*("+MC_sample_cut[i_MC] + " && " + cut_cat_JecDown + " && Sum$(" + jet_sel_JecDown + ")>0)",
				  nbin,xmin,xmax);

    h_MC_down->Scale(luminosity*XS[i_MC]/norm_MC[i_MC]);
    h_MC_down->SetNameTitle("x_"+MC_sample_name[i_MC]+"_jDown","x_"+MC_sample_name[i_MC]+"_jDown");
    makeBinContentsPositive(h_MC_down,true);
    h_MC_down->Write();

    }


  cout<<"OK JEC"<<endl;

  TH1F* h_clos_diff = (TH1F*)f_new->Get("x_ttbar_fake")->Clone();
  h_clos_diff->Add( (TH1F*)f_new->Get("x_ttbar"), -1 );

  TH1F* h_fakes_clos_Up = (TH1F*)h_fakes->Clone();
  h_fakes_clos_Up->Add(h_clos_diff,+1);
  h_fakes_clos_Up->SetNameTitle("x_data_fakes_Clos_shapeUp","x_data_fakes_Clos_shapeUp");
  makeBinContentsPositive(h_fakes_clos_Up);
  h_fakes_clos_Up->Write();

  TH1F* h_fakes_clos_Down = (TH1F*)h_fakes->Clone();
  h_fakes_clos_Down->Add(h_clos_diff,-1);
  h_fakes_clos_Down->SetNameTitle("x_data_fakes_Clos_shapeDown","x_data_fakes_Clos_shapeDown");
  makeBinContentsPositive(h_fakes_clos_Down);
  h_fakes_clos_Down->Write();


  cout<<"OK Closure"<<endl;


  //Fakes systematics

  for(unsigned int i_f=0;i_f<fake_syst_names.size();i_f++){

    for(unsigned i_MC=0;i_MC<MC_sample_name.size();i_MC++){

      if(MC_sample_name[i_MC].Contains("T1qqqq")) continue;

      TString weight = "EventWeight_FR_"+fake_syst_names[i_f]+"/(1-EventWeight_FR_"+fake_syst_names[i_f]+")*genWeight*puWeight*("+MC_sample_cut[i_MC] + " && " + cut_cat + "&& Sum$(" + jet_sel + ")==0)";

      TH1F* h_MC_FR = single_plot(filename_MC[i_MC],
				  tree_FR,
				  var,
				  weight,
				  nbin,xmin,xmax);
      
      h_MC_FR->Scale(luminosity*XS[i_MC]/norm_MC[i_MC]);
      h_MC_FR->SetNameTitle("x_"+MC_sample_name[i_MC]+"_fake_"+fake_syst_names[i_f],"x_"+MC_sample_name[i_MC]+"_fake_"+fake_syst_names[i_f]);
      makeBinContentsPositive(h_MC_FR,true);
      h_MC_FR->Write();

    }

  }

  cout<<"OK Fake syst"<<endl;


  f_new->Close();
  return;

}




void test(){

  TString var = "MET_pt";
  TString cut_cat = "nMuon_presel==0 && nElectron_presel==0 && nJet_sel>2";
  TString jet_sel = "Jet_isSel && Jet_BDT_sig1000_b>0.1";
  
  bool syst = true;
  TString cut_cat_JecUp = "nMuon_presel==0 && nElectron_presel==0 && nJetUp_sel>2";
  TString jet_sel_JecUp = "JetUp_isSel && Jet_BDT_sig1000_b>0.1";
  TString cut_cat_JecDown = "nMuon_presel==0 && nElectron_presel==0 && nJetDown_sel>2";
  TString jet_sel_JecDown = "JetDown_isSel && Jet_BDT_sig1000_b>0.1";

  datacard_maker(var, 25, 0, 500,
		 cut_cat, jet_sel, 
		 "datacard_test.root", 
		 syst,
		 var, cut_cat_JecUp, jet_sel_JecUp,
		 var, cut_cat_JecDown, jet_sel_JecDown		 
		 );

}
