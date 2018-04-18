#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TH3F.h>
#include <iostream>
#include <TLegend.h>
#include <TBranch.h>
#include <TClonesArray.h>
#include <TChain.h>
#include <TMath.h>

#include <TSystem.h>
#include <TROOT.h>
using namespace std;

//gSystem->Load("/home/llr/cms/strebler/MG5_aMC_v2_2_1/ExRootAnalysis/libExRootAnalysis.so") before compilation in Root


float binning(TH1F* h){

  float max=h->GetBinLowEdge(h->GetNbinsX()+1);
  float min=h->GetBinLowEdge(1);
  float bin=(max-min)/h->GetNbinsX();
  return bin;
}


TH1F* single_plot(TString file, TString tree_name, TString var, TString cut, int nbin, float min, float max){

  TChain * tree = new TChain(tree_name);
  tree->Add(file);

  TH1F* g=new TH1F("g","g",nbin,min,max);
  g->Sumw2();
  //g->SetBinErrorOption(TH1::kPoisson);
 
  tree->Draw(var+">>g",cut,"goff");

  delete tree;
  return (TH1F*)g->Clone();

}


TH1F* single_plot(vector<TString> file, vector<TString> tree_name, TString var, TString cut, int nbin, float min, float max){

  TChain * tree = new TChain(tree_name[0]);
  tree->Add(file[0]);
  if(tree_name.size()>1){
    for(unsigned int i=1;i<tree_name.size();i++){
      TChain* tree2 = new TChain(tree_name[i]);
      tree2->Add(file[i]);
      tree->AddFriend(tree2);
    }
  }


  TH1F* g=new TH1F("g","g",nbin,min,max);
  g->Sumw2();
  //g->SetBinErrorOption(TH1::kPoisson);
 
  tree->Draw(var+">>g",cut,"goff");

  delete tree;
  return (TH1F*)g->Clone();

}


void single_plot(TH1F* g, TString file, TString tree_name, TString var, TString cut){

  TChain * tree = new TChain(tree_name);
  tree->Add(file);
  tree->Draw(var+">>"+g->GetName(),cut,"goff");

  delete tree;
  return;

}






TH1F* single_plot(vector<TString> files, TString tree_name, TString var, TString cut, int nbin, float min, float max){

  TChain * tree = new TChain(tree_name);
  for(unsigned int i=0; i<files.size(); i++)
    tree->Add(files[i]);

  TH1F* g=new TH1F("g","g",nbin,min,max);
  g=new TH1F("g","g",nbin,min,max);

  //g->Sumw2();
  g->SetBinErrorOption(TH1::kPoisson);

  tree->Draw(var+">>g",cut,"goff");
  delete tree;

//cout<<"tree deleted"<<endl;
  return g;

}




TH1F* single_plot(vector<TString> files, TString tree_name, TString var, vector<TString> cut, int nbin, float min, float max){


  TH1F* g=new TH1F("g","g",nbin,min,max);
  //g->Sumw2();
  g->SetBinErrorOption(TH1::kPoisson);

  for(unsigned int i=0; i<files.size(); i++){
    g->Add(single_plot(files[i],tree_name,var,cut[i],nbin,min,max));
  }

  return g;

}



TH1D* single_plot_d(vector<TString> files, TString tree_name, TString var, TString cut, int nbin, float min, float max){

  TChain * tree = new TChain(tree_name);
  for(unsigned int i=0; i<files.size(); i++)
    tree->Add(files[i]);

  TH1D* g=new TH1D("g","g",nbin,min,max);
  //g->Sumw2();
  g->SetBinErrorOption(TH1::kPoisson);

  tree->Draw(var+">>g",cut,"goff");

  return g;

}




TH1F* single_plot(TString file, TString tree_name, TString var, TString cut, int nbin, double* x){

  TChain * tree = new TChain(tree_name);
  tree->Add(file);

  TH1F* g=new TH1F("g","g",nbin,x);
  //g->Sumw2();
  g->SetBinErrorOption(TH1::kPoisson);

  tree->Draw(var+">>g",cut,"goff");
  delete tree;



  return g;

}




TH1F* single_plot(vector<TString> file, vector<TString> tree_name, TString var, TString cut, int nbin, double* x){

  TChain * tree = new TChain(tree_name[0]);
  tree->Add(file[0]);
  if(tree_name.size()>1){
    for(unsigned int i=1;i<tree_name.size();i++){
      TChain* tree2 = new TChain(tree_name[i]);
      tree2->Add(file[i]);
      tree->AddFriend(tree2);
    }
  }


  TH1F* g=new TH1F("g","g",nbin,x);
  //g->Sumw2();
  g->SetBinErrorOption(TH1::kPoisson);

  tree->Draw(var+">>g",cut,"goff");
  delete tree;



  return g;

}





TH1F* single_plot(vector<TString> files, TString tree_name, TString var, TString cut, int nbin, double* x){

  TChain * tree = new TChain(tree_name);
  for(unsigned int i=0;i<files.size();i++)
    tree->Add(files[i]);

  TH1F* g=new TH1F("g","g",nbin,x);
  //g->Sumw2();
  g->SetBinErrorOption(TH1::kPoisson);

  tree->Draw(var+">>g",cut,"goff");


  return g;

}







TH2F* single_plot2D(TString file, TString tree_name, TString var1, TString var2, TString cut, int nbin1, float min1, float max1, int nbin2, float min2, float max2){

  TChain * tree = new TChain(tree_name);
  tree->Add(file);
  cout<<tree->GetEntries()<<" entries"<<endl;

  tree->Draw(var2+":"+var1+Form(">>h(%i,%f,%f,%i,%f,%f)",nbin1,min1,max1,nbin2,min2,max2),cut,"goff");
  TH2F* g=(TH2F*) ((TH2F*)gDirectory->Get("h"))->Clone();
  return g;
}





TH2F* single_plot2D(vector<TString> files, TString tree_name, TString var1, TString var2, TString cut, int nbin1, float min1, float max1, int nbin2, float min2, float max2){

  TChain * tree = new TChain(tree_name);
  for(unsigned int i=0; i<files.size(); i++){
    tree->Add(files[i]);
  }

  tree->Draw(var2+":"+var1+Form(">>h(%i,%f,%f,%i,%f,%f)",nbin1,min1,max1,nbin2,min2,max2),cut,"goff");
  TH2F* g=(TH2F*) ((TH2F*)gDirectory->Get("h"))->Clone();
  return g;
}






TH2F* single_plot2D(TString file, TString tree_name, TString var1, TString var2, TString cut, int nbinx, double* x, int nbiny, double* y){

  TChain * tree = new TChain(tree_name);
  tree->Add(file);

  TH2F* g=new TH2F("g","g",nbinx,x,nbiny,y);
  g->Sumw2();

  tree->Draw(var2+":"+var1+">>g",cut,"goff");

  return g;

}




TH2F* single_plot2D(vector<TString> file, vector<TString> tree_name, TString var1, TString var2, TString cut, int nbinx, double xmin, double xmax, int nbiny, double ymin, double ymax){

  TChain * tree = new TChain(tree_name[0]);
  tree->Add(file[0]);
  if(file.size()>1){
    for(unsigned int i=1;i<file.size();i++){
      TChain* tree2 = new TChain(tree_name[i]);
      tree2->Add(file[i]);
      tree->AddFriend(tree2);
    }
  }

  TH2F* g=new TH2F("g","g",nbinx,xmin,xmax,nbiny,ymin,ymax);
  g->Sumw2();

  tree->Draw(var2+":"+var1+">>g",cut,"goff");

  return g;

}






TH2F* single_plot2D(vector<TString> file, vector<TString> tree_name, TString var1, TString var2, TString cut, int nbinx, double* x, int nbiny, double* y){

  TChain * tree = new TChain(tree_name[0]);
  tree->Add(file[0]);
  if(file.size()>1){
    for(unsigned int i=1;i<file.size();i++){
      TChain* tree2 = new TChain(tree_name[i]);
      tree2->Add(file[i]);
      tree->AddFriend(tree2);
    }
  }
 

  TH2F* g=new TH2F("g","g",nbinx,x,nbiny,y);
  g->Sumw2();

  tree->Draw(var2+":"+var1+">>g",cut,"goff");

  return g;

}






TH3F* single_plot3D(TString file, TString tree_name, TString var1, TString var2,TString var3, TString cut, int nbin1, float min1, float max1, int nbin2, float min2, float max2, int nbin3, float min3, float max3){

  TChain * tree = new TChain(tree_name);
  tree->Add(file);

  tree->Draw(var3+":"+var2+":"+var1+Form(">>h(%i,%f,%f,%i,%f,%f,%i,%f,%f)",nbin1,min1,max1,nbin2,min2,max2,nbin3,min3,max3),cut,"goff");
  TH3F* g=(TH3F*) ((TH3F*)gDirectory->Get("h"))->Clone();
  return g;
}






TH3F* single_plot3D(vector<TString> file, vector<TString> tree_name, TString var1, TString var2, TString var3, TString cut, int nbinx, double* x, int nbiny, double* y, int nbinz, double* z){

  TChain * tree = new TChain(tree_name[0]);
  tree->Add(file[0]);
  if(file.size()>1){
    for(unsigned int i=1;i<file.size();i++)
      tree->AddFriend(tree_name[i],file[i]);
  } 

  TH3F* g=new TH3F("g","g",nbinx,x,nbiny,y,nbinz,z);
  g->Sumw2();

  tree->Draw(var3+":"+var2+":"+var1+">>g",cut,"goff");

  return g;

}






std::vector<TH1F*> sort_histo(std::vector<TH1F*> h){
  
  float max=0;
  unsigned int imax=0;
  for(unsigned int i=0;i<h.size();i++){
    float m=h[i]->GetMaximum();
    if(m>max){
      max=m;
      imax=i;
    }
  }
  vector<TH1F*> new_h;
  new_h.push_back(h[imax]);
  for(unsigned int i=0;i<h.size();i++){
    if(i!=imax)
      new_h.push_back(h[i]);
  }

  return new_h;
}







TH1F* single_plot(TString hname, vector<TString> files, TString tree_name, TString var, TString cut, int nbin, float min, float max){

  TChain * tree = new TChain(tree_name);
  for(unsigned int i=0; i<files.size(); i++)
    tree->Add(files[i]);

  TH1F* g=new TH1F(hname,hname,nbin,min,max);

  //g->Sumw2();
  g->SetBinErrorOption(TH1::kPoisson);

  tree->Draw(var+">>"+hname,cut,"goff");
  return g;

}






void makeBinContentsPositive(TH1* histogram, bool verbosity=false)
{
  if ( verbosity ) {
    std::cout << "<makeBinContentsPositive>:" << std::endl;
    std::cout << " integral(" << histogram->GetName() << ") = " << histogram->Integral() << std::endl;
  }
  double integral_original = histogram->Integral();
  if ( integral_original < 0. ) integral_original = 0.;
  /*if ( verbosity ) {
    std::cout << " integral_original = " << integral_original << std::endl;
    }*/
  int numBins = histogram->GetNbinsX();
  for ( int iBin = 0; iBin <= (numBins + 1); ++iBin ) {
    double binContent_original = histogram->GetBinContent(iBin);
    double binError2_original = pow(histogram->GetBinError(iBin),2);
    if ( binContent_original < 0. ) {
      double binContent_modified = 0.;
      double binError2_modified = binError2_original + pow(binContent_original - binContent_modified,2);
      assert(binError2_modified >= 0.);
      if ( verbosity ) {
	std::cout << "bin #" << iBin << " (x =  " << histogram->GetBinCenter(iBin) << "): binContent = " << binContent_original << " +/- " << TMath::Sqrt(binError2_original)
          << " --> setting it to binContent = " << binContent_modified << " +/- " << TMath::Sqrt(binError2_modified) << std::endl;
      }
      histogram->SetBinContent(iBin, binContent_modified);
      histogram->SetBinError(iBin, TMath::Sqrt(binError2_modified));
    }
  }
  double integral_modified = histogram->Integral();
  if ( integral_modified < 0. ) integral_modified = 0.;
  /*if ( verbosity ) {
    std::cout << " integral_modified = " << integral_modified << std::endl;
    }*/
  if ( integral_modified > 0. ) {
    double sf = integral_original/integral_modified;
    /*if ( verbosity ) {
      std::cout << "--> scaling histogram by factor = " << sf << std::endl;
      }*/
    histogram->Scale(sf);
  } else {
    for ( int iBin = 0; iBin <= (numBins + 1); ++iBin ) {
      histogram->SetBinContent(iBin, 0.);
    }
  }
  /*if ( verbosity ) {
    std::cout << " integral(" << histogram->GetName() << ") = " << histogram->Integral() << std::endl;
    }*/
}




TH1F* unfold(TH2F* histo){

  int nbins = histo->GetSize();

  TH1F* h = new TH1F("h","h",nbins+1,-0.5,nbins+0.5);
  for(int i=0;i<=nbins;i++)
    h->SetBinContent(i+1,histo->GetBinContent(i));


  return h;

}
