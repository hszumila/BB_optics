#include <TSystem.h>
#include <TString.h>
#include "TFile.h"
#include "TTree.h"
#include <TNtuple.h>
#include "TCanvas.h"
#include <iostream>
#include <fstream>
#include <TMath.h>
#include "TH1F.h"
#include <TH2.h>
#include <TCutG.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TMath.h>
#include <TLegend.h>
#include <TPaveLabel.h>
#include <TProfile.h>
#include <TPolyLine.h>
#include <TObjArray.h>
#include "TVector3.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include<math.h>
#include <iostream>
#include <fstream>
using namespace std;

void make_hist_bb_optics(Int_t nrun=1813,Bool_t CutYtarFlag=kTRUE,Bool_t CutYpFpYFpFlag=kTRUE,Bool_t CutXpFpXFpFlag=kTRUE,Int_t FileID=-2){
  gStyle->SetPalette(1,0);
  gStyle->SetOptStat(1000011);
  gStyle->SetOptFit(11);
  gStyle->SetTitleOffset(1.,"Y");
  gStyle->SetTitleOffset(.7,"X");
  gStyle->SetLabelSize(0.04,"XY");
  gStyle->SetTitleSize(0.06,"XY");
  gStyle->SetPadLeftMargin(0.12);
  //  Get info for that optics run
  TString OpticsFile = "list_of_optics_run.dat";
  ifstream file_optics(OpticsFile.Data());
  TString opticsline;
  TString OpticsID="";
  Int_t RunNum=0.;
  Double_t CentAngle=0.;
  Int_t SieveFlag=1;
  Int_t NumFoil=0;
  TString temp;
  //
  vector <Double_t> ztar_foil;
  Int_t ndelcut=-1;
  vector<Double_t > delcut;
  vector<Double_t > delwidth;
  if (file_optics.is_open()) {
    //
    cout << " Open file = " << OpticsFile << endl;
    while (RunNum!=nrun  ) {
      temp.ReadToDelim(file_optics,',');
      cout << temp << endl;
      if (temp.Atoi() == nrun) {
	RunNum = temp.Atoi();
      } else {
	temp.ReadLine(file_optics);
      }
    }
    if (RunNum==nrun) {
      temp.ReadToDelim(file_optics,',');
      OpticsID = temp;
      temp.ReadToDelim(file_optics,',');
      CentAngle = temp.Atof();
      temp.ReadToDelim(file_optics,',');
      NumFoil = temp.Atoi();
      temp.ReadToDelim(file_optics,',');
      SieveFlag = temp.Atoi();
      temp.ReadToDelim(file_optics);
      ndelcut = temp.Atoi();
      for (Int_t nf=0;nf<NumFoil-1;nf++) {
        temp.ReadToDelim(file_optics,',');
	ztar_foil.push_back(temp.Atof());
      }
      temp.ReadToDelim(file_optics);
      ztar_foil.push_back(temp.Atof());
      for (Int_t nd=0;nd<ndelcut-1;nd++) {
        temp.ReadToDelim(file_optics,',');
      	delcut.push_back(temp.Atof());
      }
      temp.ReadToDelim(file_optics);
      delcut.push_back(temp.Atof());
      for (Int_t nw=0;nw<ndelcut-1;nw++) {
	temp.ReadToDelim(file_optics,',');
	delwidth.push_back(temp.Atof());
      }
      temp.ReadToDelim(file_optics);
      delwidth.push_back(temp.Atof());
    }
  } else {
    cout << " No file = " << OpticsFile << endl;    
  }
  cout << RunNum << " " << OpticsID << " " << CentAngle << " " << NumFoil << " " << SieveFlag << endl;
  if (NumFoil==0) return;
  //
  CentAngle *= 3.14/180.0;
  
  TString inputroot;
  TString outputhist;
  inputroot=Form("../sim/gmn_4.5GeV2_opcal_job0.root");//shms_replay_matrixopt_%s_%d.root",OpticsID.Data(),FileID);
  outputhist=Form("hist/Optics_%s_%d_hist.root",OpticsID.Data(),FileID);
  cout << " input root = " << inputroot << endl;
  TObjArray HList(0);
  //
  TString YtarDeltaCutFile;
  TFile *fYtarDeltaCut;
  vector <TCutG*> ytar_delta_cut;
  if (CutYtarFlag) {
    YtarDeltaCutFile=Form("cuts/ytar_delta_%s_%d_cut.root",OpticsID.Data(),FileID);
    fYtarDeltaCut = new TFile(YtarDeltaCutFile);
    cout << "Ytar Cut file = " << YtarDeltaCutFile << endl;
    for (Int_t nc=0;nc<NumFoil;nc++) {
      fYtarDeltaCut->cd();
      TCutG* tempcut = (TCutG*)gROOT->FindObject(Form("delta_vs_ytar_cut_foil%d",nc));
      if (tempcut) {
	Int_t npt = tempcut->GetN();
	cout << "hYtarDelta_cut = " << nc << " npts = " << npt << endl;
	ytar_delta_cut.push_back(tempcut);
      } else {
	cout << " No hYtarDelta_cut = " << nc << endl;
      }
    }
  }
  //
  TString outCutFile;
  TFile *fcut;
  vector<vector<vector<TCutG*> > > ypfp_yfp_cut;
  vector<vector<vector<Int_t> > > ypfp_yfp_cut_flag;
  ypfp_yfp_cut.resize(NumFoil);
  ypfp_yfp_cut_flag.resize(NumFoil);
  for  (Int_t nf=0;nf<NumFoil;nf++) {
    ypfp_yfp_cut[nf].resize(ndelcut);
    ypfp_yfp_cut_flag[nf].resize(ndelcut);
  }
  if (CutYpFpYFpFlag) {
    outCutFile=Form("cuts/YpFpYFp_%s_%d_cut.root",OpticsID.Data(),FileID);
    fcut = new TFile(outCutFile);
    cout << " Cut file = " << outCutFile << endl;
    fcut->cd();
    for  (Int_t nf=0;nf<NumFoil;nf++) {
      for  (Int_t nd=0;nd<ndelcut;nd++) {
        for (Int_t nc=0;nc<7;nc++) {
	  TCutG* tempg  = (TCutG*)gROOT->FindObject(Form("hYpFpYFp_cut_yscol_%d_nfoil_%d_ndel_%d",nc,nf,nd));
	  if (tempg)  {
	    //cout << "hYpFpYFp_cut = " << nc << " " << nf << " " << nd << endl;
	    ypfp_yfp_cut[nf][nd].push_back(tempg);
	  } else {
	    //cout << " No hYpFpYFp_cut = " << nc << " " << nf << " " << nd << endl;
	    ypfp_yfp_cut[nf][nd].push_back(tempg);
	  }
	}}}
  }
  //
//
  TString xpfp_xfp_outCutFile;
  TFile *xpfp_xfp_fcut;
  vector<vector<vector<TCutG*> > > xpfp_xfp_cut;
  vector<vector<vector<Int_t> > > xpfp_xfp_cut_flag;
  xpfp_xfp_cut.resize(NumFoil);
  xpfp_xfp_cut_flag.resize(NumFoil);
  for  (Int_t nf=0;nf<NumFoil;nf++) {
    xpfp_xfp_cut[nf].resize(ndelcut);
    xpfp_xfp_cut_flag[nf].resize(ndelcut);
  }
  if (CutXpFpXFpFlag) {
    xpfp_xfp_outCutFile=Form("cuts/XpFpXFp_%s_%d_cut.root",OpticsID.Data(),FileID);
    xpfp_xfp_fcut = new TFile(xpfp_xfp_outCutFile);
    cout << "xpfp_xfp_ Cut file = " << xpfp_xfp_outCutFile << endl;
    xpfp_xfp_fcut->cd();
    for  (Int_t nf=0;nf<NumFoil;nf++) {
      for  (Int_t nd=0;nd<ndelcut;nd++) {
        for (Int_t nc=0;nc<13;nc++) {
	  TCutG* tempg  = (TCutG*)gROOT->FindObject(Form("hXpFpXFp_cut_yscol_%d_nfoil_%d_ndel_%d",nc,nf,nd));
	  if (tempg)  {
	    //cout << "hXpFpXFp_cut = " << nc << " " << nf << " " << nd << endl;
	    xpfp_xfp_cut[nf][nd].push_back(tempg);
	  } else {
	    //cout << " No hXpFpXFp_cut = " << nc << " " << nf << " " << nd << endl;
	    xpfp_xfp_cut[nf][nd].push_back(tempg);
	  }
	}}}
  }
  //
  TFile *fsimc = new TFile(inputroot); 
  TTree *tsimc = (TTree*) fsimc->Get("T");
  // Define branches
  
  vector<double> *yfp = 0;
  tsimc->SetBranchAddress("Earm.BBGEM.Track.Y",&yfp);
  vector<double>  *ypfp = 0;
  tsimc->SetBranchAddress("Earm.BBGEM.Track.Yp",&ypfp);
  vector<double>  *xfp = 0;
  tsimc->SetBranchAddress("Earm.BBGEM.Track.X",&xfp);
  vector<double>  *xpfp = 0;
  tsimc->SetBranchAddress("Earm.BBGEM.Track.Xp",&xpfp);
  vector<double>  *yfp_fit = 0;
  tsimc->SetBranchAddress("Earm.BBGEM.Track.Yfit",&yfp_fit);
  vector<double>  *ypfp_fit = 0;
  tsimc->SetBranchAddress("Earm.BBGEM.Track.Ypfit",&ypfp_fit);
  vector<double>  *xfp_fit = 0;
  tsimc->SetBranchAddress("Earm.BBGEM.Track.Xfit",&xfp_fit);
  vector<double>  *xpfp_fit = 0;
  tsimc->SetBranchAddress("Earm.BBGEM.Track.Xpfit",&xpfp_fit);
  int ntracks = 0;
  tsimc->SetBranchAddress("Earm.BBGEM.Track.ntracks",&ntracks);
  vector<int> *nhits= 0;
  tsimc->SetBranchAddress("Earm.BBGEM.Track.NumHits",&nhits);
  vector<double> *ndf = 0;
  tsimc->SetBranchAddress("Earm.BBGEM.Track.NDF",&ndf);
  vector<double> *chisq = 0;
  tsimc->SetBranchAddress("Earm.BBGEM.Track.Chi2fit",&chisq);

  
  //define the variables
  double vx, vy, vz, px, py, pz;
  double p, xptar, yptar, ytar, xtar;
  double p_fit, xptar_fit, yptar_fit, ytar_fit; //Fit is reconstructed using fit coefficients, no smearing for detector resolution
  double p_recon, xptar_recon, yptar_recon, ytar_recon; //recon is reconstructed using fit coefficients, fp quantities smeared by det. resolution
  double pthetabend_true;
  double pthetabend_fit, pthetabend_recon;
  double pinv_fit, pinv_recon;
  double vz_fit, vz_recon;
  double thetabend_true;
  double thetabend_fit;
  double thetabend_recon;
  double xtar_recon, xtar_fit;
  double xsieve, ysieve;
  double z0 = 1.172;//distance to face of sieve,[m]?
  TVector3 spec_xaxis_fp,spec_yaxis_fp, spec_zaxis_fp;
  TVector3 spec_xaxis_tgt,spec_yaxis_tgt, spec_zaxis_tgt;
  double tracker_pitch_angle = 10.0*3.14/180.0;//put this into the input file

  TVector3 BB_zaxis( sin(CentAngle), 0.0, cos(CentAngle) ); //BB is on beam right, global x axis points to beam left
  TVector3 BB_xaxis(0,-1,0); //X axis of transport coordinates is vertically down:
  TVector3 BB_yaxis = (BB_zaxis.Cross(BB_xaxis)).Unit();
  
  spec_xaxis_tgt = BB_xaxis;
  spec_yaxis_tgt = BB_yaxis;
  spec_zaxis_tgt = BB_zaxis;
  
  spec_zaxis_fp = BB_zaxis;
  spec_yaxis_fp = BB_yaxis;
  spec_zaxis_fp.Rotate(-tracker_pitch_angle, spec_yaxis_fp);
  
  // Define histograms
  TH1F *hytar = new TH1F("hytar",Form("Run %d ; Ytar; Counts",nrun),500,-0.2,0.2);
  HList.Add(hytar);
  TH1F *hztar = new TH1F("hztar",Form("Run %d ; Ztar; Counts",nrun),500,-0.2,0.2);
  HList.Add(hztar);
  TH2F *hXptarDelta = new TH2F("hXptarDelta",Form("Run %d ; Xptar ; Pinv x #theta_{bend}",nrun),120,-.6,.6,100,0,0.25);
  HList.Add(hXptarDelta);
  TH2F *hYptarDelta = new TH2F("hYptarDelta",Form("Run %d ; Yptar ; Pinv x #theta_{bend}",nrun),120,-.4,.4,100,0,0.25);
  HList.Add(hYptarDelta);
  TH2F *hYtarDelta = new TH2F("hYtarDelta",Form("Run %d ; Ytar ; Pinv x #theta_{bend}",nrun),100,-0.3,0.3,100,0,0.25);
  HList.Add(hYtarDelta);
  //
  TH2F *hYpFpYFp_all = new TH2F("hYpFpYFp_all",Form("Run %d ; Ypfp ; Yfp",nrun),100,-.3,.3,100,-0.3,0.3);
  HList.Add(hYpFpYFp_all);
  TH2F *hYFpXFp = new TH2F("hYFpXFp",Form("Run %d ; Yfp ; Xfp",nrun),100,-0.3,0.3,100,-0.7,0.7);
  HList.Add(hYFpXFp);
  TH2F *hXpFpXFp = new TH2F("hXpFpXFp",Form("Run %d ; Xpfp ; Xfp",nrun),100,-.7,.7,100,-0.7,0.7);
  HList.Add(hXpFpXFp);
  TH2F *hYtarYptar = new TH2F("hYtarYptar",Form("Run %d ; Yptar ; Ytar",nrun),100,-.3,.3,100,-0.15,0.15);
  HList.Add(hYtarYptar);
  TH2F *hZtarDelta = new TH2F("hZtarDelta",Form("Run %d ; Ztar ; Pinv x #theta_{bend}",nrun),100,-0.3,0.3,100,0,0.25);
  HList.Add(hZtarDelta);

  TH1F *h_p = new TH1F("h_p",Form("Run %d ; P, recon",nrun),100,0,10);
  HList.Add(h_p);//p
  TH1F *h_pinvtheta = new TH1F("h_pinvtheta",Form("Run %d ; Pinv x #theta_{bend}",nrun),100,0,0.5);
  HList.Add(h_pinvtheta);//ptheta
  TH2F *hPinvthetaVx = new TH2F("hPinvthetaVx",Form("Run %d ; Pinv x #theta_{bend}; xfp",nrun),100,0,0.25,100,-0.7,0.7);
  HList.Add(hPinvthetaVx);//ptheta vs x
  TH2F *hPinvthetaVxtar = new TH2F("hPinvthetaVxtar",Form("Run %d ; Pinv x #theta_{bend}; xTar",nrun),100,0,0.2,100,-0.05,0.05);
  HList.Add(hPinvthetaVxtar);//ptheta vs x
   TH1F *h_ptheta = new TH1F("h_ptheta",Form("Run %d ; P x #theta_{bend}",nrun),100,0,4);
  HList.Add(h_ptheta);//ptheta
  TH2F *hPthetaVx = new TH2F("hPthetaVx",Form("Run %d ; P x #theta_{bend}; xfp",nrun),100,0,4,100,-0.7,0.7);
  HList.Add(hPthetaVx);//ptheta vs x
  TH2F *hPthetaVxtar = new TH2F("hPthetaVxtar",Form("Run %d ; P x #theta_{bend}; xTar",nrun),100,0,4,100,-0.05,0.05);
  HList.Add(hPthetaVxtar);//ptheta vs x
  TH2F *hthetaVp = new TH2F("hthetaVp",Form("Run %d ; #theta_{bend}; P [GeV/c]",nrun),100,0,1,100,0,10);
  HList.Add(hthetaVp);
  TH2F *h_xsVys = new TH2F("h_xsVys",Form("Run %d ; y_{sieve}; x_{sieve}",nrun),100,-0.2,0.2,100,-0.4,0.4);
  HList.Add(h_xsVys);
  TH1F *h_zmsc = new TH1F("h_zmsc",Form("Run %d ; z_{fit} - z_{true}",nrun),100,-0.001,0.001);
  HList.Add(h_zmsc);
  TH2F *h_zmscVz = new TH2F("h_zmscVz",Form("Run %d ; z_{fit}; z_{fit} - z_{true}",nrun),100,-0.2,0.2,100,-0.001,0.001);
  HList.Add(h_zmscVz);
   TH2F *hYpFpYFp_cut0 = new TH2F("hYpFpYFp_cut0",Form("Run %d, yS=2 ; Ypfp ; Yfp",nrun),100,-.3,.3,100,-0.3,0.3);
  HList.Add(hYpFpYFp_cut0);
  TH2F *hXpFpXFp_cut0 = new TH2F("hXpFpXFp_cut0",Form("Run %d, xS=2 ; Xpfp ; Xfp",nrun),100,-.7,.7,100,-0.7,0.7);
  HList.Add(hXpFpXFp_cut0);
  
  
  //
  vector <TH2F*> hYsDelta;
  hYsDelta.resize(NumFoil);
  vector <TH2F*> hXsDelta;
  hXsDelta.resize(NumFoil);
  vector <TH2F*> hYpFpYFp;
  hYpFpYFp.resize(NumFoil);
  vector<vector<vector<TH2F*> > > hYsXs_DelCut_YpYfpCut;
  vector<vector<vector<TH2F*> > > hYsXs_DelCut_XpXfpCut;
  vector<vector<vector<TH1F*> > > hXs_DelCut_YpYfpCut;
  vector<vector<TH2F*> > hYsXs_DelCut;
  vector<vector<TH2F*> > hYpFpYFp_DelCut;
  vector<vector<TH2F*> > hXpFpXFp_DelCut;
  cout << " setup DelCut 2d" << endl;
  hYsXs_DelCut.resize(NumFoil);
  hYsXs_DelCut_YpYfpCut.resize(NumFoil);
  hYsXs_DelCut_XpXfpCut.resize(NumFoil);
  hXs_DelCut_YpYfpCut.resize(NumFoil);
  hYpFpYFp_DelCut.resize(NumFoil);
  hXpFpXFp_DelCut.resize(NumFoil);
  for  (Int_t nf=0;nf<NumFoil;nf++) {
    hYsXs_DelCut[nf].resize(ndelcut);
    hYsXs_DelCut_YpYfpCut[nf].resize(ndelcut);
    hYsXs_DelCut_XpXfpCut[nf].resize(ndelcut);
    hXs_DelCut_YpYfpCut[nf].resize(ndelcut);
    hYpFpYFp_DelCut[nf].resize(ndelcut);
    hXpFpXFp_DelCut[nf].resize(ndelcut);
    for  (Int_t nd=0;nd<ndelcut;nd++) {
      hYsXs_DelCut_YpYfpCut[nf][nd].resize(13);
      hYsXs_DelCut_XpXfpCut[nf][nd].resize(13);
      hXs_DelCut_YpYfpCut[nf][nd].resize(13);
    }
  }
  cout << " finish setup Cut 2d" << endl;
  for  (Int_t nc=0;nc<NumFoil;nc++) {
    hYsDelta[nc] = new TH2F(Form("hYsDelta_Foil_%d",nc),Form("Run %d Foil %d; Ys ; Pinv x #theta_{bend}",nc,nrun),100,-0.2,0.2,50,0.,0.25);
    HList.Add(hYsDelta[nc]);
    hXsDelta[nc] = new TH2F(Form("hXsDelta_Foil_%d",nc),Form("Run %d Foil %d; Xs ; Pinv x #theta_{bend}",nc,nrun),100,-0.4,0.4,50,0,0.25);
    HList.Add(hXsDelta[nc]);
    hYpFpYFp[nc] = new TH2F(Form("hYpFpYFp_%d",nc),Form("Run %d Foil %d; Ypfp ; Yfp",nrun,nc),100,-.3,.3,100,-.3,.3);
    HList.Add(hYpFpYFp[nc]);
    for  (Int_t nd=0;nd<ndelcut;nd++) {
      hYsXs_DelCut[nc][nd]  = new TH2F(Form("hYsXs_Foil_%d_DelCut_%d",nc,nd),Form("Run %d Foil %d Cut %3.1f; Ys ; Xs",nrun,nc,delcut[nd]),50,-0.2,0.2,100,-0.4,0.4);
      HList.Add(hYsXs_DelCut[nc][nd]);
      for  (Int_t ny=0;ny<13;ny++) {
	hYsXs_DelCut_YpYfpCut[nc][nd][ny]  = new TH2F(Form("hYsXs_Foil_%d_DelCut_%d_FpCut_%d",nc,nd,ny),Form("Run %d Foil %d Cut %3.1f Ys=%d; Ys ; Xs",nrun,nc,delcut[nd],ny),100,-0.2,0.2,100,-0.4,0.4);
	HList.Add(hYsXs_DelCut_YpYfpCut[nc][nd][ny]);
	hYsXs_DelCut_XpXfpCut[nc][nd][ny]  = new TH2F(Form("hYsXs_Foil_%d_DelCut_%d_XFpCut_%d",nc,nd,ny),Form("Run %d Foil %d Cut %3.1f Xs=%d; Ys ; Xs",nrun,nc,delcut[nd],ny),100,-0.2,0.2,100,-0.4,0.4);
	HList.Add(hYsXs_DelCut_XpXfpCut[nc][nd][ny]);
	hXs_DelCut_YpYfpCut[nc][nd][ny]  = new TH1F(Form("hXs_Foil_%d_DelCut_%d_FpCut_%d",nc,nd,ny),Form("Run %d Foil %d Cut %3.1f Ys=%d; Xs",nrun,nc,delcut[nd],ny),100,-0.4,0.4);
	HList.Add(hXs_DelCut_YpYfpCut[nc][nd][ny]);
      }
      hYpFpYFp_DelCut[nc][nd]  = new TH2F(Form("hYpFpYFp_%d_DelCut_%d",nc,nd),Form("Run %d Foil %d Cut %3.1f; Ypfp ; Yfp",nrun,nc,delcut[nd]),75,-.3,.3,150,-0.3,0.3);
      HList.Add(hYpFpYFp_DelCut[nc][nd]);
      hXpFpXFp_DelCut[nc][nd]= new TH2F(Form("hXpFpXFp_%d_DelCut_%d",nc,nd),Form("Run %d Foil %d Cut %3.1f; Xpfp ; Xfp",nrun,nc,delcut[nd]),150,-0.7,0.7,150,-0.7,0.7);
      HList.Add(hXpFpXFp_DelCut[nc][nd]);
    }
  }	  
  

  //reading the model file and storing the data in a matrix, M
  //string  modelfilename = "BB_optics_fit_GMN13.txt";
  string  modelfilename = "newfit.dat";
  ifstream modelfile(modelfilename.c_str());
  TString currentline;
  //while( currentline.ReadLine(inputfile) ){}
  
  int row_M = 0, col_M = 9;
  modelfile >> row_M;
  TMatrixD M(row_M,col_M);
  for(int row=0; row<row_M; row++){
    for(int col=0; col<col_M; col++){ 
      modelfile >> M(row,col);
      cout<<M(row,col)<<" ";
    }
    cout<<endl;
  }

  // loop over entries
  Long64_t nentries = tsimc->GetEntries();
  cout << " start loop " << nentries << endl;
  for (int i = 0; i < nentries; i++) {
    tsimc->GetEntry(i);
    if (i%500==0)cout << " Entry = " << i << endl;

    //determine if good track
    bool goodtrack = false;

    if (ntracks==1 && nhits->at(0)==5 && xfp->at(0)<0.55 && xfp->at(0)>-0.55 && chisq->at(0)/ndf->at(0)<=30.0){
      goodtrack=1;
    }

    vy=0;
    if (goodtrack){
      //reconstruct the target quantities
	xtar_fit = -vy;
	xtar_recon = -vy;
	
	for( int iter=0; iter<3; iter++ ){

	  xptar_fit = 0.0;
	  yptar_fit = 0.0;
	  ytar_fit = 0.0;
	  pthetabend_fit = 0.0;
	  pinv_fit = 0.0;
	
	  xptar_recon = 0.0;
	  yptar_recon = 0.0;
	  ytar_recon = 0.0;
	  pthetabend_recon = 0.0;
	  pinv_recon = 0.0;
	
	  for (int row=0; row<row_M; row++){
	    xptar_fit += M(row,0)*pow(xfp->at(0),M(row,4))*pow(yfp->at(0),M(row,5))*pow(xpfp->at(0),M(row,6))*pow(ypfp->at(0),M(row,7))*pow(xtar_recon,M(row,8));
	    yptar_fit += M(row,1)*pow(xfp->at(0),M(row,4))*pow(yfp->at(0),M(row,5))*pow(xpfp->at(0),M(row,6))*pow(ypfp->at(0),M(row,7))*pow(xtar_recon,M(row,8));
	    ytar_fit += M(row,2)*pow(xfp->at(0),M(row,4))*pow(yfp->at(0),M(row,5))*pow(xpfp->at(0),M(row,6))*pow(ypfp->at(0),M(row,7))*pow(xtar_recon,M(row,8));
	    pinv_fit += M(row,3)*pow(xfp->at(0),M(row,4))*pow(yfp->at(0),M(row,5))*pow(xpfp->at(0),M(row,6))*pow(ypfp->at(0),M(row,7))*pow(xtar_recon,M(row,8));
	    pthetabend_fit += M(row,3)*pow(xfp->at(0),M(row,4))*pow(yfp->at(0),M(row,5))*pow(xpfp->at(0),M(row,6))*pow(ypfp->at(0),M(row,7))*pow(xtar_recon,M(row,8));
	    
	    xptar_recon += M(row,0)*pow(xfp_fit->at(0),M(row,4))*pow(yfp_fit->at(0),M(row,5))*pow(xpfp_fit->at(0),M(row,6))*pow(ypfp_fit->at(0),M(row,7))*pow(xtar_fit,M(row,8));
	    yptar_recon += M(row,1)*pow(xfp_fit->at(0),M(row,4))*pow(yfp_fit->at(0),M(row,5))*pow(xpfp_fit->at(0),M(row,6))*pow(ypfp_fit->at(0),M(row,7))*pow(xtar_fit,M(row,8));
	    ytar_recon += M(row,2)*pow(xfp_fit->at(0),M(row,4))*pow(yfp_fit->at(0),M(row,5))*pow(xpfp_fit->at(0),M(row,6))*pow(ypfp_fit->at(0),M(row,7))*pow(xtar_fit,M(row,8));
	    pinv_recon += M(row,3)*pow(xfp_fit->at(0),M(row,4))*pow(yfp_fit->at(0),M(row,5))*pow(xpfp_fit->at(0),M(row,6))*pow(ypfp_fit->at(0),M(row,7))*pow(xtar_fit,M(row,8));
	    pthetabend_recon += M(row,3)*pow(xfp_fit->at(0),M(row,4))*pow(yfp_fit->at(0),M(row,5))*pow(xpfp_fit->at(0),M(row,6))*pow(ypfp_fit->at(0),M(row,7))*pow(xtar_fit,M(row,8));
	  }

	  //BB, beam left:
	  vz_fit = -ytar_fit / (sin(CentAngle) + cos(CentAngle)*yptar_fit);
	  vz_recon = -ytar_recon / (sin(CentAngle) + cos(CentAngle)*yptar_recon);
	   
	  xtar_recon = -vy - vz_recon * cos(CentAngle) * xptar_recon;
	  xtar_fit   = -vy - vz_fit * cos(CentAngle) * xptar_fit;
	}

	//calculate theta bend:
	TVector3 phat_tgt_recon(xptar_recon, yptar_recon, 1.0 );
	phat_tgt_recon = phat_tgt_recon.Unit();

	TVector3 phat_tgt_recon_global = phat_tgt_recon.X() * spec_xaxis_tgt +
	  phat_tgt_recon.Y() * spec_yaxis_tgt +
	  phat_tgt_recon.Z() * spec_zaxis_tgt;

	TVector3 phat_fp_recon(xpfp_fit->at(0), ypfp_fit->at(0), 1.0 );
	phat_fp_recon = phat_fp_recon.Unit();
	
	TVector3 phat_fp_recon_global = phat_fp_recon.X() * spec_xaxis_fp +
	  phat_fp_recon.Y() * spec_yaxis_fp +
	  phat_fp_recon.Z() * spec_zaxis_fp;

	thetabend_recon = acos( phat_fp_recon_global.Dot( phat_tgt_recon_global ) );

	int pexpansion_flag = 1;
	if( pexpansion_flag == 0 ){
	  p_recon = pthetabend_recon/thetabend_recon;
	} else {
	  p_recon = 1.0/pthetabend_recon;
	}
	pinv_recon = 1.0/p_recon;
	
	xsieve = xtar_recon + xptar_recon*z0;
	ysieve = ytar_recon + yptar_recon*z0;
	
	hytar->Fill(ytar_recon);
	hztar->Fill(vz_recon);
	hXptarDelta->Fill(xptar_recon,pinv_recon*thetabend_recon);
	hYptarDelta->Fill(yptar_recon,pinv_recon*thetabend_recon);
	hYtarDelta->Fill(ytar_recon,pinv_recon*thetabend_recon);
	hYtarYptar->Fill(yptar_recon,ytar_recon);
	hYpFpYFp_all->Fill(ypfp_fit->at(0),yfp_fit->at(0));
	hXpFpXFp->Fill(xpfp_fit->at(0),xfp_fit->at(0));
	hYFpXFp->Fill(yfp_fit->at(0),xfp_fit->at(0)); 
	hYtarYptar->Fill(yptar_recon,ytar_recon);
	hZtarDelta->Fill(vz_recon,pinv_recon*thetabend_recon);

	if (abs(xsieve+0.2017)<0.03){
	  hXpFpXFp_cut0->Fill(xpfp_fit->at(0),xfp_fit->at(0));
	}
	

	h_p->Fill(p_recon);
	h_pinvtheta->Fill(pinv_recon*thetabend_recon);
	hPinvthetaVx->Fill(pinv_recon*thetabend_recon,xfp_fit->at(0));
	hPinvthetaVxtar->Fill(pinv_recon*thetabend_recon,xtar_recon);
	
	h_ptheta->Fill(p_recon*thetabend_recon);
	hPthetaVx->Fill(p_recon*thetabend_recon,xfp_fit->at(0));
	hPthetaVxtar->Fill(p_recon*thetabend_recon,xtar_recon);
	hthetaVp->Fill(thetabend_recon,p_recon);
	h_xsVys->Fill(ysieve,xsieve);
	h_zmsc->Fill(vz_recon-vz_fit);
	h_zmscVz->Fill(vz_recon,vz_recon-vz_fit);

	for  (UInt_t nc=0;nc<ytar_delta_cut.size();nc++) {
	  if (ytar_delta_cut[nc]->IsInside(ytar_recon,pinv_recon*thetabend_recon))	{ 
	    hYsDelta[nc]->Fill(ysieve,pinv_recon*thetabend_recon);
	    hXsDelta[nc]->Fill(xsieve,pinv_recon*thetabend_recon);
	    hYpFpYFp[nc]->Fill(ypfp->at(0),yfp->at(0));
	    if (nc==0 && abs(ysieve-(2.*0.0381-0.0381*3))<0.03){
	      hYpFpYFp_cut0->Fill(ypfp_fit->at(0),yfp_fit->at(0));
	    }
	    
	    for  (Int_t nd=0;nd<ndelcut;nd++) {
	      if ( pinv_recon*thetabend_recon >=delcut[nd]-delwidth[nd] && pinv_recon*thetabend_recon <delcut[nd]+delwidth[nd]) {
		hYsXs_DelCut[nc][nd]->Fill(ysieve,xsieve); 
		hYpFpYFp_DelCut[nc][nd]->Fill(ypfp->at(0),yfp->at(0));
		hXpFpXFp_DelCut[nc][nd]->Fill(xpfp->at(0),xfp->at(0));
		Int_t f_ny=-1;
		for  (UInt_t ny=0;ny<7;ny++) {
		  if (CutYpFpYFpFlag && ypfp_yfp_cut[nc][nd][ny] && ypfp_yfp_cut[nc][nd][ny]->IsInside(ypfp->at(0),yfp->at(0))) {
		    hYsXs_DelCut_YpYfpCut[nc][nd][ny]->Fill(ysieve,xsieve);
		    hXs_DelCut_YpYfpCut[nc][nd][ny]->Fill(xsieve);
		    f_ny=ny;
		  }
		}
		for  (UInt_t nx=0;nx<13;nx++) {
		  if (f_ny != -1 && CutXpFpXFpFlag && xpfp_xfp_cut[nc][nd][nx] && xpfp_xfp_cut[nc][nd][nx]->IsInside(xpfp->at(0),xfp->at(0))) {
		    hYsXs_DelCut_XpXfpCut[nc][nd][nx]->Fill(ysieve,xsieve);		    
		  }
		}
	      }
	    }
	  }
	}	
    }//end if good track
        
  }//end loop entries

  //
  TFile hsimc(outputhist,"recreate");
  HList.Write();
  //
}
