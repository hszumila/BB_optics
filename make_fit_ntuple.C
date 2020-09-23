#include <TSystem.h>
#include <TString.h>
#include "TFile.h"
#include "TTree.h"
#include <TNtuple.h>
#include "TCanvas.h"
#include <iostream>
#include <fstream>
#include "TMath.h"
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
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include<math.h>
#include <iostream>
#include <fstream>
using namespace std;

void make_fit_ntuple(Int_t nrun=1814,Int_t FileID=-2){
  Double_t yMP = 0.0;
  Double_t xMP = 0.0;
  
  Bool_t CutYtarFlag=kTRUE;
  Bool_t CutYpFpYFpFlag=kTRUE;
  Bool_t CutXpFpXFpFlag=kTRUE;
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
  Int_t ndelcut;
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
  for (Int_t nf=0;nf<NumFoil;nf++) {
    cout << nf << " foil = " << ztar_foil[nf] << endl;
  }
  
  vector <Double_t> xs_cent{-(0.3+0.0492)+0.0493/cos(18.*3.14/180.),
      -(0.3+0.0492)+(0.0493+0.0492)/cos(18.*3.14/180),
      -(0.3+0.0492)+0.1493/cos(9.*3.14/180.),
      -(0.3+0.0492)+(0.1493+0.0492)/cos(9.*3.14/180.),
      -(0.3+0.0492)+(0.1493+0.0492*2.)/cos(9.*3.14/180.),
      -0.0492,
      0.0,
      0.0492,
      0.3+0.0492-(0.1493+0.0492*2.)/cos(9.*3.14/180.),
      0.3+0.0492-(0.1493+0.0492)/cos(9.*3.14/180.),
      0.3+0.0492-0.1493/cos(9.*3.14/180.),
      0.3+0.0492-(0.0493+0.0492)/cos(18.*3.14/180),
      0.3+0.0492-0.0493/cos(18.*3.14/180.)};
  
  vector <Double_t> ys_cent;
  for (Int_t nys=0;nys<7;nys++) {
    Double_t pos=nys*0.0381-0.0381*3;
    ys_cent.push_back(pos);
  }
  
  //
  //
  //
  TString inputroot;
  TString outputroot;
  inputroot=Form("../sim/gmn_4.5GeV2_opcal_job0.root");
  outputroot= Form("hist/Optics_%s_%d_fit_tree.root",OpticsID.Data(),FileID);
  //
  //
  TString YtarDeltaCutFile;
  TFile *fYtarDeltaCut;
  vector <TCutG*> ytar_delta_cut;
  if (CutYtarFlag) {
    YtarDeltaCutFile=Form("cuts/ytar_delta_%s_%d_cut.root",OpticsID.Data(),FileID);
    fYtarDeltaCut = new TFile(YtarDeltaCutFile);
    cout << " Cut file = " << YtarDeltaCutFile << endl;
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
	    Int_t npt = tempg->GetN();
	    //cout << "hYpFpYFp_cut = " << nf << " " << nd << " " << nc << " npts = " << npt << endl;
	    //cout << "hYpFpYFp_cut = " << nc << " " << nf << " " << nd << endl;
	    ypfp_yfp_cut[nf][nd].push_back(tempg);
	  } else {
	    //cout << " No hYpFpYFp_cut = " << nc << " " << nf << " " << nd << endl;
	    ypfp_yfp_cut[nf][nd].push_back(tempg);
	  }
	}}}
  }
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

  Double_t xptarT,ytarT,yptarT,ysieveT,xsieveT,ztarT,ztar,pinvtheta;
  //
  TFile hroot(outputroot,"recreate");
  TTree *otree = new TTree("TFit","FitTree");
  otree->Branch("ys",&ysieve);
  otree->Branch("ysT",&ysieveT);
  otree->Branch("xs",&xsieve);
  otree->Branch("xsT",&xsieveT);
  otree->Branch("ztarT",&ztarT);
  otree->Branch("xtar",&xtar_recon);
  otree->Branch("ztar",&ztar);
  otree->Branch("xptar",&xptar_recon);
  otree->Branch("yptar",&yptar_recon);
  otree->Branch("ytar",&ytar_recon);
  otree->Branch("xptarT",&xptarT);
  otree->Branch("yptarT",&yptarT);
  otree->Branch("ytarT",&ytarT);
  otree->Branch("xpfp",&xpfp);
  otree->Branch("ypfp",&ypfp);
  otree->Branch("xfp",&xfp);
  otree->Branch("yfp",&yfp);
  otree->Branch("pinvtheta",&pinvtheta);

  CentAngle=CentAngle*3.14159/180.;
  TVector3 BB_zaxis( sin(CentAngle), 0.0, cos(CentAngle) ); //BB is on beam right, global x axis points to beam left
  TVector3 BB_xaxis(0,-1,0); //X axis of transport coordinates is vertically down:
  TVector3 BB_yaxis = (BB_zaxis.Cross(BB_xaxis)).Unit();
  
  spec_xaxis_tgt = BB_xaxis;
  spec_yaxis_tgt = BB_yaxis;
  spec_zaxis_tgt = BB_zaxis;
  
  spec_zaxis_fp = BB_zaxis;
  spec_yaxis_fp = BB_yaxis;
  spec_zaxis_fp.Rotate(-tracker_pitch_angle, spec_yaxis_fp);


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
    }
  }

  // loop over entries
  Double_t zdis_sieve = 1.172;//front of sieve
  Long64_t nentries = tsimc->GetEntries();
  cout << " start loop " << nentries << endl;
  for (int i = 0; i < nentries; i++) {
    tsimc->GetEntry(i);
    if (i%50000==0)cout << " Entry = " << i << endl;

    //determine if good track
    bool goodtrack = false;

    //if(tsimc->Earm.BBGEM.Track.ntracks==1&&tsimc->Earm.BBGEM.Track.NumHits[0]==5&&abs(tsimc->Earm.BBGEM.Track.X)<0.55){
    if (ntracks==1 && nhits->at(0)==5 && xfp->at(0)<0.55 && xfp->at(0)>-0.55 && chisq->at(0)/ndf->at(0)<=30.0){
      goodtrack=1;
    }
 
    vy = 0;
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

	  //beam left:
	  vz_fit = -ytar_fit / (sin(CentAngle) + cos(CentAngle)*yptar_fit);
	  vz_recon = -ytar_recon / (sin(CentAngle) + cos(CentAngle)*yptar_recon);
	  //beam right:
	  //vz_fit = ytar_fit / (sin(CentAngle) - cos(CentAngle)*yptar_fit);
	  //vz_recon = ytar_recon / (sin(CentAngle) - cos(CentAngle)*yptar_recon);

	  
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
	
	xsieve = xtar_recon + xptar_recon*zdis_sieve;
	ysieve = ytar_recon + yptar_recon*zdis_sieve;

	Int_t nf_found=-1, nd_found=-1,ny_found=-1,nx_found=-1;
	for  (UInt_t nf=0;nf<ytar_delta_cut.size();nf++) {
	  if (ytar_delta_cut[nf]->IsInside(ytar_recon,pinv_recon*thetabend_recon)) nf_found=nf;
	} 
	for  (UInt_t nd=0;nd<ndelcut;nd++) {
	  if ( pinv_recon*thetabend_recon >=delcut[nd]-delwidth[nd] && pinv_recon*thetabend_recon <delcut[nd]+delwidth[nd])  nd_found=nd;
	}
	if (nf_found!=-1 && nd_found!=-1) {
	  for  (UInt_t ny=0;ny<7;ny++) {
	    if (ypfp_yfp_cut[nf_found][nd_found][ny] && ypfp_yfp_cut[nf_found][nd_found][ny]->IsInside(ypfp_fit->at(0),yfp_fit->at(0))) ny_found=ny;
	  }
	  for  (UInt_t nx=0;nx<13;nx++) {
	    if (xpfp_xfp_cut[nf_found][nd_found][nx] && xpfp_xfp_cut[nf_found][nd_found][nx]->IsInside(xpfp_fit->at(0),xfp_fit->at(0))) nx_found=nx;
	  }
	}
	if (nf_found !=-1 && nd_found!=-1 && ny_found!=-1 && nx_found!=-1) {
	  //BB beam left:
	  yptarT = (ys_cent[ny_found]+ztar_foil[nf_found]*TMath::Sin(CentAngle))/(zdis_sieve-ztar_foil[nf_found]*TMath::Cos(CentAngle));
	  xptarT = (xs_cent[nx_found])/(zdis_sieve-ztar_foil[nf_found]*TMath::Cos(CentAngle));
	  ytarT = -ztar_foil[nf_found]*(TMath::Sin(CentAngle)+yptarT*TMath::Cos(CentAngle));

	  //BB beam right:
	  //yptarT = (ys_cent[ny_found]-ztar_foil[nf_found]*TMath::Sin(CentAngle))/(zdis_sieve-ztar_foil[nf_found]*TMath::Cos(CentAngle));
	  //xptarT = (xs_cent[nx_found])/(zdis_sieve-ztar_foil[nf_found]*TMath::Cos(CentAngle));
	  //ytarT = -ztar_foil[nf_found]*(-TMath::Sin(CentAngle)+yptarT*TMath::Cos(CentAngle));

	  
	  ysieveT=ys_cent[ny_found];
	  xsieveT=xs_cent[nx_found];
	  ztarT=ztar_foil[nf_found];
	  ztar=vz_recon;
	  pinvtheta =pinv_recon*thetabend_recon;

	  
	  otree->Fill();
	}
    }//end if good track
  }//end entries
  otree->Write();
}
