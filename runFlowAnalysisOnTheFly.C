////////////////////////////////////////////////////////////////////////////////
// simple macro to test flow methods
// adapted from  AliPhysics/PWGCF/FLOW/macros/runFlowAnalysisOnTheFly.C
// to test Multi-Subset Multi-Particle Correlators
// author: Jacopo Margutti (jacopo.margutti@cern.ch)
////////////////////////////////////////////////////////////////////////////////

// Settings for the simulation of events 'on the fly':
//  a) Determine how many events you want to create;
//  b) Set random or same seed for random generator;
//  c) Determine multiplicites of events;
//  d) Parametrize the phi distribution;
//  d1) Enable/disable uniform event-wise fluctuations of v2;
//  d2) Enable/diable pt dependence of v2;
//  e) Parametrize the pt distribution;
//  f) Determine how many times each sampled particle will be taken (simulating nonflow);
//  g) Configure detector's:
//  g1) acceptance;
//  g2) efficiency;
//  h) Decide which flow analysis methods you will use;
//  i) Define simple cuts for Reference Particle (RP) selection;
//  j) Define simple cuts for Particle of Interest (POI) selection;

// a) Determine how many events you want to create:
Int_t iNevts = 1000; // total statistics

// b) Set random or same seed for random generator:
Bool_t bSameSeed = kFALSE; // if kTRUE, the created events are the same when re-doing flow analysis 'on the fly'

// c) Determine multiplicites of events:
//    Remark 1: Multiplicity M for each event is sampled uniformly from interval iMinMult <= M < iMaxMult;
//    Remark 2: For constant M of e.g. 500 for each event, set iMinMult = 500 and iMaxMult = 501.
Int_t iMinMult = 500; // uniformly sampled multiplicity is >= iMinMult
Int_t iMaxMult = 501; // uniformly sampled multiplicity is < iMaxMult

// d) Parametrize the phi distribution:
//    Remark 1: Hardwired is Fourier-like distribution f(phi) = (1/2pi)(1+sum_{n=1}^{4} 2v_n cos[n(phi-rp)]),
//              where reaction plane (rp) is sampled uniformly for each event from interval [0,2pi]
Double_t dV1 = 0.0; // constant harmonic v1
Double_t dV2 = 0.05; // constant harmonic v2
Double_t dV3 = 0.0; // constant harmonic v3
Double_t dV4 = 0.0; // constant harmonic v4
Double_t dV5 = 0.0; // constant harmonic v5
Double_t dV6 = 0.0; // constant harmonic v6
//    Remark 2: By default all harmonics are constant for each event and for each particle. However, for v2
//              the uniform event-wise fluctuations or pt dependence can be enabled:
//  d1) Enable/disable uniform event-wise fluctuations of v2:
Bool_t bUniformFluctuationsV2 = kFALSE; // enable uniform event-wise flow fluctuations (set than also dMinV2 and dMaxV2 bellow)
Double_t dMinV2 = 0.04; // lower boundary on v2, when bUniformFluctuationsV2 = kTRUE
Double_t dMaxV2 = 0.06; // upper boundary on v2, when bUniformFluctuationsV2 = kTRUE
//  d2) Enable/disable pt dependence of v2:
Bool_t bPtDependentV2 = kFALSE; // enable pt dependence of v2 (set then also dV2vsPtMax and dV2vsPtCutOff bellow)
Double_t dV2vsPtCutOff = 2.0; // up to pt = dV2vsPtCutOff v2 is growing linearly as a function of pt
Double_t dV2vsPtMax = 0.20; // for pt >= dV2vsPtCutOff, v2(pt) = dV2vsPtMax

// e) Parametrize the pt distribution:
//    Remark: Hardwired is Boltzmann distribution f(pt) = pt*exp[-sqrt(dMass^2+pt^2)/dT]
Double_t dMass = 0.13957; // mass in GeV/c^2 (e.g. m_{pions} = 0.13957)
Double_t dTemperature = 0.44; // "temperature" in GeV/c (increase this parameter to get more high pt particles)

// f) Determine how many times each sampled particle will be taken in the analysis (simulating nonflow):
Int_t nTimes = 1; // e.g. for nTimes = 2, strong 2-particle nonflow correlations are introduced

// g1) Configure detector's acceptance:
Bool_t uniformAcceptance = kTRUE; // if kTRUE: detectors has uniform azimuthal acceptance.
                                  // if kFALSE: you will simulate detector with non-uniform acceptance in one or
                                  // two sectors. For each sector you specify phiMin, phiMax and probability p.
                                  // Then all particles emitted in direction phiMin < phi < phiMax will be taken
                                  // with probability p. If p = 0, that sector is completely blocked. Set bellow
                                  // phiMin1, phiMax1, p1 for the first sector and phiMin2, phiMax2, p2 for the second
                                  // sector. If you set phiMin2 = phiMax2 = p2 = 0, only first non-uniform sector is
                                  // simulated.
// 1st non-uniform sector:
Double_t phiMin1 = 60; // first non-uniform sector starts at this azimuth (in degrees)
Double_t phiMax1 = 120; // first non-uniform sector ends at this azimuth (in degrees)
Double_t p1 = 0.5; // probablitity that particles emitted in [phiMin1,phiMax1] are taken
// 2nd non-uniform sector:
Double_t phiMin2 = 0.; // first non-uniform sector starts at this azimuth (in degrees)
Double_t phiMax2 = 0.; // first non-uniform sector ends at this azimuth (in degrees)
Double_t p2 = 0.; // probablitity that particles emitted in [phiMin2,phiMax2] are taken

// g2) Configure detector's efficiency:
Bool_t uniformEfficiency = kTRUE; // if kTRUE: detectors has uniform pT efficiency
                                  // if kFALSE: you will simulate detector with non-uniform pT efficiency.
                                  // Then all particles emitted in ptMin <= pt < ptMax will be taken
                                  // with probability p, to be specified in lines just below.
Double_t ptMin = 0.8; // non-uniform efficiency vs pT starts at pT = fPtMin
Double_t ptMax = 1.2; // non-uniform efficiency vs pT ends at pT = fPtMax
Double_t p = 0.5; // probablitity that particles emitted in [ptMin,ptMax> are taken

// h) Decide which flow analysis methods you will use:
Bool_t MCEP     = kTRUE; // Monte Carlo Event Plane
Bool_t MSMPC    = kTRUE; // Multi-Subevent Multi-Particle Correlations (NEW!)

// i) Define simple cuts for Reference Particle (RP) selection:
Double_t ptMinRP = 0.0; // in GeV
Double_t ptMaxRP = 10.0; // in GeV
Double_t etaMinRP = -1.;
Double_t etaMaxRP = 1.;
Double_t phiMinRP = 0.0; // in degrees
Double_t phiMaxRP = 360.0; // in degrees
Bool_t bUseChargeRP = kFALSE; // if kFALSE, RPs with both sign of charges are taken
Int_t chargeRP = 1; // +1 or -1

// j) Define simple cuts for Particle of Interest (POI) selection:
Double_t ptMinPOI = 0.0; // in GeV
Double_t ptMaxPOI = 10.0; // in GeV
Double_t etaMinPOI = -1.; //
Double_t etaMaxPOI = 1.;
Double_t phiMinPOI = 0.0; // in degrees
Double_t phiMaxPOI = 360.0; // in degrees
Bool_t bUseChargePOI = kFALSE; // if kFALSE, POIs with both sign of charges are taken
Int_t chargePOI = -1; // +1 or -1

enum anaModes {mLocal,mLocalSource,mLocalPAR};
// mLocal: Analyze data on your computer using aliroot
// mLocalPAR: Analyze data on your computer using root + PAR files
// mLocalSource: Analyze data on your computer using root + source files

#include "TStopwatch.h"
#include "TObjArray"
#include "Riostream.h"
#include "TFile.h"

int runFlowAnalysisOnTheFly(Int_t mode=mLocal)
{
 // Beging analysis 'on the fly'.

 // a) Formal necessities....;
 // b) Initialize the flow event maker 'on the fly';
 // c) If enabled, access particle weights from external file;
 // d) Configure the flow analysis methods;
 // e) Simple cuts for RPs;
 // f) Simple cuts for POIs;
 // g) Create and analyse events 'on the fly';
 // h) Create the output file and directory structure for the final results of all methods;
 // i) Calculate and store the final results of all methods.

 // a) Formal necessities....:
 CheckUserSettings();
 WelcomeMessage();
 TStopwatch timer;
 timer.Start();
 LoadLibraries(mode);

 // b) Initialize the flow event maker 'on the fly':
 UInt_t uiSeed = 0; // if uiSeed is 0, the seed is determined uniquely in space and time via TUUID
 if(bSameSeed){uiSeed = 44;}
 AliFlowEventSimpleMakerOnTheFly* eventMakerOnTheFly = new AliFlowEventSimpleMakerOnTheFly(uiSeed);
 eventMakerOnTheFly->SetMinMult(iMinMult);
 eventMakerOnTheFly->SetMaxMult(iMaxMult);
 eventMakerOnTheFly->SetMass(dMass);
 eventMakerOnTheFly->SetTemperature(dTemperature);
 eventMakerOnTheFly->SetV1(dV1);
 eventMakerOnTheFly->SetV2(dV2);
 eventMakerOnTheFly->SetV3(dV3);
 eventMakerOnTheFly->SetV4(dV4);
 eventMakerOnTheFly->SetV5(dV5);
 eventMakerOnTheFly->SetV6(dV6);
 if(bUniformFluctuationsV2)
 {
  eventMakerOnTheFly->SetUniformFluctuationsV2(bUniformFluctuationsV2);
  eventMakerOnTheFly->SetMinV2(dMinV2);
  eventMakerOnTheFly->SetMaxV2(dMaxV2);
 }
 if(bPtDependentV2)
 {
  eventMakerOnTheFly->SetPtDependentV2(bPtDependentV2);
  eventMakerOnTheFly->SetV2vsPtCutOff(dV2vsPtCutOff);
  eventMakerOnTheFly->SetV2vsPtMax(dV2vsPtMax);
 }
 eventMakerOnTheFly->SetNTimes(nTimes);
 if(!uniformAcceptance)
 {
  eventMakerOnTheFly->SetUniformAcceptance(kFALSE);
  eventMakerOnTheFly->SetFirstSectorPhiMin(phiMin1);
  eventMakerOnTheFly->SetFirstSectorPhiMax(phiMax1);
  eventMakerOnTheFly->SetFirstSectorProbability(p1);
  eventMakerOnTheFly->SetSecondSectorPhiMin(phiMin2);
  eventMakerOnTheFly->SetSecondSectorPhiMax(phiMax2);
  eventMakerOnTheFly->SetSecondSectorProbability(p2);
 }
 if(!uniformEfficiency)
 {
  eventMakerOnTheFly->SetUniformEfficiency(kFALSE);
  eventMakerOnTheFly->SetPtMin(ptMin);
  eventMakerOnTheFly->SetPtMax(ptMax);
  eventMakerOnTheFly->SetPtProbability(p);
 }
 eventMakerOnTheFly->Init();

 // c) If enabled, access particle weights from external file:
 // TBI

 // d) Configure the flow analysis methods:

 // MCEP = monte carlo event plane
 if(MCEP)
 {
  AliFlowAnalysisWithMCEventPlane *FlowAnalysisMCEP = new AliFlowAnalysisWithMCEventPlane();
  FlowAnalysisMCEP->SetHarmonic(2); // default is v2
  FlowAnalysisMCEP->Init();
 } // end of if(MCEP)

 if(MSMPC)
 {
   AliFlowAnalysisMSMPC* FlowAnalysisMSMPC = new AliFlowAnalysisMSMPC();
   // define subsets
   FlowAnalysisMSMPC->SetNSubsets(3);
   FlowAnalysisMSMPC->DefineSubset(1, 1, -1., 0.); // positive charge, negative rapidity
   FlowAnalysisMSMPC->DefineSubset(2, -1, 0., 1.); // negative charge, positive rapidity
   FlowAnalysisMSMPC->DefineSubset(3, 0, -1., 1.); // any charge, any rapidity
   // set overlapping subsets
   FlowAnalysisMSMPC->SetOverlappingSubsets(1,3);
   FlowAnalysisMSMPC->SetOverlappingSubsets(2,3);
   // define correlators
   FlowAnalysisMSMPC->SetNCorrelators(2);
   Int_t harmonics_cor1[] = {2,2,-2,-2};
   Int_t subevents_cor1[] = {1,1,2,2};
   FlowAnalysisMSMPC->SetInfoCorrelator(1,4,TArrayI(4,harmonics_cor1),TArrayI(4,subevents_cor1));
   Int_t harmonics_cor2[] = {2,2,-4};
   Int_t subevents_cor2[] = {1,1,3};
   FlowAnalysisMSMPC->SetInfoCorrelator(2,3,TArrayI(3,harmonics_cor2),TArrayI(3,subevents_cor2));
   // initialize flow analysis
   FlowAnalysisMSMPC->Init();
 } // end of if(MSMPC)

 // e) Simple cuts for RPs:
 AliFlowTrackSimpleCuts *cutsRP = new AliFlowTrackSimpleCuts();
 cutsRP->SetPtMax(ptMaxRP);
 cutsRP->SetPtMin(ptMinRP);
 cutsRP->SetEtaMax(etaMaxRP);
 cutsRP->SetEtaMin(etaMinRP);
 cutsRP->SetPhiMax(phiMaxRP*TMath::Pi()/180.);
 cutsRP->SetPhiMin(phiMinRP*TMath::Pi()/180.);
 if(bUseChargeRP){cutsRP->SetCharge(chargeRP);}

 // f) Simple cuts for POIs:
 AliFlowTrackSimpleCuts *cutsPOI = new AliFlowTrackSimpleCuts();
 cutsPOI->SetPtMax(ptMaxPOI);
 cutsPOI->SetPtMin(ptMinPOI);
 cutsPOI->SetEtaMax(etaMaxPOI);
 cutsPOI->SetEtaMin(etaMinPOI);
 cutsPOI->SetPhiMax(phiMaxPOI*TMath::Pi()/180.);
 cutsPOI->SetPhiMin(phiMinPOI*TMath::Pi()/180.);
 if(bUseChargePOI){cutsPOI->SetCharge(chargePOI);}

 // g) Create and analyse events 'on the fly':
 for(Int_t i=0;i<iNevts;i++)
 {
  // Creating the event 'on the fly':
  AliFlowEventSimple *event = eventMakerOnTheFly->CreateEventOnTheFly(cutsRP,cutsPOI);
  // Passing the created event to flow analysis methods:
  if(MCEP){FlowAnalysisMCEP->Make(event);}
  if(MSMPC){FlowAnalysisMSMPC->Make(event);}
  delete event;
 } // end of for(Int_t i=0;i<iNevts;i++)

 // h) Create the output file and directory structure for the final results of all methods:
 TString outputFileName = "AnalysisResults.root";
 TFile *outputFile = new TFile(outputFileName.Data(),"RECREATE");
 const Int_t nMethods = 2;
 TString method[nMethods] = {"MCEP","MSMPC"};
 TDirectoryFile *dirFileFinal[nMethods] = {NULL};
 TString fileName[nMethods];
 for(Int_t i=0;i<nMethods;i++)
 {
  fileName[i]+="output";
  fileName[i]+=method[i].Data();
  fileName[i]+="analysis";
  dirFileFinal[i] = new TDirectoryFile(fileName[i].Data(),fileName[i].Data());
 }

 // i) Calculate and store the final results of all methods:
 if(MCEP){FlowAnalysisMCEP->Finish();FlowAnalysisMCEP->WriteHistograms(dirFileFinal[0]);}
 if(MSMPC){FlowAnalysisMSMPC->Finish();FlowAnalysisMSMPC->WriteHistograms(dirFileFinal[1]);}

 outputFile->Close();
 delete outputFile;

 cout<<endl;
 cout<<endl;
 cout<<" ---- LANDED SUCCESSFULLY ---- "<<endl;
 cout<<endl;

 timer.Stop();
 cout << endl;
 timer.Print();

} // end of int runFlowAnalysisOnTheFly(Int_t mode=mLocal)

void SetupPar(char* pararchivename)
{
  //Load par files, create analysis libraries
  //For testing, if par file already decompressed and modified
  //classes then do not decompress.

  TString cdir(Form("%s", gSystem->WorkingDirectory() )) ;
  TString parpar(Form("%s.par", pararchivename)) ;
  if ( gSystem->AccessPathName(parpar.Data()) ) {
    gSystem->ChangeDirectory(gSystem->Getenv("ALICE_PHYSICS")) ;
    TString processline(Form(".! make %s", parpar.Data())) ;
    gROOT->ProcessLine(processline.Data()) ;
    gSystem->ChangeDirectory(cdir) ;
    processline = Form(".! mv /tmp/%s .", parpar.Data()) ;
    gROOT->ProcessLine(processline.Data()) ;
  }
  if ( gSystem->AccessPathName(pararchivename) ) {
    TString processline = Form(".! tar xvzf %s",parpar.Data()) ;
    gROOT->ProcessLine(processline.Data());
  }

  TString ocwd = gSystem->WorkingDirectory();
  gSystem->ChangeDirectory(pararchivename);

  // check for BUILD.sh and execute
  if (!gSystem->AccessPathName("PROOF-INF/BUILD.sh")) {
    printf("*******************************\n");
    printf("*** Building PAR archive    ***\n");
    cout<<pararchivename<<endl;
    printf("*******************************\n");

    if (gSystem->Exec("PROOF-INF/BUILD.sh")) {
      Error("runProcess","Cannot Build the PAR Archive! - Abort!");
      return -1;
    }
  }
  // check for SETUP.C and execute
  if (!gSystem->AccessPathName("PROOF-INF/SETUP.C")) {
    printf("*******************************\n");
    printf("*** Setup PAR archive       ***\n");
    cout<<pararchivename<<endl;
    printf("*******************************\n");
    gROOT->Macro("PROOF-INF/SETUP.C");
  }

  gSystem->ChangeDirectory(ocwd.Data());
  printf("Current dir: %s\n", ocwd.Data());
}

void CheckUserSettings()
{
 // Check if user settings make sense before taking off.

 if(iNevts <= 0)
 {
  printf("\n WARNING: nEvts <= 0 !!!! Please check your settings before taking off.\n\n");
  exit(0);
 }
 if(iMinMult < 0.)
 {
  printf("\n WARNING: iMinMult < 0 !!!! Please check your settings before taking off.\n\n");
  exit(0);
 }
 if(iMaxMult <= 0.)
 {
  printf("\n WARNING: iMaxMult <= 0 !!!! Please check your settings before taking off.\n\n");
  exit(0);
 }
 if(iMinMult >= iMaxMult)
 {
  printf("\n WARNING: iMinMult >= iMaxMult !!!! Please check your settings before taking off.\n\n");
  exit(0);
 }
 if(dMass < 0.)
 {
  printf("\n WARNING: dMass < 0 !!!! Please check your settings before taking off.\n\n");
  exit(0);
 }
 if(dTemperature <= 1e-44)
 {
  printf("\n WARNING: dTemperature <= 0 !!!! Please check your settings before taking off.\n\n");
  exit(0);
 }
 if(TMath::Abs(dV1) > 0.5)
 {
  printf("\n WARNING: |dV1| > 0.5 !!!! Please check your settings before taking off.\n\n");
  exit(0);
 }
 if(TMath::Abs(dV2) > 0.5)
 {
  printf("\n WARNING: |dV2| > 0.5 !!!! Please check your settings before taking off.\n\n");
  exit(0);
 }
 if(TMath::Abs(dV3) > 0.5)
 {
  printf("\n WARNING: |dV3| > 0.5 !!!! Please check your settings before taking off.\n\n");
  exit(0);
 }
 if(TMath::Abs(dV4) > 0.5)
 {
  printf("\n WARNING: |dV4| > 0.5 !!!! Please check your settings before taking off.\n\n");
  exit(0);
 }
 if(TMath::Abs(dV5) > 0.5)
 {
  printf("\n WARNING: |dV5| > 0.5 !!!! Please check your settings before taking off.\n\n");
  exit(0);
 }
 if(TMath::Abs(dV6) > 0.5)
 {
  printf("\n WARNING: |dV6| > 0.5 !!!! Please check your settings before taking off.\n\n");
  exit(0);
 }

 if(!uniformAcceptance && phiMin1 > phiMax1)
 {
  cout<<" WARNING: You must have phiMin1 < phiMax1 !!!!"<<endl;
  exit(0);
 }
 if(!uniformAcceptance && !((TMath::Abs(phiMin2) < 1.e-44) && (TMath::Abs(phiMax2) < 1.e-44) && (TMath::Abs(p2) < 1.e-44))
    && (phiMin2 < phiMax1 || phiMin2 > phiMax2))
 {
  cout<<" WARNING: You must have phiMin2 > phiMax1 and phiMin2 < phiMax2 !!!!"<<endl;
  exit(0);
 }
 if((phiMin1 < 0 || phiMin1 > 360) || (phiMax1 < 0 || phiMax1 > 360) ||
    (phiMin2 < 0 || phiMin2 > 360) || (phiMax2 < 0 || phiMax2 > 360) )
 {
  cout<<" WARNING: You must take azimuthal angles from interval [0,360] !!!!"<<endl;
  exit(0);
 }
 if((p1 < 0 || p1 > 1) || (p2 < 0 || p2 > 1))
 {
  cout<<" WARNING: you must take p1 and p2 from interval [0,1] !!!!"<<endl;
  exit(0);
 }
 if(bPtDependentV2 && bUniformFluctuationsV2)
 {
  cout<<" WARNING: Uniform fluctuations not supported for pt denependent v2 !!!!"<<endl;
  exit(0);
 }

} // end of void CheckUserSettings()

void WelcomeMessage()
{
 // Welcome.

 cout<<endl;
 cout<<endl;
 cout<<"      ---- ARE YOU READY TO FLY ? ----      "<<endl;
 cout<<endl;

 gSystem->Sleep(1544);

 cout<<endl;
 cout<<" ---- BEGIN FLOW ANALYSIS 'ON THE FLY' ---- "<<endl;
 cout<<endl;
 cout<<endl;

 gSystem->Sleep(1544);

} // end of void WelcomeMessage()

void LoadLibraries(const anaModes mode) {

  //--------------------------------------
  // Load the needed libraries most of them already loaded by aliroot
  //--------------------------------------
  //gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libXMLIO");
  gSystem->Load("libPhysics");

  //----------------------------------------------------------
  // >>>>>>>>>>> Local mode <<<<<<<<<<<<<<
  //----------------------------------------------------------
  if (mode==mLocal) {
    //--------------------------------------------------------
    // If you want to use already compiled libraries
    // in the aliroot distribution
    //--------------------------------------------------------
    gSystem->Load("libSTEERBase");
    gSystem->Load("libESD");
    gSystem->Load("libAOD");
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");
    gSystem->Load("libCORRFW");
    cerr<<"libCORRFW loaded..."<<endl;
    gSystem->Load("libPWGflowBase");
    cerr<<"libPWGflowBase loaded..."<<endl;
    gSystem->Load("libPWGflowTasks");
    cerr<<"libPWGflowTasks loaded..."<<endl;
  }

  else if (mode == mLocalPAR) {
    //--------------------------------------------------------
    //If you want to use root and par files from aliroot
    //--------------------------------------------------------
     //If you want to use root and par files from aliroot
    //--------------------------------------------------------
    SetupPar("STEERBase");
    SetupPar("ESD");
    SetupPar("AOD");
    SetupPar("ANALYSIS");
    SetupPar("ANALYSISalice");
    SetupPar("PWG2AOD");
    SetupPar("CORRFW");
    SetupPar("PWGflowBase");
    cerr<<"PWGflowBase.par loaded..."<<endl;
    SetupPar("PWGflowTasks");
    cerr<<"PWGflowTasks.par loaded..."<<endl;
  }

  //---------------------------------------------------------
  // <<<<<<<<<< Source mode >>>>>>>>>>>>
  //---------------------------------------------------------
  else if (mode==mLocalSource) {

    // In root inline compile

    // Constants
    gROOT->LoadMacro("Base/AliFlowCommonConstants.cxx+");

    // Flow event
    gROOT->LoadMacro("Base/AliFlowVector.cxx+");
    gROOT->LoadMacro("Base/AliFlowTrackSimple.cxx+");
    gROOT->LoadMacro("Base/AliFlowTrackSimpleCuts.cxx+");
    gROOT->LoadMacro("Base/AliFlowEventSimple.cxx+");

    // Output histosgrams
    gROOT->LoadMacro("Base/AliFlowCommonHist.cxx+");
    gROOT->LoadMacro("Base/AliFlowCommonHistResults.cxx+");

    // Flow Analysis code for various methods
    gROOT->LoadMacro("Base/AliFlowAnalysisWithMCEventPlane.cxx+");

    // Class to fill the FlowEvent on the fly (generate Monte Carlo events)
    gROOT->LoadMacro("Base/AliFlowEventSimpleMakerOnTheFly.cxx+");

    cout << "finished loading macros!" << endl;

  }
  // load extra classes for MSMPC
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  gROOT->ProcessLine(".include $ALICE_PHYSICS/include");
  gROOT->LoadMacro("MSMPCorrelator.cxx++g");
  gROOT->LoadMacro("AliFlowAnalysisMSMPC.cxx++g");
}
