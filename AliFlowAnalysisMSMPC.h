/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////////////
// Flow analysis with Multi-Subset Multi-Particle Correlators
// author: Jacopo Margutti (jacopo.margutti@cern.ch)
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIFLOWANALYSISMSMPC_H
#define ALIFLOWANALYSISMSMPC_H

class AliFlowTrackSimple;
class AliFlowEventSimple;
class AliFlowCommonHist;
class AliFlowCommonHistResults;

#include "TMatrixD.h"
#include "TList.h"
#include "AliFlowAnalysis.h"
#include "MSMPCorrelator.h"
#include <complex>

class TH1D;
class TProfile;
class TDirectoryFile;
class MSMPCorrelator;

class AliFlowAnalysisMSMPC : public AliFlowAnalysis {

public:

  AliFlowAnalysisMSMPC();            // default constructor
  virtual  ~AliFlowAnalysisMSMPC();  // destructor

  void Init();                                                // define output objects
  void Make(AliFlowEventSimple* anEvent);                     // process one event
  void GetOutputHistograms(TList *outputListHistos);          // copy output objects from TList
  void Finish();                                              // compute final values
  void WriteHistograms(TDirectoryFile *outputFileName) const; // writes histograms locally (for OnTheFly)
  void ResetEventByEventQuantities();

  void SetHistList(TList* const hlist) {this->fHistList = hlist;};
  TList* GetHistList() const {return this->fHistList;};

  void SetNCorrelators(Int_t num);                                                           // set number of correlators
  void SetNSubsets(Int_t num);                                                               // set number of subsets
  void SetInfoCorrelator(Int_t num, Int_t ord, TArrayI& har, TArrayI& sub);                  // set information of correlator (harmonics, subsets)
  void DefineSubset(Int_t num, Int_t charge, Double_t etamin, Double_t etamax);              // define subset
  void SetOverlappingSubsets(Int_t sub1, Int_t sub2);                                        // set if two subsets are not disjoint
  Bool_t IsTrackInSubset(AliFlowTrackSimple *pTrack, Int_t iSub);                            // check if track is in subset
  void CalculateFlowMSMPC();                                                                 // calculate multi-subset multi-particle Correlators
  std::complex<double> ucN(const Int_t n, const TArrayI& harmonics, const TArrayI& subsets); // recursive algorithm (1/2)
  std::complex<double> ucN2(const Int_t n, TArrayI& h, TArrayI& cnt, const TArrayI& se);     // recursive algorithm (2/2)

private:

  AliFlowAnalysisMSMPC(const AliFlowAnalysisMSMPC& anAnalysis);            // copy constructor
  AliFlowAnalysisMSMPC& operator=(const AliFlowAnalysisMSMPC& anAnalysis); // assignment operator
  TList* fHistList;                                                        //! base list to hold all output object

  const static Int_t fkMaxNSub = 10; // maximum number of Subsets
  TMatrixD *fReQ[fkMaxNSub];         //! Real part of Q-vector
  TMatrixD *fImQ[fkMaxNSub];         //! Imaginary part of Q-vector
  Int_t fChargeSub[fkMaxNSub];       // particle charge of Subsets
  Double_t fEtaMinSub[fkMaxNSub];    // minimum eta of subsets
  Double_t fEtaMaxSub[fkMaxNSub];    // maximum eta of subsets
  Int_t fNSubsets;                   // number of subsets
  TMatrixD *fAreSubsetsDisjoint;     //! matrix defining which subsets are disjoint
  Int_t fNCorrelators;               // number of correlators
  const static Int_t fkMaxNumCorr = 1000;     // maximum number of correlators
  TProfile* fProfileCorrelator[fkMaxNumCorr]; // correlators - TProfiles
  TH1D* fHistogramCorrelator[fkMaxNumCorr];   // correlators - THistograms
  MSMPCorrelator fInfoCorrelators[fkMaxNumCorr]; // information on correlators

  ClassDef(AliFlowAnalysisMSMPC,1);  // class version
};

#endif
