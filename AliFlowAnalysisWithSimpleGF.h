/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


#ifndef ALIFLOWANALYSISWITHSIMPLEGF_H
#define ALIFLOWANALYSISWITHSIMPLEGF_H

class AliFlowTrackSimple;
class AliFlowEventSimple;
class AliFlowCommonHist;
class AliFlowCommonHistResults;

#include "TMatrixD.h"
#include "TList.h"
#include "AliFlowAnalysis.h"
#include "AliFlowAnalysisCorrelator.h"
#include <complex>
#include <cmath>

class TH1D;
class TH1F;
class TH2D;
class TProfile;
class TDirectoryFile;
class AliFlowAnalysisCorrelator;

//////////////////////////////////////////////////////////////////////////////////
// Description: Flow analysis with Generic Framework (arXiv:1312.3572) + subevents
// author: Jacopo Margutti (jacopo.margutti@cern.ch)
//////////////////////////////////////////////////////////////////////////////////

class AliFlowAnalysisWithSimpleGF : public AliFlowAnalysis {

public:

  AliFlowAnalysisWithSimpleGF();            //default constructor
  virtual  ~AliFlowAnalysisWithSimpleGF();  //destructor

  void Init();                                       //Define output objects
  void Make(AliFlowEventSimple* anEvent);            //Main routine
  void GetOutputHistograms(TList *outputListHistos); //Copy output objects from TList
  void Finish();                                     //Fill results
  void WriteHistograms(TDirectoryFile *outputFileName) const; //writes histograms locally (for OnTheFly)
  void ResetEventByEventQuantities();

  void SetHistList(TList* const hlist) {this->fHistList = hlist;};
  TList* GetHistList() const {return this->fHistList;};

  void SetNCorrelators(Int_t num);
  void SetNSubevents(Int_t num);
  void SetInfoCorrelator(Int_t num, Int_t ord, TArrayI& har, TArrayI& sub);
  void DefineSubevent(Int_t num, Int_t charge, Double_t etamin, Double_t etamax);
  Int_t IsTrackInSubevent(AliFlowTrackSimple *pTrack);
  void CalculateFlowGF();
  std::complex<double> ucN(const Int_t n, const TArrayI& harmonics, const TArrayI& subevents);
  std::complex<double> ucN2(const Int_t n, TArrayI& h, TArrayI& cnt, const TArrayI& se);

private:

  AliFlowAnalysisWithSimpleGF(const AliFlowAnalysisWithSimpleGF& anAnalysis);            //copy constructor
  AliFlowAnalysisWithSimpleGF& operator=(const AliFlowAnalysisWithSimpleGF& anAnalysis); //assignment operator
  TList* fHistList; //! base list to hold all output object

  const static Int_t fkMaxNSub = 10; // maximum number of subevents
  TMatrixD *fReQ[fkMaxNSub]; //! Real part of Q-vector: fReQ[m][k] = sum_{i=1}^{M} w_{i}^{k} cos(m*phi_{i})
  TMatrixD *fImQ[fkMaxNSub]; //! Imaginary part of Q-vector: fImQ[m][k] = sum_{i=1}^{M} w_{i}^{k} sin(m*phi_{i})
  Int_t fChargeSub[fkMaxNSub];    // particle charge of subevents
  Double_t fEtaMinSub[fkMaxNSub]; // minimum eta of subevents
  Double_t fEtaMaxSub[fkMaxNSub]; // maximum eta of subevents
  Int_t fNSubevents;   // number of subevents

  Int_t fNCorrelators; // number of correlators
  const static Int_t fkMaxNumCorr = 100; // maximum number of correlators
  TProfile* fProfileCorrelator[fkMaxNumCorr]; // m-particle correlators - profiles
  TH1D* fHistogramCorrelator[fkMaxNumCorr]; // m-particle correlators - histograms
  AliFlowAnalysisCorrelator fInfoCorrelators[fkMaxNumCorr]; // information on correlators

  ClassDef(AliFlowAnalysisWithSimpleGF,1);  // class version
};

#endif
