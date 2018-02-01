/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


#ifndef ALIFLOWANALYSISCORRELATOR_H
#define ALIFLOWANALYSISCORRELATOR_H

#include "TMatrixD.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TRandom3.h"
#include "AliFlowCommonConstants.h"
#include "TNamed.h"
#include <complex>
#include <cmath>

class TObjArray;
class TList;
class TFile;
class TGraph;
class TH1;
class TH3;
class TProfile;
class TProfile2D;
class TProfile3D;
class TDirectoryFile;
class TRandom3;
class TNtuple;
class THnSparse;

//////////////////////////////////////////////////////////////////////////
// Description: base class for correlators with subevents
// author: Jacopo Margutti (jacopo.margutti@cern.ch)
//////////////////////////////////////////////////////////////////////////

class AliFlowAnalysisCorrelator {

public:

  AliFlowAnalysisCorrelator();                                      // empty constructor
  AliFlowAnalysisCorrelator(Int_t ord, TArrayI& har, TArrayI& sub); // default constructor
  virtual  ~AliFlowAnalysisCorrelator();                            // destructor

  void SetOrder(Int_t ord) {this->fOrder = ord;};
  void SetHarmonics(TArrayI& har) {this->fHarmonics = har;};
  void SetSubevents(TArrayI& sub) {this->fSubevents = sub;};

  Int_t GetOrder() const {return this->fOrder;};
  Int_t GetHarmonic(Int_t ind) const;
  Int_t GetSubevent(Int_t ind) const;
  TArrayI GetHarmonics() const {return this->fHarmonics;};
  TArrayI GetSubevents() const {return this->fSubevents;};

private:

  Int_t fOrder;
  TArrayI fHarmonics;
  TArrayI fSubevents;

  ClassDef(AliFlowAnalysisCorrelator,1)  // class version
};


#endif
