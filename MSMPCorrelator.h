/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////////////
// Base class of Multi-Subset Multi-Particle Correlators
// author: Jacopo Margutti (jacopo.margutti@cern.ch)
////////////////////////////////////////////////////////////////////////////////

#ifndef MSMPCORRELATOR_H
#define MSMPCORRELATOR_H
#include "TArrayI.h"
class TObjArray;

class MSMPCorrelator {

public:

  MSMPCorrelator();                                      // dummy constructor
  MSMPCorrelator(Int_t ord, TArrayI& har, TArrayI& sub); // default constructor
  virtual  ~MSMPCorrelator();                            // destructor

  void SetOrder(Int_t ord) {this->fOrder = ord;};             // set order of correlator
  void SetHarmonics(TArrayI& har) {this->fHarmonics = har;};  // set harmonics of correlator
  void SetSubsets(TArrayI& sub) {this->fSubsets = sub;};      // set subsets of correlator

  Int_t GetOrder() const {return this->fOrder;};             // get order of correlator
  Int_t GetHarmonic(Int_t ind) const;                        // get harmonic of one element of the correlator
  Int_t GetSubset(Int_t ind) const;                          // get subset of one element of the correlator
  TArrayI GetHarmonics() const {return this->fHarmonics;};   // get harmonics of correlator
  TArrayI GetSubsets() const {return this->fSubsets;};       // get subsets of correlator

private:

  Int_t fOrder;         // order (number of particles to correlate)
  TArrayI fHarmonics;   // harmonics
  TArrayI fSubsets;     // subsets

  ClassDef(MSMPCorrelator,1)  // class version
};


#endif
