/*************************************************************************
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

////////////////////////////////////////////////////////////////////////////////
// Base class of Multi-Subset Multi-Particle Correlators
// author: Jacopo Margutti (jacopo.margutti@cern.ch)
////////////////////////////////////////////////////////////////////////////////

#define MSMPCORRELATOR_CXX

#include "Riostream.h"
#include "MSMPCorrelator.h"

using std::cout;
using std::endl;

ClassImp(MSMPCorrelator)

//-----------------------------------------------------------------------
MSMPCorrelator::MSMPCorrelator():
fOrder(2)
{
  // dummy constructor
  Int_t iHar[2] = {2,2};
  fHarmonics = TArrayI(2,iHar);
  Int_t iSub[2] = {1,1};
  fSubsets = TArrayI(2,iSub);
}
//-----------------------------------------------------------------------
MSMPCorrelator::MSMPCorrelator(Int_t ord, TArrayI& har, TArrayI& sub):
fOrder(ord),
fHarmonics(har),
fSubsets(sub)
{
  // default constructor
}
//-----------------------------------------------------------------------
MSMPCorrelator::~MSMPCorrelator()
{
  // destructor
}
//-----------------------------------------------------------------------
Int_t MSMPCorrelator::GetHarmonic(Int_t ind) const {
  // get harmonic of one element of the correlator
  if(ind>fOrder || ind<0) return 0;
  else return fHarmonics[ind];
}
//-----------------------------------------------------------------------
Int_t MSMPCorrelator::GetSubset(Int_t ind) const {
  // get subset of one element of the correlator
  if(ind>fOrder || ind<0) return 0;
  else return fSubsets[ind];
}
