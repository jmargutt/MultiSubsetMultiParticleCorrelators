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

#define ALIFLOWANALYSISCORRELATOR_CXX

#include "Riostream.h"
#include "AliFlowAnalysisCorrelator.h"

using std::cout;
using std::endl;

ClassImp(AliFlowAnalysisCorrelator)

//-----------------------------------------------------------------------
AliFlowAnalysisCorrelator::AliFlowAnalysisCorrelator():
fOrder(2)
{
  Int_t iHar[2] = {2,2};
  fHarmonics = TArrayI(2,iHar);
  Int_t iSub[2] = {1,1};
  fSubevents = TArrayI(2,iSub);
}
//-----------------------------------------------------------------------
AliFlowAnalysisCorrelator::AliFlowAnalysisCorrelator(Int_t ord, TArrayI& har, TArrayI& sub):
fOrder(ord),
fHarmonics(har),
fSubevents(sub)
{
}
//-----------------------------------------------------------------------
AliFlowAnalysisCorrelator::~AliFlowAnalysisCorrelator()
{
  //destructor
}
//-----------------------------------------------------------------------
Int_t AliFlowAnalysisCorrelator::GetHarmonic(Int_t ind) const {
  if(ind>fOrder || ind<0) return 0;
  else return fHarmonics[ind];
}
//-----------------------------------------------------------------------
Int_t AliFlowAnalysisCorrelator::GetSubevent(Int_t ind) const {
  if(ind>fOrder || ind<0) return 0;
  else return fSubevents[ind];
}
