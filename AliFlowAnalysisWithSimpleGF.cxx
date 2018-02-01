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

#define ALIFLOWANALYSISWITHSIMPLEGF_CXX

#include "TFile.h"
#include "TList.h"
#include "TMath.h"
#include "TString.h"
#include "TProfile.h"
#include "TVector2.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TMatrixD.h"

#include "AliFlowCommonConstants.h"
#include "AliFlowEventSimple.h"
#include "AliFlowVector.h"
#include "AliFlowTrackSimple.h"
#include "AliFlowCommonHist.h"
#include "AliFlowCommonHistResults.h"
#include "AliFlowAnalysisWithSimpleGF.h"
#include "AliLog.h"

using std::cout;
using std::endl;


ClassImp(AliFlowAnalysisWithSimpleGF)

//-----------------------------------------------------------------------
AliFlowAnalysisWithSimpleGF::AliFlowAnalysisWithSimpleGF():
fHistList(NULL),
fNSubevents(1),
fNCorrelators(1)
{
  for(Int_t i=0; i<fkMaxNSub; i++) {
    fReQ[i] = NULL;
    fImQ[i] = NULL;
    fChargeSub[i] = 0;
    fEtaMinSub[i] = -10.;
    fEtaMaxSub[i] = 10.;
  }
  for(Int_t i=0; i<fkMaxNumCorr; i++) {
    fProfileCorrelator[i] = NULL;
    fHistogramCorrelator[i] = NULL;
    fInfoCorrelators[i] = AliFlowAnalysisCorrelator();
  }
}

//-----------------------------------------------------------------------
AliFlowAnalysisWithSimpleGF::~AliFlowAnalysisWithSimpleGF()
{
  //destructor
  delete fHistList;
}

//-----------------------------------------------------------------------
void AliFlowAnalysisWithSimpleGF::Init()
{
  printf("****************************************************************\n");
  printf("    initialize AliFlowAnalysis with simple Generic Framework    \n");

  //Define all histograms
  Bool_t oldHistAddStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  fHistList = new TList();
  fHistList->SetName("cobjGF");
  fHistList->SetOwner();

  for(Int_t i=0; i<fNSubevents; i++) {
    fReQ[i] = new TMatrixD(21,9);
    fImQ[i] = new TMatrixD(21,9);
  }

  printf("\n Subevents defined: %d \n",fNSubevents);
  for(Int_t i=0; i<fNSubevents; i++) {
    printf("%d) charge: %s, eta-range: (%.1f, %.1f) \n",i+1,(fChargeSub[i]==0?"both":(fChargeSub[i]>0?"pos":"neg")),fEtaMinSub[i],fEtaMaxSub[i]);
  }

  printf("\n Correlators defined: %d \n",fNCorrelators);
  for(Int_t i=0; i<fNCorrelators; i++) {
    TString NameCorrelator = "har(";
    Int_t ord = fInfoCorrelators[i].GetOrder();
    for(Int_t k=0; k<ord; k++) {
      if(k!=0) NameCorrelator += ",";
      NameCorrelator += Form("%d",fInfoCorrelators[i].GetHarmonic(k));
    }
    NameCorrelator += ")_sub(";
    for(Int_t k=0; k<ord; k++) {
      if(k!=0) NameCorrelator += ",";
      NameCorrelator += Form("%d",fInfoCorrelators[i].GetSubevent(k));
    }
    NameCorrelator += ")";
    fProfileCorrelator[i] = new TProfile(Form("fProfileCorrelator_%s",NameCorrelator.Data()),Form("fProfileCorrelator_%s",NameCorrelator.Data()),1,0.,1.);
    fProfileCorrelator[i]->Sumw2();
    fHistList->Add(fProfileCorrelator[i]);
    fHistogramCorrelator[i] = new TH1D(Form("fHistogramCorrelator_%s",NameCorrelator.Data()),Form("fHistogramCorrelator_%s",NameCorrelator.Data()),1,0.,1.);
    fHistList->Add(fHistogramCorrelator[i]);
    printf("%d) %s\n",i+1,NameCorrelator.Data());
  }

  printf("\n****************************************************************\n");

  TH1::AddDirectory(oldHistAddStatus);
}

//-----------------------------------------------------------------------
void AliFlowAnalysisWithSimpleGF::SetInfoCorrelator(Int_t num, Int_t ord, TArrayI& har, TArrayI& sub)
{
  Int_t ind = num-1; // position in array
  if(ind<0 || ind>fNCorrelators) {
    printf("error: correlator %d out-of-range, set correct number of correlators first \n",num);
    exit(0);
  }
  fInfoCorrelators[ind] = AliFlowAnalysisCorrelator(ord, har, sub);
}

//-----------------------------------------------------------------------
void AliFlowAnalysisWithSimpleGF::DefineSubevent(Int_t num, Int_t charge, Double_t etamin, Double_t etamax)
{
  Int_t ind = num-1; // position in array
  if(ind<0 || ind>fNSubevents) {
    printf("error: subevent %d out-of-range, set correct number of subevents first \n",num);
    exit(0);
  }
  fChargeSub[ind] = charge;
  fEtaMinSub[ind] = etamin;
  fEtaMaxSub[ind] = etamax;
}

//-----------------------------------------------------------------------
void AliFlowAnalysisWithSimpleGF::Make(AliFlowEventSimple* anEvent)
{
  // Scalar Product method
  if (!anEvent) return; // for coverity

  //loop over the tracks of the event
  AliFlowTrackSimple*   pTrack = NULL;
  Int_t iNumberOfTracks = anEvent->NumberOfTracks();
  Double_t dMq = 0;
  for (Int_t i=0;i<iNumberOfTracks;i++) {
    // so this is a track loop ...
    pTrack = anEvent->GetTrack(i) ;
    if (!pTrack) continue;
    Double_t dPhi = pTrack->Phi();
    Double_t dPt  = pTrack->Pt();
    Double_t dEta = pTrack->Eta();
    Int_t iSub = IsTrackInSubevent(pTrack);
    if(iSub < 0) continue;

    // determine track weight (TBI)
    Double_t wTrack = 1.;

    for(Int_t m=0;m<21;m++) {
      for(Int_t k=0;k<9;k++) {
        (*fReQ[iSub])(m,k) += pow(wTrack,k)*TMath::Cos(m*dPhi);
        (*fImQ[iSub])(m,k) += pow(wTrack,k)*TMath::Sin(m*dPhi);
      }
    }

  }//loop over tracks

  this->CalculateFlowGF();

  // o) Reset all event-by-event quantities
  this->ResetEventByEventQuantities();
}

//--------------------------------------------------------------------
void AliFlowAnalysisWithSimpleGF::CalculateFlowGF()
{
  for(Int_t i=0; i<fNCorrelators; i++) {

    // get info on correlator to calculate
    Int_t CorOrder = fInfoCorrelators[i].GetOrder();           // order: how many particles to correlate
    TArrayI CorHarmonics(fInfoCorrelators[i].GetHarmonics());  // harmonics
    TArrayI CorSubevents(fInfoCorrelators[i].GetSubevents());  // subevents

    // check that there are enough particles in the event
    Int_t iMult = 0;
    for(Int_t k=0; k<CorOrder; k++) {
      iMult += (*fReQ[CorSubevents[k]-1])(0,0);
    }
    if(iMult<CorOrder) continue;

    TArrayI CorNullArray(CorOrder);
    for(Int_t k=0; k<CorOrder; k++) CorNullArray[k] = 0;

    // compute correlator
    std::complex<double> N = this->ucN(CorOrder, CorHarmonics, CorSubevents); // numerator
    std::complex<double> D = this->ucN(CorOrder, CorNullArray, CorSubevents); // denominator

    // fill TProfiles, use denominator (number of combinations) as weight
    if(D.real()>0.) {
      fProfileCorrelator[i]->Fill(0.5,N.real()/D.real(),D.real());
    }

  } // end of for(Int_t i=0; i<fNCorrelators; i++)
}

//--------------------------------------------------------------------
std::complex<double> AliFlowAnalysisWithSimpleGF::ucN(const Int_t n, const TArrayI& harmonics, const TArrayI& subevents)
{
  TArrayI counts(n);
  for (Int_t i = 0; i < n; i++) {
    counts[i] = 1;
  }
  TArrayI hh(harmonics);

  return ucN2(n, hh, counts, subevents);
}

//--------------------------------------------------------------------
std::complex<double> AliFlowAnalysisWithSimpleGF::ucN2(const Int_t n, TArrayI& h, TArrayI& cnt, const TArrayI& se)
{
  Int_t j = n-1;
  std::complex<double> c;
  if(h[j] >= 0) {
    c = std::complex<double>((*fReQ[se[j]-1])(h[j],cnt[j]),(*fImQ[se[j]-1])(h[j],cnt[j]));
  } else {
    c = std::complex<double>((*fReQ[se[j]-1])(-h[j],cnt[j]),-(*fImQ[se[j]-1])(-h[j],cnt[j]));
  }

  if (n == 1) return c;

  c *= ucN2(j, h, cnt, se);

  if (cnt[j] > 1) return c;

  for (Int_t i = 0; i < (n-1); i++) {
    h[i] += h[j];
    cnt[i] = cnt[i] + 1;
    Double_t factor = 1.*(cnt[i]-1);
    if(se[i]==se[j]) {
      c -= factor * ucN2(j, h, cnt, se);
    }
    cnt[i]--;
    h[i] -= h[j];
  }
  return c;
}

//--------------------------------------------------------------------
Int_t AliFlowAnalysisWithSimpleGF::IsTrackInSubevent(AliFlowTrackSimple *pTrack)
{
  Int_t subev = -1;
  for(Int_t s=0; s<fNSubevents; s++) {
    Bool_t IsInSubevent = kTRUE;
    if((fChargeSub[s] != 0) && (pTrack->Charge() != fChargeSub[s])) IsInSubevent = kFALSE;
    if((pTrack->Eta() < fEtaMinSub[s]) || (pTrack->Eta() > fEtaMaxSub[s])) IsInSubevent = kFALSE;
    if(IsInSubevent) subev = s;
  }
  return subev;
}

//--------------------------------------------------------------------
void AliFlowAnalysisWithSimpleGF::GetOutputHistograms(TList *outputListHistos)
{
  //get pointers to all output histograms (called before Finish())
  fHistList = outputListHistos;
}

//--------------------------------------------------------------------
void AliFlowAnalysisWithSimpleGF::ResetEventByEventQuantities()
{
  for(Int_t i=0; i<fNSubevents; i++) {
    fReQ[i]->Zero();
    fImQ[i]->Zero();
  }
}

//--------------------------------------------------------------------
void AliFlowAnalysisWithSimpleGF::Finish()
{
  //calculate flow and fill the AliFlowCommonHistResults
  printf("**********************************************\n");
  printf("**********************************************\n");
  printf("    multi-particle correlators from simple         \n");
  printf("       Generic Framework with subevents           \n\n");

  for(Int_t i=0; i<fNCorrelators; i++) {
    for(Int_t bx=1; bx<=fProfileCorrelator[i]->GetNbinsX(); bx++) {
      Double_t stats[6]={0.};
      fProfileCorrelator[i]->GetXaxis()->SetRange(bx,bx);
      fProfileCorrelator[i]->GetStats(stats);
      Double_t SumWeig   = stats[0];
      Double_t SumWeigSq  = stats[1];
      Double_t SumTwo  = stats[4];
      Double_t SumTwoSq = stats[5];
      if(SumWeig>0.) {
        Double_t Corr = SumTwo/SumWeig;
        Double_t SqCorr = SumTwoSq/SumWeig;
        Double_t Weig = SumWeig;
        Double_t SqWeig = SumWeigSq;
        Double_t spread=0., termA=0., termB=0.;
        if(SqCorr-pow(Corr,2.)>=0.) { spread = pow(SqCorr-pow(Corr,2.),0.5); }
        if(TMath::Abs(Weig)>0.) { termA = (pow(SqWeig,0.5)/Weig); }
        if(1.-pow(termA,2.)>0.) { termB = 1./pow(1.-pow(termA,2.),0.5); }
        Double_t CorrErr = termA*spread*termB; // final error (unbiased estimator for standard deviation)
        if(CorrErr) {
          fHistogramCorrelator[i]->SetBinContent(bx,Corr);
          fHistogramCorrelator[i]->SetBinError(bx,CorrErr);
          printf("%s: %e +/- %e \n",fHistogramCorrelator[i]->GetName(),Corr,CorrErr);
        }
      }
    } // end of for(Int_t pt=1;pt<=fPtDiffNBins;pt++)
    fProfileCorrelator[i]->GetXaxis()->SetRange(1,fProfileCorrelator[i]->GetNbinsX());
  }

  printf("\n");
  printf("**********************************************\n");
  printf("**********************************************\n");
}

//-----------------------------------------------------------------------
void AliFlowAnalysisWithSimpleGF::SetNCorrelators(Int_t num)
{
  if(num>fkMaxNumCorr) {
    printf("Error: maximum number of correlators is %d \n",fkMaxNumCorr);
    exit(0);
  }
  fNCorrelators = num;
}

//-----------------------------------------------------------------------
void AliFlowAnalysisWithSimpleGF::SetNSubevents(Int_t num)
{
  if(num>fkMaxNSub) {
    printf("Error: maximum number of subevents is %d \n",fkMaxNSub);
    exit(0);
  }
  fNSubevents = num;
}

//-----------------------------------------------------------------------
void AliFlowAnalysisWithSimpleGF::WriteHistograms(TDirectoryFile *outputFileName) const
{
  //store the final results in output .root file
  outputFileName->Add(fHistList);
  outputFileName->Write(outputFileName->GetName(), TObject::kSingleKey);
}

//--------------------------------------------------------------------
