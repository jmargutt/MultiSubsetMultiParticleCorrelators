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
// Flow analysis with Multi-Subset Multi-Particle Correlators
// author: Jacopo Margutti (jacopo.margutti@cern.ch)
////////////////////////////////////////////////////////////////////////////////

#define ALIFLOWANALYSISMSMPC_CXX

#include "TFile.h"
#include "TList.h"
#include "TMath.h"
#include "TString.h"
#include "TProfile.h"
#include "TH1D.h"
#include "TMatrixD.h"

#include "AliFlowCommonConstants.h"
#include "AliFlowEventSimple.h"
#include "AliFlowTrackSimple.h"
#include "AliFlowCommonHist.h"
#include "AliFlowAnalysisMSMPC.h"

using std::cout;
using std::endl;


ClassImp(AliFlowAnalysisMSMPC)

//-----------------------------------------------------------------------
AliFlowAnalysisMSMPC::AliFlowAnalysisMSMPC():
fHistList(NULL),
fNSubsets(1),
fNCorrelators(1)
{
  for(Int_t i=0; i<fkMaxNSub; i++) {
    fReQ[i] = NULL;
    fImQ[i] = NULL;
    fChargeSub[i] = 0;
    fEtaMinSub[i] = -10.;
    fEtaMaxSub[i] = 10.;
  }
  fAreSubsetsDisjoint = new TMatrixD(fkMaxNSub,fkMaxNSub);
  for(Int_t i=0; i<fkMaxNSub; i++) {
    for(Int_t j=0; j<fkMaxNSub; j++) {
      (*fAreSubsetsDisjoint)(i,j) = 1.;
    }
  }
  for(Int_t i=0; i<fkMaxNumCorr; i++) {
    fProfileCorrelator[i] = NULL;
    fHistogramCorrelator[i] = NULL;
    fInfoCorrelators[i] = MSMPCorrelator();
  }
}

//-----------------------------------------------------------------------
AliFlowAnalysisMSMPC::~AliFlowAnalysisMSMPC()
{
  // destructor
  delete fHistList;
  delete fAreSubsetsDisjoint;
}

//-----------------------------------------------------------------------
void AliFlowAnalysisMSMPC::Init()
{
  printf("****************************************************************\n");
  printf(" initialize analysis with multi-subset multi-particle correlators \n");

  // define all histograms
  Bool_t oldHistAddStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  fHistList = new TList();
  fHistList->SetName("cobjGF");
  fHistList->SetOwner();

  for(Int_t i=0; i<fNSubsets; i++) {
    fReQ[i] = new TMatrixD(21,9);
    fImQ[i] = new TMatrixD(21,9);
  }

  printf("\n subsets defined: %d \n",fNSubsets);
  for(Int_t i=0; i<fNSubsets; i++) {
    printf("%d) charge: %s, eta-range: (%.1f, %.1f) \n",i+1,(fChargeSub[i]==0?"both":(fChargeSub[i]>0?"pos":"neg")),fEtaMinSub[i],fEtaMaxSub[i]);
  }
  Int_t oversub=0;
  printf("\n overlapping subsets: \n");
  for(Int_t i=0; i<fNSubsets; i++) {
    for(Int_t j=i; j<fNSubsets; j++) {
      if(!(Int_t)(*fAreSubsetsDisjoint)(i,j)) {printf(" (%d, %d)\n",i+1,j+1); oversub++;}
    }
  }
  if(!oversub) printf(" none \n");

  printf("\n correlators defined: %d \n",fNCorrelators);
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
      NameCorrelator += Form("%d",fInfoCorrelators[i].GetSubset(k));
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
void AliFlowAnalysisMSMPC::SetInfoCorrelator(Int_t num, Int_t ord, TArrayI& har, TArrayI& sub)
{
  // set information on correlator (harmonics, subsets)
  Int_t ind = num-1; // position in array
  if(ind<0 || ind>fNCorrelators) {
    printf("error: correlator %d out-of-range, set correct number of correlators first \n",num);
    exit(0);
  }
  fInfoCorrelators[ind] = MSMPCorrelator(ord, har, sub);
}

//-----------------------------------------------------------------------
void AliFlowAnalysisMSMPC::DefineSubset(Int_t num, Int_t charge, Double_t etamin, Double_t etamax)
{
  // define subset (particle charge and pseudorapidity)
  Int_t ind = num-1; // position in array
  if(ind<0 || ind>fNSubsets) {
    printf("error: subset %d out-of-range, set correct number of subsets first \n",num);
    exit(0);
  }
  fChargeSub[ind] = charge;
  fEtaMinSub[ind] = etamin;
  fEtaMaxSub[ind] = etamax;
}

//-----------------------------------------------------------------------
void AliFlowAnalysisMSMPC::SetOverlappingSubsets(Int_t sub1, Int_t sub2)
{
  // set overlapping subsets
  Int_t ind1 = sub1-1, ind2 = sub2-1; // position in array
  if((ind1<0 || ind1>fNSubsets) || (ind2<0 || ind2>fNSubsets)) {
    printf("error: subset %d and/or %d out-of-range, set correct number of subsets first \n",sub1,sub2);
    exit(0);
  }
  (*fAreSubsetsDisjoint)(ind1,ind2) = 0.;
  (*fAreSubsetsDisjoint)(ind2,ind1) = 0.;
}

//-----------------------------------------------------------------------
void AliFlowAnalysisMSMPC::Make(AliFlowEventSimple* anEvent)
{
  // process event
  if (!anEvent) return; // for coverity

  // loop over the tracks of the event
  AliFlowTrackSimple*   pTrack = NULL;
  Int_t iNumberOfTracks = anEvent->NumberOfTracks();

  for (Int_t i=0;i<iNumberOfTracks;i++) {
    pTrack = anEvent->GetTrack(i) ;
    if (!pTrack) continue;
    Double_t dPhi = pTrack->Phi();
    Double_t dPt  = pTrack->Pt();
    Double_t dEta = pTrack->Eta();
    // determine track weight (TBI)
    Double_t wTrack = 1.;

    for (Int_t iSub=0; iSub<fNSubsets; iSub++) {
      if(IsTrackInSubset(pTrack,iSub)) {
        for(Int_t m=0;m<21;m++) {
          for(Int_t k=0;k<9;k++) {
            (*fReQ[iSub])(m,k) += pow(wTrack,k)*TMath::Cos(m*dPhi);
            (*fImQ[iSub])(m,k) += pow(wTrack,k)*TMath::Sin(m*dPhi);
          }
        }
      }
    } // end of for (Int_t sub=0; sub<fNSubsets; sub++)

  }// end of loop over tracks

  // calculate correlators event-by-event
  this->CalculateFlowMSMPC();

  // reset all event-by-event quantities
  this->ResetEventByEventQuantities();
}

//--------------------------------------------------------------------
void AliFlowAnalysisMSMPC::CalculateFlowMSMPC()
{
  // calculate correlators event-by-event
  for(Int_t i=0; i<fNCorrelators; i++) {

    // get info on correlator to calculate
    Int_t CorOrder = fInfoCorrelators[i].GetOrder();           // order: how many particles to correlate
    TArrayI CorHarmonics(fInfoCorrelators[i].GetHarmonics());  // harmonics
    TArrayI CorSubsets(fInfoCorrelators[i].GetSubsets());      // subsets

    // check that there are enough particles in the event
    Int_t iMult = 0;
    for(Int_t k=0; k<CorOrder; k++) {
      iMult += (*fReQ[CorSubsets[k]-1])(0,0);
    }
    if(iMult<CorOrder) continue;

    TArrayI CorNullArray(CorOrder);
    for(Int_t k=0; k<CorOrder; k++) CorNullArray[k] = 0;

    // compute correlator
    std::complex<double> N = this->ucN(CorOrder, CorHarmonics, CorSubsets); // numerator
    std::complex<double> D = this->ucN(CorOrder, CorNullArray, CorSubsets); // denominator

    // fill fProfileCorrelator, use denominator (number of combinations) as weight
    if(D.real()>0.) {
      fProfileCorrelator[i]->Fill(0.5,N.real()/D.real(),D.real());
    }

  } // end of for(Int_t i=0; i<fNCorrelators; i++)
}

//--------------------------------------------------------------------
std::complex<double> AliFlowAnalysisMSMPC::ucN(const Int_t n, const TArrayI& harmonics, const TArrayI& subsets)
{
  // recursive algorithm to calculate correlators (1/2)
  TArrayI counts(n);
  for (Int_t i = 0; i < n; i++) {
    counts[i] = 1;
  }
  TArrayI hh(harmonics);

  return ucN2(n, hh, counts, subsets);
}

//--------------------------------------------------------------------
std::complex<double> AliFlowAnalysisMSMPC::ucN2(const Int_t n, TArrayI& h, TArrayI& cnt, const TArrayI& se)
{
  // recursive algorithm to calculate correlators (2/2)
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
    if (se[i]==se[j] || !(Int_t)((*fAreSubsetsDisjoint)(se[i]-1,se[j]-1))) {
      c -= factor * ucN2(j, h, cnt, se);
    }
    cnt[i]--;
    h[i] -= h[j];
  }
  return c;
}

//--------------------------------------------------------------------
Bool_t AliFlowAnalysisMSMPC::IsTrackInSubset(AliFlowTrackSimple *pTrack, Int_t iSub)
{
  // check if track is in subset
  Bool_t IsInSubset = kTRUE;
  if((fChargeSub[iSub] != 0) && (pTrack->Charge() != fChargeSub[iSub])) IsInSubset = kFALSE;
  if((pTrack->Eta() < fEtaMinSub[iSub]) || (pTrack->Eta() > fEtaMaxSub[iSub])) IsInSubset = kFALSE;
  return IsInSubset;
}

//--------------------------------------------------------------------
void AliFlowAnalysisMSMPC::GetOutputHistograms(TList *outputListHistos)
{
  // get pointers to all output histograms (called before Finish())
  fHistList = outputListHistos;
}

//--------------------------------------------------------------------
void AliFlowAnalysisMSMPC::ResetEventByEventQuantities()
{
  // reset Q-vectors (called at the end of Make())
  for(Int_t i=0; i<fNSubsets; i++) {
    fReQ[i]->Zero();
    fImQ[i]->Zero();
  }
}

//--------------------------------------------------------------------
void AliFlowAnalysisMSMPC::Finish()
{
  // calculate flow with correct statistical uncertainties
  printf("**********************************************\n");
  printf(" finalize analysis with multi-subset multi-particle correlators \n");

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
}

//-----------------------------------------------------------------------
void AliFlowAnalysisMSMPC::SetNCorrelators(Int_t num)
{
  // set number of correlators
  if(num>fkMaxNumCorr) {
    printf("Error: maximum number of correlators is %d \n",fkMaxNumCorr);
    exit(0);
  }
  fNCorrelators = num;
}

//-----------------------------------------------------------------------
void AliFlowAnalysisMSMPC::SetNSubsets(Int_t num)
{
  // set number of subsets
  if(num>fkMaxNSub) {
    printf("Error: maximum number of subsets is %d \n",fkMaxNSub);
    exit(0);
  }
  fNSubsets = num;
}

//-----------------------------------------------------------------------
void AliFlowAnalysisMSMPC::WriteHistograms(TDirectoryFile *outputFileName) const
{
  // store the final results in output .root file
  outputFileName->Add(fHistList);
  outputFileName->Write(outputFileName->GetName(), TObject::kSingleKey);
}

//--------------------------------------------------------------------
