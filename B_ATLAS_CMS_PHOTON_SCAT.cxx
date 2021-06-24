//------------------------------------------------------------------------------
//
// BLUE: A ROOT class implementing the Best Linear Unbiased Estimate method.
// 
// Copyright (C) 2012-2019, Richard.Nisius@mpp.mpg.de
// All rights reserved
//
// This file is part of BLUE - Version 2.2.0.
//
// BLUE is free software: you can redistribute it and/or modify it under the 
// terms of the GNU Lesser General Public License as published by the Free 
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// For the licensing terms see the file COPYING or http://www.gnu.org/licenses.
//
//------------------------------------------------------------------------------
#include "Blue.h"

void B_ATLAS_CMS_PHOTON_SCAT(Int_t Flag = 0){

  //----------------------------------------------------------------------------
  // Flag steers which of the results should be calculated
  //  0: The LHC Combination but with the correlation rho for the W+jets = 1
  //     [Shows how to use SolveInfWei(), PrintInfWei() and PrintMatrix()]
  //  1: The individual variations of the correlations
  //----------------------------------------------------------------------------

  // The number of estimates, uncertainties and observables
  static const Int_t NumEst =  2;
  static const Int_t NumUnc =  8;
  static const Int_t NumObs =  1;

  // Set the names  
  TString NamEst[NumEst] = {"  ATLAS", "    CMS"};
  TString NamUnc[NumUnc] = {"   Stat", " Trigger", " PhotonReco", " PhotonAngRes"
			    " Electron", " Bkg", " Lumi", " Theory"};
  TString NamObs[NumObs] = {"     #sigma"};
  TString NamFil= "B_ATLAS_CMS_PHOTON_SCAT";

  // Preset according to Flag
  if(Flag == 0){
    printf("... B_ATLAS_CMS_PHOTON_SCAT: ------------------------------");
    printf("--------------------------\n");
    printf("... B_ATLAS_CMS_PHOTON_SCAT: The 2014 LHC Combination with");
    printf(" Rho(W+jets) = 1, Flag = %2i \n",Flag);
  }else if(Flag == 1){
    printf("... B_ATLAS_CMS_PHOTON_SCAT: ----------------------------- \n");
    printf("... B_ATLAS_CMS_PHOTON_SCAT: The rho variations, Flag = %2i \n",Flag);
  }else{
    printf("... ATLAS_CMS_PHOTON_SCAT: Not implemented Flag = %2i \n",Flag);
    return;
  }
  
  static const Int_t LenXEst = NumEst * (NumUnc+1);
  Double_t XEst[LenXEst] = {
    //    0  1  2  3  4  5  6  7  
    120, 17, 6, 6, 2, 0, 7, 4, 2,
    91, 36, 11, 16, 0, 5, 5, 0, 10
  };

  //                           0    1    2    3    4    5    6   7
  Double_t RhoVal[NumUnc] = {0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 1.0};

  //-- Local Structures for Blue output
  // TMatrices
  TMatrixD* LocRho = new TMatrixD(NumEst,NumEst);
  //-- End

  // Define formats for Figures and Latex file
  TString FilNam = "B_ATLAS_CMS_PHOTON_SCAT";
  const TString ForVal = "%5.3f";
  const TString ForUnc = ForVal;
  const TString ForWei = "%4.3f";
  const TString ForRho = "%4.2f";
  const TString ForPul = ForRho;
  const TString ForUni = "None";

  // Construct Object
  Blue *myBlue = new Blue(NumEst, NumUnc);
  myBlue->PrintStatus();
  myBlue->SetFormat(ForVal, ForUnc, ForWei, ForRho, ForPul, ForUni);

  // Fill names
  myBlue->FillNamEst(&NamEst[0]);
  myBlue->FillNamUnc(&NamUnc[0]);
  myBlue->FillNamObs(&NamObs[0]);

  // Fill estimates
  Int_t ind = 0;
  for(Int_t i = 0; i<NumEst; i++){
    myBlue->FillEst(i,&XEst[ind]);
    ind = ind + NumUnc + 1;
  }
  
  // Fill correlations
  for(Int_t k = 0; k<NumUnc; k++)myBlue->FillCor(k,RhoVal[k]);

  // Solve according to Flag
  if(Flag == 0){

    // The central result
    myBlue->FixInp();
    myBlue->SolveInfWei();

    printf("... B_ATLAS_CMS_PHOTON_SCAT: For the estimates and the result");
    printf(" see Table 1 \n");
    myBlue->PrintEst();
    myBlue->PrintResult();
    myBlue->LatexResult(FilNam,ForVal,ForVal,ForVal,ForVal,ForVal);

    printf("... B_ATLAS_CMS_PHOTON_SCAT: The correlation of the estimates \n");
    myBlue->GetRho(LocRho);
    LocRho->operator*=(100);
    myBlue->PrintMatrix(LocRho,"%+6.1f%%");

    printf("... B_ATLAS_CMS_PHOTON_SCAT: For the information weights and pulls");
    printf(" see Table 2 \n");
    myBlue->PrintInfWei();
    myBlue->PrintPull();

  }else if(Flag == 1){

    // Individual scan of correlations
    myBlue->FixInp();
    myBlue->SolveScaRho(0);
    myBlue->PrintScaRho(NamFil);

  }
  // Delete Object
  delete myBlue; myBlue = NULL;
  LocRho->Delete(); LocRho = NULL;
  return;
}
