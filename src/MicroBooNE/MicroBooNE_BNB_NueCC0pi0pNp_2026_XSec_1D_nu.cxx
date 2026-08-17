// Copyright 2016 L. Pickering, P Stowell, R. Terri, C. Wilkinson, C. Wret

/*******************************************************************************
 *    This file is part of NUISANCE.
 *
 *    NUISANCE is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    NUISANCE is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with NUISANCE.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/

#include "MicroBooNE_BNB_NueCC0pi0pNp_2026_XSec_1D_nu.h"

#include "TH1D.h"
#include "TH2D.h"

//********************************************************************
MicroBooNE_BNB_NueCC0pi0pNp_2026_XSec_1D_nu::MicroBooNE_BNB_NueCC0pi0pNp_2026_XSec_1D_nu(nuiskey samplekey) {
  //********************************************************************
  fSettings = LoadSampleSettings(samplekey);

  std::string name = fSettings.GetS("name");
  std::string objSuffix;

  // work out which sample you need, and set axii
  if (!name.compare("MicroBooNE_BNB_NueCC0pi0pNp_2026_XSec_1DEe_nu")) {
    fDist = kEe;
    objSuffix = "Ee";
    fSettings.SetXTitle("E_{e} (GeV)");
    fSettings.SetYTitle("d#sigma/dE_{e} (cm^{2}/GeV/nucleon)");
    is0p = false;
  }
  else if (!name.compare("MicroBooNE_BNB_NueCC0pi0pNp_2026_XSec_1Dcosthe_nu")) {
    fDist = kcos_the;
    objSuffix = "cos_the";
    fSettings.SetXTitle("cos#theta_{e} ");
    fSettings.SetYTitle("d#sigma/dcos#theta_{e} (cm^{2}/nucleon)"); 
    is0p = false;
  }
  else if (!name.compare("MicroBooNE_BNB_NueCC0pi0pNp_2026_XSec_1Dcosthp_nu")) {
    fDist = kcos_thp;
    objSuffix = "cos_thp";
    fSettings.SetXTitle("cos#theta_{P^{lead}}");
    fSettings.SetYTitle("d#sigma/dcos#theta_{P^{lead}} (cm^{2}/nucleon)"); 
    is0p = false;
  }
  else if (!name.compare("MicroBooNE_BNB_NueCC0pi0pNp_2026_XSec_1Dcosthep_nu")) {
    fDist = kcos_thep;
    objSuffix = "cos_thep";
    fSettings.SetXTitle("cos#theta_{eP}");
    fSettings.SetYTitle("d#sigma/dcos#theta_{eP} (cm^{2}/nucleon)");
    is0p = false;
  }  
  else if (!name.compare("MicroBooNE_BNB_NueCC0pi0pNp_2026_XSec_1DKEp_nu")) {
    fDist = kKEp;
    objSuffix = "KEp";
    fSettings.SetXTitle("KE_{P^{lead}} (GeV)");
    fSettings.SetYTitle("d#sigma/dKE_{P^{lead}} (cm^{2}/GeV/nucleon)");
    is0p = true;
  }
  else if (!name.compare("MicroBooNE_BNB_NueCC0pi0pNp_2026_XSec_1D2Bin_nu")) {
    fDist = k2Bin;
    objSuffix = "2Bin";
    fSettings.SetXTitle("1e0#pi0p vs 1e0#piNp");
    fSettings.SetYTitle("#sigma (cm^{2}/nucleon)");
    is0p = true;
  }
  else {
    assert(false);
  }

  // Sample overview ---------------------------------------------------
  std::string descrip = name + " sample.\n"
                               "Target: Ar\n"
                               "Flux: BNB FHC Nue\n"
                               "Signal: CC0pi0pNp \n"
                               "Contact: microboone_info@fnal.gov\n"
                               "Reference: Phys. Rev. D 114, 032006\n"
                               "DOI: https://doi.org/10.1103/2zwr-h81t\n";

  fSettings.SetDescription(descrip);
  fSettings.SetTitle(name);
  fSettings.SetAllowedTypes("FULL,DIAG/FREE,SHAPE,FIX/SYSTCOV/STATCOV",
                            "FIX/FULL");
  fSettings.SetEnuRange(0.0, 6.8);
  fSettings.DefineAllowedTargets("Ar");
  fSettings.DefineAllowedSpecies("nue");
  FinaliseSampleSettings();
  
  // Load data ---------------------------------------------------------
  std::string inputFile = FitPar::GetDataBase() +
              "/MicroBooNE/BNB_NueCC0pi0pNp_2026/BNB_NueCC0pi0pNp_2026_data_release.root";
  SetDataFromRootFile(inputFile, "DataXSec_" + objSuffix);
  ScaleData(1E-39);

  // ScaleFactor for DiffXSec/cm2/Nucleon
  fScaleFactor = (GetEventHistogram()->Integral("width") *1E-38) / double(fNEvents) / TotalIntegratedFlux();

  SetCovarFromRootFile(inputFile, "NuisanceScaledCovarianceMatrix_" + objSuffix);

  // Load regularization matrix
  TFile* inputRootFile = TFile::Open(inputFile.c_str());
  assert(inputRootFile && inputRootFile->IsOpen());
  TH2D* hreg = (TH2D*)inputRootFile->Get(("RegularizationMatrix_" + objSuffix).c_str());
  assert(hreg);

  int nrows = hreg->GetNbinsX();
  int ncols = hreg->GetNbinsY();
  fRegularizationMatrix = new TMatrixD(nrows, ncols);
  for (int i=0; i<nrows; i++) {
    for (int j=0; j<ncols; j++) {
      (*fRegularizationMatrix)(i,j) = hreg->GetBinContent(i+1, j+1);
    }
  }

  inputRootFile->Close();

  // Final setup ------------------------------------------------------
  FinaliseMeasurement();
}

// check if event meets Np signal definition
bool MicroBooNE_BNB_NueCC0pi0pNp_2026_XSec_1D_nu::isNpSignal(FitEvent *event) {
  // Check CC nue
  if (!SignalDef::isCCINC(event, 12, EnuMin, EnuMax)) return false;

  // Veto events which don't have exactly 1 FS electron
  if (event->NumFSParticle(11) != 1 ) return false;

  // Veto events where the electron doesn't have KE > 30 MeV
  if (event->GetHMFSParticle(11)->KE() <= 30.0) return false;

  // Veto events with neutral pions of any energy
  if (event->NumFSParticle(111) != 0) return false;

  // Veto events with charged pions with KE > 40MeV
  if (event->NumFSParticle(211) != 0 && event->GetHMFSParticle(211)->KE() > 40.0) return false;
  if (event->NumFSParticle(-211) != 0 && event->GetHMFSParticle(-211)->KE() > 40.0) return false;

  // Np events must have at least 1 proton in FS with KE >= 50MeV and cos of opening angle > -0.9
  if (event->NumFSParticle(2212) == 0) return false;
  if (event->GetHMFSParticle(2212)->KE() < 50.0) return false;

  TLorentzVector Pe = event->GetHMFSParticle(11)->fP;
  TLorentzVector Pp = event->GetHMFSParticle(2212)->fP;
  double cos_opening_ep = cos(Pe.Vect().Angle(Pp.Vect()));
  if (cos_opening_ep <= -0.9) return false;

  // Events that pass selection are NueCC0piNp
  return true;
}

// check if event meets 0p signal definition
bool MicroBooNE_BNB_NueCC0pi0pNp_2026_XSec_1D_nu::is0pSignal(FitEvent *event) {
  // Check CC nue
  if (!SignalDef::isCCINC(event, 12, EnuMin, EnuMax)) return false;

  // Veto events which don't have exactly 1 FS electron
  if (event->NumFSParticle(11) != 1 ) return false;

  // Veto events where the electron doesn't have KE > 30 MeV
  if (event->GetHMFSParticle(11)->KE() <= 30.0) return false;

  // Veto events with neutral pions of any energy
  if (event->NumFSParticle(111) != 0) return false;

  // Veto events with charged pions with KE > 40MeV
  if (event->NumFSParticle(211) != 0 && event->GetHMFSParticle(211)->KE() > 40.0) return false;
  if (event->NumFSParticle(-211) != 0 && event->GetHMFSParticle(-211)->KE() > 40.0) return false;

  // 0p events either have no final state proton or the leading proton has KE < 50MeV
  if (event->NumFSParticle(2212) != 0 && event->GetHMFSParticle(2212)->KE() >=50.) return false;
    
  // additional phase space requirements on FS electron 
  // Electron energy > 0.5 GeV and cos th_e > 0.6
  TLorentzVector pe = event->GetHMFSParticle(11)->fP;
  if (pe.E() <= 500. || pe.CosTheta() < 0.6 ) return false;
                     
  // Events that pass selection are NueCC0pi0p
  return true;
}

// check if event is signal based on whether distribution is Np or 0p+Np
bool MicroBooNE_BNB_NueCC0pi0pNp_2026_XSec_1D_nu::isSignal(FitEvent *event) {
  if (is0p){
    if (!isNpSignal(event) && !is0pSignal(event)){ return false;}
    return true;
  }
  else { return isNpSignal(event); }

}

void MicroBooNE_BNB_NueCC0pi0pNp_2026_XSec_1D_nu::FillEventVariables(FitEvent *event) {  
  // check that event is signal
  if (!isSignal(event)) { return; }

  // Must have an electron in the final state
  if (event->NumFSParticle(11) == 0 ) { return; }

  // get highest momentum electron vector
  TLorentzVector ve = event->GetHMFSParticle(11)->fP;

  // Calculate electron variables
  double ElecEnergy = ve.E() / 1000.0; // electron energy in GeV
  double ElecCosTheta = ve.CosTheta(); // electron angle

  // initialize proton variables
  TLorentzVector vp;
  double ProtonCosTheta;
  double ProtonKE = 0.; // 0p events have proton KE set to 0 to fill lowest bin
  double OpeningAngle;

  // if Np signal, grab highest momentum proton vector
  if (isNpSignal(event)){ 
    vp = event->GetHMFSParticle(2212)->fP;

    // Calculate proton variables
    ProtonCosTheta = vp.CosTheta();
    ProtonKE = (vp.E() - vp.M()) / 1000.0; // in GeV 

    OpeningAngle = cos(ve.Vect().Angle(vp.Vect()));;

  }

    if (fDist == kEe) {
      fXVar = ElecEnergy;
    }
    else if (fDist == kcos_the) {
      fXVar = ElecCosTheta;
    }
    else if (fDist == kcos_thp) {
      fXVar = ProtonCosTheta;
    }
    else if (fDist == kcos_thep) {
      fXVar = OpeningAngle;
    }
    // full proton KE distribution
    else if (fDist == kKEp) {
      fXVar = ProtonKE;
    }
    // Np events fill 2nd bin for 2Bin KE 
    else if (fDist == k2Bin){
      fXVar = ProtonKE;
    }
}    

void MicroBooNE_BNB_NueCC0pi0pNp_2026_XSec_1D_nu::ConvertEventRates() {
  // standard conversion
  Measurement1D::ConvertEventRates();

  // Apply regularization matrix
  int nBins = fMCHist->GetNbinsX();

  // Convert to TVectorD
  TVectorD MC_beforeReg(nBins);
  for (int iBin=0; iBin<nBins; iBin++) {
      MC_beforeReg(iBin) = fMCHist->GetBinContent(iBin+1)*fMCHist->GetBinWidth(iBin+1);
  }
  
  // Apply regularization matrix
  TVectorD MC_afterReg = (*fRegularizationMatrix) * MC_beforeReg;
  
  // Copy results back to hist
  for (int iBin=0; iBin<nBins; iBin++) {
    if (fDist != k2Bin){
      fMCHist->SetBinContent(iBin+1, MC_afterReg(iBin)/fMCHist->GetBinWidth(iBin+1));
    } else{
        // 2Bin proton KE distribution doesn't get divided by bin width
        fMCHist->SetBinContent(iBin+1, MC_afterReg(iBin));
    }
  }  

}

