// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

// Skeleton derived from RS's code in ITSOffStudy

#include "ITSStudies/Efficiency.h"

#include "DataFormatsITS/TrackITS.h"
#include "Framework/Task.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "DataFormatsGlobalTracking/RecoContainer.h"
#include "ReconstructionDataFormats/V0.h"
#include "ReconstructionDataFormats/PrimaryVertex.h"
#include "CommonUtils/TreeStreamRedirector.h"
#include "DataFormatsParameters/GRPObject.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "ITSBase/GeometryTGeo.h"
#include <SimulationDataFormat/MCTrack.h>
#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/DCA.h"
#include <Steer/MCKinematicsReader.h>

#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TLine.h>
#include <numeric>
#include <TStyle.h>
#include <TDatabasePDG.h>

namespace o2
{
namespace its
{
namespace study
{
using namespace o2::framework;
using namespace o2::globaltracking;
using namespace o2::its;
using namespace std;

using GTrackID = o2::dataformats::GlobalTrackID;
using MCLabel = o2::MCCompLabel;

class EfficiencyStudy : public Task
{
 public:
  EfficiencyStudy(shared_ptr<DataRequest> dr, mask_t src, bool mc) : mDataRequest(dr), mTracksSrc(src), useMC(mc) {};
  ~EfficiencyStudy() final = default;
  void init(InitContext& ic) final;
  void run(ProcessingContext&) final;
  void endOfStream(EndOfStreamContext&) final;
  void finaliseCCDB(ConcreteDataMatcher&, void*) final;
  void process(o2::globaltracking::RecoContainer&);

 private:
  void updateTimeDependentParams(ProcessingContext& pc);

  // helper functions
  void printHistograms();
  void saveHistograms();

  // input
  shared_ptr<DataRequest> mDataRequest;
  GTrackID::mask_t mTracksSrc{};
  bool useMC;

  unique_ptr<o2::steer::MCKinematicsReader> MCReader;

  // internal
  vector<o2::MCTrack> MCTracks;

  // output
  TH2F* h2d_eta_phi_reco;
  TH2F* h2d_eta_phi_truth;
  TH2F* h2d_eta_err;
  TH2F* h2d_phi_err;

  const string outFile{"o2standalone_efficiency_study.root"};
};

void EfficiencyStudy::init(InitContext& ic)
{
  LOGP(info, "Starting track reconstruction efficiency study...");

  // read in MC truth tracks
  if (useMC) {
    MCReader = make_unique<o2::steer::MCKinematicsReader>("collisioncontext.root");
    for (int iEvent{0}; iEvent < MCReader->getNEvents(0); iEvent++) {
      auto mctrk = MCReader->getTracks(0, iEvent);
      MCTracks.insert(MCTracks.end(), mctrk.begin(), mctrk.end());
    }
  }

  // prepare output
  int etabin = 50;
  double eta_lo = -1.5;
  double eta_hi = 1.5;
  int phibin = 50;
  double phi_lo = 0;
  double phi_hi = 2*TMath::Pi();

  h2d_eta_phi_reco = new TH2F("h2d_eta_phi_reco", "h2d_eta_phi_reco", etabin, eta_lo, eta_hi, phibin, phi_lo, phi_hi);
  h2d_eta_phi_truth = new TH2F("h2d_eta_phi_truth", "h2d_eta_phi_truth", etabin, eta_lo, eta_hi, phibin, phi_lo, phi_hi);
  h2d_eta_err = new TH2F("h2d_eta_err", "h2d_eta_err", etabin, eta_lo, eta_hi, phibin, phi_lo, phi_hi);
  h2d_phi_err = new TH2F("h2d_phi_err", "h2d_phi_err", etabin, eta_lo, eta_hi, phibin, phi_lo, phi_hi);

  LOGP(info, "Efficiency study initialized.");
}

void EfficiencyStudy::run(ProcessingContext& pc)
{
  o2::globaltracking::RecoContainer recoData;
  recoData.collectData(pc, *mDataRequest.get());
  updateTimeDependentParams(pc); // Make sure this is called after recoData.collectData, which may load some conditions
  process(recoData);
}

void EfficiencyStudy::process(o2::globaltracking::RecoContainer& recoData)
{
  TDatabasePDG* pdg_db = TDatabasePDG::Instance();

  auto itsTracks = recoData.getITSTracks();
  auto itsClusters = recoData.getITSClusters();

  gsl::span<const MCLabel> mcLabels;

  mcLabels = recoData.getITSTracksMCLabels();
  MCLabel trackMCLabel;

  if (useMC) cout<<"USING MC"<<endl;
  else cout<<"NOT USING MC"<<endl;

  LOGP(info, "Found {} ITS tracks", itsTracks.size());
  LOGP(info, "Found {} ITS clusters", itsClusters.size());
  LOGP(info, "Found {} labels", mcLabels.size());

  const o2::MCTrack* truth_track;

  // comparing reco and MC tracks one-to-one phi values
  TrackITS reco_track;
  for (int i = 0; i < itsTracks.size(); i++)
  {
    if (i % 10000 == 0) LOGP(info, "Processing track {}", i);

    reco_track = itsTracks[i];

    Double_t reco_eta = reco_track.getEta();
    Double_t reco_phi = reco_track.getPhi();

    h2d_eta_phi_reco->Fill(reco_eta, reco_phi);

    trackMCLabel = mcLabels[i];
    if (trackMCLabel.isValid())
    {
      truth_track = MCReader->getTrack(trackMCLabel);

      Double_t truth_eta = truth_track->GetEta();
      Double_t truth_phi = truth_track->GetPhi();
      h2d_eta_phi_truth->Fill(truth_eta, truth_phi);

      Double_t eta_err = (reco_eta - truth_eta) / truth_eta;
      Double_t phi_err = (reco_phi - truth_phi) / truth_phi;
      h2d_eta_err->Fill(truth_eta, truth_phi, eta_err);
      h2d_phi_err->Fill(truth_eta, truth_phi, phi_err);

    }

  }

  // printing out all MC tracks (including those that don't pass reco)
  /*
  for (int i = 0; i < MCTracks.size(); i++)
  {
    mcTrack = (MCTrack*) &MCTracks[i];

    int pdg_code = mcTrack->GetPdgCode();
    Double_t charge = pdg_db->GetParticle(pdg_code)->Charge();

    // if is charged final state particle and in eta, pt acceptance range of ITS
    if (charge != 0
          && abs(mcTrack->GetEta()) < 0.9
          && mcTrack->GetPt() > 0.15)
    {
      cout<<"MC TRACK "<<i<<" with truth phi "<<mcTrack->GetPhi()<<endl;
    }

  }
  */


}

void EfficiencyStudy::printHistograms()
{
  // setup standard 2D hist canvas
  TCanvas* c = new TCanvas("temp", "temp", 0, 0, 400, 400);
  gStyle->SetOptStat(0); gStyle->SetOptTitle(0); gStyle->SetMarkerSize(1.6);
  c->cd();
  TPad* mpad = new TPad("temp","temp",0.02,0.02,0.99,0.99,0,0,0);
  mpad->SetTickx();
  mpad->SetTicky();
  mpad->SetTopMargin(0.1);
  mpad->SetBottomMargin(0.13);
  mpad->SetLeftMargin(0.12);
  mpad->SetRightMargin(0.15);
  mpad->Draw();
  mpad->cd();
  c->Modified();
  c->Update();

  // phi-eta distributions for truth and reco tracks
  Double_t plot_zrange_lo = 0;
  Double_t plot_zrange_hi = 250;

  // reco
  h2d_eta_phi_reco->GetXaxis()->SetTitle("#eta");
  h2d_eta_phi_reco->GetYaxis()->SetTitle("#phi");
  h2d_eta_phi_reco->GetZaxis()->SetRangeUser(plot_zrange_lo, plot_zrange_hi);

  h2d_eta_phi_reco->Draw("colz");
  c->Print("./h2d_eta_phi_reco.pdf");

  //truth
  h2d_eta_phi_truth->GetXaxis()->SetTitle("#eta");
  h2d_eta_phi_truth->GetYaxis()->SetTitle("#phi");
  h2d_eta_phi_truth->GetZaxis()->SetRangeUser(plot_zrange_lo, plot_zrange_hi);

  h2d_eta_phi_truth->Draw("colz");
  c->Print("./h2d_eta_phi_truth.pdf");

  // phi and eta relative error vs phi and eta (each)
  plot_zrange_lo = 0;
  plot_zrange_hi = 250;

  // eta error
  h2d_eta_err->GetXaxis()->SetTitle("#eta");
  h2d_eta_err->GetYaxis()->SetTitle("#phi");

  h2d_eta_err->Draw("colz");
  c->Print("./h2d_eta_err.pdf");

  // phi error
  h2d_phi_err->GetXaxis()->SetTitle("#eta");
  h2d_phi_err->GetYaxis()->SetTitle("#phi");

  h2d_phi_err->Draw("colz");
  c->Print("./h2d_phi_err.pdf");


}

void EfficiencyStudy::saveHistograms()
{
  TFile* fout = new TFile(outFile.c_str(), "UPDATE");

  h2d_eta_phi_reco->Write();
  h2d_eta_phi_truth->Write();
  h2d_eta_err->Write();
  h2d_phi_err->Write();

  fout->Close();
  LOGP(important, "Stored histograms into {}", outFile.c_str());
}

void EfficiencyStudy::updateTimeDependentParams(ProcessingContext& pc)
{
  static bool initOnceDone = false;
  if (!initOnceDone) { // this params need to be queried only once
    initOnceDone = true;
  }
}

void EfficiencyStudy::endOfStream(EndOfStreamContext& ec)
{
  printHistograms();
  saveHistograms();
}

void EfficiencyStudy::finaliseCCDB(ConcreteDataMatcher& matcher, void* obj)
{
}

DataProcessorSpec getEfficiencyStudy(mask_t srcTracksMask, mask_t srcClustersMask, bool useMC)
{
  vector<OutputSpec> outputs;
  auto dataRequest = make_shared<DataRequest>();
  dataRequest->requestTracks(srcTracksMask, useMC);
  dataRequest->requestClusters(srcClustersMask, useMC);
  // dataRequest->requestPrimaryVertertices(useMC);

  return DataProcessorSpec{
    "its-study-Efficiency",
    dataRequest->inputs,
    outputs,
    AlgorithmSpec{adaptFromTask<EfficiencyStudy>(dataRequest, srcTracksMask, useMC)},
    Options{}};
}
} // namespace study
} // namespace its
} // namespace o2