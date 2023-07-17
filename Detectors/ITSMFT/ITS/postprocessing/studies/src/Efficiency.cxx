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
#include <DataFormatsITSMFT/TopologyDictionary.h>
#include "DetectorsBase/GRPGeomHelper.h"
#include "DataFormatsITSMFT/ROFRecord.h"

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
using namespace std;

class EfficiencyStudy : public Task
{
 public:
  EfficiencyStudy(shared_ptr<DataRequest> dr, shared_ptr<o2::base::GRPGeomRequest> gr, bool mc) : mDataRequest(dr), mGGCCDBRequest(gr), useMC(mc) {};
  ~EfficiencyStudy() final = default;
  void init(InitContext& ic) final;
  void run(ProcessingContext&) final;
  void endOfStream(EndOfStreamContext&) final;
  void finaliseCCDB(ConcreteDataMatcher&, void*) final;
  void process(o2::globaltracking::RecoContainer&);
  void setClusterDictionary(const o2::itsmft::TopologyDictionary* d) { dict = d; }

 private:
  void updateTimeDependentParams(ProcessingContext& pc);

  // helper functions
  void loadROFrameData(gsl::span<const o2::itsmft::ROFRecord> rofs, gsl::span<const itsmft::CompClusterExt> clusters,
                        gsl::span<const unsigned char>::iterator& pattIt);
  void printHistograms();
  void saveHistograms();

  // input
  shared_ptr<DataRequest> mDataRequest;
  shared_ptr<o2::base::GRPGeomRequest> mGGCCDBRequest;
  bool useMC;

  unique_ptr<o2::steer::MCKinematicsReader> MCReader;

  // internal
  const o2::itsmft::TopologyDictionary* dict;
  vector<o2::MCTrack> MCTracks;

  // output
  static const int nlayers = 7;
  TH2F* h2d_cluster_z_phi[nlayers] = {};

  const string outFile{"o2standalone_efficiency_study.root"};
};

void EfficiencyStudy::init(InitContext& ic)
{
  o2::base::GRPGeomHelper::instance().setRequest(mGGCCDBRequest);
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
  int zbin = 25;
  double z_length[nlayers] = {20, 20, 20, 60, 60, 80, 80};
  int phibin = 25;
  double phi_lo = 0;
  double phi_hi = TMath::Pi();

  for (int ilayer = 0; ilayer < nlayers; ilayer++)
  {
    int z_bound = z_length[ilayer] / 2;
    h2d_cluster_z_phi[ilayer] = new TH2F(Form("h2d_cluster_z_phi_%d", ilayer), Form("h2d_cluster_z_phi_%d", ilayer),
                                        zbin, -z_length[ilayer], z_length[ilayer], phibin, phi_lo, phi_hi);
  }

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
  auto rofs = recoData.getITSClustersROFRecords();
  auto clusters = recoData.getITSClusters();

  auto pattIt_span = recoData.getITSClustersPatterns();
  gsl::span<const unsigned char>::iterator pattIt = pattIt_span.begin();

  auto tracks = recoData.getITSTracks();

  if (useMC) cout<<"USING MC"<<endl;
  else cout<<"NOT USING MC"<<endl;

  LOGP(info, "Found {} ITS tracks", tracks.size());
  LOGP(info, "Found {} ITS clusters", clusters.size());
  LOGP(info, "Found {} ROFs.", recoData.getITSTracksROFRecords().size());

  // fill 2d z-phi cluster distribution, INCLUSIVE all reco clusters
  loadROFrameData(rofs, clusters, pattIt);

}

// adapted from AliceO2/Detectors/ITSMFT/ITS/tracking/src/TimeFrame.cxx
void EfficiencyStudy::loadROFrameData(gsl::span<const o2::itsmft::ROFRecord> rofs,
                               gsl::span<const itsmft::CompClusterExt> clusters,
                               gsl::span<const unsigned char>::iterator& pattIt)
{
  o2::its::GeometryTGeo* geom = o2::its::GeometryTGeo::Instance();
  geom->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L, o2::math_utils::TransformType::L2G));

  for (auto& rof : rofs) {
    for (int clusterId{rof.getFirstEntry()}; clusterId < rof.getFirstEntry() + rof.getNEntries(); ++clusterId) {
      auto& c = clusters[clusterId];

      int layer = geom->getLayer(c.getSensorID());

      auto pattID = c.getPatternID();
      o2::math_utils::Point3D<float> locXYZ;
      if (pattID != itsmft::CompCluster::InvalidPatternID) {
        if (!dict->isGroup(pattID)) {
          locXYZ = dict->getClusterCoordinates(c);
        } else {
          o2::itsmft::ClusterPattern patt(pattIt);
          locXYZ = dict->getClusterCoordinates(c, patt);
        }
      } else {
        o2::itsmft::ClusterPattern patt(pattIt);
        locXYZ = dict->getClusterCoordinates(c, patt, false);
      }
      auto sensorID = c.getSensorID();
      // Inverse transformation to the local --> tracking
      auto trkXYZ = geom->getMatrixT2L(sensorID) ^ locXYZ;
      // Transformation to the local --> global
      auto gloXYZ = geom->getMatrixL2G(sensorID) * locXYZ;

      // fill histograms
      cout<<"CLUSTER ID "<<clusterId<<" layer "<<layer<<" z "<<gloXYZ.Z()<<" phi "<<gloXYZ.Phi()<<endl;
      h2d_cluster_z_phi[layer]->Fill(gloXYZ.Z(), gloXYZ.Phi());
    }

  }
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

  // cluster z-eta distributions, inclusive
  Double_t plot_zrange_lo = 0;
  Double_t plot_zrange_hi = 250;

  for (int ilayer = 0; ilayer < nlayers; ilayer++)
  {
    h2d_cluster_z_phi[ilayer]->GetXaxis()->SetTitle("z");
    h2d_cluster_z_phi[ilayer]->GetYaxis()->SetTitle("#phi");
    // h2d_cluster_z_phi->GetZaxis()->SetRangeUser(plot_zrange_lo, plot_zrange_hi);

    h2d_cluster_z_phi[ilayer]->Draw("colz");
    c->Print(Form("./h2d_cluster_z_phi_%d.pdf", ilayer));
  }

}

void EfficiencyStudy::saveHistograms()
{
  TFile* fout = new TFile(outFile.c_str(), "UPDATE");

  for (int ilayer = 0; ilayer < nlayers; ilayer++)
  {
    h2d_cluster_z_phi[ilayer]->Write();
  }

  fout->Close();
  LOGP(important, "Stored histograms into {}", outFile.c_str());
}

void EfficiencyStudy::updateTimeDependentParams(ProcessingContext& pc)
{
  o2::base::GRPGeomHelper::instance().checkUpdates(pc);
  static bool initOnceDone = false;
  if (!initOnceDone) { // this params need to be queried only once
    initOnceDone = true;
    o2::its::GeometryTGeo* geom = o2::its::GeometryTGeo::Instance();
    geom->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L, o2::math_utils::TransformType::T2GRot, o2::math_utils::TransformType::T2G));
  }
}

void EfficiencyStudy::endOfStream(EndOfStreamContext& ec)
{
  printHistograms();
  saveHistograms();
}

void EfficiencyStudy::finaliseCCDB(ConcreteDataMatcher& matcher, void* obj)
{
  if (o2::base::GRPGeomHelper::instance().finaliseCCDB(matcher, obj)) {
    return;
  }
  if (matcher == ConcreteDataMatcher("ITS", "CLUSDICT", 0)) {
    setClusterDictionary((const o2::itsmft::TopologyDictionary*)obj);
    return;
  }
}

DataProcessorSpec getEfficiencyStudy(mask_t srcTracksMask, mask_t srcClustersMask, bool useMC)
{
  vector<OutputSpec> outputs;
  auto dataRequest = make_shared<DataRequest>();
  dataRequest->requestTracks(srcTracksMask, useMC);
  dataRequest->requestClusters(srcClustersMask, useMC);
  // dataRequest->requestPrimaryVertertices(useMC);

  auto ggRequest = std::make_shared<o2::base::GRPGeomRequest>(false,                             // orbitResetTime
                                                              true,                              // GRPECS=true
                                                              false,                             // GRPLHCIF
                                                              false,                             // GRPMagField
                                                              false,                             // askMatLUT
                                                              o2::base::GRPGeomRequest::Aligned, // geometry
                                                              dataRequest->inputs,
                                                              true);

  return DataProcessorSpec{
    "its-study-Efficiency",
    dataRequest->inputs,
    outputs,
    AlgorithmSpec{adaptFromTask<EfficiencyStudy>(dataRequest, ggRequest, useMC)},
    Options{}};
}
} // namespace study
} // namespace its
} // namespace o2