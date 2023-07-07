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
#include <TFile.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TLine.h>
#include <numeric>
#include <TStyle.h>

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

  shared_ptr<DataRequest> mDataRequest;
  GTrackID::mask_t mTracksSrc{};
  bool useMC;

  unique_ptr<o2::steer::MCKinematicsReader> MCReader;
  vector<o2::MCTrack> MCTracks;
};

void EfficiencyStudy::init(InitContext& ic)
{
  LOGP(info, "Starting track reconstruction efficiency study...");

  // read in MC truth tracks

  if (useMC) { // for counting the missed K0shorts
    MCReader = make_unique<o2::steer::MCKinematicsReader>("collisioncontext.root");
    for (int iEvent{0}; iEvent < MCReader->getNEvents(0); iEvent++) {
      auto mctrk = MCReader->getTracks(0, iEvent);
      MCTracks.insert(MCTracks.end(), mctrk.begin(), mctrk.end());
    }
  }

  // prepare output

  // TODO
  

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
  auto itsTracks = recoData.getITSTracks();
  auto itsClusters = recoData.getITSClusters();

  gsl::span<const MCLabel> mcLabels;
  mcLabels = recoData.getITSTracksMCLabels();

  cout<<"number of tracks = "<<itsTracks.size()<<endl;
  cout<<"number of clusters = "<<itsClusters.size()<<endl;

  TrackITS this_track;
  for (int i = 0; i < itsTracks.size(); i++)
  {
    this_track = itsTracks[i]; 
    cout<<"TRACK "<<i<<" phi = "<<this_track.getPhi()<<endl;

  }


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