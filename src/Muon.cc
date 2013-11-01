#include "Analysis/NtupleProducer/interface/NtupleProducer.h"

void NtupleProducer::DoMuonAnalysis(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using pat::Muon;
    using pat::MuonCollection;
    using namespace std;
    using namespace reco;
    using namespace edm;
    using namespace pat;

    (m->PreSelectedMuons).clear();

    edm::Handle< pat::TriggerEvent > triggerEvent;
    iEvent.getByLabel(triggerEvent_, triggerEvent);
    // PAT trigger helper for trigger matching information
    const pat::helper::TriggerMatchHelper matchHelper;


    // Get the B-field
    ESHandle<MagneticField> B;
    iSetup.get<IdealMagneticFieldRecord > ().get(B);
    const MagneticField* magField = B.product();

    //BeamSpot
    Handle<reco::BeamSpot> theBeamSpotHandle;
    iEvent.getByLabel(std::string("offlineBeamSpot"), theBeamSpotHandle);
    const reco::BeamSpot* beamSpot = theBeamSpotHandle.product();

    // Get the geometry
    ESHandle<GlobalTrackingGeometry> geomHandle;
    iSetup.get<GlobalTrackingGeometryRecord > ().get(geomHandle);



    Handle<pat::MuonCollection> muonsHandle;
    iEvent.getByLabel(PreSelectedMuonCollection_, muonsHandle);
    const MuonCollection & muons = *(muonsHandle.product());


  


    edm::Handle<edm::ValueMap<double> > isoAllMuMap;
    iEvent.getByLabel(edm::InputTag("muPFIsoValueChargedAll04"), isoAllMuMap);

    edm::Handle<edm::ValueMap<double> > isoChargedMuMap;
    iEvent.getByLabel(edm::InputTag("muPFIsoValueCharged04"), isoChargedMuMap);

    edm::Handle<edm::ValueMap<double> > isoNeutralMuMap;
    iEvent.getByLabel(edm::InputTag("muPFIsoValueNeutral04"), isoNeutralMuMap);

    edm::Handle<edm::ValueMap<double> > isoGammaMuMap;
    iEvent.getByLabel(edm::InputTag("muPFIsoValueGamma04"), isoGammaMuMap);

    edm::Handle<edm::ValueMap<double> > isoPUMuMap;
    iEvent.getByLabel(edm::InputTag("muPFIsoValuePU04"), isoPUMuMap);

//    edm::Handle<edm::ValueMap<double> > isoPULowMuMap;
//    iEvent.getByLabel(edm::InputTag("muPFIsoValuePULow"), isoPULowMuMap);

//    edm::Handle<edm::ValueMap<double> > isoAllMuMap;
//    iEvent.getByLabel(edm::InputTag("muPFIsoValueAll"), isoAllMuMap);
//
//    edm::Handle<edm::ValueMap<double> > isoChargedMuMap;
//    iEvent.getByLabel(edm::InputTag("muPFIsoValueCharged"), isoChargedMuMap);
//
//    edm::Handle<edm::ValueMap<double> > isoNeutralMuMap;
//    iEvent.getByLabel(edm::InputTag("muPFIsoValueNeutral"), isoNeutralMuMap);
//
//    edm::Handle<edm::ValueMap<double> > isoGammaMuMap;
//    iEvent.getByLabel(edm::InputTag("muPFIsoValueGamma"), isoGammaMuMap);
//
//    edm::Handle<edm::ValueMap<double> > isoPUMuMap;
//    iEvent.getByLabel(edm::InputTag("muPFIsoValuePU"), isoPUMuMap);
//
//    edm::Handle<edm::ValueMap<double> > isoPULowMuMap;
//    iEvent.getByLabel(edm::InputTag("muPFIsoValuePULow"), isoPULowMuMap);


    // Get primary vertex collection                                                                                                                                                                                                 
    Handle<reco::VertexCollection> recoPVCollection;
    iEvent.getByLabel(vertexCollectionForLeptonIP_, recoPVCollection);

    //  useBeamSpot_  = pset.getParameter<bool>("useBeamSpot");                                                                                                                                                                      
    bool useBeamSpot_ = true;

    //                                                                                                                                                                                                                               
    reco::Vertex primVertex;
    bool pvfound = (recoPVCollection->size() != 0);

    if (pvfound) {
      PrimaryVertexSorter pvs;
      vector<reco::Vertex> sortedList = pvs.sortedList(*(recoPVCollection.product()));
      primVertex = (sortedList.front());
    } else {
      //creating a dummy PV                                                                                                                                                                                                        
      reco::Vertex::Point p(0, 0, 0);
      reco::Vertex::Error e;
      e(0, 0) = 0.0015 * 0.0015;
      e(1, 1) = 0.0015 * 0.0015;
      e(2, 2) = 15. * 15.;
      primVertex = reco::Vertex(p, e, 1, 1, 1);
    }

    MuonCollection::const_iterator amuon;
    int qq = 0;
    for (amuon = muons.begin(); amuon != muons.end(); amuon++, ++qq) {
        double DepInTracker = amuon->isolationR03().sumPt;
        double DepInEcal = amuon->isolationR03().emEt;
        double DepInHcal = amuon->isolationR03().hadEt;
        myobject muo;
        muo.pt = amuon->pt();
        muo.eta = amuon->eta();
        muo.phi = amuon->phi();
        muo.px = amuon->px();
        muo.py = amuon->py();
        muo.pz = amuon->pz();
        muo.z = amuon->vz();

        muo.E = amuon->p();
        muo.Energy = amuon->energy();
        muo.mass = amuon->mass();
        muo.mt = amuon->mt();
        muo.et = amuon->et();
        muo.charge = amuon->charge();
        muo.dB = amuon->dB();
        muo.d0 = (amuon->track().isNonnull() ? amuon->track()->dxy(primVertex.position()) : 0.);
	muo.IP3D = amuon->dB(pat::Muon::PV3D);
        muo.dxy_PV = (amuon->track().isNonnull() ? amuon->track()->dxy(primVertex.position()) : 0.);
        muo.dz_PV = (amuon->track().isNonnull() ? amuon->track()->dz(primVertex.position()) : 0.);

        muo.DepositR03Ecal = DepInEcal;
        muo.DepositR03Hcal = DepInHcal;
        muo.DepositR03TrackerOfficial = DepInTracker;
        muo.GlobalMuonPromptTight = muon::isGoodMuon(*amuon, muon::GlobalMuonPromptTight);
        muo.TMOneStationLoose = muon::isGoodMuon(*amuon, muon::TMOneStationLoose);
        muo.TM2DCompatibilityLoose = muon::isGoodMuon(*amuon, muon::TM2DCompatibilityLoose);


        muo.normalizedChi2 = (amuon->globalTrack().isNonnull() ? amuon->globalTrack()->normalizedChi2() : 0);
        muo.numberOfValidMuonHits = (amuon->globalTrack().isNonnull() ? amuon->globalTrack()->hitPattern().numberOfValidMuonHits() : 0);
        muo.numberOfHits = (amuon->globalTrack().isNonnull() ? amuon->globalTrack()->hitPattern().numberOfHits() : 0);




        muo.trkLayerMeasure =(amuon->track().isNonnull() ? amuon->track()->hitPattern().trackerLayersWithMeasurement():0);
        muo.intrkLayerMeasure =(amuon->innerTrack().isNonnull() ? amuon->innerTrack()->hitPattern().trackerLayersWithMeasurement():0);
        muo.intrkLayerpixel =(amuon->innerTrack().isNonnull() ? amuon->innerTrack()->hitPattern().numberOfValidPixelHits():0);
        muo.dxy_in =(amuon->innerTrack().isNonnull() ? amuon->innerTrack()->dxy(primVertex.position()) : 0.);
        muo.dZ_in =(amuon->innerTrack().isNonnull() ? amuon->innerTrack()->dz(primVertex.position()) : 0.);


        muo.isGlobalMuon = amuon->isGlobalMuon();
        muo.isTrackerMuon = amuon->isTrackerMuon();
        muo.isStandAloneMuon = amuon->isStandAloneMuon();





muo.isPFMuon = amuon->isPFMuon();
muo.numMatchStation= amuon->numberOfMatchedStations();
muo.normalizedChi2_innTrk = (amuon->innerTrack().isNonnull() ? amuon->innerTrack()->normalizedChi2() : 0);
        muo.numberOfValidMuonHits_innTrk = (amuon->innerTrack().isNonnull() ? amuon->innerTrack()->hitPattern().numberOfValidMuonHits() : 0);
        muo.numberOfHits_innTrk = (amuon->innerTrack().isNonnull() ? amuon->innerTrack()->hitPattern().numberOfHits() : 0);
//
//        edm::Ref<edm::View<reco::Muon> > themu(MuCandidates, indexbis);
        edm::Ref<pat::MuonCollection> muRef(muonsHandle, qq);

        muo.pfIsoAll = (*isoAllMuMap)[muRef];
        muo.pfIsoCharged = (*isoChargedMuMap)[muRef];
        muo.pfIsoNeutral = (*isoNeutralMuMap)[muRef];
        muo.pfIsoGamma = (*isoGammaMuMap)[muRef];
        muo.pfIsoPU = (*isoPUMuMap)[muRef];
//        muo.pfIsoPULow = (*isoPULowMuMap)[muRef];

        //        cout<<"IsoMu_PF"<<(*isoGammaMuMap)[muRef]<<endl;

        muo.z_expo = 0.0;
        if (amuon->innerTrack().isNonnull()) {
            reco::TransientTrack track(amuon->innerTrack(), magField, geomHandle);
            //FreeTrajectoryState state = track.impactPointTSCP().theState();
            //states.push_back(state);
            //calculate Z at beamspot
            TransverseImpactPointExtrapolator extrapolator(magField);
            TrajectoryStateOnSurface closestOnTransversePlaneState = extrapolator.extrapolate(track.impactPointState(), GlobalPoint(beamSpot->position().x(), beamSpot->position().y(), 0.0));
            muo.z_expo = closestOnTransversePlaneState.globalPosition().z();

        }
	
	const pat::TriggerObjectRef trigRef_loose(matchHelper.triggerMatchObject(muonsHandle, qq, muonMatch_Loose_, iEvent, *triggerEvent));
        muo.hasTrgObject_loose = false;
        muo.TrgObjectEta_loose = -100;
        muo.TrgObjectPt_loose = -100;
        muo.TrgObjectPhi_loose = -100;
        // finally we can fill the histograms
        if (trigRef_loose.isAvailable()) { // check references (necessary!)

	  muo.hasTrgObject_loose = true;
	  muo.TrgObjectEta_loose = trigRef_loose->eta();
	  muo.TrgObjectPt_loose = trigRef_loose->pt();
	  muo.TrgObjectPhi_loose = trigRef_loose->phi();
        }

	const pat::TriggerObjectRef trigRef_medium(matchHelper.triggerMatchObject(muonsHandle, qq, muonMatch_Medium_, iEvent, *triggerEvent));
	muo.hasTrgObject_medium = false;
        muo.TrgObjectEta_medium = -100;
        muo.TrgObjectPt_medium = -100;
        muo.TrgObjectPhi_medium = -100;
        // finally we can fill the histograms                                                                                                                                                                                        
        if (trigRef_medium.isAvailable()) { // check references (necessary!)                                                                                                                                                          

          muo.hasTrgObject_medium = true;
          muo.TrgObjectEta_medium = trigRef_medium->eta();
          muo.TrgObjectPt_medium = trigRef_medium->pt();
          muo.TrgObjectPhi_medium = trigRef_medium->phi();
        }


        (m->PreSelectedMuons).push_back(muo);


    }//for muons


}
