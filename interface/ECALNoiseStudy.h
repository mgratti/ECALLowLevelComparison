#ifndef ECALNoiseStudy_h
#define ECALNoiseStudy_h

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "RecoEcal/EgammaCoreTools/interface/PositionCalc.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexSorter.h"
#include "RecoVertex/PrimaryVertexProducer/interface/VertexHigherPtSquared.h"

// ROOT include
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include <vector>


// Less than operator for sorting EcalRecHits according to energy.
class ecalRecHitLess : public std::binary_function<EcalRecHit, EcalRecHit, bool>
{
public:
  bool operator()(EcalRecHit x, EcalRecHit y)
  {
    return (x.energy() > y.energy());
  }
};


class ECALNoiseStudy : public edm::EDAnalyzer {

      public:
         explicit ECALNoiseStudy(const edm::ParameterSet&);
	 ~ECALNoiseStudy();

      private:
	 virtual void beginJob() ;
	 virtual void analyze(const edm::Event&, const edm::EventSetup&);
	 virtual void endJob() ;

	 // ----------member data ---------------------------

	 edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
         edm::EDGetTokenT<reco::GenParticleCollection> genParticleCollection_;
         edm::EDGetTokenT<EcalRecHitCollection>                    recHitCollection_EB_;
	 edm::EDGetTokenT<EcalRecHitCollection>                    recHitCollection_EE_;
         edm::EDGetTokenT<reco::PFRecHitCollection>                PFrecHitCollection_;
         edm::EDGetTokenT<reco::PFClusterCollection>                         PFclusterCollection_;
	 edm::EDGetTokenT<reco::SuperClusterCollection>            superClusterCollection_EB_;//reco::SuperClusterCollection
	 edm::EDGetTokenT<reco::SuperClusterCollection>            superClusterCollection_EE_;
	 edm::EDGetTokenT<reco::BeamSpot>                      	  beamSpot_ ;//reco::BeamSpot

         bool SaveSrFlag_;

	 double ethrEB_;
	 double ethrEE_;
	 double scEtThrEB_;
	 double scEtThrEE_;

         std::string anaName_;

   // tree
   //TTree *outTree;

	 // ------------- HISTOGRAMS ------------------------------------
	 int naiveId_;
	 TH1D *h_nPVs;
	 TH1D *h_numberOfEvents;

   // gen particles
   TH1D *h_genP_pt;
   TH1D *h_genP_eta;
   TH1D *h_genP_phi;
   TH1D *h_genP_status;
   TH1D *h_genP_pdgid;
   TH1D *h_genP_pt_EB;
   TH1D *h_genP_pt_EEP;
   TH1D *h_genP_pt_EEM;

   std::vector<TH2D*> h_PFclusters_etaVsPhi;
   std::vector<TH2D*> h_PFclusters_genMatched_etaVsPhi;
   std::vector<TH2D*> h_superClusters_etaVsPhi;
   std::vector<TH2D*> h_superClusters_genMatched_etaVsPhi;
   std::vector<TH2D*> h_recHits_etaVsPhi;
   std::vector<TH2D*> h_genP_etaVsPhi;

   // Rechits and PfRechit vs eta
   std::vector<TString> regions={"EB", "EEM", "EEP"};
   std::map<TString, std::map<TString, TH1F*>> h_recHits_energy_etaBinned;
   std::map<TString, std::map<TString, TH1F*>> h_recHits_et_etaBinned;
   std::map<TString, std::map<TString, TH1F*>> h_PFrecHits_energy_etaBinned;

   // PFClusters vs eta and Et
   std::map<TString, std::map<TString, TH1F*>> h_PFclusters_genMatched_eOverEtrue_EtaEtBinned;
   std::map<TString, std::map<TString, TH1F*>> h_genP_nEvts_EtaEtBinned;

   // SuperClusters vs eta and Et
   std::map<TString, std::map<TString, TH1F*>> h_superClusters_genMatched_eOverEtrue_EtaEtBinned;

   std::map<TString, std::vector<TString>> eta_keys;
   std::map<TString, std::map<TString, std::pair<Float_t,Float_t>>> eta_edges;

   //std::map<TString, std::vector<TString>> Et_keys;
   //std::map<TString, std::map<TString, std::pair<Float_t,Float_t>>> Et_edges;
   std::vector<TString> Et_keys;
   std::map<TString, std::pair<Float_t,Float_t>> Et_edges;
   std::vector<TString> Eta_keys;
   std::map<TString, std::pair<Float_t,Float_t>> Eta_edges;

	 // RecHits ----------------------------------------------
   TH1D *h_recHits_EB_size;
   TH1D *h_recHits_EB_eta;
   TH1D *h_recHits_EB_maxEneEta;
   TH1D *h_recHits_EB_energy;
   TH1D *h_recHits_EB_energyMax;
   TH1D *h_recHits_EB_time;
   TH1D *h_recHits_EB_Chi2;
   TH1D *h_recHits_EB_OutOfTimeChi2;
   TH1D *h_recHits_EB_E1oE4;
   TH1D *h_recHits_EB_iPhiOccupancy;
   TH1D *h_recHits_EB_iEtaOccupancy;
   TH2D *h_recHits_EB_occupancy;
   TH2D *h_recHits_EB_energy_etaphi;
   TH2D *h_recHits_EB_energy_ietaiphi;
   TH2D *h_recHits_EB_occupancy_gt10;
   TH2D *h_recHits_EB_occupancy_lt10;
   TH1D *h_recHits_EB_energy_spike;
   TH2D *h_recHits_EB_eneVSieta;

   TH1D *h_PFrecHits_EB_eta;
   TH1D *h_PFrecHits_EB_energy;
   TH1D *h_PFrecHits_EB_time;
   TH2D *h_PFrecHits_EB_occupancy;
   TH2D *h_PFrecHits_EB_eneVSieta;

	 //... barrel ( with spike cleaning )
   TH1D *h_recHits_EB_size_cleaned;
   TH1D *h_recHits_EB_energy_cleaned;
   TH1D *h_recHits_EB_time_cleaned;
   TH1D *h_recHits_EB_Chi2_cleaned;

	 // ... endcap
   TH1D *h_recHits_EE_size;
   TH1D *h_recHits_EEP_size;
   TH1D *h_recHits_EEP_eta;
   TH1D *h_recHits_EEP_maxEneEta;
   TH1D *h_recHits_EEP_energy;
   TH1D *h_recHits_EEP_energy_gt25;
   TH1D *h_recHits_EEP_energyMax;
   TH1D *h_recHits_EEP_time;
   TH1D *h_recHits_EEP_Chi2;
   TH1D *h_recHits_EEP_OutOfTimeChi2;
   TH1D *h_recHits_EEP_E1oE4;
   TH1D *h_recHits_EEP_iXoccupancy;
   TH1D *h_recHits_EEP_iYoccupancy;
   TH2D *h_recHits_EEP_occupancy;
   TH2D *h_recHits_EEP_occupancy_etaphi;
   TH2D *h_recHits_EEP_occupancy_gt10;
   TH2D *h_recHits_EEP_occupancy_lt10;
   TH2D *h_recHits_EEP_energy_etaphi;
   TH2D *h_recHits_EEP_energy_ixiy;

   TH1D *h_PFrecHits_EEP_eta;
   TH1D *h_PFrecHits_EEP_energy;
   TH1D *h_PFrecHits_EEP_time;
   TH2D *h_PFrecHits_EEP_occupancy;
   TH2D *h_PFrecHits_EEP_eneVSieta;


   TH1D *h_recHits_EEM_size;
   TH1D *h_recHits_EEM_eta;
   TH1D *h_recHits_EEM_maxEneEta;
   TH1D *h_recHits_EEM_energy;
   TH1D *h_recHits_EEM_energy_gt25;
   TH1D *h_recHits_EEM_energyMax;
   TH1D *h_recHits_EEM_time;
   TH1D *h_recHits_EEM_Chi2;
   TH1D *h_recHits_EEM_OutOfTimeChi2;
   TH1D *h_recHits_EEM_E1oE4;
   TH1D *h_recHits_EEM_iXoccupancy;
   TH1D *h_recHits_EEM_iYoccupancy;
   TH2D *h_recHits_EEM_occupancy;
   TH2D *h_recHits_EEM_occupancy_gt10;
   TH2D *h_recHits_EEM_occupancy_lt10;

   TH1D *h_PFrecHits_EEM_eta;
   TH1D *h_PFrecHits_EEM_energy;
   TH1D *h_PFrecHits_EEM_time;
   TH2D *h_PFrecHits_EEM_occupancy;
   TH2D *h_PFrecHits_EEM_eneVSieta;

   TH1D *h_recHits_eta;  // all

   TH2D *h_recHits_EEP_neighbourEt_eta20;
   TH2D *h_recHits_EEP_neighbourEt_eta24;
   TH1D *h_recHits_EEP_sumneighbourEt_eta20;
   TH1D *h_recHits_EEP_sumneighbourEt_eta24;

   // PF clusters ----------------------------------------------
   TH1D *h_PFclusters_EB_size;
   TH1D *h_PFclusters_EB_nXtals;
   TH1D *h_PFclusters_EB_energy;
   TH1D *h_PFclusters_EB_et;
   TH1D *h_PFclusters_EB_eta;
   TH1D *h_PFclusters_EB_phi;

   TH1D *h_PFclusters_EEP_size;
   TH1D *h_PFclusters_EEP_nXtals;
   TH1D *h_PFclusters_EEP_energy;
   TH1D *h_PFclusters_EEP_et;
   TH1D *h_PFclusters_EEP_eta;
   TH1D *h_PFclusters_EEP_phi;

   TH1D *h_PFclusters_EEM_size;
   TH1D *h_PFclusters_EEM_nXtals;
   TH1D *h_PFclusters_EEM_energy;
   TH1D *h_PFclusters_EEM_et;
   TH1D *h_PFclusters_EEM_eta;
   TH1D *h_PFclusters_EEM_phi;

   TH1D *h_PFclusters_eta;
   TH1D *h_PFclusters_deltaR_gen;
   TH1D *h_PFclusters_deltaR_gen_EB;
   TH1D *h_PFclusters_deltaR_gen_EEP;
   TH1D *h_PFclusters_deltaR_gen_EEM;
   TH1D *h_PFclusters500_deltaR_gen;
   TH1D *h_PFclusters500_deltaR_gen_EB;
   TH1D *h_PFclusters500_deltaR_gen_EEP;
   TH1D *h_PFclusters500_deltaR_gen_EEM;
   TH1D *h_PFclusters1000_deltaR_gen;
   TH1D *h_PFclusters1000_deltaR_gen_EB;
   TH1D *h_PFclusters1000_deltaR_gen_EEP;
   TH1D *h_PFclusters1000_deltaR_gen_EEM;

   TH1D *h_PFclusters_genMatched_EB_size;
   TH1D *h_PFclusters_genMatched_EB_nXtals;
   TH1D *h_PFclusters_genMatched_EB_energy;
   TH1D *h_PFclusters_genMatched_EB_et;
   TH1D *h_PFclusters_genMatched_EB_eta;
   TH1D *h_PFclusters_genMatched_EB_phi;
   TH1D *h_PFclusters_genMatched_EB_eOverEtrue;

   TH1D *h_PFclusters_genMatched_EEP_size;
   TH1D *h_PFclusters_genMatched_EEP_nXtals;
   TH1D *h_PFclusters_genMatched_EEP_energy;
   TH1D *h_PFclusters_genMatched_EEP_et;
   TH1D *h_PFclusters_genMatched_EEP_eta;
   TH1D *h_PFclusters_genMatched_EEP_phi;
   TH1D *h_PFclusters_genMatched_EEP_eOverEtrue;

   TH1D *h_PFclusters_genMatched_EEM_size;
   TH1D *h_PFclusters_genMatched_EEM_nXtals;
   TH1D *h_PFclusters_genMatched_EEM_energy;
   TH1D *h_PFclusters_genMatched_EEM_et;
   TH1D *h_PFclusters_genMatched_EEM_eta;
   TH1D *h_PFclusters_genMatched_EEM_phi;
   TH1D *h_PFclusters_genMatched_EEM_eOverEtrue;

   // Super Clusters ----------------------------------------------
   // ... barrel
   TH1D *h_superClusters_EB_size;
   TH1D *h_superClusters_EB_nXtals;
   TH1D *h_superClusters_EB_nBC;
   TH1D *h_superClusters_EB_energy;
   TH1D *h_superClusters_EB_rawEnergy;
   TH1D *h_superClusters_EB_et;

   // ... endcap
   TH1D *h_superClusters_EEP_size;
   TH1D *h_superClusters_EEP_nXtals;
   TH1D *h_superClusters_EEP_nBC;
   TH1D *h_superClusters_EEP_energy;
   TH1D *h_superClusters_EEP_rawEnergy;
   TH1D *h_superClusters_EEP_et;

   TH1D *h_superClusters_EEM_size;
   TH1D *h_superClusters_EEM_nXtals;
   TH1D *h_superClusters_EEM_nBC;
   TH1D *h_superClusters_EEM_energy;
   TH1D *h_superClusters_EEM_rawEnergy;
   TH1D *h_superClusters_EEM_et;

   TH1D *h_superClusters_eta;
   TH1D *h_superClusters_EB_eta;
   TH1D *h_superClusters_EE_eta;

   TH2D *h_superClusters_nBCvsEta;
   TH1D *h_superClusters_nBC_0to1;
   TH1D *h_superClusters_nBC_1to1d5;
   TH1D *h_superClusters_nBC_1d5to1d8;
   TH1D *h_superClusters_nBC_1d8to2d1;
   TH1D *h_superClusters_nBC_2d1to2d5;
   TH1D *h_superClusters_nBC_2d5to3;

   TH1D *h_superClusters_deltaR_gen;
   TH1D *h_superClusters_deltaR_gen_EB;
   TH1D *h_superClusters_deltaR_gen_EEP;
   TH1D *h_superClusters_deltaR_gen_EEM;

   TH1D *h_superClusters500_deltaR_gen;
   TH1D *h_superClusters500_deltaR_gen_EB;
   TH1D *h_superClusters500_deltaR_gen_EEP;
   TH1D *h_superClusters500_deltaR_gen_EEM;

   TH1D *h_superClusters_genMatched_EB_size;
   TH1D *h_superClusters_genMatched_EB_nXtals;
   TH1D *h_superClusters_genMatched_EB_energy;
   TH1D *h_superClusters_genMatched_EB_rawEnergy;
   TH1D *h_superClusters_genMatched_EB_et;
   TH1D *h_superClusters_genMatched_EB_eta;
   TH1D *h_superClusters_genMatched_EB_phi;
   TH1D *h_superClusters_genMatched_EB_eOverEtrue;

   TH1D *h_superClusters_genMatched_EEP_size;
   TH1D *h_superClusters_genMatched_EEP_nXtals;
   TH1D *h_superClusters_genMatched_EEP_energy;
   TH1D *h_superClusters_genMatched_EEP_rawEnergy;
   TH1D *h_superClusters_genMatched_EEP_et;
   TH1D *h_superClusters_genMatched_EEP_eta;
   TH1D *h_superClusters_genMatched_EEP_phi;
   TH1D *h_superClusters_genMatched_EEP_eOverEtrue;

   TH1D *h_superClusters_genMatched_EEM_size;
   TH1D *h_superClusters_genMatched_EEM_nXtals;
   TH1D *h_superClusters_genMatched_EEM_energy;
   TH1D *h_superClusters_genMatched_EEM_rawEnergy;
   TH1D *h_superClusters_genMatched_EEM_et;
   TH1D *h_superClusters_genMatched_EEM_eta;
   TH1D *h_superClusters_genMatched_EEM_phi;
   TH1D *h_superClusters_genMatched_EEM_eOverEtrue;


};


#endif
