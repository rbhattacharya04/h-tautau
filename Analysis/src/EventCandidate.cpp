/*! Definition of class EventCandidate.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "h-tautau/Analysis/include/EventCandidate.h"

namespace analysis {

EventCandidate::EventCandidate(const ntuple::Event& _event, UncertaintySource _unc_source,
                               UncertaintyScale _unc_scale) :
    event(&_event), unc_source(_unc_source), unc_scale(_unc_scale), event_id(_event), same_as_central(true),
    tuple_met(*event, MetType::PF), met(tuple_met, tuple_met.cov())
{
    CreateLeptons();
    CreateJets();
    CreateFatJets();
}

void EventCandidate::InitializeUncertainties(Period period, bool is_full, const std::string& working_path,
                                             TauIdDiscriminator tau_id_discriminator)
{
    static const std::map<analysis::Period,std::string> file_unc_sources = {
        { analysis::Period::Run2016,
          "h-tautau/McCorrections/data/2016/JES/Summer16_23Sep2016V4_MC_UncertaintySources_AK4PFchs.txt" },
        { analysis::Period::Run2017,
          "h-tautau/McCorrections/data/2017/JES/Fall17_17Nov2017_V32_MC_UncertaintySources_AK4PFchs.txt" },
        { analysis::Period::Run2018,
          "h-tautau/McCorrections/data/2018/JES/Autumn18_V8_MC_UncertaintySources_AK4PFchs.txt" }
    };
    //taken from: https://twiki.cern.ch/twiki/bin/viewauth/CMS/JECUncertaintySources
    static const std::map<analysis::Period,std::string> file_reduced_unc_sources = {
        { analysis::Period::Run2016,
          "h-tautau/McCorrections/data/2016/JES/RegroupedV2_Summer16_07Aug2017_V11_MC_UncertaintySources_AK4PFchs.txt" },
        { analysis::Period::Run2017,
          "h-tautau/McCorrections/data/2017/JES/RegroupedV2_Fall17_17Nov2017_V32_MC_UncertaintySources_AK4PFchs.txt" },
        { analysis::Period::Run2018,
          "h-tautau/McCorrections/data/2018/JES/RegroupedV2_Autumn18_V19_MC_UncertaintySources_AK4PFchs.txt" }
    };

    std::string full_path_source;
    if(is_full){
        if(!file_unc_sources.count(period))
            throw exception("Period not found in file uncertainty source.");
        full_path_source = tools::FullPath({working_path, file_unc_sources.at(period)});
    }
    else{
        if(!file_reduced_unc_sources.count(period))
            throw exception("Period not found in file reduced uncertainty source.");
        full_path_source = tools::FullPath({working_path, file_reduced_unc_sources.at(period)});
    }

    *jecUncertainties = jec::JECUncertaintiesWrapper(full_path_source,is_full,period);

    const std::map<analysis::Period, std::map<TauIdDiscriminator, std::vector<std::string>>> file_tes = {
        { analysis::Period::Run2016, { { TauIdDiscriminator::byDeepTau2017v2p1VSjet,
                                       { "TauPOG/TauIDSFs/data/TauES_dm_DeepTau2017v2p1VSjet_2016Legacy.root",
                                         "TauPOG/TauIDSFs/data/TauES_dm_DeepTau2017v2p1VSjet_2016Legacy_ptgt100.root" }},
                                       { TauIdDiscriminator::byIsolationMVArun2017v2DBoldDMwLT2017,
                                       { "TauPOG/TauIDSFs/data/TauES_dm_MVAoldDM2017v2_2016Legacy.root",
                                         "TauPOG/TauIDSFs/data/TauES_dm_MVAoldDM2017v2_2016Legacy_ptgt100.root" }}} },
        { analysis::Period::Run2017, { { TauIdDiscriminator::byDeepTau2017v2p1VSjet,
                                       { "TauPOG/TauIDSFs/data/TauES_dm_DeepTau2017v2p1VSjet_2017ReReco.root",
                                         "TauPOG/TauIDSFs/data/TauES_dm_DeepTau2017v2p1VSjet_2017ReReco_ptgt100.root" }},
                                       { TauIdDiscriminator::byIsolationMVArun2017v2DBoldDMwLT2017,
                                       { "TauPOG/TauIDSFs/data/TauES_dm_MVAoldDM2017v2_2017ReReco.root",
                                         "TauPOG/TauIDSFs/data/TauES_dm_MVAoldDM2017v2_2017ReReco_ptgt100.root" }}} },
        { analysis::Period::Run2018, { { TauIdDiscriminator::byDeepTau2017v2p1VSjet,
                                       { "TauPOG/TauIDSFs/data/TauES_dm_DeepTau2017v2p1VSjet_2018ReReco.root",
                                         "TauPOG/TauIDSFs/data/TauES_dm_DeepTau2017v2p1VSjet_2018ReReco_ptgt100.root" }},
                                       { TauIdDiscriminator::byIsolationMVArun2017v2DBoldDMwLT2017,
                                       { "TauPOG/TauIDSFs/data/TauES_dm_MVAoldDM2017v2_2018ReReco.root",
                                         "TauPOG/TauIDSFs/data/TauES_dm_MVAoldDM2017v2_2018ReReco_ptgt100.root" }}} }
    };

    if(!file_tes.count(period) || !file_tes.at(period).count(tau_id_discriminator))
        throw exception("Period not found in files for TauES.");

    const std::map<analysis::Period, std::string> file_ele_faking_tau = {
        { analysis::Period::Run2016, "TauPOG/TauIDSFs/data/TauFES_eta-dm_DeepTau2017v2p1VSe_2016Legacy.root" },
        { analysis::Period::Run2017, "TauPOG/TauIDSFs/data/TauFES_eta-dm_DeepTau2017v2p1VSe_2017ReReco.root" },
        { analysis::Period::Run2018, "TauPOG/TauIDSFs/data/TauFES_eta-dm_DeepTau2017v2p1VSe_2018ReReco.root" },
    };
    if(!file_ele_faking_tau.count(period))
        throw exception("Period not found in files for electron faking tau.");

    *tauESUncertainties = TauESUncertainties(
            tools::FullPath({ working_path, file_tes.at(period).at(tau_id_discriminator).at(0) }),
            tools::FullPath({ working_path, file_tes.at(period).at(tau_id_discriminator).at(1) }),
            tools::FullPath({ working_path, file_ele_faking_tau.at(period) }));
}

const jec::JECUncertaintiesWrapper& EventCandidate::GetJecUncertainties()
{
    if(!(*jecUncertainties))
        throw exception("JEC uncertainties are not initialized.");
    return **jecUncertainties;
}

const TauESUncertainties& EventCandidate::GetTauESUncertainties()
{
    if(!(*tauESUncertainties))
        throw exception("TauES uncertainties are not initialized.");
    return **tauESUncertainties;
}

const LepCollection& EventCandidate::GetLeptons() const { return lepton_candidates; }
const JetCollection& EventCandidate::GetJets() const { return jet_candidates; }
const FatJetCollection& EventCandidate::GetFatJets() const { return fatJets; }
const MET& EventCandidate::GetMET() const { return met; }
bool EventCandidate::IsSameAsCentral() const { return same_as_central; }
const ntuple::Event& EventCandidate::GetEvent() const { return *event; }
UncertaintySource EventCandidate::GetUncSource() const { return unc_source; }
UncertaintyScale EventCandidate::GetUncScale() const { return unc_scale; }
const EventIdentifier& EventCandidate::GetEventId() const { return event_id; }
Channel EventCandidate::GetChannel() const { return static_cast<Channel>(event->channelId); }
Period EventCandidate::GetPeriod() const { return static_cast<Period>(event->period); }

UncertaintySource EventCandidate::GetCacheUncSource() const
{
    return IsSameAsCentral() ? UncertaintySource::None : unc_source;
}

UncertaintyScale EventCandidate::GetCacheUncScale() const
{
    return IsSameAsCentral() ? UncertaintyScale::Central : unc_scale;
}

const EventCacheProvider& EventCandidate::GetCacheProvider()
{
    Lock lock(mutex);
    if(!cache_provider)
        cache_provider = std::make_shared<EventCacheProvider>(*event);
    return *cache_provider;
}

void EventCandidate::SetCacheProvider(const std::shared_ptr<EventCacheProvider>& _cache_provider)
{
    Lock lock(mutex);
    cache_provider = _cache_provider;
}

void EventCandidate::SetHHTagScores(size_t htt_index)
{
    Lock lock(mutex);
    const auto& cache = GetCacheProvider();
    for(size_t jet_index = 0; jet_index < tuple_jets.size(); ++jet_index) {
        const auto cached_score = cache.TryGetHHbtag(htt_index, jet_index, GetCacheUncSource(), GetCacheUncScale());
        float score = -1;
        if(cached_score.is_initialized())
            score = *cached_score;
        tuple_jets.at(jet_index).set_hh_btag(score);
    }
}

void EventCandidate::CreateLeptons()
{
    static constexpr bool preserve_dm0_mass = false;
    double shifted_met_px = 0;
    double shifted_met_py = 0;

    for(size_t n = 0; n < event->lep_type.size(); ++n)
        tuple_leptons.emplace_back(*event, n);
    for(size_t n = 0; n < tuple_leptons.size(); ++n)
        lepton_candidates.emplace_back(tuple_leptons.at(n), tuple_leptons.at(n).iso());
    for(size_t n = 0; n < tuple_leptons.size(); ++n) {
        const auto& tuple_lepton = tuple_leptons.at(n);
        LorentzVectorM lepton_p4(tuple_lepton.p4());
        LorentzVectorM corrected_lepton_p4(tuple_lepton.p4());

        if(!event->isData && tuple_lepton.leg_type() == analysis::LegType::tau) {
            bool tau_same_as_central = true;
            const double sf = GetTauESUncertainties().GetCorrectionFactor(
                    tuple_lepton.decayMode(), tuple_lepton.gen_match(), unc_source, unc_scale,
                    tuple_lepton.p4().pt(), tuple_lepton.p4().eta(), &tau_same_as_central);
            same_as_central = same_as_central && tau_same_as_central;

            if(preserve_dm0_mass && tuple_lepton.decayMode() == 0) {
                const double shifted_pt = lepton_p4.pt() * sf;
                corrected_lepton_p4 = LorentzVectorM(shifted_pt, lepton_p4.eta(), lepton_p4.phi(),lepton_p4.M());
            } else {
                corrected_lepton_p4 = lepton_p4 * sf;
            }

            shifted_met_px += tuple_lepton.p4().px() - corrected_lepton_p4.px();
            shifted_met_py += tuple_lepton.p4().py() - corrected_lepton_p4.py();
        }
        lepton_candidates.at(n).SetMomentum(corrected_lepton_p4);
    }

    shifted_met_px += met.GetMomentum().px();
    shifted_met_py += met.GetMomentum().py();
    analysis::LorentzVectorXYZ shifted_met;
    double E = std::hypot(shifted_met_px,shifted_met_py);
    shifted_met.SetPxPyPzE(shifted_met_px,shifted_met_py,0,E);
    met.SetMomentum(shifted_met);
}

void EventCandidate::CreateJets()
{
    for(size_t n = 0; n < event->jets_p4.size(); ++n)
        tuple_jets.emplace_back(*event, n);
    for(size_t n = 0; n < tuple_jets.size(); ++n)
        jet_candidates.emplace_back(tuple_jets.at(n));
  
    bool m_enabled = true; //Option to switch on or JER smearing. Ultimately have to pass from outside
    bool m_genJetMatcher = true; //Option to switch on gen matching. Ultimately have to pass from outside



    if(m_enabled && !event->isData){
	m_enabled = false;
	m_genJetMatcher = false;
    }
    double MIN_JET_Energy = 1e-2;
    for(const auto& jet : jet_candidates){
	if((!m_enabled) || (jet.GetMomentum().pt() == 0)) continue;

        //get resolution and scale factor
        double jet_resolution;
        double jer_sf;

        LorentzVectorE matched_genJet = 0;

        //match gen jet and reco jet
        //
        if (m_genJetMatcher){
            auto genJets_p4 = event->genJets_p4;
            double min_dR = std::numeric_limits<double>::infinity();
            double m_dR_max = 0.2; //R/2 (R =0.4 for AK4 jets)
            double m_dPt_max_factor; //have to intialise and pass from outside

            for(const auto& genJet : genJets_p4){
                double dR = genJet.DeltaR(genJet,jet);
                if (dR > min_dR) continue;

                if (dR < m_dR_max) {
                    double dPt = std::abs(genJet.pt() - jet.pt());
                    if (dPt > m_dPt_max_factor * jet_resolution) continue;

                    min_dR = dR;
                    matched_genJet = genJet;
                }

            }
        }


        double smearFactor = 1.;
  
        if(matched_genJet){
            /*
	     * Case 1: we have a "good" gen jet matched to the reco jet
	     */
                 
             double dPt = jet.GetMomentum().pt() - genJet.GetMomentum().pt();
      	     smearFactor = 1 + (jer_sf - 1.) * dPt / jet.GetMomentum().pt(); 
        }      
        else if(jer_sf > 1) {
             /*
	      * Case 2: we don't have a gen jet. Smear jet pt using random gaussian variation
	      */ 


	    double sigma = jet_resolution * std::sqrt(jer_sf*jer_sf - 1);
	    std::normal_distribution<> d(0,sigma);
	    smearFactor = 1. + d(m_random_generator) 
        }
        else 
	    std::cout<<"Impossible to smear this jet"<<std::endl;
        
        if(jet.GetMomentum().energy() * smearFactor < MIN_JET_Energy){
	    //Negative or too small smearFactor. We weould change direction of the jet
	    //and this is not what we want.
	    //Recompute the smearing factor in order to have jet energy == MIN_JET_ENERGY	

            double newSmearFactor = MIN_JET_ENERGY/ jet.GetMomentum().energy();
            smearFactor = newSmearFactor;     

	}
        jet.scaleEnergy(smearFactor);
        corrected_jet_candidates.push_back(jet);
    }
   
    

    if(!event->isData && jec::JECUncertaintiesWrapper::IsJetUncertainties(unc_source)) {
        same_as_central = false;
        const auto& other_jets_p4 = event->other_jets_p4;
        auto shifted_met_p4(met.GetMomentum());
        jet_candidates = GetJecUncertainties().ApplyShift(jet_candidates, unc_source, unc_scale, &other_jets_p4,
                                                          &shifted_met_p4);
        met.SetMomentum(shifted_met_p4);
    }
}

void EventCandidate::CreateFatJets()
{
    const size_t N = std::min(event->fatJets_p4.size(), event->fatJets_m_softDrop.size()); // workaround
    for(size_t n = 0; n < N; ++n)
        tuple_fatJets.emplace_back(*event, n);
    for(size_t n = 0; n < tuple_fatJets.size(); ++n)
        fatJets.emplace_back(tuple_fatJets.at(n));
}

const std::unique_ptr<boost::optional<jec::JECUncertaintiesWrapper>> EventCandidate::jecUncertainties
    = std::make_unique<boost::optional<jec::JECUncertaintiesWrapper>>();
const std::unique_ptr<boost::optional<TauESUncertainties>> EventCandidate::tauESUncertainties
    = std::make_unique<boost::optional<TauESUncertainties>>();

} // namespace analysis
