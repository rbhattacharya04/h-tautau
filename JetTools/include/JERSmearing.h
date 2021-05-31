/*! Apply jet uncertainties to the event.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include <string>
#include <vector>
#include <random>
#include "AnalysisTools/Core/include/TextIO.h"
#include "h-tautau/Core/include/AnalysisTypes.h"
#include "h-tautau/Core/include/Candidate.h"
#include "h-tautau/Core/include/TupleObjects.h"
#include "JetCorrectorParameters.h"   // CondFormats/JetMETObjects/interface
#include "JetCorrectionUncertainty.h" // CondFormats/JetMETObjects/interface
#include "h-tautau/JetTools/include/JetResolution.h" // CondFormats/JetMETObjects/interface/ &  JetMETCorrections/Modules/interface/

namespace jer {

using JetCandidate = analysis::Candidate<ntuple::TupleJet>;
using JetCollection = std::vector<JetCandidate>;
using LorentzVectorE = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float>>;
using analysis::UncertaintySource;
using analysis::UncertaintyScale;

class JERSmearing
{
public:
    //static const std::set<UncertaintySource>& JetFullUncertainties();
    //static const std::set<UncertaintySource>& JetReducedUncertainties();

    JERSmearing(const std::string& Resolution_File_name, const std::string& SF_file_name) { //, analysis::Period& period){
        m_resolution_from_file.reset(new JME::JetResolution(Resolution_File_name));
        m_scale_factor_from_file.reset(new JME::JetResolutionScaleFactor(SF_file_name));

    }

    template<typename JetCollection, typename LorentzVector1 = analysis::LorentzVector,
             typename LorentzVector2 = analysis::LorentzVector>
    JetCollection ApplyShift(const JetCollection& jet_candidates,
                             const ntuple::Event* event,
                             analysis::UncertaintySource uncertainty_source,
                             analysis::UncertaintyScale scale,
                             LorentzVector2* met = nullptr) const
    {
        //static const std::map<UncertaintyScale, bool> scales = {
        //    { UncertaintyScale::Up, true }, { UncertaintyScale::Down, false }
        //};

        static const std::map<UncertaintyScale, Variation> scales_variation = {
            {UncertaintyScale::Central, Variation::NOMINAL}, { UncertaintyScale::Up, Variation::UP },
            { UncertaintyScale::Down, Variation::DOWN }
        };

        JetCollection corrected_jets;
        if(uncertainty_source != UncertaintySource::JetResolution){
            uncertainty_source = UncertaintySource::JetResolution;
            scale = UncertaintyScale::Central;
        }

        double shifted_met_px = 0;
        double shifted_met_py = 0;

        JME::JetResolution resolution = *m_resolution_from_file;
        JME::JetResolutionScaleFactor resolution_sf = *m_scale_factor_from_file;

        for(const auto& jet : jet_candidates){

            if(jet.GetMomentum().pt() == 0) {
                corrected_jets.push_back(jet);
                continue;
            }

            double smearFactor =  GetSmearFactor(jet, event, resolution, resolution_sf, scales_variation.at(scale));

            JetCandidate corr_jet(jet);
            //smearedJet.scaleEnergy(smearFactor);
            corr_jet.GetMomentum().SetE(jet.GetMomentum().energy()*smearFactor);
            corrected_jets.push_back(corr_jet);
            shifted_met_px += jet.GetMomentum().px() - corr_jet.GetMomentum().px();
            shifted_met_py += jet.GetMomentum().py() - corr_jet.GetMomentum().py();
         }

         if(met){
             shifted_met_px += met->px();
             shifted_met_py += met->py();
             double E = std::hypot(shifted_met_px,shifted_met_py);
             met->SetPxPyPzE(shifted_met_px,shifted_met_py,0,E);
         }

        return corrected_jets;
    }

    double GetSmearFactor(JetCandidate jet, const ntuple::Event* event,
                          JME::JetResolution resolution, JME::JetResolutionScaleFactor resolution_sf,
                          Variation variation){

        //get resolution and scale factor
        double jet_resolution = resolution.getResolution({{JME::Binning::JetPt, jet.GetMomentum().pt()},
                                                          {JME::Binning::JetEta, jet.GetMomentum().eta()}});
                                                          //,{JME::Binning::Rho, *rho}});
        double jer_sf = resolution_sf.getScaleFactor({{JME::Binning::JetEta, jet.GetMomentum().eta()}}, variation);

        size_t seed = event->evt + size_t(jet.GetMomentum().pt() * 100)
                       + size_t(std::abs(jet.GetMomentum().eta()) * 100) * 100
                       + size_t(std::abs(jet.GetMomentum().Phi()) * 100) * 10000;

        std::mt19937_64 m_random_generator = std::mt19937_64(seed);

        LorentzVectorE matched_genJet; //= 0;
        bool match_success = false;
        //match gen jet and reco jet
        //

        auto genJets_p4 = event->genJets_p4;
        double min_dR = std::numeric_limits<double>::infinity();
        double m_dR_max = 0.2; //R/2 (R =0.4 for AK4 jets)
        double m_dPt_max_factor; //have to intialise and pass from outside

        for(const auto& genJet : genJets_p4){
            double dR = ROOT::Math::VectorUtil::DeltaR(genJet, jet.GetMomentum());
            if (dR > min_dR) continue;

            if (dR < m_dR_max) {
                double dPt = std::abs(genJet.Pt() - jet.GetMomentum().pt());
                if (dPt > m_dPt_max_factor * jet_resolution) continue;

                min_dR = dR;
                matched_genJet = genJet;
                match_success = true;
            }
        }

        double smearFactor = 1.;

        if(match_success){
            /*
            * Case 1: we have a "good" gen jet matched to the reco jet
            */

            double dPt = jet.GetMomentum().pt() - matched_genJet.Pt();
            smearFactor = 1 + (jer_sf - 1.) * dPt / jet.GetMomentum().pt();
        }
        else if(jer_sf > 1) {
            /*
            * Case 2: we don't have a gen jet. Smear jet pt using random gaussian variation
            */


            double sigma = jet_resolution * std::sqrt(jer_sf*jer_sf - 1);
            std::normal_distribution<> d(0,sigma);
            smearFactor = 1. + d(m_random_generator);
        }
        else
            std::cout<<"Impossible to smear this jet"<<std::endl;

        if(jet.GetMomentum().energy() * smearFactor < MIN_JET_ENERGY){
            //Negative or too small smearFactor. We weould change direction of the jet
            //and this is not what we want.
            //Recompute the smearing factor in order to have jet energy == MIN_JET_ENERGY

            double newSmearFactor = MIN_JET_ENERGY/ jet.GetMomentum().energy();
            smearFactor = newSmearFactor;
        }

        return smearFactor;
    }


private:
    std::unique_ptr<JME::JetResolution> m_resolution_from_file;
    std::unique_ptr<JME::JetResolutionScaleFactor> m_scale_factor_from_file;
    static constexpr const double MIN_JET_ENERGY = 1e-2;

};

} // namespace jec
