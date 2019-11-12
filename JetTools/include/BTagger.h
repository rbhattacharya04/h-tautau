/*! b-jet tagging.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "h-tautau/Core/include/AnalysisTypes.h"
#include "h-tautau/Core/include/EventTuple.h"
#include "h-tautau/Core/include/TupleObjects.h"

namespace analysis {

enum class JetOrdering { NoOrdering, Pt, CSV, DeepCSV, DeepFlavour, HHJetTag };
ENUM_NAMES(JetOrdering) = {
    { JetOrdering::NoOrdering, "NoOrdering" },
    { JetOrdering::Pt, "Pt" },
    { JetOrdering::CSV, "CSV" },
    { JetOrdering::DeepCSV, "DeepCSV" },
    { JetOrdering::DeepFlavour, "DeepFlavour" },
    { JetOrdering::HHJetTag, "HHJetTag" },
};

struct BTagger {
public:
    BTagger(Period _period, JetOrdering _ordering);

    double BTag(const ntuple::Event& event, size_t jet_index,
        analysis::UncertaintySource unc_source,analysis::UncertaintyScale unc_scale,
        bool use_base_ordering) const;
    double BTag(const ntuple::TupleJet& jet, analysis::UncertaintySource unc_source,
        analysis::UncertaintyScale unc_scale, bool use_base_ordering) const;
    bool Pass(const ntuple::Event& event, size_t jet_index,
        analysis::UncertaintySource unc_source,analysis::UncertaintyScale unc_scale,
        DiscriminatorWP wp = DiscriminatorWP::Medium) const;
    bool Pass(const ntuple::TupleJet& jet, analysis::UncertaintySource unc_source,
        analysis::UncertaintyScale unc_scale, DiscriminatorWP wp = DiscriminatorWP::Medium) const;

    double PtCut() const;
    double EtaCut() const;

private:
    Period period;
    JetOrdering ordering;
    JetOrdering base_ordering;
    const std::map<DiscriminatorWP,double>* cut;
};

}
