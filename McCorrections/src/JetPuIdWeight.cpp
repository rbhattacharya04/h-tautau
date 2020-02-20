/*! b tag weight.
Original code: https://github.com/ajgilbert/ICHiggsTauTau/blob/master/Analysis/Utilities/src/BTagWeight.cc
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "h-tautau/McCorrections/include/JetPUIdWeight.h"

namespace analysis {
namespace mc_corrections {

namespace detail {

JetInfo::JetInfo(EventInfoBase& eventInfo, size_t jetIndex) :
    eff(0.), SF(0.)
{
    const ntuple::Event& event = *eventInfo;
    pt  = event.jets_p4.at(jetIndex).pt();
    eta = event.jets_p4.at(jetIndex).eta();

}

JetPUIdReaderInfo::JetPUIdReaderInfo(FileNamesMap Eff_files, FileNamesMap Sf_files, int _isHardJet, Period _period, DiscriminatorWP wp) :
 isHardJet(_isHardJet), period(_period)
{

    static const std::map<DiscriminatorWP, std::string> wp_prefixes = {
        { DiscriminatorWP::Loose, "L" }, { DiscriminatorWP::Medium, "M" },
        { DiscriminatorWP::Tight, "T" }
    };

    if(!wp_prefixes.count(wp))
        throw exception("Jet pu Id working point %1% not supported.") % wp;

    static const std::map<int, std::string> isHardJet_prefixes = {
        {0 , "mistag"}, {1 , "eff"},
    };
    FilePtr eff_file = root_ext::OpenRootFile(Eff_files.at(period).at(wp).at(isHardJet));
    FilePtr sf_file = root_ext::OpenRootFile(Sf_files.at(period).at(wp).at(isHardJet));

    const std::string eff_name = boost::str(boost::format("h2_%1%_mc%3%_%4%")
                                        % isHardJet_prefixes.at(isHardJet) % period % wp_prefixes.at(wp));
    const std::string sf_name = boost::str(boost::format("h2_%1%_sf%3%_%4%")
                                       % isHardJet_prefixes.at(isHardJet) % period % wp_prefixes.at(wp));
    eff_hist = HistPtr(root_ext::ReadCloneObject<TH2D>(*eff_file, eff_name, "", true));
    sf_hist = HistPtr(root_ext::ReadCloneObject<TH2D>(*sf_file, sf_name, "", true));
}

void JetPUIdReaderInfo::Eval(JetInfo& jetInfo, const std::string& unc_name)
{
    jetInfo.SF  = GetSF(jetInfo.pt, std::abs(jetInfo.eta));
    jetInfo.eff = GetEfficiency(jetInfo.pt, std::abs(jetInfo.eta));
}

double JetPUIdReaderInfo::GetEfficiency(double pt, double eta) const
{
    int xBin = eff_hist->GetXaxis()->FindFixBin(pt);
    xBin = std::min(eff_hist->GetXaxis()->GetNbins(), std::max(1, xBin));
    int yBin = eff_hist->GetYaxis()->FindFixBin(eta);
    yBin = std::min(eff_hist->GetYaxis()->GetNbins(), std::max(1, yBin));
    return eff_hist->GetBinContent(xBin, yBin);
}

double JetPUIdReaderInfo::GetSF(double pt, double eta) const
{
    int xBin = sf_hist->GetXaxis()->FindFixBin(pt);
    xBin = std::min(sf_hist->GetXaxis()->GetNbins(), std::max(1, xBin));
    int yBin = sf_hist->GetYaxis()->FindFixBin(eta);
    yBin = std::min(sf_hist->GetYaxis()->GetNbins(), std::max(1, yBin));
    return sf_hist->GetBinContent(xBin, yBin);
}


} // namespace detail

JetPUIdWeight::JetPUIdWeight(FileNamesMap Eff_files, FileNamesMap Sf_files,
                             Period period, DiscriminatorWP _wp) :
    wp(_wp)
{

    static const std::vector<DiscriminatorWP> op_vec = {
        DiscriminatorWP::Loose, DiscriminatorWP::Medium,
        DiscriminatorWP::Tight
    };

    static const std::map<std::string, int> jet_configurations {
        { "HardJet", 1 }, { "Not_HardJet", 0 },
    };

    if(std::find(op_vec.begin(), op_vec.end(), wp) != op_vec.end())
        throw exception("Working point %1% is not supported.") % wp;

    auto jetPuIdEffFile = root_ext::OpenRootFile(jetPuIdEffFileName);
    auto jetPuIdSFFile = root_ext::OpenRootFile(jetPuIdSFFileName);
    auto jetPuIdMisTagFile = root_ext::OpenRootFile(jetPuIdMisTagFileName);

    readerInfos[jet_configurations.at("HardJet")] =
            ReaderInfoPtr(new ReaderInfo(Eff_files, Sf_files, jet_configurations.at("HardJet"),period,wp));
    readerInfos[jet_configurations.at("Not_HardJet")] =
            ReaderInfoPtr(new ReaderInfo(Eff_files, Sf_files, jet_configurations.at("Not_HardJet"),period,wp));



}

double JetPUIdWeight::Get(EventInfoBase& eventInfo) const
{
    return GetEx(eventInfo, UncertaintyScale::Central);
}

double JetPUIdWeight::Get(const ntuple::ExpressEvent& /*event*/) const
{
    throw exception("ExpressEvent is not supported in BTagWeight::Get.");
}

double JetPUIdWeight::GetEx(EventInfoBase& eventInfo, UncertaintyScale unc) const
{
    const std::string unc_name = GetUncertantyName(unc);

    JetInfoVector jetInfos;
    const ntuple::Event& event = *eventInfo;
    for (size_t jetIndex = 0; jetIndex < event.jets_p4.size(); ++jetIndex) {
        JetInfo jetInfo(eventInfo, jetIndex);
        if(jetInfo.pt > cuts::hh_bbtautau_2017::jetID::max_pt_veto) continue;
        jetInfo.isHardJet = 0;
        LorentzVectorE jet_p4 = event.jets_p4.at(jetIndex);
        double min_dr = 999.9;
        for(size_t genJetIndex = 0; genJetIndex < event.genJets_p4.size(); ++ genJetIndex){
            LorentzVectorE genJet_p4 = event.jets_p4.at(genJetIndex);
            double dr = ROOT::Math::VectorUtil::DeltaR(genJet_p4, jet_p4);
            if (dr < min_dr) min_dr = dr;
        }
        if(min_dr<0.4) jetInfo.isHardJet = 1;
        GetReader(jetInfo.isHardJet).Eval(jetInfo, unc_name);
        analysis::DiscriminatorIdResults jet_pu_id = jet->GetPuId();
        jetInfo.jetPuIdOutcome = (jet_pu_id.Passed(wp));
        jetInfos.push_back(jetInfo);
    }

    return GetJetPUIdWeight(jetInfos);
}

std::string BTagWeight::GetUncertantyName(UncertaintyScale unc)
{
    std::string unc_name = ToString(unc);
    std::transform(unc_name.begin(), unc_name.end(), unc_name.begin(), ::tolower);
    return unc_name;
}

double JetPUIdWeight::GetJetPUIdWeight(const JetInfoVector& jetInfos)
{
    double MC = 1;
    double Data = 1;
    for (const auto& jetInfo : jetInfos) {
        MC *= jetInfo.JetPuIdOutcome ? jetInfo.eff : 1 - jetInfo.eff;
        Data *= jetInfo.JetPuIdOutcome ? jetInfo.eff * jetInfo.SF : 1 - jetInfo.eff * jetInfo.SF;
    }
    return MC != 0 ? Data/MC : 0;
}

JetPUIdWeight::ReaderInfo& JetPUIdWeight::GetReader(int isHardJet) const
{
    int default_isHardJet = 0;
    int hardJet = std::abs(isHardJet);
    if(!readerInfos.count(isHardJet))
        hardJet = default_isHardJet;
    return *readerInfos.at(hardJet);
}

} // namespace mc_corrections
} // namespace analysis
