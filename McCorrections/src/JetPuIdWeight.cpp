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

JetPUIdReaderInfo::JetPUIdReaderInfo(FilePtr eff_file, FilePtr sf_file, DiscriminatorWP wp) :
{

    static const std::map<DiscriminatorWP, std::string> wp_prefixes = {
        { DiscriminatorWP::Loose, "L" }, { DiscriminatorWP::Medium, "M" },
        { DiscriminatorWP::Tight, "T" }
    };

    if(!wp_prefixes.count(wp))
        throw exception("Jet pu Id working point %1% not supported.") % wp;

    const std::string eff_name = boost::str(boost::format("All/Efficiency/%1%_%2%_all")
                                        % flavor_prefixes.at(flavor) % wp_prefixes.at(wp));
    const std::string sf_name = boost::str(boost::format("All/Efficiency/%1%_%2%_all")
                                        % flavor_prefixes.at(flavor) % wp_prefixes.at(wp));
    eff_hist = HistPtr(root_ext::ReadCloneObject<TH2D>(*eff_file, name, "", true));
    sf_hist = HistPtr(root_ext::ReadCloneObject<TH2D>(*sf_file, name, "", true));
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

JetPUIdWeight::JetPUIdWeight(const std::string& jetPuIdEffFileName, const std::string& jetPuIdSFFileName,Period period,
                        DiscriminatorWP wp) :
{

    static const std::map<DiscriminatorWP, OperatingPoint> op_map = {
        {DiscriminatorWP::Loose, BTagEntry::OP_LOOSE }, {DiscriminatorWP::Medium, BTagEntry::OP_MEDIUM} ,
        {DiscriminatorWP::Tight, BTagEntry::OP_TIGHT}
    };

    if(!op_map.count(wp))
        throw exception("Working point %1% is not supported.") % wp;

    auto jetPuIdEffFile = root_ext::OpenRootFile(jetPuIdEffFileName);
    auto jetPuIdSFFile = root_ext::OpenRootFile(jetPuIdSFFileName);

    //ReaderPtr reader_b(new BTagCalibrationReader(op_map.at(wp), "central", {"up", "down"}));
    //reader_b->load(calib, BTagEntry::FLAV_B, "comb");
    readerInfo = new ReaderInfo(jetPuIdEffFile,jetPuIdSFFile, wp);

    //readerInfos[jet_flavors.at(BTagEntry::FLAV_B)] =
    //        ReaderInfoPtr(new ReaderInfo(reader_b, BTagEntry::FLAV_B, bTagEffFile, wp));

    //ReaderPtr reader_c(new BTagCalibrationReader(op_map.at(wp), "central", {"up", "down"}));
    //reader_c->load(calib, BTagEntry::FLAV_C, "comb");
    //readerInfos[jet_flavors.at(BTagEntry::FLAV_C)] =
    //        ReaderInfoPtr(new ReaderInfo(reader_c, BTagEntry::FLAV_C, bTagEffFile, wp));

    //ReaderPtr reader_light(new BTagCalibrationReader(op_map.at(wp), "central", {"up", "down"}));
    //reader_light->load(calib, BTagEntry::FLAV_UDSG, "incl");
    //readerInfos[jet_flavors.at(BTagEntry::FLAV_UDSG)] =
    //        ReaderInfoPtr(new ReaderInfo(reader_light, BTagEntry::FLAV_UDSG, bTagEffFile, wp));*/
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
        //if(std::abs(jetInfo.eta) >= bTagger.EtaCut()) continue;
        GetReader().Eval(jetInfo, unc_name);
        jetInfo.jetPuIdOutcome = jet.GetMomentum().pt() > cuts::hh_bbtautau_2017::jetID::max_pt_veto || (jet_pu_id.Passed(analysis::DiscriminatorWP::Loose));
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

JetPUIdWeight::ReaderInfo& JetPUIdWeight::GetReader() const
{
    return *readerInfo;
}

} // namespace mc_corrections
} // namespace analysis
