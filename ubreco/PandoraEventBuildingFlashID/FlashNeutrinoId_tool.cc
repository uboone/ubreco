/**
 *  @file   larpandora/LArPandoraEventBuilding/SliceIdTools/FlashNeutrinoId_tool.cc
 *
 *  @brief  implementation of the flash based neutrino id tool
 */

#include "FlashNeutrinoId_tool.h"

namespace lar_pandora
{

//------------------------------------------------------------------------------------------------------------------------------------------

void FlashNeutrinoId::ClassifySlices(SliceVector &slices, const art::Event &evt)
{
    std::cout << "[FlashNeutrinoId::ClassifySlices] Start SliceID" << std::endl;
    // Reset the output addresses in case we are writing monitoring details to an output file
    m_outputEvent.Reset(evt);

    FlashCandidateVector flashCandidates;
    SliceCandidateVector sliceCandidates;
    FlashNeutrinoId::FlashCandidate beamFlash;
    unsigned int bestSliceIndex;

    try
    {
        // Find the flash, if any, in time with the beam with the largest number of photoelectrons that is sufficiently bright
        this->GetFlashCandidates(evt, flashCandidates);
        beamFlash = this->GetBeamFlash(flashCandidates);
    }
    catch (const FailureMode &)
    {
    }
    try
    {
        this->GetSliceCandidates(evt, slices, sliceCandidates);
        this->IdentifySliceWithBestTopologicalScore(sliceCandidates);
        if (m_outputEvent.m_hasBeamFlash)
        {
            //// WOUTER: Order decides if cosmci mtching always runs or only runs if event is selected.
            // Find the slice - if any that matches best with the beamFlash
            bestSliceIndex = this->GetBestSliceIndex(beamFlash, sliceCandidates);

            // Obvious-Cosmic Beam-Flash Matching
            GetBestObviousCosmicMatch(evt, beamFlash);
            m_outputEvent.m_bestCosmicMatchRatio = sliceCandidates.at(bestSliceIndex).m_flashMatchScore / m_outputEvent.m_bestCosmicMatch;
            std::cout << "[FlashNeutrinoId::ClassifySlices] Obvious Cosmic Rejection ratio: " << m_outputEvent.m_bestCosmicMatchRatio << std::endl;
            if (m_obviousMatchingCut < m_outputEvent.m_bestCosmicMatchRatio)
            {
                m_outputEvent.m_foundATargetSlice = false;
                sliceCandidates.at(bestSliceIndex).m_isTaggedAsTarget = false;
                throw FailureMode("An obvious cosmic matches a lot better with the flash");
            }
            slices.at(bestSliceIndex).TagAsTarget();
        }
        else
        {

            for (auto &slice : sliceCandidates)
            {
                slice.m_isConsideredByFlashId = 0;
            }
        }
    }
    catch (const FailureMode &)
    {
    }

    if (!m_shouldWriteToFile)
        return;

    this->FillFlashTree(flashCandidates);
    this->FillSliceTree(evt, slices, sliceCandidates);
    this->FillEventTree();
    std::cout << "[FlashNeutrinoId::ClassifySlices] Finished SliceID" << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void FlashNeutrinoId::GetFlashCandidates(const art::Event &event, FlashCandidateVector &flashCandidates)
{
    // Collect all flashes from the event
    art::InputTag flashTag(m_flashLabel);
    const auto flashes(*event.getValidHandle<FlashVector>(flashTag));

    for (const auto &flash : flashes)
        flashCandidates.emplace_back(event, flash, m_PMTch_correction);

    m_outputEvent.m_nFlashes = flashCandidates.size();
}

//------------------------------------------------------------------------------------------------------------------------------------------

FlashNeutrinoId::FlashCandidate &FlashNeutrinoId::GetBeamFlash(FlashCandidateVector &flashCandidates)
{
    bool foundFlashInBeamWindow(false);
    unsigned int brightestFlashIndex(std::numeric_limits<unsigned int>::max());
    float maxTotalPE(-std::numeric_limits<float>::max());
    m_outputEvent.m_nFlashesInBeamWindow = 0;

    // Find the brightest flash in the beam window
    for (unsigned int flashIndex = 0; flashIndex < flashCandidates.size(); ++flashIndex)
    {
        // ATTN non const reference is required since monitoring variables are stored in the slice candidate
        auto &flashCandidate(flashCandidates.at(flashIndex));

        if (!flashCandidate.IsInBeamWindow(m_beamWindowStart, m_beamWindowEnd))
            continue;

        m_outputEvent.m_nFlashesInBeamWindow++;

        const auto totalPE(flashCandidate.m_totalPE);
        if (totalPE < maxTotalPE)
            continue;

        foundFlashInBeamWindow = true;
        maxTotalPE = totalPE;
        brightestFlashIndex = flashIndex;
    }

    if (!foundFlashInBeamWindow)
        throw FailureMode("There were no flashes in the beam window");

    // Ensure it is sufficiently bright
    auto &brightestFlash(flashCandidates.at(brightestFlashIndex));
    brightestFlash.m_isBrightestInWindow = true;

    if (!brightestFlash.PassesPEThreshold(m_minBeamFlashPE))
        throw FailureMode("No flashes in the beam window passed the PE threshold");

    // Save the monitoring information
    brightestFlash.m_isBeamFlash = true;
    m_outputEvent.m_hasBeamFlash = true;

    return brightestFlash;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void FlashNeutrinoId::GetSliceCandidates(const art::Event &event, SliceVector &slices, SliceCandidateVector &sliceCandidates)
{
    m_outputEvent.m_nSlices = slices.size();
    if (slices.empty())
        throw FailureMode("No slices to choose from");

    // Collect the PFParticles and their associations to SpacePoints and Hits
    PFParticleVector pfParticles;
    SpacePointVector spacePoints;
    TrackVector pftracks;
    SpacePointsToHits spacePointToHitMap;
    PFParticleMap pfParticleMap;
    PFParticlesToTracks particlesToTracks;

    PFParticlesToSpacePoints pfParticleToSpacePointMap;
    LArPandoraHelper::CollectPFParticles(event, m_pandoraLabel, pfParticles, pfParticleToSpacePointMap);
    LArPandoraHelper::CollectSpacePoints(event, m_pandoraLabel, spacePoints, spacePointToHitMap);
    LArPandoraHelper::BuildPFParticleMap(pfParticles, pfParticleMap);
    LArPandoraHelper::CollectTracks(event, m_pandoraTrackLabel, pftracks, particlesToTracks);
    art::Handle<std::vector<recob::Track>> track_h;
    event.getByLabel(m_pandoraTrackLabel, track_h);

    for (unsigned int sliceIndex = 0; sliceIndex < slices.size(); ++sliceIndex)
    {
        const auto &slice = slices[sliceIndex];
        if (m_hasCRT)
        {
            const art::FindMany<anab::T0> trk_t0_assn_v(track_h, event, m_crtTrackMatchLabel);
            sliceCandidates.emplace_back(event, slice, pfParticleMap, pfParticleToSpacePointMap, spacePointToHitMap, particlesToTracks,
                                         trk_t0_assn_v, m_mcsfitter,
                                         m_pandoraLabel,
                                         m_cosmictagmanager, m_chargeToNPhotonsTrack, m_chargeToNPhotonsShower, m_xclCoef, sliceIndex + 1,
                                         m_verbose, m_ophitLabel, m_UP, m_DOWN, m_anodeTime, m_cathodeTime, m_driftVel, m_ophitPE,
                                         m_nOphit, m_ophit_time_res, m_ophit_pos_res, m_min_track_length, m_dt_resolution_ophit);
        }
        else
        {

            sliceCandidates.emplace_back(event, slice, pfParticleMap, pfParticleToSpacePointMap, spacePointToHitMap, particlesToTracks,
                                         m_mcsfitter,
                                         m_pandoraLabel,
                                         m_cosmictagmanager,
                                         m_chargeToNPhotonsTrack, m_chargeToNPhotonsShower, m_xclCoef, sliceIndex + 1,
                                         m_verbose, m_ophitLabel, m_UP, m_DOWN, m_anodeTime, m_cathodeTime, m_driftVel, m_ophitPE,
                                         m_nOphit, m_ophit_time_res, m_ophit_pos_res, m_min_track_length, m_dt_resolution_ophit);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int FlashNeutrinoId::GetBestSliceIndex(const FlashCandidate &beamFlash, SliceCandidateVector &sliceCandidates)
{
    bool foundViableSlice(false);
    bool foundHighestTopoligicalScore(false);
    unsigned int bestFlashMatchSliceIndex(std::numeric_limits<unsigned int>::max());
    unsigned int bestCombinedSliceIndex(std::numeric_limits<unsigned int>::max());
    float minScore(std::numeric_limits<float>::max());
    m_outputEvent.m_nSlicesAfterPrecuts = 0;

    for (unsigned int sliceIndex = 0; sliceIndex < sliceCandidates.size(); ++sliceIndex)
    {
        auto &sliceCandidate(sliceCandidates.at(sliceIndex));

        // Apply the pre-selection cuts to ensure that the slice is compatible with the beam flash
        if (!sliceCandidate.IsCompatibleWithBeamFlash(beamFlash, m_maxDeltaY, m_maxDeltaZ, m_maxDeltaYSigma, m_maxDeltaZSigma,
                                                      m_minChargeToLightRatio, m_maxChargeToLightRatio))
        {
            //// WOUTER: This line guarantees that the score is availible for every slice!
            sliceCandidate.GetFlashMatchScore(beamFlash, m_flashMatchManager);
            continue;
        }

        foundViableSlice = true;
        m_outputEvent.m_nSlicesAfterPrecuts++;
        if (sliceCandidate.m_hasBestTopologicalScore == true)
        {
            foundHighestTopoligicalScore = true;
            bestCombinedSliceIndex = sliceIndex;
        }
        // ATTN if there is only one slice that passes the pre-selection cuts, then the score won't be used
        const auto &score(sliceCandidate.GetFlashMatchScore(beamFlash, m_flashMatchManager));
        if (score > minScore)
            continue;

        bestFlashMatchSliceIndex = sliceIndex;
        minScore = score;
    }
    if (!foundViableSlice)
        throw FailureMode("None of the slices passed the pre-selection cuts");

    sliceCandidates.at(bestFlashMatchSliceIndex).m_hasBestFlashMatchScore = true;

    if (!foundHighestTopoligicalScore)
        bestCombinedSliceIndex = bestFlashMatchSliceIndex;
    sliceCandidates.at(bestCombinedSliceIndex).m_isTaggedAsTarget = true;
    if (m_outputEvent.m_nSlicesAfterPrecuts < 2)
        sliceCandidates.at(bestCombinedSliceIndex).m_targetMethod = 0;
    else if (foundHighestTopoligicalScore)
        sliceCandidates.at(bestCombinedSliceIndex).m_targetMethod = 1;
    else
        sliceCandidates.at(bestCombinedSliceIndex).m_targetMethod = 2;
    m_outputEvent.m_foundATargetSlice = true;
    m_outputEvent.m_targetSliceMethod = sliceCandidates.at(bestCombinedSliceIndex).m_targetMethod;
    std::cout << "[FlashNeutrinoId] Slice with index " << bestCombinedSliceIndex << " is tagged as neutrino. Method: " << m_outputEvent.m_targetSliceMethod << ", ch2: " << sliceCandidates.at(bestCombinedSliceIndex).m_flashMatchScore << std::endl;
    return bestCombinedSliceIndex;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void FlashNeutrinoId::FillEventTree()
{
    if (!m_pEventTree)
        throw cet::exception("FlashNeutrinoId") << "Trying to fill the event tree which hasn't been configured" << std::endl;

    m_pEventTree->Fill();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void FlashNeutrinoId::FillFlashTree(const FlashCandidateVector &flashCandidates)
{
    if (!m_pFlashTree)
        throw cet::exception("FlashNeutrinoId") << "Trying to fill the flash tree which hasn't been configured" << std::endl;

    for (const auto &flashCandidate : flashCandidates)
    {
        m_outputFlash = flashCandidate;
        m_pFlashTree->Fill();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void FlashNeutrinoId::FillSliceTree(const art::Event &evt, const SliceVector &slices, const SliceCandidateVector &sliceCandidates)
{
    if (!m_pSliceTree)
        throw cet::exception("FlashNeutrinoId") << "Trying to fill the slice tree which hasn't been configured" << std::endl;

    SliceCandidateVector allSliceCandidates(sliceCandidates);

    if (slices.size() != allSliceCandidates.size())
        throw cet::exception("FlashNeutrinoId") << "The number of slice candidates doesn't match the number of slices" << std::endl;

    // If available, get the information from the MC neutrino
    LArPandoraSliceIdHelper::SliceMetadataVector sliceMetadata;
    if (m_hasMCNeutrino)
    {
        simb::MCNeutrino mcNeutrino;
        LArPandoraSliceIdHelper::GetSliceMetadata(slices, evt, m_truthLabel, m_mcParticleLabel, m_hitLabel, m_backtrackLabel,
                                                  m_pandoraLabel, sliceMetadata, mcNeutrino);

        m_nuMode = mcNeutrino.Mode();
        m_nuCCNC = mcNeutrino.CCNC();
        m_nuX = mcNeutrino.X();
        m_nuW = mcNeutrino.W();
        m_nuPt = mcNeutrino.Pt();
        m_nuTheta = mcNeutrino.Theta();

        const auto nuMCParticle(mcNeutrino.Nu());
        const auto leptonMCParticle(mcNeutrino.Lepton());

        m_nuEnergy = nuMCParticle.E();
        m_leptonEnergy = leptonMCParticle.E();
        m_nuVertexX = nuMCParticle.Vx();
        m_nuVertexY = nuMCParticle.Vy();
        m_nuVertexZ = nuMCParticle.Vz();
        m_nuTime = nuMCParticle.T();
        m_nuPdgCode = nuMCParticle.PdgCode();

        if (slices.size() != sliceMetadata.size())
            throw cet::exception("FlashNeutrinoId") << "The number of slice metadata doesn't match the number of slices" << std::endl;
    }

    // Output the info for each slice
    for (unsigned int sliceIndex = 0; sliceIndex < slices.size(); ++sliceIndex)
    {
        m_outputSlice = allSliceCandidates.at(sliceIndex);
        if (m_hasMCNeutrino)
            m_outputSliceMetadata = sliceMetadata.at(sliceIndex);

        m_pSliceTree->Fill();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void FlashNeutrinoId::IdentifySliceWithBestTopologicalScore(SliceCandidateVector &sliceCandidates) const
{
    if (sliceCandidates.empty())
        return;

    float bestTopologicalScore(-std::numeric_limits<float>::max());
    unsigned int bestSliceIndex(std::numeric_limits<unsigned int>::max());

    for (unsigned int sliceIndex = 0; sliceIndex < sliceCandidates.size(); ++sliceIndex)
    {
        const float topologicalScore(sliceCandidates.at(sliceIndex).m_topologicalNeutrinoScore);
        if (topologicalScore < bestTopologicalScore)
            continue;

        bestTopologicalScore = topologicalScore;
        bestSliceIndex = sliceIndex;
    }

    if (bestSliceIndex > sliceCandidates.size())
        throw cet::exception("FlashNeutrinoId") << "Couldn't find slice the best topological score" << std::endl;
    sliceCandidates.at(bestSliceIndex).m_hasBestTopologicalScore = true;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

FlashNeutrinoId::FailureMode::FailureMode(const std::string &reason) : m_reason(reason)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

FlashNeutrinoId::FailureMode::~FailureMode()
{
    std::cout << "[Flash neutrino ID] Failed to find neutrino slice: ";
    std::cout << m_reason << std::endl
              << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

FlashNeutrinoId::SliceCandidate::Deposition::Deposition(const float x, const float y, const float z, const float charge, const float nPhotons) : m_x(x),
                                                                                                                                                 m_y(y),
                                                                                                                                                 m_z(z),
                                                                                                                                                 m_charge(charge),
                                                                                                                                                 m_nPhotons(nPhotons)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

void FlashNeutrinoId::OutputEvent::Reset(const art::Event &event)
{
    m_run = event.run();
    m_subRun = event.subRun();
    m_event = event.event();
    m_timeHigh = event.time().timeHigh();
    m_timeLow = event.time().timeLow();
    m_nFlashes = -std::numeric_limits<int>::max();
    m_nFlashesInBeamWindow = -std::numeric_limits<int>::max();
    m_hasBeamFlash = false;
    m_nSlices = -std::numeric_limits<int>::max();
    m_nSlicesAfterPrecuts = -std::numeric_limits<int>::max();
    m_foundATargetSlice = false;
    m_targetSliceMethod = -1;
    m_bestCosmicMatchRatio = -1;
    m_bestCosmicMatch = -1;
    m_cosmicMatchHypothesis.clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

FlashNeutrinoId::FlashCandidate::FlashCandidate() : m_run(-std::numeric_limits<int>::max()),
                                                    m_subRun(-std::numeric_limits<int>::max()),
                                                    m_event(-std::numeric_limits<int>::max()),
                                                    m_timeHigh(0),
                                                    m_timeLow(0),
                                                    m_time(-std::numeric_limits<float>::max()),
                                                    m_totalPE(-std::numeric_limits<float>::max()),
                                                    m_centerY(-std::numeric_limits<float>::max()),
                                                    m_centerZ(-std::numeric_limits<float>::max()),
                                                    m_widthY(-std::numeric_limits<float>::max()),
                                                    m_widthZ(-std::numeric_limits<float>::max()),
                                                    m_inBeamWindow(false),
                                                    m_isBrightestInWindow(false),
                                                    m_isBeamFlash(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

FlashNeutrinoId::FlashCandidate::FlashCandidate(const art::Event &event, const recob::OpFlash &flash, const std::vector<float> PMTch_correction) : m_run(event.run()),
                                                                                                                                                   m_subRun(event.subRun()),
                                                                                                                                                   m_event(event.event()),
                                                                                                                                                   m_timeHigh(event.time().timeHigh()),
                                                                                                                                                   m_timeLow(event.time().timeLow()),
                                                                                                                                                   m_time(flash.Time()),
                                                                                                                                                   m_inBeamWindow(false),
                                                                                                                                                   m_isBrightestInWindow(false),
                                                                                                                                                   m_isBeamFlash(false)
{
    // Correct the PE values using the gain database and save the values in the OpDet order
    // gain service
    const ::lariov::PmtGainProvider &gain_provider = art::ServiceHandle<lariov::PmtGainService>()->GetProvider();
    // pmt remapping service
    const ::util::PMTRemapProvider &pmtremap_provider = art::ServiceHandle<util::PMTRemapService>()->GetProvider();
    const art::ServiceHandle<geo::Geometry> geometry;
    uint nOpDets(geometry->NOpDets());
    m_peSpectrum.resize(nOpDets);

    for (uint OpChannel = 0; OpChannel < nOpDets; ++OpChannel)
    {
        auto oldch = pmtremap_provider.OriginalOpChannel(OpChannel);
        auto gain = gain_provider.Gain(oldch);
        uint opdet = geometry->OpDetFromOpChannel(OpChannel);

        if (gain != 0)
        {
            m_peSpectrum[opdet] = flash.PEs().at(OpChannel) * 120 / gain * PMTch_correction.at(OpChannel);
        }
        else
        {
            m_peSpectrum[opdet] = 0;
        }
    }

    GetFlashLocation();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void FlashNeutrinoId::FlashCandidate::GetFlashLocation()
{
    // Reset variables
    m_centerY = m_centerZ = 0.;
    m_widthY = m_widthZ = -999.;
    m_totalPE = 0.;
    float sumy = 0., sumz = 0., sumy2 = 0., sumz2 = 0.;

    art::ServiceHandle<geo::Geometry> geometry;
    for (unsigned int opdet = 0; opdet < m_peSpectrum.size(); opdet++)
    {
        double PMTxyz[3];
        geometry->OpDetGeoFromOpDet(opdet).GetCenter(PMTxyz);
        // Add up the position, weighting with PEs
        sumy += m_peSpectrum[opdet] * PMTxyz[1];
        sumy2 += m_peSpectrum[opdet] * PMTxyz[1] * PMTxyz[1];
        sumz += m_peSpectrum[opdet] * PMTxyz[2];
        sumz2 += m_peSpectrum[opdet] * PMTxyz[2] * PMTxyz[2];
        m_totalPE += m_peSpectrum[opdet];
    }
    m_centerY = sumy / m_totalPE;
    m_centerZ = sumz / m_totalPE;
    // This is just sqrt(<x^2> - <x>^2)
    if ((sumy2 * m_totalPE - sumy * sumy) > 0.)
        m_widthY = std::sqrt(sumy2 * m_totalPE - sumy * sumy) / m_totalPE;

    if ((sumz2 * m_totalPE - sumz * sumz) > 0.)
        m_widthZ = std::sqrt(sumz2 * m_totalPE - sumz * sumz) / m_totalPE;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool FlashNeutrinoId::FlashCandidate::IsInBeamWindow(const float beamWindowStart, const float beamWindowEnd)
{
    m_inBeamWindow = (m_time > beamWindowStart && m_time < beamWindowEnd);
    return m_inBeamWindow;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool FlashNeutrinoId::FlashCandidate::PassesPEThreshold(const float minBeamFlashPE) const
{
    return (m_totalPE > minBeamFlashPE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

flashana::Flash_t FlashNeutrinoId::FlashCandidate::ConvertFlashFormat() const
{
    // Ensure the input flash is valid
    const art::ServiceHandle<geo::Geometry> geometry;
    uint nOpDets(geometry->NOpDets());
    if (m_peSpectrum.size() != nOpDets)
        throw cet::exception("FlashNeutrinoId") << "Number of channels in beam flash doesn't match the number of OpDets!" << std::endl;

    // Set the flash properties
    flashana::Flash_t flash;
    flash.x = 0;
    flash.x_err = 0;
    flash.y = m_centerY;
    flash.y_err = m_widthY;
    flash.z = m_centerZ;
    flash.z_err = m_widthZ;
    flash.time = m_time;
    flash.pe_v.resize(nOpDets);
    flash.pe_err_v.resize(nOpDets);

    // Fill the flash with the PE spectrum
    for (unsigned int i = 0; i < nOpDets; ++i)
    {
        const auto PE(m_peSpectrum.at(i));
        flash.pe_v.at(i) = PE;
        flash.pe_err_v.at(i) = std::sqrt(PE);
    }

    return flash;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

FlashNeutrinoId::SliceCandidate::SliceCandidate()
    : m_sliceId(-std::numeric_limits<int>::max()),
      m_run(-std::numeric_limits<int>::max()),
      m_subRun(-std::numeric_limits<int>::max()),
      m_event(-std::numeric_limits<int>::max()),
      m_timeHigh(0),
      m_timeLow(0),
      m_hasDeposition(false),
      m_totalCharge(-std::numeric_limits<float>::max()),
      m_centerX(-std::numeric_limits<float>::max()),
      m_centerY(-std::numeric_limits<float>::max()),
      m_centerZ(-std::numeric_limits<float>::max()),
      m_minCRTdist(std::numeric_limits<float>::max()),
      m_CRTtime(-std::numeric_limits<float>::max()),
      m_deltaY(-std::numeric_limits<float>::max()),
      m_deltaZ(-std::numeric_limits<float>::max()),
      m_deltaYSigma(-std::numeric_limits<float>::max()),
      m_deltaZSigma(-std::numeric_limits<float>::max()),
      m_chargeToLightRatio(-std::numeric_limits<float>::max()),
      m_xChargeLightVariable(-std::numeric_limits<float>::max()),
      m_passesPrecuts(false),
      m_flashMatchScore(-std::numeric_limits<float>::max()),
      m_totalPEHypothesis(-std::numeric_limits<float>::max()),
      m_isTaggedAsTarget(false),
      m_targetMethod(-std::numeric_limits<int>::max()),
      m_isConsideredByFlashId(false),
      m_topologicalNeutrinoScore(-std::numeric_limits<float>::max()),
      m_hasBestTopologicalScore(false),
      m_hasBestFlashMatchScore(false),
      m_chargeToNPhotonsTrack(-std::numeric_limits<float>::max()),
      m_chargeToNPhotonsShower(-std::numeric_limits<float>::max()),
      m_xclCoef(-std::numeric_limits<float>::max()),
      m_maxDeltaLLMCS(-std::numeric_limits<float>::max()),
      m_lengthDeltaLLMCS(-std::numeric_limits<float>::max()),
      m_ct_result_michel_plane0(false),
      m_ct_result_michel_plane1(false),
      m_ct_result_michel_plane2(false),
      m_ct_result_bragg_plane0(false),
      m_ct_result_bragg_plane1(false),
      m_ct_result_bragg_plane2(false),
      m_dqds_startend_percdiff_plane0(-std::numeric_limits<float>::max()),
      m_dqds_startend_percdiff_plane1(-std::numeric_limits<float>::max()),
      m_dqds_startend_percdiff_plane2(-std::numeric_limits<float>::max()),
      m_bragg_local_lin_plane0(-std::numeric_limits<float>::max()),
      m_bragg_local_lin_plane1(-std::numeric_limits<float>::max()),
      m_bragg_local_lin_plane2(-std::numeric_limits<float>::max()),
      m_n_michel_hits_plane0(-std::numeric_limits<int>::max()),
      m_n_michel_hits_plane1(-std::numeric_limits<int>::max()),
      m_n_michel_hits_plane2(-std::numeric_limits<int>::max()),
      m_min_lin_braggalgonly_plane0(-std::numeric_limits<int>::max()),
      m_min_lin_braggalgonly_plane1(-std::numeric_limits<int>::max()),
      m_min_lin_braggalgonly_plane2(-std::numeric_limits<int>::max())
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

FlashNeutrinoId::SliceCandidate::SliceCandidate(const art::Event &event, const Slice &slice)
    : m_sliceId(-std::numeric_limits<int>::max()),
      m_run(event.run()),
      m_subRun(event.subRun()),
      m_event(event.event()),
      m_timeHigh(event.time().timeHigh()),
      m_timeLow(event.time().timeLow()),
      m_hasDeposition(false),
      m_totalCharge(-std::numeric_limits<float>::max()),
      m_centerX(-std::numeric_limits<float>::max()),
      m_centerY(-std::numeric_limits<float>::max()),
      m_centerZ(-std::numeric_limits<float>::max()),
      m_minCRTdist(std::numeric_limits<float>::max()),
      m_CRTtime(-std::numeric_limits<float>::max()),
      m_deltaY(-std::numeric_limits<float>::max()),
      m_deltaZ(-std::numeric_limits<float>::max()),
      m_deltaYSigma(-std::numeric_limits<float>::max()),
      m_deltaZSigma(-std::numeric_limits<float>::max()),
      m_chargeToLightRatio(-std::numeric_limits<float>::max()),
      m_xChargeLightVariable(-std::numeric_limits<float>::max()),
      m_passesPrecuts(false),
      m_flashMatchScore(-std::numeric_limits<float>::max()),
      m_totalPEHypothesis(-std::numeric_limits<float>::max()),
      m_isTaggedAsTarget(false),
      m_targetMethod(-std::numeric_limits<int>::max()),
      m_isConsideredByFlashId(false),
      m_topologicalNeutrinoScore(slice.GetTopologicalScore()),
      m_hasBestTopologicalScore(false),
      m_hasBestFlashMatchScore(false),
      m_chargeToNPhotonsTrack(-std::numeric_limits<float>::max()),
      m_chargeToNPhotonsShower(-std::numeric_limits<float>::max()),
      m_xclCoef(-std::numeric_limits<float>::max()),
      m_maxDeltaLLMCS(-std::numeric_limits<float>::max()),
      m_lengthDeltaLLMCS(-std::numeric_limits<float>::max()),
      m_ct_result_michel_plane0(false),
      m_ct_result_michel_plane1(false),
      m_ct_result_michel_plane2(false),
      m_ct_result_bragg_plane0(false),
      m_ct_result_bragg_plane1(false),
      m_ct_result_bragg_plane2(false),
      m_dqds_startend_percdiff_plane0(-std::numeric_limits<float>::max()),
      m_dqds_startend_percdiff_plane1(-std::numeric_limits<float>::max()),
      m_dqds_startend_percdiff_plane2(-std::numeric_limits<float>::max()),
      m_bragg_local_lin_plane0(-std::numeric_limits<float>::max()),
      m_bragg_local_lin_plane1(-std::numeric_limits<float>::max()),
      m_bragg_local_lin_plane2(-std::numeric_limits<float>::max()),
      m_n_michel_hits_plane0(-std::numeric_limits<int>::max()),
      m_n_michel_hits_plane1(-std::numeric_limits<int>::max()),
      m_n_michel_hits_plane2(-std::numeric_limits<int>::max()),
      m_min_lin_braggalgonly_plane0(-std::numeric_limits<int>::max()),
      m_min_lin_braggalgonly_plane1(-std::numeric_limits<int>::max()),
      m_min_lin_braggalgonly_plane2(-std::numeric_limits<int>::max())
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

FlashNeutrinoId::SliceCandidate::SliceCandidate(const art::Event &event, const Slice &slice, const PFParticleMap &pfParticleMap,
                                                const PFParticlesToSpacePoints &pfParticleToSpacePointMap, const SpacePointsToHits &spacePointToHitMap,
                                                const PFParticlesToTracks &particlesToTracks,
                                                const trkf::TrajectoryMCSFitter &mcsfitter, std::string &pandoraLabel, fhicl::ParameterSet &cosmictagmanager,
                                                const float chargeToNPhotonsTrack, const float chargeToNPhotonsShower, const float xclCoef, const int sliceId,
                                                bool m_verbose, std::string m_ophitLabel, float m_UP, float m_DOWN, float m_anodeTime, float m_cathodeTime,
                                                float m_driftVel, float m_ophitPE, int m_nOphit, float m_ophit_time_res, float m_ophit_pos_res, float m_min_track_length,
                                                float m_dt_resolution_ophit)
    : m_sliceId(sliceId),
      m_run(event.run()),
      m_subRun(event.subRun()),
      m_event(event.event()),
      m_timeHigh(event.time().timeHigh()),
      m_timeLow(event.time().timeLow()),
      m_hasDeposition(false),
      m_totalCharge(-std::numeric_limits<float>::max()),
      m_centerX(-std::numeric_limits<float>::max()),
      m_centerY(-std::numeric_limits<float>::max()),
      m_centerZ(-std::numeric_limits<float>::max()),
      m_minCRTdist(std::numeric_limits<float>::max()),
      m_CRTtime(-std::numeric_limits<float>::max()),
      m_deltaY(-std::numeric_limits<float>::max()),
      m_deltaZ(-std::numeric_limits<float>::max()),
      m_deltaYSigma(-std::numeric_limits<float>::max()),
      m_deltaZSigma(-std::numeric_limits<float>::max()),
      m_chargeToLightRatio(-std::numeric_limits<float>::max()),
      m_xChargeLightVariable(-std::numeric_limits<float>::max()),
      m_passesPrecuts(false),
      m_flashMatchScore(-std::numeric_limits<float>::max()),
      m_totalPEHypothesis(-std::numeric_limits<float>::max()),
      m_isTaggedAsTarget(false),
      m_targetMethod(-std::numeric_limits<int>::max()),
      m_isConsideredByFlashId(true),
      m_topologicalNeutrinoScore(slice.GetTopologicalScore()),
      m_hasBestTopologicalScore(false),
      m_hasBestFlashMatchScore(false),
      m_chargeToNPhotonsTrack(chargeToNPhotonsTrack),
      m_chargeToNPhotonsShower(chargeToNPhotonsShower),
      m_xclCoef(xclCoef),
      m_maxDeltaLLMCS(-std::numeric_limits<float>::max()),
      m_lengthDeltaLLMCS(-std::numeric_limits<float>::max()),
      m_ct_result_michel_plane0(false),
      m_ct_result_michel_plane1(false),
      m_ct_result_michel_plane2(false),
      m_ct_result_bragg_plane0(false),
      m_ct_result_bragg_plane1(false),
      m_ct_result_bragg_plane2(false),
      m_dqds_startend_percdiff_plane0(-std::numeric_limits<float>::max()),
      m_dqds_startend_percdiff_plane1(-std::numeric_limits<float>::max()),
      m_dqds_startend_percdiff_plane2(-std::numeric_limits<float>::max()),
      m_bragg_local_lin_plane0(-std::numeric_limits<float>::max()),
      m_bragg_local_lin_plane1(-std::numeric_limits<float>::max()),
      m_bragg_local_lin_plane2(-std::numeric_limits<float>::max()),
      m_n_michel_hits_plane0(-std::numeric_limits<int>::max()),
      m_n_michel_hits_plane1(-std::numeric_limits<int>::max()),
      m_n_michel_hits_plane2(-std::numeric_limits<int>::max()),
      m_min_lin_braggalgonly_plane0(-std::numeric_limits<int>::max()),
      m_min_lin_braggalgonly_plane1(-std::numeric_limits<int>::max()),
      m_min_lin_braggalgonly_plane2(-std::numeric_limits<int>::max())
{
    const auto chargeDeposition(this->GetDepositionVector(pfParticleMap, pfParticleToSpacePointMap, spacePointToHitMap, slice));
    m_lightCluster = this->GetLightCluster(chargeDeposition);
    m_totalCharge = this->GetTotalCharge(chargeDeposition);
    m_hasDeposition = (m_totalCharge > std::numeric_limits<float>::epsilon());

    if (!m_hasDeposition)
        return;

    const auto chargeCenter(this->GetChargeWeightedCenter(chargeDeposition));
    m_centerX = chargeCenter.GetX();
    m_centerY = chargeCenter.GetY();
    m_centerZ = chargeCenter.GetZ();

    this->RejectStopMuByDirMCS(slice.GetCosmicRayHypothesis(), event, particlesToTracks, mcsfitter);

    this->RejectStopMuByCalo(slice.GetCosmicRayHypothesis(), event, particlesToTracks, pfParticleToSpacePointMap, pandoraLabel, cosmictagmanager);

    // copying variables from one class to another
    mm_verbose = m_verbose;
    mm_ophitLabel = m_ophitLabel;
    mm_UP = m_UP;
    mm_DOWN = m_DOWN;
    mm_anodeTime = m_anodeTime;
    mm_cathodeTime = m_cathodeTime;
    mm_driftVel = m_driftVel;
    mm_ophitPE = m_ophitPE;
    mm_nOphit = m_nOphit;
    mm_ophit_pos_res = m_ophit_pos_res;
    mm_ophit_time_res = m_ophit_time_res;
    mm_min_track_length = m_min_track_length;
    mm_dt_resolution_ophit = m_dt_resolution_ophit;

    // ACPT tagger
    this->ACPTtagger(slice.GetCosmicRayHypothesis(), event, particlesToTracks);
}

//------------------------------------------------------------------------------------------------------------------------------------------

FlashNeutrinoId::SliceCandidate::SliceCandidate(const art::Event &event, const Slice &slice, const PFParticleMap &pfParticleMap,
                                                const PFParticlesToSpacePoints &pfParticleToSpacePointMap, const SpacePointsToHits &spacePointToHitMap,
                                                const PFParticlesToTracks &particlesToTracks, const art::FindMany<anab::T0> &trk_t0_assn_v,
                                                const trkf::TrajectoryMCSFitter &mcsfitter, std::string &pandoraLabel, fhicl::ParameterSet &cosmictagmanager,
                                                const float chargeToNPhotonsTrack, const float chargeToNPhotonsShower, const float xclCoef, const int sliceId,
                                                bool m_verbose, std::string m_ophitLabel, float m_UP, float m_DOWN, float m_anodeTime, float m_cathodeTime,
                                                float m_driftVel, float m_ophitPE, int m_nOphit, float m_ophit_time_res, float m_ophit_pos_res, float m_min_track_length, float m_dt_resolution_ophit)
    : FlashNeutrinoId::SliceCandidate::SliceCandidate(event, slice, pfParticleMap,
                                                      pfParticleToSpacePointMap, spacePointToHitMap,
                                                      particlesToTracks, mcsfitter, pandoraLabel, cosmictagmanager, chargeToNPhotonsTrack, chargeToNPhotonsShower, xclCoef,
                                                      sliceId, m_verbose, m_ophitLabel, m_UP, m_DOWN, m_anodeTime, m_cathodeTime,
                                                      m_driftVel, m_ophitPE, m_nOphit, m_ophit_time_res, m_ophit_pos_res, m_min_track_length, m_dt_resolution_ophit)
{

    if (trk_t0_assn_v.size() != 0)
    {
        this->GetClosestCRTCosmic(slice.GetCosmicRayHypothesis(), event, particlesToTracks, trk_t0_assn_v);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

FlashNeutrinoId::SliceCandidate::DepositionVector FlashNeutrinoId::SliceCandidate::GetDepositionVector(const PFParticleMap &pfParticleMap,
                                                                                                       const PFParticlesToSpacePoints &pfParticleToSpacePointMap, const SpacePointsToHits &spacePointToHitMap, const Slice &slice) const
{
    // Collect all PFParticles in the slice, including those downstream of the primaries
    // ATTN here we only use the neutrino hypothesis, in theory this should work with either (or indeed both with some thought)
    PFParticleVector allParticlesInSlice;
    this->CollectDownstreamPFParticles(pfParticleMap, slice.GetTargetHypothesis(), allParticlesInSlice);

    DepositionVector depositionVector;
    for (const auto &particle : allParticlesInSlice)
    {
        // Get the associated spacepoints
        const auto &partToSpacePointIter(pfParticleToSpacePointMap.find(particle));
        if (partToSpacePointIter == pfParticleToSpacePointMap.end())
            continue;

        for (const auto &spacePoint : partToSpacePointIter->second)
        {
            // Get the associated hit
            const auto &spacePointToHitIter(spacePointToHitMap.find(spacePoint));
            if (spacePointToHitIter == spacePointToHitMap.end())
                continue;

            // Only use hits from the collection plane
            const auto &hit(spacePointToHitIter->second);
            if (hit->View() != geo::kZ)
                continue;

            // Add the charged point to the vector
            const auto &position(spacePoint->XYZ());
            const auto charge(hit->Integral());

            depositionVector.emplace_back(position[0], position[1], position[2], charge, this->GetNPhotons(charge, particle));
        }
    }

    return depositionVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void FlashNeutrinoId::SliceCandidate::CollectDownstreamPFParticles(const PFParticleMap &pfParticleMap, const PFParticleVector &parentPFParticles,
                                                                   PFParticleVector &downstreamPFParticles) const
{
    for (const auto &particle : parentPFParticles)
        this->CollectDownstreamPFParticles(pfParticleMap, particle, downstreamPFParticles);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void FlashNeutrinoId::SliceCandidate::CollectDownstreamPFParticles(const PFParticleMap &pfParticleMap, const art::Ptr<recob::PFParticle> &particle,
                                                                   PFParticleVector &downstreamPFParticles) const
{
    if (std::find(downstreamPFParticles.begin(), downstreamPFParticles.end(), particle) == downstreamPFParticles.end())
        downstreamPFParticles.push_back(particle);

    for (const auto &daughterId : particle->Daughters())
    {
        const auto iter(pfParticleMap.find(daughterId));
        if (iter == pfParticleMap.end())
            throw cet::exception("FlashNeutrinoId") << "Scrambled PFParticle IDs" << std::endl;

        this->CollectDownstreamPFParticles(pfParticleMap, iter->second, downstreamPFParticles);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float FlashNeutrinoId::SliceCandidate::GetNPhotons(const float charge, const art::Ptr<recob::PFParticle> &particle) const
{
    return charge * (LArPandoraHelper::IsTrack(particle) ? m_chargeToNPhotonsTrack : m_chargeToNPhotonsShower);
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::CartesianVector FlashNeutrinoId::SliceCandidate::GetChargeWeightedCenter(const DepositionVector &depositionVector) const
{
    pandora::CartesianVector center(0.f, 0.f, 0.f);
    float totalCharge(0.f);

    for (const auto &chargePoint : depositionVector)
    {
        center += pandora::CartesianVector(chargePoint.m_x, chargePoint.m_y, chargePoint.m_z) * chargePoint.m_charge;
        totalCharge += chargePoint.m_charge;
    }

    if (totalCharge <= std::numeric_limits<float>::epsilon())
        throw cet::exception("FlashNeutrinoId") << "Can't find charge weighted center of slice with zero total charge" << std::endl;

    center *= (1.f / totalCharge);

    return center;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float FlashNeutrinoId::SliceCandidate::GetTotalCharge(const DepositionVector &depositionVector) const
{
    float totalCharge(0.f);

    for (const auto &chargePoint : depositionVector)
        totalCharge += chargePoint.m_charge;

    return totalCharge;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void FlashNeutrinoId::SliceCandidate::GetClosestCRTCosmic(const PFParticleVector &parentPFParticles, const art::Event &event, const PFParticlesToTracks &particlesToTracks,
                                                          const art::FindMany<anab::T0> &trk_t0_assn_v)
{
    m_numcosmictrack = 0;
    m_minCRTdist = 100; //Initialise on a value higher than we would call a match

    for (const art::Ptr<recob::PFParticle> pfp : parentPFParticles)
    {
        if (LArPandoraHelper::IsTrack(pfp))
        {
            if (particlesToTracks.count(pfp))
            {
                const art::Ptr<recob::Track> this_track = particlesToTracks.at(pfp).front();
                if (this_track->Length() > 20)
                {
                    m_numcosmictrack++;
                }
                const std::vector<const anab::T0 *> &T0_v = trk_t0_assn_v.at(this_track.key());
                if (T0_v.size() == 1)
                {
                    std::cout << "[GetClosestCRTCosmic] CRT-track-match found with dca: " << T0_v.front()->TriggerConfidence();
                    std::cout << "\tTime: " << T0_v.front()->Time();
                    std::cout << "\tPlane: " << T0_v.front()->TriggerBits() << std::endl;
                    if (T0_v.front()->TriggerConfidence() < m_minCRTdist)
                    {
                        m_minCRTdist = T0_v.front()->TriggerConfidence();
                        m_CRTplane = T0_v.front()->TriggerBits();
                        m_CRTtime = T0_v.front()->Time();
                        m_CRTtracklength = this_track->Length();
                    }
                }
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------DAVIDC

void FlashNeutrinoId::SliceCandidate::ACPTtagger(const PFParticleVector &parentPFParticles,
                                                 const art::Event &event,
                                                 const PFParticlesToTracks &particlesToTracks)
{

    mm_y_up = -9999.;
    mm_y_dn = -9999.;
    mm_x_up = -9999.;
    mm_x_dn = -9999.;
    mm_z_up = -9999.;
    mm_z_dn = -9999.;

    mm_ACPTdt = -9999.;

    mm_z_center = -9999.;

    mm_flash_timeanode_u = -9999.;
    mm_flash_timeanode_d = -9999.;
    mm_flash_timecathode_u = -9999.;
    mm_flash_timecathode_d = -9999.;

    // save OpHits in the event
    if (mm_verbose)
    {
        std::cout << "[ACPTTagger] \t Loading ophits from producer " << mm_ophitLabel << std::endl;
    }
    art::Handle<std::vector<recob::OpHit>> ophit_h;
    event.getByLabel(mm_ophitLabel, ophit_h);
    mm_ophit_v.clear();
    art::fill_ptr_vector(mm_ophit_v, ophit_h);

    // is the slice cosmic?
    bool isCosmic = false;

    for (const art::Ptr<recob::PFParticle> pfp : parentPFParticles)
    {
        if (LArPandoraHelper::IsTrack(pfp))
        {
            if (particlesToTracks.count(pfp))
            {
                const art::Ptr<recob::Track> this_track = particlesToTracks.at(pfp).front();

                float y_up = -9999;
                float y_dn = -9999;

                // Will store sorted points for the object [assuming downwards going]
                // The first vector is to consider end points estimated via different methods
                // (spacepoints, hits, tracks). The second vector has length==2, and
                // contains the start and end points of the track
                std::vector<TVector3> sorted_pts;
                sorted_pts.clear();

                if (this_track->Length() < mm_min_track_length)
                    continue;
                this->SortTrackPoints(*this_track, sorted_pts);
                if (sorted_pts.size() >= 2)
                {
                    //_trk_len.emplace_back(this_track->Length());
                    //_trk_x_up.emplace_back(pts[0].X());
                    //_trk_x_down.emplace_back(pts[sorted_pts.size()-1].X());
                    mm_z_center = sorted_pts[0].Z();
                    mm_z_center += sorted_pts[sorted_pts.size() - 1].Z();
                    mm_z_center /= 2.;
                    //_trk_z_center.emplace_back(z_center);

                    TVector3 start = sorted_pts[0];
                    TVector3 end = sorted_pts[sorted_pts.size() - 1];
                    sorted_pts.clear();
                    sorted_pts.resize(2);
                    sorted_pts.at(0) = start;
                    sorted_pts.at(1) = end;
                    y_up = sorted_pts.at(0).Y();
                    y_dn = sorted_pts.at(1).Y();
                } // if more then two points

                mm_x_up = sorted_pts.at(0).X();
                mm_x_dn = sorted_pts.at(1).X();
                mm_y_up = y_up;
                mm_y_dn = y_dn;
                mm_z_up = sorted_pts.at(0).Z();
                mm_z_dn = sorted_pts.at(1).Z();

                // can this track be a viable ACPT candidate?
                // reconstruct the track time assuming it hits the anode/cathode

                // check if track is compatible with OpHits
                if ((y_up != -9999 && y_dn != -9999) && ((y_up - y_dn) > 1.0) && (abs(y_up) > 0.001))
                {
                    float flash_zcenter, flash_time;
                    mm_ACPTdt = this->GetClosestDt_OpHits(sorted_pts, y_up, y_dn, mm_ophit_v, flash_zcenter, flash_time);
                    mm_flashTime = flash_time;
                    mm_flashZCenter = flash_zcenter;
                    if (mm_verbose)
                        std::cout << "[ACPTTagger] \t dt_ophits is " << mm_ACPTdt << std::endl;
                    if (mm_ACPTdt != -9999 && fabs(mm_ACPTdt) < mm_dt_resolution_ophit)
                    {
                        isCosmic = true;
                        if (mm_verbose)
                            std::cout << "[ACPTTagger] \t ===> Tagged! (ophit)" << std::endl;
                    }
                }

                float cosmicScore = 0;
                if (isCosmic)
                {
                    cosmicScore = 1;
                }

                if (mm_verbose)
                    std::cout << "Cosmic score is: " << cosmicScore << std::endl
                              << std::endl;

            } // if there is  track associated to this PFP
        }     // if this is a track
    }         // for all PFParticles

    if ((mm_y_up - mm_y_dn) < 1.)
        mm_ACPTdt = -9999.;

    return;
} // ACPT tagging

//------------------------------------------------------------------------------------------------------------------------------------------DAVIDC

void FlashNeutrinoId::SliceCandidate::SortTrackPoints(const recob::Track &track, std::vector<TVector3> &sorted_points)
{

    // vector to store 3D coordinates of
    // ordered track
    sorted_points.clear();

    // take the reconstructed 3D track
    // and assuming it is downwards
    // going, sort points so that
    // the track starts at the top
    // which point is further up in Y coord?
    // start or end?
    auto const &N = track.CountValidPoints();
    auto const &start = track.Vertex();
    auto const &end = track.End();

    if (mm_verbose)
    {
        std::cout << "[ACPTTagger] \t Track start " << start.X() << " " << start.Y() << " " << start.Z() << std::endl;
        std::cout << "[ACPTTagger] \t Track end   " << end.X() << " " << end.Y() << " " << end.Z() << std::endl;
    }

    // if points are ordered correctly
    if (start.Y() > end.Y())
    {
        for (size_t i = 0; i < N; i++)
            sorted_points.push_back(track.LocationAtPoint<TVector3>(track.NextValidPoint(i)));
    }

    // otherwise flip order
    else
    {
        if (mm_verbose)
            std::cout << "[ACPTTagger] \t\t These two points will be flipped" << std::endl;
        for (size_t i = 0; i < N; i++)
            sorted_points.push_back(track.LocationAtPoint<TVector3>(track.NextValidPoint(N - i - 1)));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------DAVIDC

float FlashNeutrinoId::SliceCandidate::GetClosestDt_OpHits(std::vector<TVector3> &sorted_points, double y_up, double y_down, const std::vector<art::Ptr<recob::OpHit>> mm_ophit_v,
                                                           float &flash_zcenter, float &flash_time)
{

    mm_flash_timeanode_u = sorted_points.at(0).X() / mm_driftVel - mm_anodeTime;
    mm_flash_timeanode_d = sorted_points.at(sorted_points.size() - 1).X() / mm_driftVel - mm_anodeTime;
    mm_flash_timecathode_u = sorted_points.at(0).X() / mm_driftVel - mm_cathodeTime;
    mm_flash_timecathode_d = sorted_points.at(sorted_points.size() - 1).X() / mm_driftVel - mm_cathodeTime;

    if (mm_verbose)
    {
        std::cout << "[ACPTTagger] >>> Using OpHits." << std::endl;
        std::cout << "[ACPTTagger] \t Estimated times we will be looking for: "
                  << "\t flash_time_anode_u: " << mm_flash_timeanode_u
                  << "\t flash_time_anode_d: " << mm_flash_timeanode_d
                  << "\t flash_time_cathode_u: " << mm_flash_timecathode_u
                  << "\t flash_time_cathode_d: " << mm_flash_timecathode_d << std::endl;
    }

    double trk_z_start = sorted_points.at(0).Z();
    double trk_z_end = sorted_points.at(sorted_points.size() - 1).Z();

    if (trk_z_start > trk_z_end)
        std::swap(trk_z_start, trk_z_end);

    bool sign = this->GetSign(sorted_points);
    if (mm_verbose)
        std::cout << "[ACPTTagger] \t Sign is " << (sign ? "positive." : "negative.") << std::endl;

    if (mm_verbose)
        std::cout << "[ACPTTagger] \t y_up = " << y_up << ", y_down = " << y_down << std::endl;
    bool upper_det = this->IsInUpperDet(y_up);
    bool lower_det = this->IsInLowerDet(y_down);

    std::vector<double> dt_v;
    dt_v.clear();
    std::vector<double> zpos_average_v;
    std::vector<double> time_average_v;
    float zpos_average, time_average;

    if (sign && upper_det)
    {
        if (mm_verbose)
            std::cout << "[ACPTTagger] \t Looking at cathode-down" << std::endl;
        dt_v.emplace_back(this->RunOpHitFinder(mm_flash_timecathode_d, trk_z_start, trk_z_end, mm_ophit_v, time_average, zpos_average));
        time_average_v.push_back(time_average);
        zpos_average_v.push_back(zpos_average);
    }
    if (sign && lower_det)
    {
        if (mm_verbose)
            std::cout << "[ACPTTagger] \t Looking at anode-up" << std::endl;
        dt_v.emplace_back(this->RunOpHitFinder(mm_flash_timeanode_u, trk_z_start, trk_z_end, mm_ophit_v, time_average, zpos_average));
        time_average_v.push_back(time_average);
        zpos_average_v.push_back(zpos_average);
    }
    if (!sign && upper_det)
    {
        if (mm_verbose)
            std::cout << "[ACPTTagger] \t Looking at anode-down" << std::endl;
        dt_v.emplace_back(this->RunOpHitFinder(mm_flash_timeanode_d, trk_z_start, trk_z_end, mm_ophit_v, time_average, zpos_average));
        time_average_v.push_back(time_average);
        zpos_average_v.push_back(zpos_average);
    }
    if (!sign && lower_det)
    {
        if (mm_verbose)
            std::cout << "[ACPTTagger] \t Looking at cathode-up" << std::endl;
        dt_v.emplace_back(this->RunOpHitFinder(mm_flash_timecathode_u, trk_z_start, trk_z_end, mm_ophit_v, time_average, zpos_average));
        time_average_v.push_back(time_average);
        zpos_average_v.push_back(zpos_average);
    }

    double min_dt = 1e9;
    double min_dta = 1e9; // absolute value of minimum dt
    bool min_dt_found = false;
    for (size_t t = 0; t < dt_v.size(); t++)
    {
        auto dt = dt_v[t];
        auto dta = fabs(dt);
        if (dta == -9999)
            continue;
        if (dta < min_dta)
        {
            min_dt = dt;
            min_dta = dta;
            min_dt_found = true;
            flash_zcenter = zpos_average_v[t];
            flash_time = time_average_v[t];
        }
    }

    if (mm_verbose && min_dt_found)
        std::cout << "[ACPTTagger] \t Found dt min from OpHits, dt min is " << min_dt << std::endl;

    if (min_dt_found)
        return min_dt;
    else
        return -9999;
}

//------------------------------------------------------------------------------------------------------------------------------------------DAVIDC

bool FlashNeutrinoId::SliceCandidate::GetSign(std::vector<TVector3> sorted_points)
{

    double t_down = sorted_points[sorted_points.size() - 1].X();
    double t_up = sorted_points[0].X();

    bool is_positive = (t_down - t_up) > 0.;

    return is_positive;
}

//------------------------------------------------------------------------------------------------------------------------------------------DAVIDC

bool FlashNeutrinoId::SliceCandidate::IsInUpperDet(double y_up)
{

    const art::ServiceHandle<geo::Geometry> geometry;
    if (y_up > geometry->DetHalfHeight() - mm_UP)
    {
        return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------DAVIDC

bool FlashNeutrinoId::SliceCandidate::IsInLowerDet(double y_down)
{

    const art::ServiceHandle<geo::Geometry> geometry;
    if (y_down < -geometry->DetHalfHeight() + mm_DOWN)
    {
        return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------DAVIDC

float FlashNeutrinoId::SliceCandidate::RunOpHitFinder(double the_time, double trk_z_start, double trk_z_end, const std::vector<art::Ptr<recob::OpHit>> mm_ophit_v,
                                                      float &time_average, float &zpos_average)
{

    ::art::ServiceHandle<geo::Geometry> geo;

    zpos_average = 0.;
    time_average = 0.;

    if (mm_verbose)
    {
        std::cout << "[ACPTTagger] \t _ophit_time_res is " << mm_ophit_time_res << std::endl;
        std::cout << "[ACPTTagger] \t _ophit_pos_res is " << mm_ophit_pos_res << std::endl;
        std::cout << "[ACPTTagger] \t _n_ophit is " << mm_nOphit << std::endl;
        std::cout << "[ACPTTagger] \t _ophit_pe is " << mm_ophitPE << std::endl;
    }

    std::vector<double> ophit_sel_time;
    std::vector<double> ophit_sel_zpos;
    std::vector<double> ophit_sel_pe;

    ophit_sel_time.clear();
    ophit_sel_pe.clear();

    for (auto oh : mm_ophit_v)
    {

        double time_diff = std::abs(oh->PeakTime() - the_time);

        if (!geo->IsValidOpChannel(oh->OpChannel()))
            continue;
        if (oh->OpChannel() < 200 || oh->OpChannel() > 231)
            continue;

        //size_t opdet = geo->OpDetFromOpChannel(oh->OpChannel());

        double pmt_xyz[3] = {-9999, -9999, -9999};
        geo->OpDetGeoFromOpChannel(oh->OpChannel()).GetCenter(pmt_xyz);
        double pmt_z = pmt_xyz[2];

        double dz = 1e9;
        if (pmt_z > trk_z_start && pmt_z < trk_z_end)
        {
            dz = 0.;
        }
        else
        {
            if (pmt_z < trk_z_start)
                dz = std::abs(pmt_z - trk_z_start);
            if (pmt_z > trk_z_end)
                dz = std::abs(pmt_z - trk_z_end);
        }

        auto ophitPE = oh->Area();

        if (time_diff < mm_ophit_time_res && dz < mm_ophit_pos_res)
        {
            if (mm_verbose)
                std::cout << "[ACPTTagger] \t\t Found ophit, time is " << oh->PeakTime()
                          << ", pmt_z is " << pmt_z
                          << ", dz is " << dz
                          << ", opchannel is " << oh->OpChannel()
                          << ", PE is " << ophitPE << std::endl;

            ophit_sel_time.emplace_back(oh->PeakTime());
            ophit_sel_zpos.emplace_back(pmt_z);
            ophit_sel_pe.emplace_back(ophitPE);
        }
    }

    if (ophit_sel_time.size() < mm_nOphit)
    {
        if (mm_verbose)
            std::cout << "[ACPTTagger] \t Not enough ophits" << std::endl;
        return -9999;
    }

    double total_pe = std::accumulate(ophit_sel_pe.begin(), ophit_sel_pe.end(), 0.);

    if (total_pe < mm_ophitPE)
    {
        if (mm_verbose)
            std::cout << "[ACPTTagger] \t Not enough pe" << std::endl;
        return -9999;
    }

    // Calculate and return average time
    for (size_t i = 0; i < ophit_sel_time.size(); i++)
    {
        time_average += ophit_sel_time.at(i) * ophit_sel_pe.at(i);
        zpos_average += ophit_sel_zpos.at(i) * ophit_sel_pe.at(i);
    }

    time_average /= total_pe;
    zpos_average /= total_pe;

    if (mm_verbose)
        std::cout << "[ACPTTagger] \t time_average: " << time_average << ", the_time: " << the_time << std::endl;

    return time_average - the_time;
}

//------------------------------------------------------------------------------------------------------------------------------------------

flashana::QCluster_t FlashNeutrinoId::SliceCandidate::GetLightCluster(const DepositionVector &depositionVector) const
{
    flashana::QCluster_t lightCluster;

    for (const auto &point : depositionVector)
        lightCluster.emplace_back(point.m_x, point.m_y, point.m_z, point.m_nPhotons);

    return lightCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool FlashNeutrinoId::SliceCandidate::IsCompatibleWithBeamFlash(const FlashCandidate &beamFlash, const float maxDeltaY,
                                                                const float maxDeltaZ, const float maxDeltaYSigma, const float maxDeltaZSigma, const float minChargeToLightRatio,
                                                                const float maxChargeToLightRatio)
{
    // Check the flash is usable
    if (beamFlash.m_totalPE <= std::numeric_limits<float>::epsilon())
    {
        return false;
    }

    if (beamFlash.m_widthY <= std::numeric_limits<float>::epsilon())
    {
        return false;
    }

    if (beamFlash.m_widthZ <= std::numeric_limits<float>::epsilon())
    {
        return false;
    }

    if (m_totalCharge <= std::numeric_limits<float>::epsilon())
    {
        return false;
    }

    // Calculate the pre-selection variables
    m_deltaY = (m_centerY - beamFlash.m_centerY);
    m_deltaZ = (m_centerZ - beamFlash.m_centerZ);
    m_deltaYSigma = m_deltaY / beamFlash.m_widthY;
    m_deltaZSigma = m_deltaZ / beamFlash.m_widthZ;
    m_chargeToLightRatio = m_totalCharge / beamFlash.m_totalPE;
    m_xChargeLightVariable = m_xclCoef * log10(m_chargeToLightRatio) - m_centerX;

    // Check if the slice passes the pre-selection cuts
    m_passesPrecuts = (std::abs(m_deltaY) < maxDeltaY &&
                       std::abs(m_deltaZ) < maxDeltaZ &&
                       std::abs(m_deltaYSigma) < maxDeltaYSigma &&
                       std::abs(m_deltaZSigma) < maxDeltaZSigma &&
                       m_xChargeLightVariable > minChargeToLightRatio &&
                       m_xChargeLightVariable < maxChargeToLightRatio &&
                       (fabs(mm_ACPTdt) > mm_dt_resolution_ophit)); // DAVIDC

    return m_passesPrecuts;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float FlashNeutrinoId::SliceCandidate::GetFlashMatchScore(const FlashCandidate &beamFlash, flashana::FlashMatchManager &flashMatchManager)
{
    flashMatchManager.Reset();

    // Convert the flash and the charge cluster into the required format for flash matching
    auto flash(beamFlash.ConvertFlashFormat());

    // Perform the match
    flashMatchManager.Emplace(std::move(flash));
    flashMatchManager.Emplace(std::move(m_lightCluster));
    const auto matches(flashMatchManager.Match());

    // Unable to match
    if (matches.empty())
        return -1.f;

    if (matches.size() != 1)
        throw cet::exception("FlashNeutrinoId") << "Flash matching returned multiple matches!" << std::endl;

    // Fill the slice candidate with the details of the matching
    const auto match(matches.front());

    m_flashMatchScore = match.score;
    m_totalPEHypothesis = std::accumulate(match.hypothesis.begin(), match.hypothesis.end(), 0.f);

    //std::cout << "DAVIDC Slice has " << m_lightCluster.size() << " spacepoints and flash-match score of " << m_flashMatchScore << std::endl;

    // Fill the slice with the hypothesized PE spectrum
    if (!m_peHypothesisSpectrum.empty())
        throw cet::exception("FlashNeutrinoId") << "Hypothesized PE spectrum already set for this flash" << std::endl;

    for (auto hypo_pe : match.hypothesis)
        m_peHypothesisSpectrum.push_back(static_cast<float>(hypo_pe));

    return m_flashMatchScore;
}

void FlashNeutrinoId::GetBestObviousCosmicMatch(const art::Event &event, const FlashCandidate &beamFlash)
{
    float bestCosmicMatch = -1;
    std::vector<float> cosmicMatchHypothesis = {};
    bool foundCosmic = false;

    PFParticleVector pfParticles;
    PFParticlesToMetadata particlesToMetadata;
    SpacePointVector spacePoints;
    SpacePointsToHits spacePointToHitMap;
    PFParticleMap pfParticleMap;
    PFParticlesToSpacePoints pfParticleToSpacePointMap;
    LArPandoraHelper::CollectPFParticles(event, m_pandoraLabel, pfParticles, pfParticleToSpacePointMap);
    LArPandoraHelper::CollectSpacePoints(event, m_pandoraLabel, spacePoints, spacePointToHitMap);
    LArPandoraHelper::BuildPFParticleMap(pfParticles, pfParticleMap);
    LArPandoraHelper::CollectPFParticleMetadata(event, m_pandoraLabel, pfParticles, particlesToMetadata);

    m_flashMatchManager.Reset();
    // Convert the flash and the charge cluster into the required format for flash matching
    auto flash(beamFlash.ConvertFlashFormat());
    // Perform the match
    m_flashMatchManager.Emplace(std::move(flash));

    if (pfParticles.size() == 0)
    {
        std::cout << "There were no PFParticles in the event!" << std::endl;
        return;
    }

    for (const art::Ptr<recob::PFParticle> &pfp : pfParticles)
    {

        MetadataVector pfp_metadata_vec = particlesToMetadata.at(pfp);
        const larpandoraobj::PFParticleMetadata::PropertiesMap &pfp_properties = pfp_metadata_vec.front()->GetPropertiesMap();

        if (pfp_properties.count("IsClearCosmic"))
        {
            if (pfp_properties.at("IsClearCosmic") && pfp->IsPrimary())
            {

                PFParticleVector downstreamPFParticles;
                CollectDownstreamPFParticles(pfParticleMap, pfp, downstreamPFParticles);
                flashana::QCluster_t lightCluster;
                for (const auto &particle : downstreamPFParticles)
                {

                    // Get the associated spacepoints
                    const auto &partToSpacePointIter(pfParticleToSpacePointMap.find(particle));
                    if (partToSpacePointIter == pfParticleToSpacePointMap.end())
                        continue;

                    for (const auto &spacePoint : partToSpacePointIter->second)
                    {

                        // Get the associated hit
                        const auto &spacePointToHitIter(spacePointToHitMap.find(spacePoint));
                        if (spacePointToHitIter == spacePointToHitMap.end())
                            continue;

                        // Only use hits from the collection plane
                        const auto &hit(spacePointToHitIter->second);
                        if (hit->View() != geo::kZ)
                            continue;

                        // Add the charged point to the vector
                        const auto &position(spacePoint->XYZ());
                        const auto charge(hit->Integral());
                        lightCluster.emplace_back(position[0], position[1], position[2], charge * (LArPandoraHelper::IsTrack(particle) ? m_chargeToNPhotonsTrack : m_chargeToNPhotonsShower));
                    }
                }
                m_flashMatchManager.Emplace(std::move(lightCluster));
                foundCosmic = true;
            }
        }
    }

    if (foundCosmic)
    {
        const auto matches(m_flashMatchManager.Match());
        if (!matches.empty())
        {
            const auto match(matches.back());
            bestCosmicMatch = match.score;
            for (auto hypo_pe : match.hypothesis)
                cosmicMatchHypothesis.push_back(static_cast<float>(hypo_pe));
        }
        std::cout << "[FlashNeutrinoId] Chi2 best cosmic (out of " << matches.size() << " matches): " << matches.back().score << "\tworst match: " << matches.front().score << std::endl;
        m_outputEvent.m_bestCosmicMatch = bestCosmicMatch;
        m_outputEvent.m_cosmicMatchHypothesis = cosmicMatchHypothesis;
    }
    else
    {
        std::cout << "There were no Obvious cosmics in the event!" << std::endl;
        return;
    }
}

void FlashNeutrinoId::CollectDownstreamPFParticles(const PFParticleMap &pfParticleMap, const art::Ptr<recob::PFParticle> &particle,
                                                   PFParticleVector &downstreamPFParticles) const
{
    if (std::find(downstreamPFParticles.begin(), downstreamPFParticles.end(), particle) == downstreamPFParticles.end())
        downstreamPFParticles.push_back(particle);

    for (const auto &daughterId : particle->Daughters())
    {
        const auto iter(pfParticleMap.find(daughterId));
        if (iter == pfParticleMap.end())
            throw cet::exception("FlashNeutrinoId") << "Scrambled PFParticle IDs" << std::endl;

        this->CollectDownstreamPFParticles(pfParticleMap, iter->second, downstreamPFParticles);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void FlashNeutrinoId::SliceCandidate::RejectStopMuByDirMCS(const PFParticleVector &parentPFParticles, const art::Event &event,
                                                           const PFParticlesToTracks &particlesToTracks,
                                                           const trkf::TrajectoryMCSFitter &mcsfitter)
{

    if (mm_verbose)
        std::cout << "[RejectStopMuByDirMCS] Slice with N pfps = " << parentPFParticles.size() << std::endl;

    ::art::ServiceHandle<geo::Geometry> geo;
    float bnd = 20.;
    for (const art::Ptr<recob::PFParticle> pfp : parentPFParticles)
    {
        if (LArPandoraHelper::IsTrack(pfp))
        {
            if (particlesToTracks.count(pfp))
            {
                const art::Ptr<recob::Track> this_track = particlesToTracks.at(pfp).front();

                auto InFV = [&geo, bnd](recob::tracking::Point_t p) -> bool { return (p.X() > bnd && p.X() < (2. * geo->DetHalfWidth() - bnd) &&
                                                                                      p.Y() > (-geo->DetHalfHeight() + bnd) && p.Y() < (geo->DetHalfHeight() - bnd) &&
                                                                                      p.Z() > bnd && p.Z() < (geo->DetLength() - bnd)); };

                bool vtx_contained = InFV(this_track->Vertex());
                bool end_contained = InFV(this_track->End());

                // If fully contained, exit
                if (vtx_contained && end_contained)
                {
                    continue;
                }
                // If fully uncontained, exit
                if (!vtx_contained && !end_contained)
                {
                    continue;
                }

                auto _result = mcsfitter.fitMcs(*this_track);
                double fwd_ll = _result.fwdLogLikelihood();
                double bwd_ll = _result.bwdLogLikelihood();

                if (mm_verbose)
                {
                    std::cout << "[RejectStopMuByDirMCS] FWD " << fwd_ll
                              << ", BWD " << bwd_ll
                              << ", length=" << this_track->Length()
                              << ", vtx=" << this_track->Vertex()
                              << ", end=" << this_track->End()
                              << std::endl;
                }

                if (vtx_contained && !end_contained)
                {
                    m_maxDeltaLLMCS = std::max(float(fwd_ll - bwd_ll), m_maxDeltaLLMCS);
                    m_lengthDeltaLLMCS = std::max(float(this_track->Length()), m_lengthDeltaLLMCS);
                }
                else if (!vtx_contained && end_contained)
                {
                    m_maxDeltaLLMCS = std::max(float(bwd_ll - fwd_ll), m_maxDeltaLLMCS);
                    m_lengthDeltaLLMCS = std::max(float(this_track->Length()), m_lengthDeltaLLMCS);
                }
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void FlashNeutrinoId::SliceCandidate::RejectStopMuByCalo(const PFParticleVector &pfp_v, const art::Event &event, const PFParticlesToTracks &pfps_to_tracks, const PFParticlesToSpacePoints &pfps_to_spacepoints, std::string &pandoraLabel, fhicl::ParameterSet &cosmictagmanager)
{

    if (mm_verbose)
        std::cout << "[RejectStopMuByCalo] Slice with N pfps = " << pfp_v.size() << std::endl;

    ::art::ServiceHandle<geo::Geometry> geo;
    float bnd = 20.;

    // Declare fiducial volume - need this for later (copied from RejectStopMuByDirMCS above)
    auto InFV = [&geo, bnd](double p[3]) -> bool { return (p[0] > bnd && p[0] < (2. * geo->DetHalfWidth() - bnd) && p[1] > (-geo->DetHalfHeight() + bnd) && p[1] < (geo->DetHalfHeight() - bnd) && p[2] > bnd && p[2] < (geo->DetLength() - bnd)); };

    // Configure cosmic tag manager
    ::cosmictag::CosmicTagManager _ct_manager;
    _ct_manager.Configure(cosmictagmanager);

    // Detector properties
    ::detinfo::DetectorProperties const *fDetectorProperties;
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();

    // These three are needed for later
    art::Ptr<recob::PFParticle> primary_pfp;
    std::vector<art::Ptr<recob::Track>> primary_track_v;
    bool ignore_this = false;

    // Set results to false (not tagged as stopping muon) by default
    m_ct_result_bragg_plane0 = false;
    m_ct_result_bragg_plane1 = false;
    m_ct_result_bragg_plane2 = false;
    m_ct_result_michel_plane0 = false;
    m_ct_result_michel_plane1 = false;
    m_ct_result_michel_plane2 = false;

    // Map pfp->cluster
    lar_pandora::PFParticleVector pfp_v_copy = pfp_v;
    lar_pandora::PFParticlesToClusters pfps_to_clusters;
    lar_pandora::LArPandoraHelper::CollectPFParticles(event, pandoraLabel, pfp_v_copy, pfps_to_clusters);

    // Collect clusters and map cluster->hit
    lar_pandora::ClusterVector cluster_v;
    lar_pandora::ClustersToHits clusters_to_hits;
    lar_pandora::LArPandoraHelper::CollectClusters(event, pandoraLabel, cluster_v, clusters_to_hits);

    // Get all the SpacePoints from the PFPs (this will be used to look for
    // the highest point in the TPCObject). Also get all the clusters for this
    // TPCObject, and use the clusters to go to the hits. Collect all
    // these hits in a vector.
    std::vector<art::Ptr<recob::SpacePoint>> sp_v;
    sp_v.clear();

    std::vector<art::Ptr<recob::Hit>> hit_v;
    hit_v.clear();

    if (pfp_v.size() == 0)
    {
        ignore_this = true;
    }

    for (auto p : pfp_v)
    {
        auto iter = pfps_to_spacepoints.find(p);
        if (iter == pfps_to_spacepoints.end())
        {
            continue;
        }
        sp_v.reserve(sp_v.size() + iter->second.size());
        sp_v.insert(sp_v.end(), iter->second.begin(), iter->second.end());

        // Find clusters first ...
        auto iter2 = pfps_to_clusters.find(p);
        if (iter2 == pfps_to_clusters.end())
        {
            continue;
        }

        // ... then find hits
        for (auto c : iter2->second)
        {
            auto iter3 = clusters_to_hits.find(c);
            if (iter3 == clusters_to_hits.end())
            {
                if (mm_verbose)
                    std::cout << "[StoppingMuonTagger] Cluster not found by pandora !?" << std::endl;
                throw std::exception();
            }

            hit_v.reserve(hit_v.size() + iter3->second.size());
            hit_v.insert(hit_v.end(), iter3->second.begin(), iter3->second.end());
        }

        if (p->IsPrimary() && !lar_pandora::LArPandoraHelper::IsNeutrino(p))
        {

            auto iter4 = pfps_to_tracks.find(p);
            if (iter4 == pfps_to_tracks.end())
            {
                // /*if (_debug)*/ std::cout << "[StoppingMuonTagger] PFParticle not found by pandora !?" << std::endl;
                // throw cet::exception("FlashNeutrinoId") << "[StoppingMuonTagger] PFParticle not found by pandora !?" << std::endl;
                return;
            }

            primary_pfp = p;
            primary_track_v = iter4->second;
        }
    }

    if (hit_v.size() == 0)
    {
        ignore_this = true;
    }

    if (ignore_this || !primary_pfp)
    {
        return;
    }

    if (mm_verbose)
        std::cout << "[StoppingMuonTagger] Primary PFP is " << primary_pfp->Self() << std::endl;

    //
    // First exclude spacepoints outside the tpc
    //
    std::vector<art::Ptr<recob::SpacePoint>> temp;
    ::geoalgo::AABox tpcvol(0., (-1.) * geo->DetHalfHeight(),
                            0., geo->DetHalfWidth() * 2,
                            geo->DetHalfHeight(), geo->DetLength());

    for (auto s : sp_v)
    {
        const double *xyz = s->XYZ();
        ::geoalgo::Vector point(xyz[0], xyz[1], xyz[2]);
        if (tpcvol.Contain(point))
        {
            temp.push_back(s);
        }
    }
    sp_v = temp;

    //
    // Now get the highest point
    //
    std::sort(sp_v.begin(), sp_v.end(),
              [](art::Ptr<recob::SpacePoint> a, art::Ptr<recob::SpacePoint> b) -> bool {
                  const double *xyz_a = a->XYZ();
                  const double *xyz_b = b->XYZ();
                  return xyz_a[1] > xyz_b[1];
              });

    if (sp_v.size() == 0)
    {
        if (mm_verbose)
            std::cout << "[StoppingMuonTagger] Not enough spacepoints." << std::endl;
        return;
    }

    const double *highest_point_c = sp_v.at(0)->XYZ();
    double highest_point[3] = {highest_point_c[0], highest_point_c[1], highest_point_c[2]};

    // Find highest point in the detector (also called "containing" the point)
    double x = highest_point[0];
    double y = highest_point[1];
    double z = highest_point[2];
    double e = std::numeric_limits<double>::epsilon();

    if (x < 0. + e)
        highest_point[0] = 0. + e;
    if (x > 2. * geo->DetHalfWidth() - e)
        highest_point[0] = 2. * geo->DetHalfWidth() - e;

    if (y < -geo->DetHalfWidth() + e)
        highest_point[1] = -geo->DetHalfWidth() + e;
    if (y > geo->DetHalfWidth() - e)
        highest_point[1] = geo->DetHalfWidth() - e;

    if (z < 0. + e)
        highest_point[2] = 0. + e;
    if (z > geo->DetLength() - e)
        highest_point[2] = geo->DetLength() - e;

    // Create an approximate start hit on plane 0
    int highest_w = geo->NearestWire(highest_point, 0);
    double highest_t = fDetectorProperties->ConvertXToTicks(highest_point[0], geo::PlaneID(0, 0, 0)) / 4.;
    if (mm_verbose)
        std::cout << "[StoppingMuonTagger] Highest point: wire: " << geo->NearestWire(highest_point, 0)
                  << ", time: " << fDetectorProperties->ConvertXToTicks(highest_point[0], geo::PlaneID(0, 0, 0))
                  << std::endl;

    cosmictag::SimpleHit start_highest_plane0;
    start_highest_plane0.time = highest_t;
    start_highest_plane0.wire = highest_w;
    start_highest_plane0.plane = 0;

    // Create an approximate start hit on plane 1
    highest_w = geo->NearestWire(highest_point, 1);
    highest_t = fDetectorProperties->ConvertXToTicks(highest_point[0], geo::PlaneID(0, 0, 1)) / 4.;
    if (mm_verbose)
        std::cout << "[StoppingMuonTagger] Highest point: wire: " << geo->NearestWire(highest_point, 1)
                  << ", time: " << fDetectorProperties->ConvertXToTicks(highest_point[0], geo::PlaneID(0, 0, 1))
                  << std::endl;

    cosmictag::SimpleHit start_highest_plane1;
    start_highest_plane1.time = highest_t;
    start_highest_plane1.wire = highest_w;
    start_highest_plane1.plane = 1;

    // Create an approximate start hit on plane 2
    highest_w = geo->NearestWire(highest_point, 2);
    highest_t = fDetectorProperties->ConvertXToTicks(highest_point[0], geo::PlaneID(0, 0, 2)) / 4.;
    if (mm_verbose)
        std::cout << "[StoppingMuonTagger] Highest point: wire: " << geo->NearestWire(highest_point, 2)
                  << ", time: " << fDetectorProperties->ConvertXToTicks(highest_point[0], geo::PlaneID(0, 0, 2))
                  << std::endl;

    cosmictag::SimpleHit start_highest_plane2;
    start_highest_plane2.time = highest_t;
    start_highest_plane2.wire = highest_w;
    start_highest_plane2.plane = 2;

    if (mm_verbose)
        std::cout << "[StoppingMuonTagger] Now create simple hit vector, size " << hit_v.size() << std::endl;

    //
    // Create SimpleHit vector for hits in each plane
    //
    std::vector<cosmictag::SimpleHit> simple_hit_v_plane0;
    std::vector<cosmictag::SimpleHit> simple_hit_v_plane1;
    std::vector<cosmictag::SimpleHit> simple_hit_v_plane2;
    for (auto h : hit_v)
    {

        cosmictag::SimpleHit sh;

        sh.t = fDetectorProperties->ConvertTicksToX(h->PeakTime(), geo::PlaneID(0, 0, h->View()));
        sh.w = h->WireID().Wire * geo->WirePitch(geo::PlaneID(0, 0, h->View()));

        sh.plane = h->View();
        sh.integral = h->Integral();

        sh.time = h->PeakTime() / 4;
        sh.wire = h->WireID().Wire;

        if (h->View() == 0)
        {
            simple_hit_v_plane0.emplace_back(sh);
        }
        else if (h->View() == 1)
        {
            simple_hit_v_plane1.emplace_back(sh);
        }
        else if (h->View() == 2)
        {
            simple_hit_v_plane2.emplace_back(sh);
        }
    }

    if (mm_verbose)
        std::cout << "[StoppingMuonTagger] Simple hit vector size " << simple_hit_v_plane0.size() << ", " << simple_hit_v_plane1.size() << ", " << simple_hit_v_plane2.size() << std::endl;

    //
    // Running with highest point as start hit
    //

    // --- Plane 0 ---
    _ct_manager.Reset();
    // Emplacing simple hits to the manager
    cosmictag::SimpleCluster sc_plane0(simple_hit_v_plane0);
    _ct_manager.Emplace(std::move(sc_plane0));
    _ct_manager.SetStartHit(std::move(start_highest_plane0));
    // Running the cluster analyser
    bool passed = _ct_manager.Run();

    if (passed)
    {
        cosmictag::SimpleCluster processed_cluster = _ct_manager.GetCluster();

        m_ct_result_michel_plane0 = ((cosmictag::StopMuMichel *)(_ct_manager.GetCustomAlgo("StopMuMichel")))->IsStopMuMichel(processed_cluster, m_dqds_startend_percdiff_plane0, m_bragg_local_lin_plane0, m_n_michel_hits_plane0);

        bool vtx_in_fv = InFV(highest_point);
        m_ct_result_bragg_plane0 = ((cosmictag::StopMuBragg *)(_ct_manager.GetCustomAlgo("StopMuBragg")))->IsStopMuBragg(processed_cluster, m_min_lin_braggalgonly_plane0) && !vtx_in_fv;
    }

    // --- Plane 1 ---
    _ct_manager.Reset();
    // Emplacing simple hits to the manager
    cosmictag::SimpleCluster sc_plane1(simple_hit_v_plane1);
    _ct_manager.Emplace(std::move(sc_plane1));
    _ct_manager.SetStartHit(std::move(start_highest_plane1));
    // Running the cluster analyser
    passed = _ct_manager.Run();

    if (passed)
    {
        cosmictag::SimpleCluster processed_cluster = _ct_manager.GetCluster();

        m_ct_result_michel_plane1 = ((cosmictag::StopMuMichel *)(_ct_manager.GetCustomAlgo("StopMuMichel")))->IsStopMuMichel(processed_cluster, m_dqds_startend_percdiff_plane1, m_bragg_local_lin_plane1, m_n_michel_hits_plane1);

        bool vtx_in_fv = InFV(highest_point);
        m_ct_result_bragg_plane1 = ((cosmictag::StopMuBragg *)(_ct_manager.GetCustomAlgo("StopMuBragg")))->IsStopMuBragg(processed_cluster, m_min_lin_braggalgonly_plane1) && !vtx_in_fv;
    }

    // --- Plane 2 ---
    _ct_manager.Reset();
    // Emplacing simple hits to the manager
    cosmictag::SimpleCluster sc_plane2(simple_hit_v_plane2);
    _ct_manager.Emplace(std::move(sc_plane2));
    _ct_manager.SetStartHit(std::move(start_highest_plane2));
    // Running the cluster analyser
    passed = _ct_manager.Run();

    if (passed)
    {
        cosmictag::SimpleCluster processed_cluster = _ct_manager.GetCluster();

        m_ct_result_michel_plane2 = ((cosmictag::StopMuMichel *)(_ct_manager.GetCustomAlgo("StopMuMichel")))->IsStopMuMichel(processed_cluster, m_dqds_startend_percdiff_plane2, m_bragg_local_lin_plane2, m_n_michel_hits_plane2);

        bool vtx_in_fv = InFV(highest_point);
        m_ct_result_bragg_plane2 = ((cosmictag::StopMuBragg *)(_ct_manager.GetCustomAlgo("StopMuBragg")))->IsStopMuBragg(processed_cluster, m_min_lin_braggalgonly_plane2) && !vtx_in_fv;
    }

    if (mm_verbose)
    {
        std::cout << "[RejectStopMuByCalo] result_michel = " << m_ct_result_michel_plane0 << ", " << m_ct_result_michel_plane1 << ", " << m_ct_result_michel_plane2 << std::endl;
        std::cout << "[RejectStopMuByCalo] result_bragg = " << m_ct_result_bragg_plane0 << ", " << m_ct_result_bragg_plane1 << ", " << m_ct_result_bragg_plane2 << std::endl;
    }

} // void FlashNeutrinoId::SliceCandidate::RejectStopMuByCalo

} // namespace lar_pandora
