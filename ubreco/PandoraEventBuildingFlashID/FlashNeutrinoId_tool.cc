/**
 *  @file   larpandora/LArPandoraEventBuilding/SliceIdTools/FlashNeutrinoId_tool.cc
 *
 *  @brief  implementation of the flash based neutrino id tool
 */

#include "FlashNeutrinoId_tool.h"

//------------------------------------------------------------------------------------------------------------------------------------------

void FlashNeutrinoId::ClassifySlices(SliceVector &slices, const art::Event &evt)
{
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
            //// WOUTER: If in this sequence, cosmci mtching will always run, if turned around, cosmic matching will only run if event is selected.
            // Obvious-Cosmic Beam-Flash Matching
            GetBestObviousCosmicMatch(evt, beamFlash);
            // Find the slice - if any that matches best with the beamFlash
            bestSliceIndex = this->GetBestSliceIndex(beamFlash, sliceCandidates);
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

    std::cout << "Flash PE: " << brightestFlash.m_totalPE << std::endl;

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
    LArPandoraHelper::CollectTracks(event, "pandoraAllOutcomesTrack", pftracks, particlesToTracks);
    art::Handle<std::vector<recob::Track>> track_h;
    event.getByLabel("pandoraAllOutcomesTrack", track_h);

    for (unsigned int sliceIndex = 0; sliceIndex < slices.size(); ++sliceIndex)
    {
        const auto &slice = slices[sliceIndex];
        if (m_hasCRT)
        {
            const art::FindMany<anab::T0> trk_t0_assn_v(track_h, event, "trackmatch");
            sliceCandidates.emplace_back(event, slice, pfParticleMap, pfParticleToSpacePointMap, spacePointToHitMap, particlesToTracks,
                                         trk_t0_assn_v, m_chargeToNPhotonsTrack, m_chargeToNPhotonsShower, m_xclCoef, sliceIndex + 1);
        }
        else
        {
            sliceCandidates.emplace_back(event, slice, pfParticleMap, pfParticleToSpacePointMap, spacePointToHitMap, particlesToTracks,
                                         m_chargeToNPhotonsTrack, m_chargeToNPhotonsShower, m_xclCoef, sliceIndex + 1);
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

        m_nuInteractionType = mcNeutrino.InteractionType();
        m_nuCCNC = mcNeutrino.CCNC();
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

FlashNeutrinoId::SliceCandidate::SliceCandidate() : m_sliceId(-std::numeric_limits<int>::max()),
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
                                                    m_xclCoef(-std::numeric_limits<float>::max())
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

FlashNeutrinoId::SliceCandidate::SliceCandidate(const art::Event &event, const Slice &slice) : m_sliceId(-std::numeric_limits<int>::max()),
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
                                                                                               m_xclCoef(-std::numeric_limits<float>::max())
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

FlashNeutrinoId::SliceCandidate::SliceCandidate(const art::Event &event, const Slice &slice, const PFParticleMap &pfParticleMap,
                                                const PFParticlesToSpacePoints &pfParticleToSpacePointMap, const SpacePointsToHits &spacePointToHitMap,
                                                const PFParticlesToTracks &particlesToTracks,
                                                const float chargeToNPhotonsTrack, const float chargeToNPhotonsShower, const float xclCoef, const int sliceId) : m_sliceId(sliceId),
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
                                                                                                                                                                 m_xclCoef(xclCoef)

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
}

//------------------------------------------------------------------------------------------------------------------------------------------

FlashNeutrinoId::SliceCandidate::SliceCandidate(const art::Event &event, const Slice &slice, const PFParticleMap &pfParticleMap,
                                                const PFParticlesToSpacePoints &pfParticleToSpacePointMap, const SpacePointsToHits &spacePointToHitMap,
                                                const PFParticlesToTracks &particlesToTracks, const art::FindMany<anab::T0> &trk_t0_assn_v,
                                                const float chargeToNPhotonsTrack, const float chargeToNPhotonsShower, const float xclCoef, const int sliceId) : FlashNeutrinoId::SliceCandidate::SliceCandidate(event, slice, pfParticleMap,
                                                                                                                                                                                                                 pfParticleToSpacePointMap, spacePointToHitMap,
                                                                                                                                                                                                                 particlesToTracks, chargeToNPhotonsTrack, chargeToNPhotonsShower, xclCoef, sliceId)

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

void FlashNeutrinoId::SliceCandidate::GetClosestCRTCosmic(const PFParticleVector &parentPFParticles, const art::Event &event,
                                                          const PFParticlesToTracks &particlesToTracks, const art::FindMany<anab::T0> &trk_t0_assn_v)
{
    for (const art::Ptr<recob::PFParticle> pfp : parentPFParticles)
    {
        if (LArPandoraHelper::IsTrack(pfp))
        {
            if (particlesToTracks.count(pfp))
            {
                const art::Ptr<recob::Track> this_track = particlesToTracks.at(pfp).front();
                const std::vector<const anab::T0 *> &T0_v = trk_t0_assn_v.at(this_track.key());
                if (T0_v.size() == 1)
                {
                    std::cout << "CRT-track-match found with dca: " << T0_v.front()->TriggerConfidence();
                    std::cout << "\tTime: " << T0_v.front()->Time();
                    std::cout << "\tPlane: " << T0_v.front()->TriggerBits() << std::endl;
                    if (T0_v.front()->TriggerConfidence() < m_minCRTdist)
                    {
                        m_minCRTdist = T0_v.front()->TriggerConfidence();
                        m_CRTplane = T0_v.front()->TriggerBits();
                        m_CRTtime = T0_v.front()->Time();
                    }
                }
            }
        }
    }
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
                       m_xChargeLightVariable < maxChargeToLightRatio);

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
                        m_flashMatchManager.Emplace(std::move(lightCluster));
                    }
                }
            }
        }
    }

    const auto matches(m_flashMatchManager.Match());
    if (!matches.empty())
    {
        const auto match(matches.back());
        bestCosmicMatch = match.score;
        for (auto hypo_pe : match.hypothesis)
            cosmicMatchHypothesis.push_back(static_cast<float>(hypo_pe));
    }
    std::cout << "[FlashNeutrinoId] Chi2 best matching cosmic: " << matches.back().score << "\tworst match: " << matches.front().score << std::endl;
    m_outputEvent.m_bestCosmicMatch = bestCosmicMatch;
    m_outputEvent.m_cosmicMatchHypothesis = cosmicMatchHypothesis;
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

} // namespace lar_pandora
