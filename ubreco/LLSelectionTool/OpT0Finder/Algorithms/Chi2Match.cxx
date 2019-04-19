#ifndef QHI2MATCH_CXX
#define QHI2MATCH_CXX

#include "Chi2Match.h"

namespace flashana
{

static Chi2MatchFactory __global_Chi2MatchFactory__;

void MIN_vtx_qll(Int_t &, Double_t *, Double_t &, Double_t *, Int_t);

Chi2Match::Chi2Match(const std::string name)
    : BaseFlashMatch(name), _normalize(false)
{
  _current_llhd = _current_chi2 = -1.0;
}

Chi2Match::Chi2Match()
{
  throw OpT0FinderException("Use the other ctor");
}

void Chi2Match::_Configure_(const Config_t &pset)
{
  _penalty_threshold_v = pset.get<std::vector<double>>("PEPenaltyThreshold");
}

FlashMatch_t Chi2Match::Match(const QCluster_t &pt_v, const Flash_t &flash)
{

  //
  // Prepare TPC
  //
  _raw_trk.resize(pt_v.size());

  for (size_t i = 0; i < pt_v.size(); ++i)
  {
    auto const &pt = pt_v[i];
    _raw_trk[i] = pt;
  }
  //
  // Prepare Flash
  //
  if (_hypothesis.pe_v.empty())
    _hypothesis.pe_v.resize(NOpDets(), 0.);
  if (_hypothesis.pe_v.size() != NOpDets())
  {
    throw OpT0FinderException("Hypothesis vector length != PMT count");
  }

  for (auto &v : _hypothesis.pe_v)
    v = 0;

  FillEstimate(_raw_trk, _hypothesis);

  FlashMatch_t res;

  res.tpc_point.x = res.tpc_point.y = res.tpc_point.z = 0;

  double weight = 0;

  for (size_t pmt_index = 0; pmt_index < NOpDets(); ++pmt_index)
  {

    res.tpc_point.y += OpDetY(pmt_index) * _hypothesis.pe_v[pmt_index];
    res.tpc_point.z += OpDetZ(pmt_index) * _hypothesis.pe_v[pmt_index];

    weight += _hypothesis.pe_v[pmt_index];
  }

  res.tpc_point.y /= weight;
  res.tpc_point.z /= weight;

  res.tpc_point.x = 0;
  res.tpc_point_err.x = 0;
  res.hypothesis = _hypothesis.pe_v;

  // Compute score
  res.score = Chi2(_hypothesis, flash);


  return res;
}

double Chi2Match::Chi2(const Flash_t &hypothesis,
                       const Flash_t &measurement)
{

  std::cout << "Hypothesis PE: " << std::accumulate(hypothesis.pe_v.begin(), hypothesis.pe_v.end(), 0.f) << std::endl;
  std::cout << "Flash PE: " << std::accumulate(measurement.pe_v.begin(), measurement.pe_v.end(), 0.f) << std::endl;
  double nvalid_pmt = 0;

  double PEtot_Hyp = 0;
  for (auto const &pe : hypothesis.pe_v)
    PEtot_Hyp += pe;
  double PEtot_Obs = 0;
  for (auto const &pe : measurement.pe_v)
    PEtot_Obs += pe;

  _current_chi2 = 0.;

  if (measurement.pe_v.size() != hypothesis.pe_v.size())
    throw OpT0FinderException("Cannot compute Chi2 for unmatched length!");

  double O, H, Error;
  for (size_t pmt_index = 0; pmt_index < hypothesis.pe_v.size(); ++pmt_index)
  {

    O = measurement.pe_v[pmt_index]; // observation
    H = hypothesis.pe_v[pmt_index];  // hypothesis

    std::cout << "\t\t O : " << O << "\t H : " << H << std::endl;

    if (H < 0)
      throw OpT0FinderException("Cannot have hypothesis value < 0!");

    // David C -> ignore PMTs where observation and hypothesis are below threshold
    if ((O < _penalty_threshold_v[pmt_index]) && (H < _penalty_threshold_v[pmt_index] * 1.5))
      continue;

    nvalid_pmt += 1;

    Error = O + H;
    if (Error < 1.0)
      Error = 1.0;

    auto chilocal = std::pow((O - H), 2) / (Error);

    //std::cout << "\t\t Chi2 += " << chilocal << std::endl;

    _current_chi2 += chilocal;

  } // for all PMTs

  _current_chi2 /= nvalid_pmt;

  std::cout << "Chi2 final : " << _current_chi2 << " with " << nvalid_pmt << std::endl;

  return _current_chi2;
}

} // namespace flashana
#endif
