//! @file Photolysis.h

#ifndef CT_PHOTOLYSIS_H
#define CT_PHOTOLYSIS_H

#include "cantera/base/ct_defs.h"
#include "cantera/base/Units.h"
#include "cantera/kinetics/ReactionData.h"
#include "ReactionRate.h"
#include "MultiRate.h"

namespace Cantera
{

int locate(double const *xx, double x, int n);
void interpn(double *val, double const *coor, double const *data,
             double const *axis, size_t const *len, int ndim, int nval);

class ThermoPhase;
class Kinetics;

//! Data container holding shared data specific to photolysis reactions
/**
 * The data container `PhotolysisData` holds photolysis cross-section data
 * @ingroup reactionGroup
 */
struct PhotolysisData : public ReactionData {
  bool check() const;

  bool update(const ThermoPhase& thermo, const Kinetics& kin) override;
  using ReactionData::update;

  //! \brief wavelength grid
  //!
  //! The wavelength grid is a vector of size nwave.
  //! Default units are nanometers.
  vector<double> wavelength;

  //! \brief actinic flux
  //!
  //! The actinic flux is a vector of size nwave.
  //! Default units are photons cm^-2 s^-1 nm^-1.
  vector<double> actinicFlux;
};

class PhotolysisBase : public ReactionRate {
 public:
  //! Default constructor
  PhotolysisBase() {}

  //! Constructor.
  /*!
   * @param temp Temperature grid
   * @param wavelength Wavelength grid
   * @param branches Branch strings of the photolysis reaction
   * @param xsection Cross-section data
   */
  PhotolysisBase(vector<double> const& temp, vector<double> const& wavelength,
                 vector<std::string> const& branches,
                 vector<double> const& xsection);

  //! Constructor based on AnyValue content
  explicit PhotolysisBase(AnyMap const& node, UnitStack const& rate_units={});

  void setParameters(AnyMap const& node, UnitStack const& rate_units) override;

  void setRateParameters(const AnyValue& rate, vector<string> const& branches);

  void loadCrossSectionVulcan(vector<string> files, string const& branch);

  void loadCrossSectionKinetics7(vector<string> files, string const& branch);

  void getParameters(AnyMap& node) const override;

  void getRateParameters(AnyMap& node) const;

  void check(string const& equation) override;

  void validate(const string& equation, const Kinetics& kin) override;

 protected:
  //! composition of branches
  vector<Composition> m_branch;

  //! number of temperature grid points
  size_t m_ntemp;

  //! number of wavelength grid points
  size_t m_nwave;

  //! temperature grid followed by wavelength grid
  vector<double> m_temp_wave_grid;

  //! \brief photolysis cross-section data
  //!
  //! The cross-section data is a three dimensional table of size (ntemp, nwave, nbranch).
  //! The first dimension is the number of temperature grid points, the second dimension
  //! is the number of wavelength grid points, and the third dimension is the number of
  //! branches of the photolysis reaction.
  //! Default units are nanometers, cm^2, cm^2, and cm^2, respectively.
  vector<double> m_crossSection;
};

//! Photolysis reaction rate type depends on temperature and the actinic flux
/*! 
 * A reaction rate coefficient of the following form.
 *
 * \f[
 *    k(T) = \int_{\lambda_1}^{\lambda_2} \sigma(\lambda) \phi(\lambda) d\lambda
 * \f]
 *
 * where \f$ \sigma(\lambda) \f$ is the cross-section and \f$ \phi(\lambda) \f$
 * is the actinic flux. \f$ \lambda_1 \f$ and \f$ \lambda_2 \f$ are the lower
 * and upper bounds of the wavelength grid.
 */
class PhotolysisRate : public PhotolysisBase {
 public:
  using PhotolysisBase::PhotolysisBase;  // inherit constructor

  unique_ptr<MultiRateBase> newMultiRate() const override {
    return make_unique<MultiRate<PhotolysisRate, PhotolysisData>>();
  }

  const string type() const override {
    return "Photolysis";
  }

  double evalFromStruct(PhotolysisData const& data) {
    if (m_crossSection.empty()) {
      return 0.;
    }

    double wmin = m_temp_wave_grid[m_ntemp];
    double wmax = m_temp_wave_grid.back();

    if (wmin > data.wavelength.front()) {
      throw CanteraError("PhotolysisRate::evalFromStruct",
                         "Wavelength out of range: {} nm < {} nm",
                         wmin, data.wavelength.front());
    }

    if (wmax < data.wavelength.back()) {
      throw CanteraError("PhotolysisRate::evalFromStruct",
                         "Wavelength out of range: {} nm > {} nm",
                         wmax, data.wavelength.back());
    }

    int iwmin = locate(data.wavelength.data(), wmin, data.wavelength.size());
    int iwmax = locate(data.wavelength.data(), wmax, data.wavelength.size());

    double* cross1 = new double [m_branch.size()];
    double* cross2 = new double [m_branch.size()];

    double coord[2] = {data.temperature, data.wavelength[iwmin]};
    size_t len[2] = {m_ntemp, m_nwave};

    interpn(cross1, coord, m_crossSection.data(), m_temp_wave_grid.data(),
        len, 2, m_branch.size());

    double total_rate = 0.0;
    for (auto& [name, stoich] : m_net_products)
      stoich = 0.0;

    for (int i = iwmin; i < iwmax; i++) {
      coord[1] = data.wavelength[i+1];
      interpn(cross2, coord, m_crossSection.data(), m_temp_wave_grid.data(),
          len, 2, m_branch.size());

      for (size_t n = 0; n < m_branch.size(); n++) {
        double rate = 0.5 * (data.wavelength[i+1] - data.wavelength[i])
          * (cross1[n] * data.actinicFlux[i] + cross2[n] * data.actinicFlux[i+1]);
        for (auto const& [name, stoich] : m_branch[n])
          m_net_products[name] += rate * stoich;
        total_rate += rate;
        cross1[n] = cross2[n];
      }
    }

    for (auto& [name, stoich] : m_net_products)
      stoich /= total_rate;

    delete [] cross1;
    delete [] cross2;

    return total_rate;
  }

 protected:
  //! net stoichiometric coefficients of products
  Composition m_net_products;
};

}

#endif  // CT_PHOTOLYSIS_H
