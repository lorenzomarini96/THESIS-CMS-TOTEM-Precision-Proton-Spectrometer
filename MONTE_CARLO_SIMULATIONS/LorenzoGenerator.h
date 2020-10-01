/****************************************************************************
 *
 * This is a part of CTPPS offline software.
 * Authors:
 *   Jan Kašpar
 *
 ****************************************************************************/

#ifndef LorenzoGenerator_H
#define LorenzoGenerator_H

#include "FWCore/Framework/interface/one/EDProducer.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"

#include "TH1D.h"
#include "TFile.h"

namespace edm {
  
class LorenzoGenerator : public one::EDProducer<>
{

  public:
    LorenzoGenerator(const ParameterSet &);

    virtual ~LorenzoGenerator();

  private:
    virtual void produce(Event & e, const EventSetup& es) override;

    double parameterExample;

    TH1D *h_p_T_X;
    TH1D *h_p_z_X;
    TH1D *h_p_tot_X;
    TH1D *h_Theta_X;
    TH1D *h_Eta_X;

    TH1D *h_p_T_Z;
    TH1D *h_p_z_Z;
    TH1D *h_p_tot_Z;
    TH1D *h_Theta_Z;
    TH1D *h_Eta_Z;
};

} 

#endif
