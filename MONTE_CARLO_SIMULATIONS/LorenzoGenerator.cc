#include "Lorenzo/Test/interface/LorenzoGenerator.h" 

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"

#include "HepPDT/ParticleDataTable.hh"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

#include <CLHEP/Random/RandExponential.h>


using namespace edm;
using namespace std;

//----------------------------------------------------------------------------------------------------

LorenzoGenerator::LorenzoGenerator(const edm::ParameterSet& pset) : 
  parameterExample(pset.getParameter<double>("parameterExample"))
{
  printf("Hello world");
  printf(">> LorenzoGenerator::LorenzoGenerator > parameterExample = %.0f\n", parameterExample);

  h_p_T_X = new TH1D("h_p_T_X", "", 100, 0., 0.);
  h_p_z_X = new TH1D("h_p_z_X", "", 100, 0., 0.);
  h_p_tot_X = new TH1D("h_p_tot_X", "", 100, 0., 0.);

  h_p_T_Z = new TH1D("h_p_T_Z", "", 100, 0., 0.);
  h_p_z_Z = new TH1D("h_p_z_Z", "", 100, 0., 0.);
  h_p_tot_Z = new TH1D("h_p_tot_Z", "", 100, 0., 0.);

  produces<HepMCProduct>("unsmeared");
}

//----------------------------------------------------------------------------------------------------

void LorenzoGenerator::produce(edm::Event &e, const edm::EventSetup& es) 
{
  printf("\n>> LorenzoGenerator::produce > event %llu\n", e.id().event());

  // get conditions
  edm::Service<edm::RandomNumberGenerator> rng;
  CLHEP::HepRandomEngine* engine = &rng->getEngine(e.streamID());

  ESHandle<HepPDT::ParticleDataTable> pdgTable;
  es.getData(pdgTable);

  // prepare HepMC event
  HepMC::GenEvent *fEvt = new HepMC::GenEvent();
  fEvt->set_event_number(e.id().event());
   
  // generate vertex position
  HepMC::GenVertex *vtx = new HepMC::GenVertex(HepMC::FourVector(0., 0., 0., 0.));
  fEvt->add_vertex(vtx);

  // particle ids
  unsigned int particleId_Z = 23;
  unsigned int particleId_X = 999999;
  unsigned int particleId_p = 2212;

  // particle masses
  const double m_X = 1200.; // GeV
  const double m_Z = 91.;   // GeV

  //const HepPDT::ParticleData *pData = pdgTable->particle(HepPDT::ParticleID(particleId));
  //double mass_1 = pData->mass().value();
  //double mass_2 = pData->mass().value();

  // beam momentum
  const double p_beam = 6500.; // GeV

  // generate invariant mass of the X-Z system
  const double c = 0.04;
  const double c_mean = 1. / c;
  const double m_XZ = 1300. + CLHEP::RandExponential::shoot(engine, c_mean);

  // generate p_z of the 2-proton system in the LAB frame
  const double medlab = 180.0, sigmalab = 450.0;
  const double p_z_LAB_2p = CLHEP::RandGauss::shoot(engine, medlab, sigmalab);

  // generate spherical angles in the CMS frame of the X-Z system
  const double theta_c = CLHEP::RandFlat::shoot(engine) * M_PI;
  const double phi_c = CLHEP::RandFlat::shoot(engine) * 2. * M_PI;

  // determine xi's of the protons
  // proton 1: positive z momentum component
  const double xi2 = (p_z_LAB_2p + sqrt(p_z_LAB_2p*p_z_LAB_2p + m_XZ*m_XZ)) / (2. * p_beam);
  const double xi1 = m_XZ * m_XZ / (4. * p_beam * p_beam * xi2);

  printf("  m_XZ = %.1f\n", m_XZ);
  printf("  p_z_LAB_2p = %.1f\n", p_z_LAB_2p);
  printf("  xi1 = %.3f, xi2 = %.3f\n", xi1, xi2);
  printf("  p_beam * (xi2 - xi1) = %.1f\n", p_beam * (xi2 - xi1));

  // determine momenta of the X and Z particles in the CMS frame of the X-Z system
  const double p_c = sqrt( pow(m_XZ*m_XZ - m_X*m_X - m_Z*m_Z, 2.) / (4. * m_XZ * m_XZ) - m_X*m_X * m_Z*m_Z / (m_XZ*m_XZ) );

  const double E_X_CMS = sqrt(p_c*p_c + m_X*m_X);
  const double p_X_x_CMS = p_c * sin(theta_c) * cos(phi_c);
  const double p_X_y_CMS = p_c * sin(theta_c) * sin(phi_c);
  const double p_X_z_CMS = p_c * sin(theta_c);

  const double E_Z_CMS = sqrt(p_c*p_c + m_Z*m_Z);
  const double p_Z_x_CMS = - p_c * sin(theta_c) * cos(phi_c);
  const double p_Z_y_CMS = - p_c * sin(theta_c) * sin(phi_c);
  const double p_Z_z_CMS = - p_c * sin(theta_c);

  const double E_XZ_CMS = E_X_CMS + E_Z_CMS;

  // determine boost from X-Z CMS frame to the LAB frame
  const double E_XZ_LAB = p_beam * (xi1 + xi2);
  const double beta = sqrt(1. - E_XZ_CMS / E_XZ_LAB);
  const double gamma = 1. / sqrt(1. - beta*beta);

  // determine four-momenta of the outgoing particles in the LAB frame

  HepMC::FourVector momentum_X(
    p_X_x_CMS,
    p_X_y_CMS,
    gamma * (p_X_z_CMS - beta * E_X_CMS),
    gamma * (E_X_CMS - beta * p_X_z_CMS)
  );
  
  HepMC::FourVector momentum_Z(
    p_Z_x_CMS,
    p_Z_y_CMS,
    gamma * (p_Z_z_CMS - beta * E_Z_CMS),
    gamma * (E_Z_CMS - beta * p_Z_z_CMS)
  );

  HepMC::FourVector momentum_p1(0., 0., +p_beam * (1. - xi1), p_beam * (1. - xi1));

  HepMC::FourVector momentum_p2(0., 0., -p_beam * (1. - xi2), p_beam * (1. - xi2));

  // test printout
  //printf("  m_XZ = %.1f\n", m_XZ);
  //printf("  xi1 = %.3f, xi2 = %.3f\n", xi1, xi2);

    // fill in the HepMC record
    unsigned int barcode = 0;
    
    HepMC::GenParticle* particle_Z = new HepMC::GenParticle(momentum_Z, particleId_Z, 1);
    particle_Z->suggest_barcode(++barcode);
    vtx->add_particle_out(particle_Z);
    
    HepMC::GenParticle* particle_X = new HepMC::GenParticle(momentum_X, particleId_X, 1);
    particle_X->suggest_barcode(++barcode);
    vtx->add_particle_out(particle_X);
    
    HepMC::GenParticle* particle_p1 = new HepMC::GenParticle(momentum_p1, particleId_p, 1);
    particle_p1->suggest_barcode(++barcode);
    vtx->add_particle_out(particle_p1);
    
    HepMC::GenParticle* particle_p2 = new HepMC::GenParticle(momentum_p2, particleId_p, 1);
    particle_p2->suggest_barcode(++barcode);
    vtx->add_particle_out(particle_p2);

  // validation
  const double p_T_X = sqrt(momentum_X.x()*momentum_X.x() + momentum_X.y()*momentum_X.y());
  const double p_z_X = momentum_X.z();
  const double p_tot_X = momentum_X.rho();

  const double p_T_Z = sqrt(momentum_Z.x()*momentum_Z.x() + momentum_Z.y()*momentum_Z.y());
  const double p_z_Z = momentum_Z.z();
  const double p_tot_Z = momentum_Z.rho();

  h_p_T_X->Fill(p_T_X);
  h_p_z_X->Fill(p_z_X);
  h_p_tot_X->Fill(p_tot_X);

  h_p_T_Z->Fill(p_T_Z);
  h_p_z_Z->Fill(p_z_Z);
  h_p_tot_Z->Fill(p_tot_Z);

  // save output
  std::unique_ptr<HepMCProduct> output(new HepMCProduct()) ;
  output->addHepMCData(fEvt);
  e.put(std::move(output), "unsmeared");
}

//----------------------------------------------------------------------------------------------------

LorenzoGenerator::~LorenzoGenerator()
{
  TFile *f_out = TFile::Open("LorenzoGenerator_debug.root", "recreate");

  h_p_T_X->Write();
  h_p_z_X->Write();
  h_p_tot_X->Write();

  h_p_T_Z->Write();
  h_p_z_Z->Write();
  h_p_tot_Z->Write();

  delete f_out;
}
