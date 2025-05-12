#include <iostream>
#include <sstream>
#include <string>
#include <memory>
#include <cstdint>
#include <vector>

#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"

#include "Vincia/Vincia.h"
#include "Dire/Dire.h"

using namespace Pythia8;

#include "GeneratorInterface/Pythia8Interface/interface/Py8InterfaceBase.h"

#include "GeneratorInterface/Pythia8Interface/plugins/ReweightUserHooks.h"

// PS matchning prototype
//
#include "GeneratorInterface/Pythia8Interface/plugins/JetMatchingHook.h"
#include "Pythia8Plugins/JetMatching.h"
#include "Pythia8Plugins/aMCatNLOHooks.h"

#include "GeneratorInterface/Pythia8Interface/interface/MultiUserHook.h"

// Emission Veto Hooks
//
#include "Pythia8Plugins/PowhegHooks.h"
#include "GeneratorInterface/Pythia8Interface/plugins/EmissionVetoHook1.h"

// Resonance scale hook
#include "GeneratorInterface/Pythia8Interface/plugins/PowhegResHook.h"
#include "GeneratorInterface/Pythia8Interface/plugins/PowhegHooksBB4L.h"

//decay filter hook
#include "GeneratorInterface/Pythia8Interface/interface/ResonanceDecayFilterHook.h"

//decay filter hook
#include "GeneratorInterface/Pythia8Interface/interface/PTFilterHook.h"

// EvtGen plugin
//
#include "Pythia8Plugins/EvtGen.h"

#include "FWCore/Concurrency/interface/SharedResourceNames.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

#include "GeneratorInterface/Core/interface/GeneratorFilter.h"
#include "GeneratorInterface/Core/interface/HadronizerFilter.h"

#include "GeneratorInterface/Pythia8Interface/plugins/LHAupLesHouches.h"

#include "HepPID/ParticleIDTranslations.hh"

#include "GeneratorInterface/ExternalDecays/interface/ExternalDecayDriver.h"

#include "TLorentzVector.h"
namespace CLHEP {
  class HepRandomEngine;
}

using namespace gen;


class Pythia8Hadronizer : public Py8InterfaceBase {

  public:

    Pythia8Hadronizer(const edm::ParameterSet &params);
   ~Pythia8Hadronizer() override;
 
    bool initializeForInternalPartons() override;
    bool initializeForExternalPartons();

    bool TwoMuMassFilter(Event &ev);	
    bool ThreeMuMassFilter(Event &ev);	
    bool generatePartonsAndHadronize() override;
    bool hadronize();




    virtual bool residualDecay();

    void finalizeEvent() override;

    void statistics() override;

    const char *classname() const override { return "Pythia8Hadronizer"; }
    
    GenLumiInfoHeader *getGenLumiInfoHeader() const override;
    
  private:

    std::auto_ptr<Vincia::VinciaPlugin> fvincia;
    std::auto_ptr<Pythia8::Dire> fDire;

    void doSetRandomEngine(CLHEP::HepRandomEngine* v) override { p8SetRandomEngine(v); }
    std::vector<std::string> const& doSharedResources() const override { return p8SharedResources; }

    /// Center-of-Mass energy
    double       comEnergy;

    /// Number of re-decays

    int nRepeat;
    std::string ReDecayConditions;
    std::vector<int> ParticlesIDtoRedecay;

    std::string LHEInputFileName;
    std::auto_ptr<LHAupLesHouches>  lhaUP;

    enum { PP, PPbar, ElectronPositron };
    int  fInitialState ; // pp, ppbar, or e-e+

    double fBeam1PZ;
    double fBeam2PZ;

    //helper class to allow multiple user hooks simultaneously
    std::auto_ptr<MultiUserHook> fMultiUserHook;
    
    // Reweight user hooks
    //
    std::auto_ptr<UserHooks> fReweightUserHook;
    std::auto_ptr<UserHooks> fReweightEmpUserHook;
    std::auto_ptr<UserHooks> fReweightRapUserHook;  
    std::auto_ptr<UserHooks> fReweightPtHatRapUserHook;
        
    // PS matching prototype
    //
    std::auto_ptr<JetMatchingHook> fJetMatchingHook;
    std::auto_ptr<Pythia8::JetMatchingMadgraph> fJetMatchingPy8InternalHook;
    std::auto_ptr<Pythia8::amcnlo_unitarised_interface> fMergingHook;
    
    // Emission Veto Hooks
    //
    std::auto_ptr<PowhegHooks> fEmissionVetoHook;
    std::auto_ptr<EmissionVetoHook1> fEmissionVetoHook1;
    
    // Resonance scale hook
    std::auto_ptr<PowhegResHook> fPowhegResHook;
    std::auto_ptr<PowhegHooksBB4L> fPowhegHooksBB4L;
    
    //resonance decay filter hook
    std::auto_ptr<ResonanceDecayFilterHook> fResonanceDecayFilterHook;
 
    //PT filter hook
    std::auto_ptr<PTFilterHook> fPTFilterHook;
   
    int  EV1_nFinal;
    bool EV1_vetoOn;
    int  EV1_maxVetoCount;
    int  EV1_pThardMode;
    int  EV1_pTempMode;
    int  EV1_emittedMode;
    int  EV1_pTdefMode;
    bool EV1_MPIvetoOn;   
    int  EV1_QEDvetoMode;
    int  EV1_nFinalMode;

    static const std::vector<std::string> p8SharedResources;
    
    vector<float> DJR;
    int nME;
    int nMEFiltered;

    int nISRveto;
    int nFSRveto;
    
};

const std::vector<std::string> Pythia8Hadronizer::p8SharedResources = { edm::SharedResourceNames::kPythia8 };

Pythia8Hadronizer::Pythia8Hadronizer(const edm::ParameterSet &params) :
  Py8InterfaceBase(params),
  comEnergy(params.getParameter<double>("comEnergy")),
  nRepeat(params.getParameter<int>("nRepeat")),
  ReDecayConditions(params.getParameter<std::string>("ReDecayConditions")),
  ParticlesIDtoRedecay(params.getParameter< std::vector<int> >("ParticlesIDtoRedecay")),
  LHEInputFileName(params.getUntrackedParameter<std::string>("LHEInputFileName","")),
  fInitialState(PP),
  nME(-1), nMEFiltered(-1), nISRveto(0), nFSRveto(0)
{

  // J.Y.: the following 3 parameters are hacked "for a reason"
  //
  if ( params.exists( "PPbarInitialState" ) )
  {
    if ( fInitialState == PP )
    {
      fInitialState = PPbar;
      edm::LogImportant("GeneratorInterface|Pythia8Interface")
      << "Pythia8 will be initialized for PROTON-ANTIPROTON INITIAL STATE. "
      << "This is a user-request change from the DEFAULT PROTON-PROTON initial state.";
    }
    else
    {   
      // probably need to throw on attempt to override ?
    }
  }   
  else if ( params.exists( "ElectronPositronInitialState" ) )
  {
    if ( fInitialState == PP )
    {
      fInitialState = ElectronPositron;
      edm::LogInfo("GeneratorInterface|Pythia8Interface")
      << "Pythia8 will be initialized for ELECTRON-POSITRON INITIAL STATE. "
      << "This is a user-request change from the DEFAULT PROTON-PROTON initial state.";
    }
    else
    {   
       // probably need to throw on attempt to override ?
    }
  }
  else if ( params.exists( "ElectronProtonInitialState" ) || params.exists( "PositronProtonInitialState" ) )
  {
    // throw on unknown initial state !
    throw edm::Exception(edm::errors::Configuration,"Pythia8Interface")
      <<" UNKNOWN INITIAL STATE. \n The allowed initial states are: PP, PPbar, ElectronPositron \n";
  }

  // Reweight user hook
  //
  if( params.exists( "reweightGen" ) )
  {
    edm::LogInfo("Pythia8Interface") << "Start setup for reweightGen";
    edm::ParameterSet rgParams =
       params.getParameter<edm::ParameterSet>("reweightGen");
    fReweightUserHook.reset(
       new PtHatReweightUserHook(rgParams.getParameter<double>("pTRef"),
                                 rgParams.getParameter<double>("power"))
       );
    edm::LogInfo("Pythia8Interface") << "End setup for reweightGen";
  }
  if( params.exists( "reweightGenEmp" ) )
  {
    edm::LogInfo("Pythia8Interface") << "Start setup for reweightGenEmp";
    edm::ParameterSet rgeParams =
       params.getParameter<edm::ParameterSet>("reweightGenEmp");

    std::string tuneName = "";
    if(rgeParams.exists("tune"))
        tuneName =  rgeParams.getParameter<std::string>("tune");
    fReweightEmpUserHook.reset(new PtHatEmpReweightUserHook(tuneName));
    edm::LogInfo("Pythia8Interface") << "End setup for reweightGenEmp";
  }
  if( params.exists( "reweightGenRap" ) )
  {
    edm::LogInfo("Pythia8Interface") << "Start setup for reweightGenRap";
    edm::ParameterSet rgrParams =
      params.getParameter<edm::ParameterSet>("reweightGenRap");
    fReweightRapUserHook.reset(
      new RapReweightUserHook(rgrParams.getParameter<std::string>("yLabSigmaFunc"),
                              rgrParams.getParameter<double>("yLabPower"),
                              rgrParams.getParameter<std::string>("yCMSigmaFunc"),
                              rgrParams.getParameter<double>("yCMPower"),
                              rgrParams.getParameter<double>("pTHatMin"),
                              rgrParams.getParameter<double>("pTHatMax"))
                         );
    edm::LogInfo("Pythia8Interface") << "End setup for reweightGenRap";
  }
  if( params.exists( "reweightGenPtHatRap" ) )
  {
    edm::LogInfo("Pythia8Interface") << "Start setup for reweightGenPtHatRap";
    edm::ParameterSet rgrParams =
      params.getParameter<edm::ParameterSet>("reweightGenPtHatRap");
    fReweightPtHatRapUserHook.reset(
      new PtHatRapReweightUserHook(rgrParams.getParameter<std::string>("yLabSigmaFunc"),
                                   rgrParams.getParameter<double>("yLabPower"),
                                   rgrParams.getParameter<std::string>("yCMSigmaFunc"),
                                   rgrParams.getParameter<double>("yCMPower"),
                                   rgrParams.getParameter<double>("pTHatMin"),
                                   rgrParams.getParameter<double>("pTHatMax"))
                              );
    edm::LogInfo("Pythia8Interface") << "End setup for reweightGenPtHatRap";
  }

  if( params.exists( "useUserHook" ) )
    throw edm::Exception(edm::errors::Configuration,"Pythia8Interface")
      <<" Obsolete parameter: useUserHook \n Please use the actual one instead \n";

  // PS matching prototype
  //
  if ( params.exists("jetMatching") )
  {
    edm::ParameterSet jmParams =
      params.getUntrackedParameter<edm::ParameterSet>("jetMatching");
      std::string scheme = jmParams.getParameter<std::string>("scheme");
      if ( scheme == "Madgraph" || scheme == "MadgraphFastJet" )
      {
         fJetMatchingHook.reset(new JetMatchingHook( jmParams, &fMasterGen->info ));
      }
  }

  // Pythia8Interface emission veto
  //
  if ( params.exists("emissionVeto1") )
  {
    EV1_nFinal = -1;
    if(params.exists("EV1_nFinal")) EV1_nFinal = params.getParameter<int>("EV1_nFinal");
    EV1_vetoOn = true;
    if(params.exists("EV1_vetoOn")) EV1_vetoOn = params.getParameter<bool>("EV1_vetoOn");
    EV1_maxVetoCount = 10;
    if(params.exists("EV1_maxVetoCount")) EV1_maxVetoCount = params.getParameter<int>("EV1_maxVetoCount");
    EV1_pThardMode = 1;
    if(params.exists("EV1_pThardMode")) EV1_pThardMode = params.getParameter<int>("EV1_pThardMode");
    EV1_pTempMode = 0;
    if(params.exists("EV1_pTempMode")) EV1_pTempMode = params.getParameter<int>("EV1_pTempMode");
    if(EV1_pTempMode > 2 || EV1_pTempMode < 0)
      throw edm::Exception(edm::errors::Configuration,"Pythia8Interface")
        <<" Wrong value for EV1_pTempMode code\n";
    EV1_emittedMode = 0;
    if(params.exists("EV1_emittedMode")) EV1_emittedMode = params.getParameter<int>("EV1_emittedMode");
    EV1_pTdefMode = 1;
    if(params.exists("EV1_pTdefMode")) EV1_pTdefMode = params.getParameter<int>("EV1_pTdefMode");
    EV1_MPIvetoOn = false;
    if(params.exists("EV1_MPIvetoOn")) EV1_MPIvetoOn = params.getParameter<bool>("EV1_MPIvetoOn");
    EV1_QEDvetoMode = 0;
    if(params.exists("EV1_QEDvetoMode")) EV1_QEDvetoMode = params.getParameter<int>("EV1_QEDvetoMode");
    EV1_nFinalMode = 0;
    if(params.exists("EV1_nFinalMode")) EV1_nFinalMode = params.getParameter<int>("EV1_nFinalMode");
    fEmissionVetoHook1.reset(new EmissionVetoHook1(EV1_nFinal, EV1_vetoOn,
                               EV1_maxVetoCount, EV1_pThardMode, EV1_pTempMode,
                               EV1_emittedMode, EV1_pTdefMode, 
			       EV1_MPIvetoOn, EV1_QEDvetoMode, EV1_nFinalMode, 0));
  }
  
  if( params.exists( "VinciaPlugin" ) ) {
    fMasterGen.reset(new Pythia);
    fvincia.reset(new Vincia::VinciaPlugin(fMasterGen.get()));
  }
  if( params.exists( "DirePlugin" ) ) {
    fMasterGen.reset(new Pythia);
    fDire.reset(new Pythia8::Dire());
    fDire->initSettings(*fMasterGen.get());
    fDire->initShowersAndWeights(*fMasterGen.get(), nullptr, nullptr);
  }

}


Pythia8Hadronizer::~Pythia8Hadronizer()
{
  
}

bool Pythia8Hadronizer::initializeForInternalPartons()
{
  
  bool status = false, status1 = false;
  
  if (lheFile_.empty()) {
    if ( fInitialState == PP ) // default
    {
      fMasterGen->settings.mode("Beams:idA", 2212);
      fMasterGen->settings.mode("Beams:idB", 2212);
    }
    else if ( fInitialState == PPbar )
    {
      fMasterGen->settings.mode("Beams:idA", 2212);
      fMasterGen->settings.mode("Beams:idB", -2212);
    }
    else if ( fInitialState == ElectronPositron )
    {
      fMasterGen->settings.mode("Beams:idA", 11);
      fMasterGen->settings.mode("Beams:idB", -11);
    }    
    else 
    {
      // throw on unknown initial state !
      throw edm::Exception(edm::errors::Configuration,"Pythia8Interface")
        <<" UNKNOWN INITIAL STATE. \n The allowed initial states are: PP, PPbar, ElectronPositron \n";
    }
    fMasterGen->settings.parm("Beams:eCM", comEnergy);
  }
  else {
    fMasterGen->settings.mode("Beams:frameType", 4);
    fMasterGen->settings.word("Beams:LHEF", lheFile_);
  }
  
  fMultiUserHook.reset(new MultiUserHook);
  
  if(fReweightUserHook.get()) fMultiUserHook->addHook(fReweightUserHook.get());
  if(fReweightEmpUserHook.get()) fMultiUserHook->addHook(fReweightEmpUserHook.get());
  if(fReweightRapUserHook.get()) fMultiUserHook->addHook(fReweightRapUserHook.get());
  if(fReweightPtHatRapUserHook.get()) fMultiUserHook->addHook(fReweightPtHatRapUserHook.get());
  if(fJetMatchingHook.get()) fMultiUserHook->addHook(fJetMatchingHook.get());
  if(fEmissionVetoHook1.get()) { 
    edm::LogInfo("Pythia8Interface") << "Turning on Emission Veto Hook 1 from CMSSW Pythia8Interface";
    fMultiUserHook->addHook(fEmissionVetoHook1.get());
  }
  
  if (fMasterGen->settings.mode("POWHEG:veto") > 0 || fMasterGen->settings.mode("POWHEG:MPIveto") > 0) {

    if(fJetMatchingHook.get() || fEmissionVetoHook1.get())
      throw edm::Exception(edm::errors::Configuration,"Pythia8Interface")
      <<" Attempt to turn on PowhegHooks by pythia8 settings but there are incompatible hooks on \n Incompatible are : jetMatching, emissionVeto1 \n";

    fEmissionVetoHook.reset(new PowhegHooks());

    edm::LogInfo("Pythia8Interface") << "Turning on Emission Veto Hook from pythia8 code";
    fMultiUserHook->addHook(fEmissionVetoHook.get());
  }
  
  bool PowhegRes = fMasterGen->settings.flag("POWHEGres:calcScales");
  if (PowhegRes) {
    edm::LogInfo("Pythia8Interface") << "Turning on resonance scale setting from CMSSW Pythia8Interface";
    fPowhegResHook.reset(new PowhegResHook());
    fMultiUserHook->addHook(fPowhegResHook.get());
  }
  
  bool PowhegBB4L = fMasterGen->settings.flag("POWHEG:bb4l");
  if (PowhegBB4L) {
    edm::LogInfo("Pythia8Interface") << "Turning on BB4l hook from CMSSW Pythia8Interface";
    fPowhegHooksBB4L.reset(new PowhegHooksBB4L());
    fMultiUserHook->addHook(fPowhegHooksBB4L.get());
  }
  
  //adapted from main89.cc in pythia8 examples
  bool internalMatching = fMasterGen->settings.flag("JetMatching:merge");
  bool internalMerging = !(fMasterGen->settings.word("Merging:Process")=="void");
  
  if (internalMatching && internalMerging) {
    throw edm::Exception(edm::errors::Configuration,"Pythia8Interface")
      <<" Only one jet matching/merging scheme can be used at a time. \n";
  }
  
  if (internalMatching) {
    fJetMatchingPy8InternalHook.reset(new Pythia8::JetMatchingMadgraph);
    fMultiUserHook->addHook(fJetMatchingPy8InternalHook.get());
  }
  
  if (internalMerging) {
    int scheme = ( fMasterGen->settings.flag("Merging:doUMEPSTree")
                || fMasterGen->settings.flag("Merging:doUMEPSSubt")) ?
                1 :
                 ( ( fMasterGen->settings.flag("Merging:doUNLOPSTree")
                || fMasterGen->settings.flag("Merging:doUNLOPSSubt")
                || fMasterGen->settings.flag("Merging:doUNLOPSLoop")
                || fMasterGen->settings.flag("Merging:doUNLOPSSubtNLO")) ?
                2 :
                0 );
    fMergingHook.reset(new Pythia8::amcnlo_unitarised_interface(scheme));
    fMultiUserHook->addHook(fMergingHook.get());
  }
  
  bool resonanceDecayFilter = fMasterGen->settings.flag("ResonanceDecayFilter:filter");
  if (resonanceDecayFilter) {
    fResonanceDecayFilterHook.reset(new ResonanceDecayFilterHook);
    fMultiUserHook->addHook(fResonanceDecayFilterHook.get());
  }
 
  bool PTFilter = fMasterGen->settings.flag("PTFilter:filter");
  if (PTFilter) {
    fPTFilterHook.reset(new PTFilterHook);
    fMultiUserHook->addHook(fPTFilterHook.get());
  }
 
  if (fMultiUserHook->nHooks()>0) {
    fMasterGen->setUserHooksPtr(fMultiUserHook.get());
  }

  edm::LogInfo("Pythia8Interface") << "Initializing MasterGen";
  if( fvincia.get() ) {
    fvincia->init(); status = true;
  }
  else if( fDire.get() ) {
    //fDire->initTune(*fMasterGen.get());
    fDire->weightsPtr->setup();
    fMasterGen->init();
    fDire->setup(*fMasterGen.get());
    status = true;
  }
  else {
    status = fMasterGen->init();
  }
  
  //clean up temp file
  if (!slhafile_.empty()) {
    std::remove(slhafile_.c_str());
  }  

  if ( pythiaPylistVerbosity > 10 )
  {
    if ( pythiaPylistVerbosity == 11 || pythiaPylistVerbosity == 13 )
           fMasterGen->settings.listAll();
    if ( pythiaPylistVerbosity == 12 || pythiaPylistVerbosity == 13 )
           fMasterGen->particleData.listAll();
  }

  // init decayer
  fDecayer->settings.flag("ProcessLevel:all", false ); // trick
  fDecayer->settings.flag("ProcessLevel:resonanceDecays", true );
  edm::LogInfo("Pythia8Interface") << "Initializing Decayer";
  status1 = fDecayer->init();

  if (useEvtGen) {
    edm::LogInfo("Pythia8Hadronizer") << "Creating and initializing pythia8 EvtGen plugin";
    evtgenDecays.reset(new EvtGenDecays(fMasterGen.get(), evtgenDecFile, evtgenPdlFile));
    for (unsigned int i=0; i<evtgenUserFiles.size(); i++) evtgenDecays->readDecayFile(evtgenUserFiles.at(i));
  }

  return (status&&status1);
}


bool Pythia8Hadronizer::initializeForExternalPartons()
{

  edm::LogInfo("Pythia8Interface") << "Initializing for external partons";

  bool status = false, status1 = false;
  
  fMultiUserHook.reset(new MultiUserHook);
  
  if(fReweightUserHook.get()) fMultiUserHook->addHook(fReweightUserHook.get());
  if(fReweightEmpUserHook.get()) fMultiUserHook->addHook(fReweightEmpUserHook.get());
  if(fReweightRapUserHook.get()) fMultiUserHook->addHook(fReweightRapUserHook.get());
  if(fReweightPtHatRapUserHook.get()) fMultiUserHook->addHook(fReweightPtHatRapUserHook.get());
  if(fJetMatchingHook.get()) fMultiUserHook->addHook(fJetMatchingHook.get());
  if(fEmissionVetoHook1.get()) { 
    edm::LogInfo("Pythia8Interface") << "Turning on Emission Veto Hook 1 from CMSSW Pythia8Interface";
    fMultiUserHook->addHook(fEmissionVetoHook1.get());
  }
  
  if (fMasterGen->settings.mode("POWHEG:veto") > 0 || fMasterGen->settings.mode("POWHEG:MPIveto") > 0) {

    if(fJetMatchingHook.get() || fEmissionVetoHook1.get())
      throw edm::Exception(edm::errors::Configuration,"Pythia8Interface")
      <<" Attempt to turn on PowhegHooks by pythia8 settings but there are incompatible hooks on \n Incompatible are : jetMatching, emissionVeto1 \n";

    fEmissionVetoHook.reset(new PowhegHooks());

    edm::LogInfo("Pythia8Interface") << "Turning on Emission Veto Hook from pythia8 code";
    fMultiUserHook->addHook(fEmissionVetoHook.get());
  }
  
  bool PowhegRes = fMasterGen->settings.flag("POWHEGres:calcScales");
  if (PowhegRes) {
    edm::LogInfo("Pythia8Interface") << "Turning on resonance scale setting from CMSSW Pythia8Interface";
    fPowhegResHook.reset(new PowhegResHook());
    fMultiUserHook->addHook(fPowhegResHook.get());
  }

  bool PowhegBB4L = fMasterGen->settings.flag("POWHEG:bb4l");
  if (PowhegBB4L) {
    edm::LogInfo("Pythia8Interface") << "Turning on BB4l hook from CMSSW Pythia8Interface";
    fPowhegHooksBB4L.reset(new PowhegHooksBB4L());
    fMultiUserHook->addHook(fPowhegHooksBB4L.get());
  }
  
  //adapted from main89.cc in pythia8 examples
  bool internalMatching = fMasterGen->settings.flag("JetMatching:merge");
  bool internalMerging = !(fMasterGen->settings.word("Merging:Process")=="void");
  
  if (internalMatching && internalMerging) {
    throw edm::Exception(edm::errors::Configuration,"Pythia8Interface")
      <<" Only one jet matching/merging scheme can be used at a time. \n";
  }
  
  if (internalMatching) {
    fJetMatchingPy8InternalHook.reset(new Pythia8::JetMatchingMadgraph);
    fMultiUserHook->addHook(fJetMatchingPy8InternalHook.get());
  }
  
  if (internalMerging) {
    int scheme = ( fMasterGen->settings.flag("Merging:doUMEPSTree")
                || fMasterGen->settings.flag("Merging:doUMEPSSubt")) ?
                1 :
                 ( ( fMasterGen->settings.flag("Merging:doUNLOPSTree")
                || fMasterGen->settings.flag("Merging:doUNLOPSSubt")
                || fMasterGen->settings.flag("Merging:doUNLOPSLoop")
                || fMasterGen->settings.flag("Merging:doUNLOPSSubtNLO")) ?
                2 :
                0 );
    fMergingHook.reset(new Pythia8::amcnlo_unitarised_interface(scheme));
    fMultiUserHook->addHook(fMergingHook.get());
  }
  
  bool resonanceDecayFilter = fMasterGen->settings.flag("ResonanceDecayFilter:filter");
  if (resonanceDecayFilter) {
    fResonanceDecayFilterHook.reset(new ResonanceDecayFilterHook);
    fMultiUserHook->addHook(fResonanceDecayFilterHook.get());
  }
 
  bool PTFilter = fMasterGen->settings.flag("PTFilter:filter");
  if (PTFilter) {
    fPTFilterHook.reset(new PTFilterHook);
    fMultiUserHook->addHook(fPTFilterHook.get());
  }
 
  if (fMultiUserHook->nHooks()>0) {
    fMasterGen->setUserHooksPtr(fMultiUserHook.get());
  }  
  
  if(!LHEInputFileName.empty()) {

    edm::LogInfo("Pythia8Interface") << "Initialize direct pythia8 reading from LHE file "
                                     << LHEInputFileName;
    edm::LogInfo("Pythia8Interface") << "Some LHE information can be not stored";
    fMasterGen->settings.mode("Beams:frameType", 4);
    fMasterGen->settings.word("Beams:LHEF", LHEInputFileName);
    status = fMasterGen->init();

  } else {

    lhaUP.reset(new LHAupLesHouches());
    lhaUP->setScalesFromLHEF(fMasterGen->settings.flag("Beams:setProductionScalesFromLHEF"));
    lhaUP->loadRunInfo(lheRunInfo());
    
    if ( fJetMatchingHook.get() )
    {
       fJetMatchingHook->init ( lheRunInfo() );
    }
    
    fMasterGen->settings.mode("Beams:frameType", 5);
    fMasterGen->setLHAupPtr(lhaUP.get());
    edm::LogInfo("Pythia8Interface") << "Initializing MasterGen";
    status = fMasterGen->init();
  }
  
  //clean up temp file
  if (!slhafile_.empty()) {
    std::remove(slhafile_.c_str());
  }  
  
  if ( pythiaPylistVerbosity > 10 )
  {
    if ( pythiaPylistVerbosity == 11 || pythiaPylistVerbosity == 13 )
           fMasterGen->settings.listAll();
    if ( pythiaPylistVerbosity == 12 || pythiaPylistVerbosity == 13 )
           fMasterGen->particleData.listAll();
  }

  // init decayer
  fDecayer->settings.flag("ProcessLevel:all", false ); // trick
  fDecayer->settings.flag("ProcessLevel:resonanceDecays", true );
  edm::LogInfo("Pythia8Interface") << "Initializing Decayer";
  status1 = fDecayer->init();

  if (useEvtGen) {
    edm::LogInfo("Pythia8Hadronizer") << "Creating and initializing pythia8 EvtGen plugin";
    evtgenDecays.reset(new EvtGenDecays(fMasterGen.get(), evtgenDecFile, evtgenPdlFile));
    for (unsigned int i=0; i<evtgenUserFiles.size(); i++) evtgenDecays->readDecayFile(evtgenUserFiles.at(i));
  }

  return (status&&status1);
}


void Pythia8Hadronizer::statistics()
{
  fMasterGen->stat();

  if(fEmissionVetoHook.get()) {
    edm::LogPrint("Pythia8Interface") << "\n"
      << "Number of ISR vetoed = " << nISRveto;
    edm::LogPrint("Pythia8Interface")
      << "Number of FSR vetoed = " << nFSRveto;
  }

  double xsec = fMasterGen->info.sigmaGen(); // cross section in mb
  xsec *= 1.0e9; // translate to pb (CMS/Gen "convention" as of May 2009)
  double err  = fMasterGen->info.sigmaErr(); // cross section err in mb
  err  *= 1.0e9; // translate to pb (CMS/Gen "convention" as of May 2009)
  runInfo().setInternalXSec(GenRunInfoProduct::XSec(xsec,err));
}

bool Pythia8Hadronizer::TwoMuMassFilter( Event &ev){
  vector<int> negMuons;
  vector<int> posMuons;

  for (int i = 0; i < ev.size(); ++i)
    {
      if(ev.at(i).id() == 13)
	{
	  if(ev.at(i).status() > 0)
	    {
	      if(ev.at(i).pT() > 1.0)
		{
		  if(abs(ev.at(i).y())< 5.0)
		    {
		      negMuons.push_back(i);
		    }
		}
	    }
	}
      
      
      
      if(ev.at(i).id() == -13)
	{
	  if(ev.at(i).status() > 0)
	    {
	      if(ev.at(i).pT() > 1.0)
		{
		  if(abs(ev.at(i).y())< 5.0)
		    {
		      posMuons.push_back(i);
		    }
		}
	    }
	}
      

    }


  for(unsigned int iN=0; iN<negMuons.size(); iN++)
    for(unsigned  int iP=0; iP<posMuons.size(); iP++)
      if( (ev.at(negMuons.at(iN)).p() + ev.at(posMuons.at(iP)).p() ).mCalc() > 0.204 &&
	  (ev.at(negMuons.at(iN)).p() + ev.at(posMuons.at(iP)).p() ).mCalc() < 1.8) return true; //  search muons pair with 2mu < mass < mtau



  for(unsigned int iN1=0; iN1<negMuons.size(); iN1++)
    for(unsigned int iN2=iN1+1; iN2<negMuons.size(); iN2++)
      if( (ev.at(negMuons.at(iN1)).p() + ev.at(negMuons.at(iN1)).p() ).mCalc() > 0.204 &&
	  (ev.at(negMuons.at(iN2)).p() + ev.at(negMuons.at(iN2)).p() ).mCalc() < 1.8) return true; //  search muons pair with 2mu < mass < mtau


  for(unsigned int iP1=0; iP1<posMuons.size(); iP1++)
    for(unsigned int iP2=iP1+1; iP2<posMuons.size(); iP2++)
      if( (ev.at(posMuons.at(iP1)).p() + ev.at(posMuons.at(iP1)).p() ).mCalc() > 0.204 &&
	  (ev.at(posMuons.at(iP2)).p() + ev.at(posMuons.at(iP2)).p() ).mCalc() < 1.8) return true; //  search muons pair with 2mu < mass < mtau




  return false;
  
}



bool Pythia8Hadronizer::ThreeMuMassFilter( Event &ev){

  float TripleMass(0.);
  
  /*
  vector<int> negMuons;
  vector<int> posMuons;

  for (int i = 0; i < ev.size(); ++i) {

    if(ev.at(i).id() == 13)
      {
	if(ev.at(i).status() > 0)
	  {
	    if(ev.at(i).pT() > 1.0)
	      {
		if(abs(ev.at(i).y())< 2.9)
		  {
		    negMuons.push_back(i);
		  }
	      }
	  }
      }
    

    if(ev.at(i).id() == -13)
      {
	if(ev.at(i).status() > 0)
	  {
	    if(ev.at(i).pT() > 1.0)
	      {
		if(abs(ev.at(i).y())< 2.9)
		  {
		    posMuons.push_back(i);
		  }
	      }
	  }
      }
  }


  if(posMuons.size()  == 0 || negMuons.size()==0) return false;


  
  //search for Ptriplet
  for(unsigned int iN = 0; iN < negMuons.size(); iN++){
    for(unsigned int iP1 = 0; iP1 < posMuons.size()-1; iP1++){
      for(unsigned int iP2 = iP1 + 1; iP2 < posMuons.size(); iP2++){

	if((ev.at(negMuons.at(iN)).p() + ev.at(posMuons.at(iP1)).p() + ev.at(posMuons.at(iP2)).p()).mCalc() > 1.39 &&
	   (ev.at(negMuons.at(iN)).p() + ev.at(posMuons.at(iP1)).p() + ev.at(posMuons.at(iP2)).p()).mCalc() < 2.11) return true;
      }
    }
  }


  //search for Mtriplet
  for(unsigned int iP = 0; iP < posMuons.size(); iP++){
    for(unsigned int iN1 = 0; iN1 < negMuons.size()-1; iN1++){
      for(unsigned int iN2 = iN1 + 1; iN2 < negMuons.size(); iN2++){
	if((ev.at(posMuons.at(iP)).p() + ev.at(negMuons.at(iN1)).p() + ev.at(negMuons.at(iN2)).p()).mCalc() > 1.39 &&
	   (ev.at(posMuons.at(iP)).p() + ev.at(negMuons.at(iN1)).p() + ev.at(negMuons.at(iN2)).p()).mCalc() < 2.11) return true;
      }
    }
  }
  */
  
  vector<int> negMuons;
  vector<int> posMuons;
  
  vector<int> negPionsOrKaons;
  vector<int> posPionsOrKaons;

  for (int i = 0; i < ev.size(); ++i) {

    if(ev.at(i).id() == 13)
      {
	if(ev.at(i).status() > 0)
	  {
	    if(ev.at(i).pT() > 1.0)
	      {
		if(abs(ev.at(i).eta())< 2.9)
		  {
		    negMuons.push_back(i);
		  }
	      }
	  }
      }
    

    if(ev.at(i).id() == -13)
      {
	if(ev.at(i).status() > 0)
	  {
	    if(ev.at(i).pT() > 1.0)
	      {
		if(abs(ev.at(i).eta())< 2.9)
		  {
		    posMuons.push_back(i);
		  }
	      }
	  }
      }
  }
  
  for (int i = 0; i < ev.size(); ++i) {

    if(ev.at(i).id() == 211 || ev.at(i).id() == 321)
      {
	if(ev.at(i).status() > 0)
	  {
	    if(ev.at(i).pT() > 1.0)
	      {
		if(abs(ev.at(i).eta())< 2.9)
		  {
		    negPionsOrKaons.push_back(i);
		  }
	      }
	  }
      }
    

    if(ev.at(i).id() == -211 || ev.at(i).id() == -321)
      {
	if(ev.at(i).status() > 0)
	  {
	    if(ev.at(i).pT() > 1.0)
	      {
		if(abs(ev.at(i).eta())< 2.9)
		  {
		    posPionsOrKaons.push_back(i);
		  }
	      }
	  }
      }
  }
  
  vector<int> OppositeSide;//Fake jets, taus, muons, electrons
  for (int i = 0; i < ev.size(); ++i) {

    int partID = abs(ev.at(i).id());
    if(partID == 15 || partID == 13 || partID == 11 || partID == 1 || partID == 2 || partID == 3 || partID == 4 || partID == 5 || partID == 21)
      {
	if(ev.at(i).status() > 0)
	  {
	    if(ev.at(i).pT() > 1.0)
	      {
		if(abs(ev.at(i).eta())< 3.1)
		  {
		    OppositeSide.push_back(i);
		  }
	      }
	  }
      }
  }
  
  
  
  bool WhetherThreeMuonPosPass(false);
  bool WhetherThreeMuonNegPass(false);
  bool WhetherThreeMuonPosPassWithOS(false);
  bool WhetherThreeMuonNegPassWithOS(false);
  bool WhetherThreeMuonOSHadron(false);//can include muons
  bool WhetherThreeMuonOSeOrHadron(false);//can include muons
  bool WhetherThreeMuonOSmuOnly(false);
  bool WhetherTwoMuonandPionPosPass(false);
  bool WhetherTwoMuonandPionNegPass(false);
  


  //With just muons
  
  if(posMuons.size()  == 0 || negMuons.size()==0) return false;
  
  //search for Ptriplet
  for(unsigned int iN = 0; iN < negMuons.size(); iN++){
    for(unsigned int iP1 = 0; iP1 < posMuons.size(); iP1++){
      for(unsigned int iP2 = 0; iP2 < iP1; iP2++){
	auto part1 = ev.at(negMuons.at(iN)).p();
        auto part2 = ev.at(posMuons.at(iP1)).p();
        auto part3 = ev.at(posMuons.at(iP2)).p();
        auto massTriplet = part1 + part2 + part3;
        if((massTriplet).mCalc() > 1.39 &&
	   (massTriplet).mCalc() < 2.11) {
                   
                   if(massTriplet.pT()>14.5){
                           //std::cout<<" Muon Positive " <<std::endl;
                           WhetherThreeMuonPosPass = true;
                           //std::cout<<" Particle 1, x: "<< part1.px() <<" , y: "<< part1.py() <<" , z: "<< part1.pz() <<" , pT: "<< part1.pT() <<std::endl;
                           //std::cout<<" Particle 2, x: "<< part2.px() <<" , y: "<< part2.py() <<" , z: "<< part2.pz() <<" , pT: "<< part2.pT() <<std::endl;
                           //std::cout<<" Particle 3, x: "<< part3.px() <<" , y: "<< part3.py() <<" , z: "<< part3.pz() <<" , pT: "<< part3.pT() <<std::endl;
                           
                           //double dR12 = std::sqrt(std::pow(part1.eta() - part2.eta(), 2) + std::pow(std::fabs(part1.phi() - part2.phi()) > M_PI ? 2 * M_PI - std::fabs(part1.phi() - part2.phi()) : std::fabs(part1.phi() - part2.phi()), 2));
                           //double dR13 = std::sqrt(std::pow(part1.eta() - part3.eta(), 2) + std::pow(std::fabs(part1.phi() - part3.phi()) > M_PI ? 2 * M_PI - std::fabs(part1.phi() - part3.phi()) : std::fabs(part1.phi() - part3.phi()), 2));
                           //double dR32 = std::sqrt(std::pow(part3.eta() - part2.eta(), 2) + std::pow(std::fabs(part3.phi() - part2.phi()) > M_PI ? 2 * M_PI - std::fabs(part3.phi() - part2.phi()) : std::fabs(part3.phi() - part2.phi()), 2));
                           //std::cout<<" dR12: "<< dR12 <<" , dR13: "<< dR13 <<" , dR32: "<< dR32 << " , triplet pT: "<< massTriplet.pT() << " , triplet eta: " << massTriplet.eta() <<std::endl;
                           
                           
                           if(OppositeSide.size()>0){
                           for(unsigned int iO = 0; iO < OppositeSide.size(); iO++){
                                   if(OppositeSide.at(iO)!=posMuons.at(iP1)&&OppositeSide.at(iO)!=negMuons.at(iN)&&OppositeSide.at(iO)!=posMuons.at(iP2)){
                                           auto oppVect = ev.at(OppositeSide.at(iO)).p();
                                           double dRtoOpp = std::sqrt(std::pow(massTriplet.eta() - oppVect.eta(), 2) + std::pow(std::fabs(massTriplet.phi() - oppVect.phi()) > M_PI ? 2 * M_PI - std::fabs(massTriplet.phi() - oppVect.phi()) : std::fabs(massTriplet.phi() - oppVect.phi()), 2));
                                           
                                           if(dRtoOpp>0.5){
                                                   //std::cout<<" dR to opposite side: "<< dRtoOpp <<" , inv mass: "<< (massTriplet+oppVect).mCalc() <<" , particleID: "<< ev.at(OppositeSide.at(iO)).id() <<std::endl;
                                                   WhetherThreeMuonPosPassWithOS = true;
                                                   if(abs(ev.at(OppositeSide.at(iO)).id())==15||abs(ev.at(OppositeSide.at(iO)).id())==11||abs(ev.at(OppositeSide.at(iO)).id())==1||abs(ev.at(OppositeSide.at(iO)).id())==2||abs(ev.at(OppositeSide.at(iO)).id())==3||abs(ev.at(OppositeSide.at(iO)).id())==4||abs(ev.at(OppositeSide.at(iO)).id())==5||abs(ev.at(OppositeSide.at(iO)).id())==21){
                                                           WhetherThreeMuonOSeOrHadron = true;
                                                           if(abs(ev.at(OppositeSide.at(iO)).id())==15||abs(ev.at(OppositeSide.at(iO)).id())==1||abs(ev.at(OppositeSide.at(iO)).id())==2||abs(ev.at(OppositeSide.at(iO)).id())==3||abs(ev.at(OppositeSide.at(iO)).id())==4||abs(ev.at(OppositeSide.at(iO)).id())==5||abs(ev.at(OppositeSide.at(iO)).id())==21){
                                                                   WhetherThreeMuonOSHadron = true;
                                                           }
                                                   }
                                                   else{
                                                           WhetherThreeMuonOSmuOnly = true;
                                                   }
                                                   
                                           }
                                   }
                           }
                           }
                           
                   }
           }
      }
    }
  }


  
  //search for Mtriplet
  for(unsigned int iP = 0; iP < posMuons.size(); iP++){
    for(unsigned int iN1 = 0; iN1 < negMuons.size(); iN1++){
      for(unsigned int iN2 = 0; iN2 < iN1; iN2++){
	auto part1 = ev.at(posMuons.at(iP)).p();
        auto part2 = ev.at(negMuons.at(iN1)).p();
        auto part3 = ev.at(negMuons.at(iN2)).p();
        auto massTriplet = part1 + part2 + part3;
        if((massTriplet).mCalc() > 1.39 &&
	   (massTriplet).mCalc() < 2.11) {
                   
                   if(massTriplet.pT()>14.5){
                           
                           //std::cout<<" Muon Negative " <<std::endl;
                           WhetherThreeMuonNegPass = true;
                           //std::cout<<" Particle 1, x: "<< part1.px() <<" , y: "<< part1.py() <<" , z: "<< part1.pz() <<" , pT: "<< part1.pT() <<std::endl;
                           //std::cout<<" Particle 2, x: "<< part2.px() <<" , y: "<< part2.py() <<" , z: "<< part2.pz() <<" , pT: "<< part2.pT() <<std::endl;
                           //std::cout<<" Particle 3, x: "<< part3.px() <<" , y: "<< part3.py() <<" , z: "<< part3.pz() <<" , pT: "<< part3.pT() <<std::endl;
                           
                           //double dR12 = std::sqrt(std::pow(part1.eta() - part2.eta(), 2) + std::pow(std::fabs(part1.phi() - part2.phi()) > M_PI ? 2 * M_PI - std::fabs(part1.phi() - part2.phi()) : std::fabs(part1.phi() - part2.phi()), 2));
                           //double dR13 = std::sqrt(std::pow(part1.eta() - part3.eta(), 2) + std::pow(std::fabs(part1.phi() - part3.phi()) > M_PI ? 2 * M_PI - std::fabs(part1.phi() - part3.phi()) : std::fabs(part1.phi() - part3.phi()), 2));
                           //double dR32 = std::sqrt(std::pow(part3.eta() - part2.eta(), 2) + std::pow(std::fabs(part3.phi() - part2.phi()) > M_PI ? 2 * M_PI - std::fabs(part3.phi() - part2.phi()) : std::fabs(part3.phi() - part2.phi()), 2));
                           //std::cout<<" dR12: "<< dR12 <<" , dR13: "<< dR13 <<" , dR32: "<< dR32 << " , triplet pT: "<< massTriplet.pT() << " , triplet eta: " << massTriplet.eta() <<std::endl;
                           
                           
                           if(OppositeSide.size()>0){
                           for(unsigned int iO = 0; iO < OppositeSide.size(); iO++){
                                   if(OppositeSide.at(iO)!=posMuons.at(iP)&&OppositeSide.at(iO)!=negMuons.at(iN1)&&OppositeSide.at(iO)!=negMuons.at(iN2)){
                                           auto oppVect = ev.at(OppositeSide.at(iO)).p();
                                           double dRtoOpp = std::sqrt(std::pow(massTriplet.eta() - oppVect.eta(), 2) + std::pow(std::fabs(massTriplet.phi() - oppVect.phi()) > M_PI ? 2 * M_PI - std::fabs(massTriplet.phi() - oppVect.phi()) : std::fabs(massTriplet.phi() - oppVect.phi()), 2));
                                           
                                           if(dRtoOpp>0.5){
                                                   //std::cout<<" dR to opposite side: "<< dRtoOpp <<" , inv mass: "<< (massTriplet+oppVect).mCalc() <<" , particleID: "<< ev.at(OppositeSide.at(iO)).id() <<std::endl;
                                                   WhetherThreeMuonNegPassWithOS = true;
                                                   if(abs(ev.at(OppositeSide.at(iO)).id())==15||abs(ev.at(OppositeSide.at(iO)).id())==11||abs(ev.at(OppositeSide.at(iO)).id())==1||abs(ev.at(OppositeSide.at(iO)).id())==2||abs(ev.at(OppositeSide.at(iO)).id())==3||abs(ev.at(OppositeSide.at(iO)).id())==4||abs(ev.at(OppositeSide.at(iO)).id())==5||abs(ev.at(OppositeSide.at(iO)).id())==21){
                                                           WhetherThreeMuonOSeOrHadron = true;
                                                           if(abs(ev.at(OppositeSide.at(iO)).id())==15||abs(ev.at(OppositeSide.at(iO)).id())==1||abs(ev.at(OppositeSide.at(iO)).id())==2||abs(ev.at(OppositeSide.at(iO)).id())==3||abs(ev.at(OppositeSide.at(iO)).id())==4||abs(ev.at(OppositeSide.at(iO)).id())==5||abs(ev.at(OppositeSide.at(iO)).id())==21){
                                                                   WhetherThreeMuonOSHadron = true;
                                                           }
                                                   }
                                                   else{
                                                           WhetherThreeMuonOSmuOnly = true;
                                                   }
                                           }
                                   }
                           }
                           }
                           
                           
                   }
           }
      }
    }
  }
  
  
  //With Pions or Kaons
  
  //search for pospion+2muons
  if(posPionsOrKaons.size()>0){
  for(unsigned int iN = 0; iN < negMuons.size(); iN++){
    for(unsigned int iP1 = 0; iP1 < posMuons.size(); iP1++){
      for(unsigned int iP2 = 0; iP2 < posPionsOrKaons.size(); iP2++){
	auto part1 = ev.at(negMuons.at(iN)).p();
        auto part2 = ev.at(posMuons.at(iP1)).p();
        auto part3 = ev.at(posPionsOrKaons.at(iP2)).p();
        auto massTriplet = part1 + part2 + part3;
        if((massTriplet).mCalc() > 1.39 &&
	   (massTriplet).mCalc() < 2.11) {
                   
                   if(massTriplet.pT()>14.5){
                           //std::cout<<" Pion Positive " <<std::endl;
                           WhetherTwoMuonandPionPosPass = true;
                           //std::cout<<" Particle 1, x: "<< part1.px() <<" , y: "<< part1.py() <<" , z: "<< part1.pz() <<std::endl;
                           //std::cout<<" Particle 2, x: "<< part2.px() <<" , y: "<< part2.py() <<" , z: "<< part2.pz() <<std::endl;
                           //std::cout<<" Particle 3, x: "<< part3.px() <<" , y: "<< part3.py() <<" , z: "<< part3.pz() <<std::endl;
                           
                           //double dR12 = std::sqrt(std::pow(part1.eta() - part2.eta(), 2) + std::pow(std::fabs(part1.phi() - part2.phi()) > M_PI ? 2 * M_PI - std::fabs(part1.phi() - part2.phi()) : std::fabs(part1.phi() - part2.phi()), 2));
                           //double dR13 = std::sqrt(std::pow(part1.eta() - part3.eta(), 2) + std::pow(std::fabs(part1.phi() - part3.phi()) > M_PI ? 2 * M_PI - std::fabs(part1.phi() - part3.phi()) : std::fabs(part1.phi() - part3.phi()), 2));
                           //double dR32 = std::sqrt(std::pow(part3.eta() - part2.eta(), 2) + std::pow(std::fabs(part3.phi() - part2.phi()) > M_PI ? 2 * M_PI - std::fabs(part3.phi() - part2.phi()) : std::fabs(part3.phi() - part2.phi()), 2));
                           //std::cout<<" dR12: "<< dR12 <<" , dR13: "<< dR13 <<" , dR32: "<< dR32 << " , triplet pT: "<< massTriplet.pT() << " , triplet eta: " << massTriplet.eta() <<std::endl;
                           
                           /*
                           if(OppositeSide.size()>0){
                           for(unsigned int iO = 0; iO < OppositeSide.size(); iO++){
                                   if(OppositeSide.at(iO)!=posMuons.at(iP1)&&OppositeSide.at(iO)!=negMuons.at(iN)){
                                           auto oppVect = ev.at(OppositeSide.at(iO)).p();
                                           double dRtoOpp = std::sqrt(std::pow(massTriplet.eta() - oppVect.eta(), 2) + std::pow(std::fabs(massTriplet.phi() - oppVect.phi()) > M_PI ? 2 * M_PI - std::fabs(massTriplet.phi() - oppVect.phi()) : std::fabs(massTriplet.phi() - oppVect.phi()), 2));
                                           std::cout<<" dR to opposite side: "<< dRtoOpp <<" , inv mass: "<< (massTriplet+oppVect).mCalc() <<" , particleID: "<< ev.at(OppositeSide.at(iO)).id() <<std::endl;
                                           if(dRtoOpp>0.5){
                                                   return true;
                                           }
                                   }
                           }
                           }
                           */
                           
                           
                   }
           }
      }
    }
  }
  }
  
  //search for negpion+2muons
  if(negPionsOrKaons.size()>0){
  for(unsigned int iP = 0; iP < posMuons.size(); iP++){
    for(unsigned int iN1 = 0; iN1 < negMuons.size(); iN1++){
      for(unsigned int iN2 = 0; iN2 < negPionsOrKaons.size(); iN2++){
	auto part1 = ev.at(posMuons.at(iP)).p();
        auto part2 = ev.at(negMuons.at(iN1)).p();
        auto part3 = ev.at(negPionsOrKaons.at(iN2)).p();
        auto massTriplet = part1 + part2 + part3;
        if((massTriplet).mCalc() > 1.39 &&
	   (massTriplet).mCalc() < 2.11) {
                   
                   if(massTriplet.pT()>14.5){
                           //std::cout<<" Pion Negative " <<std::endl;
                           WhetherTwoMuonandPionNegPass = true;
                           //std::cout<<" Particle 1, x: "<< part1.px() <<" , y: "<< part1.py() <<" , z: "<< part1.pz() <<std::endl;
                           //std::cout<<" Particle 2, x: "<< part2.px() <<" , y: "<< part2.py() <<" , z: "<< part2.pz() <<std::endl;
                           //std::cout<<" Particle 3, x: "<< part3.px() <<" , y: "<< part3.py() <<" , z: "<< part3.pz() <<std::endl;
                           
                           //double dR12 = std::sqrt(std::pow(part1.eta() - part2.eta(), 2) + std::pow(std::fabs(part1.phi() - part2.phi()) > M_PI ? 2 * M_PI - std::fabs(part1.phi() - part2.phi()) : std::fabs(part1.phi() - part2.phi()), 2));
                           //double dR13 = std::sqrt(std::pow(part1.eta() - part3.eta(), 2) + std::pow(std::fabs(part1.phi() - part3.phi()) > M_PI ? 2 * M_PI - std::fabs(part1.phi() - part3.phi()) : std::fabs(part1.phi() - part3.phi()), 2));
                           //double dR32 = std::sqrt(std::pow(part3.eta() - part2.eta(), 2) + std::pow(std::fabs(part3.phi() - part2.phi()) > M_PI ? 2 * M_PI - std::fabs(part3.phi() - part2.phi()) : std::fabs(part3.phi() - part2.phi()), 2));
                           //std::cout<<" dR12: "<< dR12 <<" , dR13: "<< dR13 <<" , dR32: "<< dR32 << " , triplet pT: "<< massTriplet.pT() << " , triplet eta: " << massTriplet.eta() <<std::endl;
                           
                           /*
                           if(OppositeSide.size()>0){
                           for(unsigned int iO = 0; iO < OppositeSide.size(); iO++){
                                   if(OppositeSide.at(iO)!=posMuons.at(iP)&&OppositeSide.at(iO)!=negMuons.at(iN1)){
                                           auto oppVect = ev.at(OppositeSide.at(iO)).p();
                                           double dRtoOpp = std::sqrt(std::pow(massTriplet.eta() - oppVect.eta(), 2) + std::pow(std::fabs(massTriplet.phi() - oppVect.phi()) > M_PI ? 2 * M_PI - std::fabs(massTriplet.phi() - oppVect.phi()) : std::fabs(massTriplet.phi() - oppVect.phi()), 2));
                                           std::cout<<" dR to opposite side: "<< dRtoOpp <<" , inv mass: "<< (massTriplet+oppVect).mCalc() <<" , particleID: "<< ev.at(OppositeSide.at(iO)).id() <<std::endl;
                                           if(dRtoOpp>0.5){
                                                   return true;
                                           }
                                   }
                           }
                           }
                           */
                           
                           
                   }
           }
      }
    }
  }
  }
  
  //if((WhetherThreeMuonOSeOrHadron||WhetherThreeMuonOSmuOnly)&&!(WhetherThreeMuonOSeOrHadron)){
  //if(WhetherThreeMuonOSeOrHadron&&!(WhetherThreeMuonOSHadron)){
  if(WhetherThreeMuonOSHadron){
  
  
  //if((WhetherTwoMuonandPionPosPass||WhetherTwoMuonandPionNegPass)&&!(WhetherThreeMuonPosPass||WhetherThreeMuonNegPass)){
  //if(WhetherThreeMuonPosPass||WhetherThreeMuonNegPass){
  //if((WhetherThreeMuonPosPass||WhetherThreeMuonNegPass)&&!(WhetherThreeMuonPosPassWithOS||WhetherThreeMuonNegPassWithOS)){
  //if(WhetherThreeMuonOSeOrHadron){
  //if((WhetherThreeMuonPosPassWithOS||WhetherThreeMuonNegPassWithOS)&&!WhetherThreeMuonOSeOrHadron){
  
          //std::cout<<" Something passed. " <<std::endl;
          //std::cout<<" WhetherThreeMuonPosPass: "<< WhetherThreeMuonPosPass <<" , WhetherThreeMuonNegPass: "<< WhetherThreeMuonNegPass <<" , WhetherTwoMuonandPionPosPass: "<< WhetherTwoMuonandPionPosPass <<" , WhetherTwoMuonandPionNegPass: "<< WhetherTwoMuonandPionNegPass <<std::endl;
          std::cout<<" A: "<< WhetherThreeMuonPosPass <<" , B: "<< WhetherThreeMuonNegPass <<" , C: "<< WhetherTwoMuonandPionPosPass <<" , D: "<< WhetherTwoMuonandPionNegPass <<" , OS1: "<< WhetherThreeMuonPosPassWithOS <<" , OS2: "<< WhetherThreeMuonNegPassWithOS <<" , Hadrons: "<< WhetherThreeMuonOSHadron <<" , eOrHadrons: "<< WhetherThreeMuonOSeOrHadron <<" , onlyMu: "<< WhetherThreeMuonOSmuOnly <<std::endl;
          return true;
  }
  

  return false;
  
}








bool Pythia8Hadronizer::generatePartonsAndHadronize()
{

  DJR.resize(0);
  nME = -1;
  nMEFiltered = -1;
  
  if ( fJetMatchingHook.get() ) 
  {
    fJetMatchingHook->resetMatchingStatus(); 
    fJetMatchingHook->beforeHadronization( lheEvent() );
  }

  //------------------------- redecay B/D's  according to nRepeat
  for(int iC : ParticlesIDtoRedecay)fMasterGen->particleData.mayDecay( iC, false);


 if (!fMasterGen->next())
    return false;

  vector<int> iBHad;
  int nBHad = 0;

  Event *pythiaEvent = &(fMasterGen->event);
  //  Event *savedEvent;
  int nBquark = 0;
  int stat;

  for (int i = 0; i < pythiaEvent->size(); ++i) {
    stat = abs(pythiaEvent->at(i).status());
    if ( abs(pythiaEvent->at(i).id()) == 5 && (stat == 62 || stat == 63)) ++nBquark;
  }

  iBHad.resize(0);
  for (int i = 0; i < pythiaEvent->size(); ++i) {
    int idAbs = abs(pythiaEvent->at(i).id());
    for(int iC : ParticlesIDtoRedecay)
      if (idAbs == iC) {
	iBHad.push_back(i);
	break;
      }
  }
  
  nBHad = iBHad.size();
   // if (nBquark != nBHad) cout << " Warning: " << nBquark
   // 			     << " b quarks but " << nBHad << " B hadrons" << endl;
   

  pythiaEvent->saveSize();
  for(int iC : ParticlesIDtoRedecay)      fMasterGen->particleData.mayDecay( iC, true);



  int final_repetition(1); //the actual no where the loop is broken
  bool whether_evt_pass(false);
  //long long final_repetition = 1000000000000000000LL; //the actual no where the loop is broken. Events that don't pass should have low weights.
  for (int iRepeat = 0; iRepeat < nRepeat; ++iRepeat) {

    if (iRepeat > 0) {
      pythiaEvent->restoreSize();
      // Repeated decays: mark decayed B hadrons as undecayed.
      for (int iB = 0; iB < nBHad; ++iB) pythiaEvent->at(iBHad.at(iB)).statusPos();
    }
    
    if (!fMasterGen->moreDecays()) continue;
    if(ReDecayConditions=="TwoMuMass")
      {
	if(TwoMuMassFilter(*pythiaEvent)) { std::cout<<"nRepeat  "<< iRepeat <<std::endl; final_repetition = iRepeat+1; whether_evt_pass = true; break;}
      }
    if(ReDecayConditions=="ThreeMuMass")
      {
	if(ThreeMuMassFilter(*pythiaEvent)){ std::cout<<"nRepeat  "<< iRepeat <<std::endl; final_repetition = iRepeat+1; whether_evt_pass = true; break;}
      }
  }

  //  return false if gluon with status > 0
  for (int i = 0; i < pythiaEvent->size(); ++i) {
    if ( abs(pythiaEvent->at(i).id()) == 21 && pythiaEvent->at(i).status() > 0 ) return false;
  }


//------------------------- redecay B/D's 


  double mergeweight = fMasterGen.get()->info.mergingWeightNLO();
  if (fMergingHook.get()) {
    mergeweight *= fMergingHook->getNormFactor();
  }
  
  //protect against 0-weight from ckkw or similar
  if (std::abs(mergeweight)==0.)
  {
    event().reset();
    return false;
  }
  
  if (fJetMatchingPy8InternalHook.get()) {
    const std::vector<double> djrmatch = fJetMatchingPy8InternalHook->getDJR();
    //cap size of djr vector to save storage space (keep only up to first 6 elements)
    unsigned int ndjr = std::min(djrmatch.size(), std::vector<double>::size_type(6));
    for (unsigned int idjr=0; idjr<ndjr; ++idjr) {
      DJR.push_back(djrmatch[idjr]);
    }
    
    nME=fJetMatchingPy8InternalHook->nMEpartons().first;
    nMEFiltered=fJetMatchingPy8InternalHook->nMEpartons().second;
  }
  
  if (evtgenDecays.get()) evtgenDecays->decay();

  event().reset(new HepMC::GenEvent);
  bool py8hepmc =  toHepMC.fill_next_event( *(fMasterGen.get()), event().get());

  if (!py8hepmc) {
    return false;
  }
  
  // apply 1/nRepeat to all weights
  double invN = 1.0 / final_repetition;
  for (auto& w : event()->weights()){
          if(whether_evt_pass){
                  w *= invN;
          }
          else{
                  w *= 0.0;
          }
  }
  
  //add ckkw/umeps/unlops merging weight
  if (mergeweight!=1.) {
    event()->weights()[0] *= mergeweight;
  }
  
  if (fEmissionVetoHook.get()) {
    nISRveto += fEmissionVetoHook->getNISRveto();
    nFSRveto += fEmissionVetoHook->getNFSRveto();  
  }
  
  //fill additional weights for systematic uncertainties
  if (fMasterGen->info.getWeightsDetailedSize() > 0) {
    for (const string &key : fMasterGen->info.initrwgt->weightsKeys) {
      double wgt = (*fMasterGen->info.weights_detailed)[key];
      event()->weights().push_back(wgt);
    }
  }
  else if (fMasterGen->info.getWeightsCompressedSize() > 0) {
    for (unsigned int i = 0; i < fMasterGen->info.getWeightsCompressedSize(); i++) {
      double wgt = fMasterGen->info.getWeightsCompressedValue(i);
      event()->weights().push_back(wgt);
    }
  }

  // fill shower weights 
  // http://home.thep.lu.se/~torbjorn/pythia82html/Variations.html
  if( fMasterGen->info.nWeights() > 1 ){
    for(int i = 0; i < fMasterGen->info.nWeights(); ++i) {
      double wgt = fMasterGen->info.weight(i);
      event()->weights().push_back(wgt);
    }
  }
  
  // VINCIA shower weights
  // http://vincia.hepforge.org/current/share/Vincia/htmldoc/VinciaUncertainties.html
  if( fvincia.get() ) {
    event()->weights()[0] *= fvincia->weight(0);
    for (int iVar=1; iVar < fvincia->nWeights(); iVar++) {
      event()->weights().push_back(fvincia->weight(iVar));
    }
  }
  
  // Retrieve Dire shower weights
  if( fDire.get() ) {
    fDire->weightsPtr->calcWeight(0.);
    fDire->weightsPtr->reset();
    
    //Make sure the base weight comes first
    event()->weights()[0] *= fDire->weightsPtr->getShowerWeight("base");
    
    map<string, double>::iterator it;
    for ( it = fDire->weightsPtr->getShowerWeights()->begin(); it != fDire->weightsPtr->getShowerWeights()->end(); it++ ) {
      if (it->first == "base") continue;
      event()->weights().push_back(it->second);
    }
  }

  return true;
  
}


bool Pythia8Hadronizer::hadronize()
{
  DJR.resize(0);
  nME = -1;
  nMEFiltered = -1;
  if(LHEInputFileName.empty()) lhaUP->loadEvent(lheEvent());

  if ( fJetMatchingHook.get() ) 
  {
    fJetMatchingHook->resetMatchingStatus(); 
    fJetMatchingHook->beforeHadronization( lheEvent() );
  }

  bool py8next = fMasterGen->next();

  double mergeweight = fMasterGen.get()->info.mergingWeightNLO();
  if (fMergingHook.get()) {
    mergeweight *= fMergingHook->getNormFactor();
  }
  
  
  //protect against 0-weight from ckkw or similar
  if (!py8next || std::abs(mergeweight)==0.)
  {
    lheEvent()->count( lhef::LHERunInfo::kSelected, 1.0, mergeweight );
    event().reset();
    return false;
  }
  
  if (fJetMatchingPy8InternalHook.get()) {
    const std::vector<double> djrmatch = fJetMatchingPy8InternalHook->getDJR();
    //cap size of djr vector to save storage space (keep only up to first 6 elements)
    unsigned int ndjr = std::min(djrmatch.size(), std::vector<double>::size_type(6));
    for (unsigned int idjr=0; idjr<ndjr; ++idjr) {
      DJR.push_back(djrmatch[idjr]);
    }
    
    nME=fJetMatchingPy8InternalHook->nMEpartons().first;
    nMEFiltered=fJetMatchingPy8InternalHook->nMEpartons().second;
  }
  
  // update LHE matching statistics
  //
  lheEvent()->count( lhef::LHERunInfo::kAccepted, 1.0, mergeweight );

  if (evtgenDecays.get()) evtgenDecays->decay();

  event().reset(new HepMC::GenEvent);
  bool py8hepmc =  toHepMC.fill_next_event( *(fMasterGen.get()), event().get());
  if (!py8hepmc) {
    return false;
  }
  
  //add ckkw/umeps/unlops merging weight
  if (mergeweight!=1.) {
    event()->weights()[0] *= mergeweight;
  }

  if (fEmissionVetoHook.get()) {
    nISRveto += fEmissionVetoHook->getNISRveto();
    nFSRveto += fEmissionVetoHook->getNFSRveto();  
  }

  // fill shower weights
  // http://home.thep.lu.se/~torbjorn/pythia82html/Variations.html
  if( fMasterGen->info.nWeights() > 1 ){
    for(int i = 0; i < fMasterGen->info.nWeights(); ++i) {
      double wgt = fMasterGen->info.weight(i);
      event()->weights().push_back(wgt);
    }
  }

  return true;

}


bool Pythia8Hadronizer::residualDecay()
{

  Event* pythiaEvent = &(fMasterGen->event);

  int NPartsBeforeDecays = pythiaEvent->size();
  int NPartsAfterDecays = event().get()->particles_size();

  if(NPartsAfterDecays == NPartsBeforeDecays) return true;

  bool result = true;

  for ( int ipart=NPartsAfterDecays; ipart>NPartsBeforeDecays; ipart-- )
  {

    HepMC::GenParticle* part = event().get()->barcode_to_particle( ipart );

    if ( part->status() == 1 && (fDecayer->particleData).canDecay(part->pdg_id()) )
    {
      fDecayer->event.reset();
      Particle py8part(  part->pdg_id(), 93, 0, 0, 0, 0, 0, 0,
                         part->momentum().x(),
                         part->momentum().y(),
                         part->momentum().z(),
                         part->momentum().t(),
                         part->generated_mass() );
      HepMC::GenVertex* ProdVtx = part->production_vertex();
      py8part.vProd( ProdVtx->position().x(), ProdVtx->position().y(),
                     ProdVtx->position().z(), ProdVtx->position().t() );
      py8part.tau( (fDecayer->particleData).tau0( part->pdg_id() ) );
      fDecayer->event.append( py8part );
      int nentries = fDecayer->event.size();
      if ( !fDecayer->event[nentries-1].mayDecay() ) continue;
      fDecayer->next();
      int nentries1 = fDecayer->event.size();
      if ( nentries1 <= nentries ) continue; //same number of particles, no decays...

      part->set_status(2);

      result = toHepMC.fill_next_event( *(fDecayer.get()), event().get(), -1, true, part);

    }
  }

  return result;

}




//bool Pythia8Hadronizer::finalizeEvent()


void Pythia8Hadronizer::finalizeEvent()
{
  bool lhe = lheEvent() != nullptr;

  // now create the GenEventInfo product from the GenEvent and fill
  // the missing pieces
  eventInfo().reset( new GenEventInfoProduct( event().get() ) );

  // in pythia pthat is used to subdivide samples into different bins
  // in LHE mode the binning is done by the external ME generator
  // which is likely not pthat, so only filling it for Py6 internal mode
  if (!lhe) {
    eventInfo()->setBinningValues(std::vector<double>(1, fMasterGen->info.pTHat()));
  }
  
  eventInfo()->setDJR(DJR);
  eventInfo()->setNMEPartons(nME);
  eventInfo()->setNMEPartonsFiltered(nMEFiltered);

  //******** Verbosity ********

  if (maxEventsToPrint > 0 &&
      (pythiaPylistVerbosity || pythiaHepMCVerbosity ||
                                pythiaHepMCVerbosityParticles) ) {
    maxEventsToPrint--;
    if (pythiaPylistVerbosity) {
      fMasterGen->info.list(); 
      fMasterGen->event.list();
    } 

    if (pythiaHepMCVerbosity) {
      std::cout << "Event process = "
                << fMasterGen->info.code() << "\n"
                << "----------------------" << std::endl;
      event()->print();
    }
    if (pythiaHepMCVerbosityParticles) {
      std::cout << "Event process = "
                << fMasterGen->info.code() << "\n"
                << "----------------------" << std::endl;
      ascii_io->write_event(event().get());
    }
  }
}

GenLumiInfoHeader *Pythia8Hadronizer::getGenLumiInfoHeader() const {
  GenLumiInfoHeader *genLumiInfoHeader = BaseHadronizer::getGenLumiInfoHeader();
  
  //fill lhe headers
  //*FIXME* initrwgt header is corrupt due to pythia bug
  for (const std::string &key : fMasterGen->info.headerKeys()) {
    genLumiInfoHeader->lheHeaders().emplace_back(key,fMasterGen->info.header(key));
  }

  //check, if it is not only nominal weight
  int weights_number = fMasterGen->info.nWeights();
  if (fMasterGen->info.initrwgt) weights_number += fMasterGen->info.initrwgt->weightsKeys.size();
  if(weights_number > 1){
    genLumiInfoHeader->weightNames().reserve(weights_number + 1);
    genLumiInfoHeader->weightNames().push_back("nominal");
  }

  //fill weight names
  if (fMasterGen->info.initrwgt) {
    for (const std::string &key : fMasterGen->info.initrwgt->weightsKeys) {
      std::string weightgroupname;
      for (const auto &wgtgrp : fMasterGen->info.initrwgt->weightgroups) {
        const auto &wgtgrpwgt = wgtgrp.second.weights.find(key);
        if (wgtgrpwgt != wgtgrp.second.weights.end()) {
          weightgroupname = wgtgrp.first;
        }
      }
      
      std::ostringstream weightname;
      weightname << "LHE, id = " << key << ", ";
      if (!weightgroupname.empty()) {
        weightname << "group = " << weightgroupname << ", ";
      }
      weightname<< fMasterGen->info.initrwgt->weights[key].contents;
      genLumiInfoHeader->weightNames().push_back(weightname.str());    
    }
  }

  //fill shower labels
  // http://home.thep.lu.se/~torbjorn/pythia82html/Variations.html
  // http://home.thep.lu.se/~torbjorn/doxygen/classPythia8_1_1Info.html
  if( fMasterGen->info.nWeights() > 1 ){
    for(int i = 0; i < fMasterGen->info.nWeights(); ++i) {
      genLumiInfoHeader->weightNames().push_back( fMasterGen->info.weightLabel(i) );
    }
  }
  
  // VINCIA shower weights
  // http://vincia.hepforge.org/current/share/Vincia/htmldoc/VinciaUncertainties.html
  if( fvincia.get() ) {
    for (int iVar=0; iVar < fvincia->nWeights(); iVar++) {
      genLumiInfoHeader->weightNames().push_back( fvincia->weightLabel(iVar) );
    }
  }
  
  if( fDire.get() ) {
    //Make sure the base weight comes first
    genLumiInfoHeader->weightNames().push_back("base");
    
    map<string, double>::iterator it;
    for ( it = fDire->weightsPtr->getShowerWeights()->begin(); it != fDire->weightsPtr->getShowerWeights()->end(); it++ ) {
      if (it->first == "base") continue;
      genLumiInfoHeader->weightNames().push_back(it->first);
    }
  }

  return genLumiInfoHeader;
}







typedef edm::GeneratorFilter<Pythia8Hadronizer, ExternalDecayDriver> Pythia8GeneratorFilter;
DEFINE_FWK_MODULE(Pythia8GeneratorFilter);


typedef edm::HadronizerFilter<Pythia8Hadronizer, ExternalDecayDriver> Pythia8HadronizerFilter;
DEFINE_FWK_MODULE(Pythia8HadronizerFilter);
