/*-------------------------------------------+
|   Main event-interface skeleton code       |
|      for OPTIMASS model classes            |
+--------------------+-----------------------+
          < Model information >
+--------------------------------------------+
#1) Model Name : 'MainTemplate'
+--------------------------------------------+
#2) Full Decay Topology :
+--------------------------------------------+
    {process_syntax_}
+--------------------------------------------+
#3) Full Invisible Leaf Particles :
+--------------------------------------------+
    {particle_invisibles_}
+--------------------------------------------+
#4) Target Mother Particles :
+--------------------------------------------+
    {target_particles_}
+--------------------------------------------+
#5) Effective Invisible Particles :
+--------------------------------------------+
    {effective_invisible_particles_}
+--------------------------------------------+
#6) Constraint Function List :
+--------------------------------------------+
    {constraint_function_list_}
+--------------------------------------------*/
#include <cmath>
#include <string>
#include <iostream>
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

#include "MassMinimizer.h"
#include "ProcessTree.h"
#include "ALMController.h"

#include "ExRootClasses.h"
#include "ExRootTreeReader.h"

#include "MainTemplate.h"

int main(int argc, char** argv){

   /*--------------------------+
    | ROOT Error ignore level  |
    +--------------------------*/
    gErrorIgnoreLevel = 1001;

   /*--------------------------------------------------------------+
    | Load user's event files (USER-DEFINED ROUTINE)               |   
    | Example #1) ExRootAnalysis                                   |
    | Ref1) http://madgraph.hep.uiuc.edu/Downloads/ExRootAnalysis/ |
    | Ref2) http://pdg.lbl.gov/2002/montecarlorpp.pdf              |
    +--------------------------------------------------------------*/

    TChain chain("LHEF");
    // Add your event '<user-events>.root' file in LHEF format
    chain.Add("<user-events>.root");
    ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
    long int NEvent = treeReader->GetEntries();
    TClonesArray *branchParticle = treeReader->UseBranch("Particle");
        
    // Loop over all events
    for(int entry = 0; entry < NEvent; entry++){
        treeReader->ReadEntry(entry);

        TLorentzVector bquark1,bquark2,lepton1,lepton2;

        // Loop over particle entries in an event
        for(int iptl=0; iptl<branchParticle->GetEntries(); iptl++){
            TRootLHEFParticle *particle = (TRootLHEFParticle*) branchParticle->At(iptl);
            // For leptons 2
            if(particle->PID==11||particle->PID==13||particle->PID==15){      
                lepton2.SetXYZM(particle->Px,particle->Py,particle->Pz,particle->M);
            }
            // For leptons 1
            if(particle->PID==-11||particle->PID==-13||particle->PID==-15){      
                lepton1.SetXYZM(particle->Px,particle->Py,particle->Pz,particle->M);
            }
            // For bquarks 1
            if(particle->PID==5){
                bquark1.SetXYZM(particle->Px,particle->Py,particle->Pz,particle->M);
            }
            // For bquarks 2
            if(particle->PID==-5){
                bquark2.SetXYZM(particle->Px,particle->Py,particle->Pz,particle->M);
            }
        } // End of for-loop over entries in an event

     

       /*---------------------------------------------------+
        | Initialization of an OPTIMASS Calculation Object  |
        +---------------------------------------------------*/
        OptiMass::MassMinimizer *optm = 0;
        optm = static_cast<OptiMass::MassMinimizer *>(new OptiMass::MainTemplate);
        optm->InitContainers();
        OptiMass::ProcessTree& process_tree = optm->GetProcessTree();
        OptiMass::ALMController& alm_controller = optm->GetALMController();

       /*----------------------------+
        | Event Input Function Call  |    
        +----------------------------*/
        /*
        TVector2 PTmiss_tot;
        PTmiss_tot.Set(,);
        */

        // MET Input : 
        // Method (1) for the system with single MET condition 
        //optm->SetMissingET(PTmiss_tot.Px(),PTmiss_tot.Py());  // 
        // Method (2) for general systems with multiple MET conditions
        {met_input_}
 
        // Visible Momentum Input :
        {visible_momentum_input_}

        // Invisible Particle (Trial) Mass Input : 
        {invisible_particle_mass_input_}    
    
       /*-------------------------------------------------------------+
        | Constraint Function Switching :                             |
        | >> Check the definition of the functions in the list above  |
        +-------------------------------------------------------------*/
        {constraint_function_switch_}
 
       /*-----------------------+
        | OptiMass Calculation  |
        +-----------------------*/
        optm->Calc();

        //---------------------------------------+
        //  Print output of the current result   |
        //---------------------------------------+
        std::cout << "==================================================" << std::endl;
        std::cout << "Event ID    = "<< entry << std::endl;
        std::cout << "OptiMass    = "<< sqrt(optm->GetOptiMass()) << std::endl;
        std::cout << "==================================================" << std::endl;
        std::cout << "<Reconstructed Particle Masses of Interest> " << std::endl;
        {reconstructed_particle_masses_}
        std::cout << "==================================================" << std::endl;
        std::cout << "modConst    = "<< alm_controller.GetSumSquaredConstraints() << std::endl; // norm of the constraint function vector (only for turned-on)
        std::cout << "nIter       = "<< alm_controller.GetNumberIteration() << std::endl;
        std::cout << "nALM-Phase1 = "<< alm_controller.GetNumberPhase1() << std::endl;
        std::cout << "nALM-Phase2 = "<< alm_controller.GetNumberPhase2() << std::endl;
        std::cout << "==================================================" << std::endl;

        //-----------------------------------------+
        // Reset the dynamically allocated object  |
        //-----------------------------------------+
        delete optm;
 
    } // End of forloop over event entries
 
}
