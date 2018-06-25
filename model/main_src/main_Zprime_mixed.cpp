/*-------------------------------------------+
|   Main event-interface skeleton code       |
|      for OPTIMASS model classes            |
+--------------------+-----------------------+
          < Model information >
+--------------------------------------------+
#1) Model Name : 'Zprime_mixed'
+--------------------------------------------+
#2) Full Decay Topology :
+--------------------------------------------+
    chain-(1) : Zp - n1 n2 , ( n1 - l1 l2) , ( n2 - v1 v2 )
    chain-(2) : Zpp - n1p n2p , ( n1p - l1p l2p ) , ( n2p - v1p v2p )
+--------------------------------------------+
#3) Full Invisible Leaf Particles :
+--------------------------------------------+
    ['v1', 'v1p', 'v2', 'v2p']
+--------------------------------------------+
#4) Target Mother Particles :
+--------------------------------------------+
    ['Zp', 'Zpp']
+--------------------------------------------+
#5) Effective Invisible Particles :
+--------------------------------------------+
    ['v1', 'v2', 'v1p', 'v2p']
+--------------------------------------------+
#6) Constraint Function List :
+--------------------------------------------+
    constraint-(1) : n1.M() - n2.M() = 0 
    constraint-(2) : n1p.M() - n2p.M() = 0 
    constraint-(3) : n1.M() - n1p.M() = 0 
    constraint-(4) : Zp.M() - Zpp.M() = 0
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

#include "Zprime_mixed.h"

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
        optm = static_cast<OptiMass::MassMinimizer *>(new OptiMass::Zprime_mixed);
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
        optm->SetInvisibleSubsystemMomenta(0,  ,  );// (i-1, PTmiss_tot.Px, PTmiss_tot.Py) for the (i=1)-th PT-conservation system
     optm->SetInvisibleSubsystemMomenta(1,  ,  );// (i-1, PTmiss_tot.Px, PTmiss_tot.Py) for the (i=2)-th PT-conservation system
 
        // Visible Momentum Input :
        optm->SetMomentumValue("l1_x",   );
        optm->SetMomentumValue("l1_y",   );
        optm->SetMomentumValue("l1_z",   );
        optm->SetMomentumValue("l1_m",   );
        optm->SetMomentumValue("l2_x",   );
        optm->SetMomentumValue("l2_y",   );
        optm->SetMomentumValue("l2_z",   );
        optm->SetMomentumValue("l2_m",   );
        optm->SetMomentumValue("l1p_x",   );
        optm->SetMomentumValue("l1p_y",   );
        optm->SetMomentumValue("l1p_z",   );
        optm->SetMomentumValue("l1p_m",   );
        optm->SetMomentumValue("l2p_x",   );
        optm->SetMomentumValue("l2p_y",   );
        optm->SetMomentumValue("l2p_z",   );
        optm->SetMomentumValue("l2p_m",   );

        // Invisible Particle (Trial) Mass Input : 
        optm->SetInitInvisibleMomentum("v1_m",   );
        optm->SetInitInvisibleMomentum("v2_m",   );
        optm->SetInitInvisibleMomentum("v1p_m",   );
        optm->SetInitInvisibleMomentum("v2p_m",   );    
    
       /*-------------------------------------------------------------+
        | Constraint Function Switching :                             |
        | >> Check the definition of the functions in the list above  |
        +-------------------------------------------------------------*/
        alm_controller.UseConstraint(0, true);// (true or false) for the constraint : n1.M() - n2.M() = 0.
        alm_controller.UseConstraint(1, true);// (true or false) for the constraint : n1p.M() - n2p.M() = 0.
        alm_controller.UseConstraint(2, true);// (true or false) for the constraint : n1.M() - n1p.M() = 0.
        alm_controller.UseConstraint(3, true);// (true or false) for the constraint : Zp.M() - Zpp.M() = 0.
 
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
        std::cout << "Zp-mass = " << process_tree.GetSubsystemMass("Zp") << std::endl;
        std::cout << "n1-mass = " << process_tree.GetSubsystemMass("n1") << std::endl;
        std::cout << "n2-mass = " << process_tree.GetSubsystemMass("n2") << std::endl;
        std::cout << "v1-mass = " << process_tree.GetSubsystemMass("v1") << std::endl;
        std::cout << "v2-mass = " << process_tree.GetSubsystemMass("v2") << std::endl;
        std::cout << "Zpp-mass = " << process_tree.GetSubsystemMass("Zpp") << std::endl;
        std::cout << "n1p-mass = " << process_tree.GetSubsystemMass("n1p") << std::endl;
        std::cout << "n2p-mass = " << process_tree.GetSubsystemMass("n2p") << std::endl;
        std::cout << "v1p-mass = " << process_tree.GetSubsystemMass("v1p") << std::endl;
        std::cout << "v2p-mass = " << process_tree.GetSubsystemMass("v2p") << std::endl;
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
