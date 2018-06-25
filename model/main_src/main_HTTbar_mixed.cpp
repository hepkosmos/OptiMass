/*-------------------------------------------+
|   Main event-interface skeleton code       |
|      for OPTIMASS model classes            |
+--------------------+-----------------------+
          < Model information >
+--------------------------------------------+
#1) Model Name : 'HTTbar_mixed'
+--------------------------------------------+
#2) Full Decay Topology :
+--------------------------------------------+
    chain-(1) : H - t1 t2 , ( t1 - b1 w1 , w1 - e1 v1 ) , ( t2 - b2 w2 , w2 - e2 v2 )
    chain-(2) : Hp - t1p t2p , ( t1p - b1p w1p , w1p - e1p v1p ) , ( t2p - b2p w2p , w2p - e2p v2p )
+--------------------------------------------+
#3) Full Invisible Leaf Particles :
+--------------------------------------------+
    ['v1', 'v2', 'v1p', 'v2p']
+--------------------------------------------+
#4) Target Mother Particles :
+--------------------------------------------+
    ['H', 'Hp']
+--------------------------------------------+
#5) Effective Invisible Particles :
+--------------------------------------------+
    ['v1', 'v2', 'v1p', 'v2p']
+--------------------------------------------+
#6) Constraint Function List :
+--------------------------------------------+
    constraint-(1) : t1.M() - t2.M() = 0 
    constraint-(2) : w1.M() - w2.M() = 0 
    constraint-(3) : t1.M() - 173. = 0 
    constraint-(4) : w1.M() - 80. = 0 
    constraint-(5) : t1p.M() - t2p.M() = 0 
    constraint-(6) : w1p.M() - w2p.M() = 0 
    constraint-(7) : t1p.M() - 173. = 0 
    constraint-(8) : w1p.M() - 80. = 0 
    constraint-(9) : H.M() - Hp.M() = 0
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

#include "HTTbar_mixed.h"

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
        optm = static_cast<OptiMass::MassMinimizer *>(new OptiMass::HTTbar_mixed);
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
        optm->SetMomentumValue("b1_x",   );
        optm->SetMomentumValue("b1_y",   );
        optm->SetMomentumValue("b1_z",   );
        optm->SetMomentumValue("b1_m",   );
        optm->SetMomentumValue("e1_x",   );
        optm->SetMomentumValue("e1_y",   );
        optm->SetMomentumValue("e1_z",   );
        optm->SetMomentumValue("e1_m",   );
        optm->SetMomentumValue("b2_x",   );
        optm->SetMomentumValue("b2_y",   );
        optm->SetMomentumValue("b2_z",   );
        optm->SetMomentumValue("b2_m",   );
        optm->SetMomentumValue("e2_x",   );
        optm->SetMomentumValue("e2_y",   );
        optm->SetMomentumValue("e2_z",   );
        optm->SetMomentumValue("e2_m",   );
        optm->SetMomentumValue("b1p_x",   );
        optm->SetMomentumValue("b1p_y",   );
        optm->SetMomentumValue("b1p_z",   );
        optm->SetMomentumValue("b1p_m",   );
        optm->SetMomentumValue("e1p_x",   );
        optm->SetMomentumValue("e1p_y",   );
        optm->SetMomentumValue("e1p_z",   );
        optm->SetMomentumValue("e1p_m",   );
        optm->SetMomentumValue("b2p_x",   );
        optm->SetMomentumValue("b2p_y",   );
        optm->SetMomentumValue("b2p_z",   );
        optm->SetMomentumValue("b2p_m",   );
        optm->SetMomentumValue("e2p_x",   );
        optm->SetMomentumValue("e2p_y",   );
        optm->SetMomentumValue("e2p_z",   );
        optm->SetMomentumValue("e2p_m",   );

        // Invisible Particle (Trial) Mass Input : 
        optm->SetInitInvisibleMomentum("v1_m",   );
        optm->SetInitInvisibleMomentum("v2_m",   );
        optm->SetInitInvisibleMomentum("v1p_m",   );
        optm->SetInitInvisibleMomentum("v2p_m",   );    
    
       /*-------------------------------------------------------------+
        | Constraint Function Switching :                             |
        | >> Check the definition of the functions in the list above  |
        +-------------------------------------------------------------*/
        alm_controller.UseConstraint(0, true);// (true or false) for the constraint : t1.M() - t2.M() = 0.
        alm_controller.UseConstraint(1, true);// (true or false) for the constraint : w1.M() - w2.M() = 0.
        alm_controller.UseConstraint(2, true);// (true or false) for the constraint : t1.M() - 173. = 0.
        alm_controller.UseConstraint(3, true);// (true or false) for the constraint : w1.M() - 80. = 0.
        alm_controller.UseConstraint(4, true);// (true or false) for the constraint : t1p.M() - t2p.M() = 0.
        alm_controller.UseConstraint(5, true);// (true or false) for the constraint : w1p.M() - w2p.M() = 0.
        alm_controller.UseConstraint(6, true);// (true or false) for the constraint : t1p.M() - 173. = 0.
        alm_controller.UseConstraint(7, true);// (true or false) for the constraint : w1p.M() - 80. = 0.
        alm_controller.UseConstraint(8, true);// (true or false) for the constraint : H.M() - Hp.M() = 0.
 
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
        std::cout << "H-mass = " << process_tree.GetSubsystemMass("H") << std::endl;
        std::cout << "t1-mass = " << process_tree.GetSubsystemMass("t1") << std::endl;
        std::cout << "w1-mass = " << process_tree.GetSubsystemMass("w1") << std::endl;
        std::cout << "v1-mass = " << process_tree.GetSubsystemMass("v1") << std::endl;
        std::cout << "t2-mass = " << process_tree.GetSubsystemMass("t2") << std::endl;
        std::cout << "w2-mass = " << process_tree.GetSubsystemMass("w2") << std::endl;
        std::cout << "v2-mass = " << process_tree.GetSubsystemMass("v2") << std::endl;
        std::cout << "Hp-mass = " << process_tree.GetSubsystemMass("Hp") << std::endl;
        std::cout << "t1p-mass = " << process_tree.GetSubsystemMass("t1p") << std::endl;
        std::cout << "w1p-mass = " << process_tree.GetSubsystemMass("w1p") << std::endl;
        std::cout << "v1p-mass = " << process_tree.GetSubsystemMass("v1p") << std::endl;
        std::cout << "t2p-mass = " << process_tree.GetSubsystemMass("t2p") << std::endl;
        std::cout << "w2p-mass = " << process_tree.GetSubsystemMass("w2p") << std::endl;
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
