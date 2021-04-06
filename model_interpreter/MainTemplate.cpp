{proc_info_}
#include <cmath>
#include <string>
#include <iostream>
//#include "TH1.h"
//#include "TFile.h"
//#include "TTree.h"
//#include "TChain.h"
//#include "TCanvas.h"
//#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TVector2.h"

//#include "ExRootClasses.h"
//#include "ExRootTreeReader.h"

#include "MassMinimizer.h"
#include "ProcessTree.h"
#include "ALMController.h"

{header_}

int main(){

   /*--------------------------+
    | ROOT Error ignore level  |
    +--------------------------*/
    gErrorIgnoreLevel = 1001;

    // * Declaration of 
	//   1) four momenta of visible particles 2) missing ET vector(s)
	//   for each PT-conserving subsystem for each Process
	// --------------------------------------------------------------
	// => check the label of the particles in the '2) Full Decay Topolgy' of the summary above.
{mom_decl_}
	// * Turning-on/off switches of the constraints 
	// --------------------------------------------
	// => Check the definition and order of the constraints in the '6) Constraint Function List' above. 
{constraint_switch_}
    // * Output containers
	// -------------------
{outputs_}
	// * Declaration of instances of OptiMass processes
	//-------------------------------------------------
{om_instances_}
	// * Invisible particle masses (!) Modify to your true/trial value. (!)
	// ----------------------------
{inv_masses_}
	// * Loading Events and Running OptiMass
	// -------------------------------------
// for(int entry = 0; entry < NEvent; entry++){

	// (1) Assignment of visible four momenta (!) Modify to your event(s). (!)
{mom_assign_}
	// (2) Assignment of missing ET vectors
{met_assign_}	
	// (3) Running OptiMass
{run_optimass_}
	// * Checking The Results
	// -----------------------
{results_}
//} (!) Event Loop (!)

	// * Removing the instances of OptiMass processes 
	//-----------------------------------------------
{rm_om_instances_}

}

