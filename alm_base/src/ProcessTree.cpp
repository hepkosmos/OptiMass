#include "ProcessTree.h"

#include <cstdlib>
#include <map>
#include <set>
#include <unordered_map>
#include <iostream>

#include "Minuit2/MnUserParameters.h"

#include "StringUtils.h"
#include "Types.h"

// ========================================================================
// class ProcessTree: Constructor and Destructor
// ========================================================================
OptiMass::ProcessTree::ProcessTree(bool debugInput) : debug(debugInput){
}

OptiMass::ProcessTree::~ProcessTree(){
    for(auto it_node = processChain.begin(),
            it_node_end = processChain.end();
            it_node != it_node_end;
            ++it_node
    ){
        delete *it_node;
    }
}
// ========================================================================
// class ProcessTree: Member Functions
// ========================================================================
void OptiMass::ProcessTree::PrintDict(){
/*
    std::unordered_map< std::string, std::unordered_map< std::string, std::string > > dic;
    for(auto it = dic.begin(); it != dic.end(); ++it){
        std::cout << it->first << " contains" << std::endl;
        std::cout << "---------------------------------" << std::endl;
        for(auto it2 = it->second.begin(); it2 != it->second.end(); ++it2){
            std::cout << it2->first << " | " << it2->second << std:: endl;
        }
        std::cout << "=================================" << std::endl;
    }
*/
    printPtlDict(map_ptl_component_index_);
    printPtlDict(map_ptl_component_effective_index_);
}

void OptiMass::ProcessTree::ParseProcess(){
    if(debug)
        std::cout << "Interpreting process syntax" << std::endl;
    for(auto it_cmd = processCmd.begin(),
             it_cmd_end = processCmd.end();
            it_cmd != it_cmd_end;
            ++it_cmd
    ){
        processChain.push_back(processParser(*it_cmd));
    }
    if(debug)
        std::cout << "update final states" << std::endl;
    for(auto it_node = processChain.begin(),
             it_node_end = processChain.end();
            it_node != it_node_end;
            ++it_node
    ){
        auto vec_final_states = (*it_node)->getFinalStates();
        for(auto it_name = vec_final_states.begin(),
                 it_name_end = vec_final_states.end();
                it_name != it_name_end;
                ++it_name
        ){
            vec_ptl_final_state_.push_back(*it_name);
        }
    }
    if(debug)
        std::cout << "update final states visible and invisible" << std::endl;
    for(auto it_s = vec_ptl_final_state_.begin(),
             it_s_end = vec_ptl_final_state_.end();
            it_s != it_s_end;
            ++it_s
    ){
        bool is_visible = true;
        for(auto it_label = vec_ptl_invisible_state_.begin(),
                 it_label_end = vec_ptl_invisible_state_.end();
                it_label != it_label_end;
                ++it_label
        ){
            if(it_s->compare(*it_label) == 0){
                vec_ptl_invisible_label_.push_back(*it_s);
                is_visible = false;
            }
        }
        if(is_visible){
            vec_ptl_visible_label_.push_back(*it_s);
        }
    }
    
    GeneratePtlIndexDict();
    GenerateEffectiveDict();
    GenerateMomentumContainers();
    GenerateMomentumDict();
}

OptiMass::ParticleNode* OptiMass::ProcessTree::singleProcessParser(const std::string& arg){
    std::string whitespace = " \n\t";
    std::vector< std::string > processPtls = split(arg,'-');
    std::string initPtl = trim(processPtls[0],whitespace);
    std::vector< std::string > siblingPtls = split(trim(processPtls[1],whitespace),' ');

    ParticleNode* outputNode = new ParticleNode(initPtl,debug);
    vec_state_node_.push_back(outputNode);
    for(auto it = siblingPtls.begin(); it != siblingPtls.end(); ++it){
        ParticleNode* siblingNode = new ParticleNode(*it,debug);
        outputNode->AddSibling(siblingNode);
        vec_state_node_.push_back(siblingNode);
    }

    return outputNode;
}

OptiMass::ParticleNode* OptiMass::ProcessTree::processParser(const std::string& arg){
    std::vector< std::string > processBuf = splitKeepParenthesis(arg);
    std::string firstProcessBuf = processBuf[0];
    std::string argBuf = arg;
    std::string delim = "() ";
    if((trim(firstProcessBuf,delim)).compare(trim(argBuf,delim))==0){
        return singleProcessParser(firstProcessBuf);
    }else{
        std::map< std::string, ParticleNode* > nodeMapBuf;
        for(auto it = processBuf.begin()+1; it != processBuf.end(); ++it){
            ParticleNode* nodeBuf = processParser(*it);
            nodeMapBuf[nodeBuf->GetLabel()]=nodeBuf;
        }

        std::string whitespace = " \n\t";
        std::vector< std::string > processPtls = split(firstProcessBuf,'-');
        std::string initPtl = trim(processPtls[0],whitespace);
        std::vector< std::string > siblingPtls = split(trim(processPtls[1],whitespace),' ');

        ParticleNode* outputNode = new ParticleNode(initPtl,debug);
        vec_state_node_.push_back(outputNode);
        for(auto it = siblingPtls.begin(); it != siblingPtls.end(); ++it){
            auto itFound = nodeMapBuf.find(*it);
            if(itFound != nodeMapBuf.end()){
                outputNode->AddSibling(itFound->second);
                nodeMapBuf.erase(*it);
            }else{
                ParticleNode* siblingNode = new ParticleNode(*it,debug);
                outputNode->AddSibling(siblingNode);
                vec_state_node_.push_back(siblingNode);
            }
        }
        if(!nodeMapBuf.empty()){
            std::cout << "Warning! decay chain is not well-defined" << std::endl;
            for(auto it = nodeMapBuf.begin(); it != nodeMapBuf.end(); ++it){
                std::cout << "isolated particle: " << it->first << std::endl;
                delete (it->second);
            }    
        }
        return outputNode;
    }
}

void OptiMass::ProcessTree::GeneratePtlIndexDict(){
    //std::unordered_map< std::string , ParticleSubelements<int> > output;
    for(auto it_node = vec_state_node_.begin(),
             it_node_end = vec_state_node_.end();
            it_node != it_node_end;
            ++it_node
    ){
        std::vector<int> indexVisibles;
        std::vector<int> indexInvisibles;
        std::vector<int> indexEff;
        auto vec_final_states = (*it_node)->getFinalStates();
        for(auto it_ptl = vec_final_states.begin(),
                 it_ptl_end = vec_final_states.end();
                it_ptl != it_ptl_end;
                ++it_ptl
        ){
            bool updated = false;
            for(unsigned int i=0; (!updated) && i<vec_ptl_visible_label_.size(); ++i) {
                if(it_ptl->compare(vec_ptl_visible_label_[i]) == 0){
                    indexVisibles.push_back(i);
                    updated = true;
                }
            }
            for(unsigned int i=0; (!updated) && i<vec_ptl_invisible_label_.size(); ++i){
                if(it_ptl->compare(vec_ptl_invisible_label_[i]) == 0){
                    indexInvisibles.push_back(i);
                    updated = true;
                }
            }
            if(!updated){
                std::cout << "Warning! Final state particle index is missing for particle: " << (*it_node)->GetLabel() << std::endl;
            }
        }
        ParticleSubelements<int> outputIndex;
        outputIndex.visible = (indexVisibles);
        outputIndex.invisible = (indexInvisibles);
        outputIndex.optimize = (indexEff);
        map_ptl_component_index_[(*it_node)->GetLabel()] = outputIndex;        
    }
}

void OptiMass::ProcessTree::GenerateEffectiveDict(){
	//PtlIndexDict output = original_dict;

	ParticleSubelements<int> buf_ptl_index;
    for(auto it_ptl_element = map_ptl_component_index_.begin(),
             it_ptl_element_end = map_ptl_component_index_.end();
            it_ptl_element != it_ptl_element_end;
            ++it_ptl_element
    ){
    	buf_ptl_index = it_ptl_element->second;
    	for(unsigned int i = 0; i < vec_ptl_optimize_label_.size();++i){
    		if(ptlIn(buf_ptl_index,map_ptl_component_index_[vec_ptl_optimize_label_[i]])){
    			buf_ptl_index = ptlRemove(buf_ptl_index,map_ptl_component_index_[vec_ptl_optimize_label_[i]],(int) i);
    		}
    	}
    	map_ptl_component_effective_index_[it_ptl_element->first] = buf_ptl_index;
    }
}

void OptiMass::ProcessTree::GenerateMomentumContainers(){
    // Find optimize label which is not in final states
    for(auto it_optimize_name = vec_ptl_optimize_label_.begin(),
             it_optimize_name_end = vec_ptl_optimize_label_.end();
            it_optimize_name != it_optimize_name_end;
            ++it_optimize_name
    ){
        bool found = false;
        for(auto it_invisible_name = vec_ptl_invisible_label_.begin(),
             it_invisible_name_end = vec_ptl_invisible_label_.end();
            it_invisible_name != it_invisible_name_end && !found;
            ++it_invisible_name 
        ){
            if(it_optimize_name->compare(*it_invisible_name) == 0){
                found = true;
            }
        }
        if(!found){
            vec_ptl_optimize_not_final_label_.push_back(*it_optimize_name);
        }
    }    
    // Labels for visible momentum
    for(auto it_name = vec_ptl_visible_label_.begin(),
             it_name_end = vec_ptl_visible_label_.end();
            it_name != it_name_end;
            ++it_name ){
        vec_ptl_visible_param_name_.push_back({*(it_name)+"_x",*(it_name)+"_y",*(it_name)+"_z",*(it_name)+"_m"});
    }
    // Labels for invisible momentum
    for(auto it_name = vec_ptl_invisible_label_.begin(),
             it_name_end = vec_ptl_invisible_label_.end();
            it_name != it_name_end;
            ++it_name ){
        vec_ptl_invisible_param_name_.push_back({*(it_name)+"_x",*(it_name)+"_y",*(it_name)+"_z",*(it_name)+"_m"});
    }
    // Labels for momentum to be optimized
    for(auto it_name = vec_ptl_optimize_label_.begin(),
             it_name_end = vec_ptl_optimize_label_.end();
            it_name != it_name_end;
            ++it_name ){
        vec_ptl_optimize_param_name_.push_back({*(it_name)+"_x",*(it_name)+"_y",*(it_name)+"_z",*(it_name)+"_m"});
    }
    // initialize invisible subsystem index
    /////////////////////////////////////////////////
    vec_ptl_invisible_subsystem_calc_auto_.assign(vec_ptl_invisible_subsystem_label_.size(), true);
    // Loop through invisible subsystems
    for(auto it_ptl_invisible_subsystem_label = vec_ptl_invisible_subsystem_label_.begin(),
             it_ptl_invisible_subsystem_label_end = vec_ptl_invisible_subsystem_label_.end();
            it_ptl_invisible_subsystem_label != it_ptl_invisible_subsystem_label_end;
            ++it_ptl_invisible_subsystem_label ) {
        // Loop through particles in a invisible subsystem, 
        std::set< int > subsystem_index;
        for(auto it_ptl_label = it_ptl_invisible_subsystem_label->begin(),
                 it_ptl_label_end = it_ptl_invisible_subsystem_label->end();
                it_ptl_label != it_ptl_label_end;
                ++it_ptl_label
            ) {
            subsystem_index.insert(map_ptl_component_effective_index_[*it_ptl_label].optimize.begin(),map_ptl_component_effective_index_[*it_ptl_label].optimize.end());
        } 
        std::vector< int > vec_output(subsystem_index.begin(),subsystem_index.end());
  
        vec_ptl_invisible_subsystem_index_.push_back(vec_output);
    }
    // Initialize momentum containers
    for(std::vector< std::vector< std::string > >::iterator it = vec_ptl_visible_param_name_.begin(); it != vec_ptl_visible_param_name_.end(); ++it) {
        for(std::vector< std::string >::iterator it2 = it->begin(); it2 != it->end(); ++it2) {
            parMomentum[*it2] = 0.;
        }
    }
    for(std::vector< std::vector< std::string > >::iterator it = vec_ptl_invisible_param_name_.begin(); it != vec_ptl_invisible_param_name_.end(); ++it) {
        for(std::vector< std::string >::iterator it2 = it->begin(); it2 != it->end(); ++it2) {
            parMomentum[*it2] = 0.;
        }
    }

    // initial mass parameters 
    for(auto it_name = vec_ptl_visible_label_.begin(),
             it_name_end = vec_ptl_visible_label_.end();
            it_name != it_name_end;
            ++it_name
        ){
        parMomentum[*it_name + "_m"] = map_ptl_mass_[*it_name];
    }
    for(auto it_name = vec_ptl_invisible_label_.begin(),
             it_name_end = vec_ptl_invisible_label_.end();
            it_name != it_name_end;
            ++it_name
        ){
        parMomentum[*it_name + "_m"] = map_ptl_mass_[*it_name];
    }

    // Initialize TLorentzVector list of particles
    vec_ptl_visible_momenta_.assign(vec_ptl_visible_param_name_.size(),TLorentzVector());
    vec_ptl_invisible_momenta_.assign(vec_ptl_invisible_param_name_.size(),TLorentzVector());
    vec_ptl_optimize_momenta_.assign(vec_ptl_optimize_param_name_.size(),TLorentzVector());
    vec_ptl_invisible_subsystem_momenta_.assign(vec_ptl_invisible_subsystem_label_.size(),TVector2());
    vec_ptl_invisible_subsystem_momenta_eff_.assign(vec_ptl_invisible_subsystem_label_.size(),TVector2());

    // Initialize buffer for evluating MET
    buf_vec_ptl_last_momenta_.assign(vec_ptl_invisible_subsystem_momenta_eff_.size(),TVector2());
}

ROOT::Minuit2::MnUserParameters OptiMass::ProcessTree::GetMnUserParameters(){
    ROOT::Minuit2::MnUserParameters parFitting;
    // Initialize momentum containers
    for(std::vector< std::vector< std::string > >::iterator it = vec_ptl_optimize_param_name_.begin(); it != vec_ptl_optimize_param_name_.end(); ++it) {
        for(std::vector< std::string >::iterator it2 = it->begin(); it2 != it->end(); ++it2) {
            parFitting.Add(*it2,0.,10);
        }
    }
    // Fix invisible parameter mass and last invisible particle's transverse momentum.
    // These parameters will not be used while in Minuit routine
    for(unsigned int i = 0; i < vec_ptl_optimize_param_name_.size();++i) {
         parFitting.Fix(vec_ptl_optimize_param_name_[i][3]);
    } 
    for(auto it_ptl_invisible_subsystem_index = vec_ptl_invisible_subsystem_index_.begin(),
             it_ptl_invisible_subsystem_index_end = vec_ptl_invisible_subsystem_index_.end();
            it_ptl_invisible_subsystem_index != it_ptl_invisible_subsystem_index_end;
            ++it_ptl_invisible_subsystem_index
         ) {
        auto it_last = it_ptl_invisible_subsystem_index->end()-1;
        parFitting.Fix(vec_ptl_optimize_param_name_[*(it_last)][0]);
        parFitting.Fix(vec_ptl_optimize_param_name_[*(it_last)][1]);
    }
    // Put initial mass parameters
    for(auto it_name = vec_ptl_optimize_label_.begin(),
             it_name_end = vec_ptl_optimize_label_.end();
            it_name != it_name_end;
            ++it_name
        ){
        //parFitting.SetValue(*it_name + "_m",map_ptl_mass_[*it_name]); //rev
        parFitting.SetValue(*it_name + "_m",map_ptl_mass_[*it_name]);
	
		//rev
		//std::cout << " map_ptl = " << *it_name << " mass = " << map_ptl_mass_[*it_name] << std::endl;
    }
    return parFitting;
}

void OptiMass::ProcessTree::GenerateMomentumDict() {
    for(auto it_pair_ptl_index = map_ptl_component_index_.begin(),
             it_pair_ptl_index_end = map_ptl_component_index_.end();
            it_pair_ptl_index != it_pair_ptl_index_end;
            ++it_pair_ptl_index
    ){
        ParticleSubelements< TLorentzVector* > vec_momentum_group;

        std::vector< TLorentzVector * > vec_momentum_visible;
        for(auto it_index = it_pair_ptl_index->second.visible.begin(),
                 it_index_end = it_pair_ptl_index->second.visible.end();
                it_index != it_index_end;
                ++it_index
        ){
            vec_momentum_visible.push_back(&(vec_ptl_visible_momenta_[*it_index]));
        }
        std::vector< TLorentzVector * > vec_momentum_invisible;
        for(auto it_index = it_pair_ptl_index->second.invisible.begin(),
                 it_index_end = it_pair_ptl_index->second.invisible.end();
                it_index != it_index_end;
                ++it_index
        ){
            vec_momentum_invisible.push_back(&(vec_ptl_invisible_momenta_[*it_index]));
        }        
        std::vector< TLorentzVector * > vec_momentum_optimize;
        for(auto it_index = it_pair_ptl_index->second.optimize.begin(),
                 it_index_end = it_pair_ptl_index->second.optimize.end();
                it_index != it_index_end;
                ++it_index
        ){
            vec_momentum_optimize.push_back(&(vec_ptl_optimize_momenta_[*it_index]));
        }
        vec_momentum_group.visible = (vec_momentum_visible);
        vec_momentum_group.invisible = (vec_momentum_invisible);
        vec_momentum_group.optimize = (vec_momentum_optimize);

        map_ptl_component_momentum_[it_pair_ptl_index->first] = vec_momentum_group;
    }
    for(auto it_pair_ptl_index = map_ptl_component_effective_index_.begin(),
             it_pair_ptl_index_end = map_ptl_component_effective_index_.end();
            it_pair_ptl_index != it_pair_ptl_index_end;
            ++it_pair_ptl_index
    ){
        ParticleSubelements< TLorentzVector* > vec_momentum_group;

        std::vector< TLorentzVector* > vec_momentum_visible;
        for(auto it_index = it_pair_ptl_index->second.visible.begin(),
                 it_index_end = it_pair_ptl_index->second.visible.end();
                it_index != it_index_end;
                ++it_index
        ){
            vec_momentum_visible.push_back(&(vec_ptl_visible_momenta_[*it_index]));
        }
        std::vector< TLorentzVector* > vec_momentum_invisible;
        for(auto it_index = it_pair_ptl_index->second.invisible.begin(),
                 it_index_end = it_pair_ptl_index->second.invisible.end();
                it_index != it_index_end;
                ++it_index
        ){
            vec_momentum_invisible.push_back(&(vec_ptl_invisible_momenta_[*it_index]));
        }        
        std::vector< TLorentzVector* > vec_momentum_optimize;
        for(auto it_index = it_pair_ptl_index->second.optimize.begin(),
                 it_index_end = it_pair_ptl_index->second.optimize.end();
                it_index != it_index_end;
                ++it_index
        ){
            vec_momentum_optimize.push_back(&(vec_ptl_optimize_momenta_[*it_index]));
        }
        vec_momentum_group.visible = (vec_momentum_visible);
        vec_momentum_group.invisible = (vec_momentum_invisible);
        vec_momentum_group.optimize = (vec_momentum_optimize);

        map_ptl_component_effective_momentum_[it_pair_ptl_index->first] = vec_momentum_group;
    }

    // Indices for subgroup of particles
    //for(auto it_name = vec_ptl_target_label_.begin(),
    //         it_name_end = vec_ptl_target_label_.end();
    //        it_name != it_name_end;
    //        ++it_name
    //){
    //    vec_ptl_target_component_momenta_.push_back(map_ptl_component_effective_momentum_.at(*it_name));
    //}

}

////////////////////////////////////////////////////////
void OptiMass::ProcessTree::CalcMissingET(){
    //TVector2 buf_vec(0,0);
    // Check missing ET calculation automatically or not.
    for(int i = 0, size = vec_ptl_invisible_subsystem_calc_auto_.size(); i < size; ++i ){

        if(vec_ptl_invisible_subsystem_calc_auto_[i]){
            vec_ptl_invisible_subsystem_momenta_[i].Set(0.,0.);
            for(auto it_ptl_label = vec_ptl_invisible_subsystem_label_[i].begin(),
                     it_ptl_label_end = vec_ptl_invisible_subsystem_label_[i].end();
                    it_ptl_label != it_ptl_label_end;
                    ++it_ptl_label
            ) {
                for(auto it_vis_index = map_ptl_component_index_[*it_ptl_label].visible.begin(),
                         it_vis_index_end = map_ptl_component_index_[*it_ptl_label].visible.end();
                        it_vis_index != it_vis_index_end;
                        ++it_vis_index
                ) {
                    //buf_vec.Set(parMomentum.at(vec_ptl_visible_label_[*it_vis_index]+"_x"),parMomentum.at(vec_ptl_visible_label_[*it_vis_index]+"_y"));
                    //vec_ptl_invisible_subsystem_momenta_[i] -= buf_vec;
                    vec_ptl_invisible_subsystem_momenta_[i].Set(
                        vec_ptl_invisible_subsystem_momenta_[i].X() - parMomentum.at(vec_ptl_visible_label_[*it_vis_index]+"_x"),
                        vec_ptl_invisible_subsystem_momenta_[i].Y() - parMomentum.at(vec_ptl_visible_label_[*it_vis_index]+"_y")
                    );
                }
            }
        }
    }
}
void OptiMass::ProcessTree::CalcMissingETEffective(){
    // Copy Missing ET vector
    vec_ptl_invisible_subsystem_momenta_eff_ = vec_ptl_invisible_subsystem_momenta_;
    // If we need to update optimize invisible momenta
    if(!vec_ptl_optimize_not_final_label_.empty()){
        // Loop through optimization targets
        for(auto it_name = vec_ptl_optimize_label_.begin(),
                 it_name_end = vec_ptl_optimize_label_.end();
                it_name != it_name_end;
                ++it_name
        ) {
            // Loop through each invisible subsystems
            int i = 0;
            for(auto it_ptl_invisible_subsystem_label = vec_ptl_invisible_subsystem_label_.begin(),
                     it_ptl_invisible_subsystem_label_end = vec_ptl_invisible_subsystem_label_.end();
                    it_ptl_invisible_subsystem_label != it_ptl_invisible_subsystem_label_end;
                    ++it_ptl_invisible_subsystem_label
            ) {
                // Loop through particles in one invisible subsystem
                for(auto it_ptl_label = it_ptl_invisible_subsystem_label->begin(),
                         it_ptl_label_end = it_ptl_invisible_subsystem_label->end();
                        it_ptl_label != it_ptl_label_end;
                        ++it_ptl_label
                ) {
                    // if optimization target contains an invisible final state,
                    //std::cout << *it_name << " " << map_ptl_component_index_[*it_name].invisible[0] << " | " << *it_ptl_label << " " << map_ptl_component_index_[*it_ptl_label].invisible[0] << std::endl;
                    //if( map_ptl_component_index_[*it_name].invisible[0] == map_ptl_component_index_[*it_ptl_label].invisible[0] ) {
                    if( std::find(map_ptl_component_index_[*it_ptl_label].invisible.begin(),map_ptl_component_index_[*it_ptl_label].invisible.end(), map_ptl_component_index_[*it_name].invisible[0]) != map_ptl_component_index_[*it_ptl_label].invisible.end()  ) {
                        // Loop through visible final states in optimization target 
                        for(auto it_vis_index = map_ptl_component_index_[*it_name].visible.begin(),
                                 it_vis_index_end = map_ptl_component_index_[*it_name].visible.end();
                                it_vis_index != it_vis_index_end;
                                ++it_vis_index
                        ) {
                            // Add visible momenta found to effective MET
                            vec_ptl_invisible_subsystem_momenta_eff_[i].Set(
                                vec_ptl_invisible_subsystem_momenta_eff_[i].X() + parMomentum.at(vec_ptl_visible_label_[*it_vis_index]+"_x"),
                                vec_ptl_invisible_subsystem_momenta_eff_[i].Y() + parMomentum.at(vec_ptl_visible_label_[*it_vis_index]+"_y")
                            );
                            //std::cout << vec_ptl_invisible_subsystem_momenta_eff_[i].X() << " " << vec_ptl_invisible_subsystem_momenta_eff_[i].Y() << std::endl;
                        }
                    }
                }
                ++i;
            }
        }
    }
}

void OptiMass::ProcessTree::UpdateInvisibleMomenta(const std::vector<double>& par){

    for(unsigned int i = 0, size = buf_vec_ptl_last_momenta_.size(); i< size; ++i){
        buf_vec_ptl_last_momenta_[i].Set(vec_ptl_invisible_subsystem_momenta_eff_[i]);
    }
    //buf_vec_ptl_last_momenta_.assign(vec_ptl_invisible_subsystem_momenta_eff_.begin(),vec_ptl_invisible_subsystem_momenta_eff_.end());

    auto it_ptl_last_momenta = buf_vec_ptl_last_momenta_.begin();

    auto itPar = par.begin();
    std::vector< TLorentzVector >::iterator itOptimize = vec_ptl_optimize_momenta_.begin();
    unsigned int offset = 0;
    //int index = 0;
    for( auto it_ptl_invisible_subsystem_index = vec_ptl_invisible_subsystem_index_.begin(),
              it_ptl_invisible_subsystem_index_end = vec_ptl_invisible_subsystem_index_.end();
            it_ptl_invisible_subsystem_index != it_ptl_invisible_subsystem_index_end;
            ++it_ptl_invisible_subsystem_index ) {
        auto it_last = it_ptl_invisible_subsystem_index->end()-1;
        for(
            auto it = it_ptl_invisible_subsystem_index->begin();
            it != it_last;
            ++it
        ){
            offset = ( (*it) * 4);
            itPar += offset;

            it_ptl_last_momenta->Set(
                it_ptl_last_momenta->X() - *( itPar  ), 
                it_ptl_last_momenta->Y() - *( itPar + 1)
            );
            (itOptimize+*(it))->SetXYZM(
                *(itPar ),     //X
                *(itPar + 1),  //Y
                *(itPar + 2),  //Z
                *(itPar + 3)   //M
            );
            //std::cout << *it << " | " << index << std::endl;
            //std::cout << (itOptimize+*(it))->Px() << " | " << (itOptimize+*(it))->Py() << " | " << index  << std::endl;
            itPar -= offset;
        }
        offset = ( (*it_last) * 4);
        itPar += offset;
        (itOptimize+*(it_last))->SetXYZM(
            it_ptl_last_momenta->X(), //X
            it_ptl_last_momenta->Y(), //Y
            *(itPar + 2),                     //Z
            *(itPar + 3)                      //M
        );
        //std::cout << *it << " | " << index << " | last" << std::endl;
        //std::cout << (itOptimize+*(it))->Px() << " | " << (itOptimize+*(it))->Py() << " | " << index  << std::endl;
        itPar -= offset;
        //++index;
        ++it_ptl_last_momenta;
    }

}

void OptiMass::ProcessTree::RefreshMomentum(){
    if(!(vec_ptl_optimize_not_final_label_.empty())){
        for(auto it_name = vec_ptl_optimize_not_final_label_.begin(),
                 it_name_end = vec_ptl_optimize_not_final_label_.end();
                it_name != it_name_end;
                ++it_name
        ){
            auto& ptl_component_momentum = map_ptl_component_momentum_[*it_name];
            TLorentzVector& buf = *(ptl_component_momentum.invisible.at(0));
            buf = *(map_ptl_component_effective_momentum_[*it_name].optimize.at(0));

            auto& ptl_visible_component_momentum = ptl_component_momentum.visible;
            if(!ptl_visible_component_momentum.empty() ){
                for(auto it_vis_momentum = ptl_visible_component_momentum.begin(),
                         it_vis_momentum_end = ptl_visible_component_momentum.end();
                        it_vis_momentum != it_vis_momentum_end;
                        ++it_vis_momentum
                ){
                    buf -= *(*it_vis_momentum);
                }
            }
        }
    }
}

void OptiMass::ProcessTree::InitializeVisibleMomenta(){
    // Mass scale initialization
    for(unsigned int i=0;i<vec_ptl_visible_param_name_.size();++i){
        vec_ptl_visible_momenta_[i].SetXYZM(parMomentum.at(vec_ptl_visible_param_name_[i][0]),parMomentum.at(vec_ptl_visible_param_name_[i][1]),parMomentum.at(vec_ptl_visible_param_name_[i][2]),parMomentum.at(vec_ptl_visible_param_name_[i][3]));
    }
}

double OptiMass::ProcessTree::GetEffectiveScale(){
    double output = 0;
    for(auto it_visible = vec_ptl_visible_momenta_.begin(),
             it_visible_end = vec_ptl_visible_momenta_.end();
            it_visible != it_visible_end;
            ++it_visible
    ){
         output += it_visible->E();
    }
    double MissPtX = vec_ptl_invisible_subsystem_momenta_[0].X();
    double MissPtY = vec_ptl_invisible_subsystem_momenta_[0].Y();
    output += sqrt(MissPtX*MissPtX+MissPtY*MissPtY);
    return 1. * output / 6.;
}
// ========================================================================
// class ParticleNode: Constructor and Destructor
// ========================================================================
// Construct a particle node with label taken
OptiMass::ParticleNode::ParticleNode(const std::string& labelInput, bool debugInput) : label(labelInput), debug(debugInput){
    if(debug)
        std::cout << "creating " <<  label << std::endl;
}
// Destruct this node and its sibling nodes
OptiMass::ParticleNode::~ParticleNode(){
    if(debug)
        std::cout << "deleting " <<  label << std::endl;
    for(auto it = siblingPtls.begin(); it != siblingPtls.end(); ++it){
        delete *it;
    }
}

// ========================================================================
// class ParticleNode: Member Functions
// ========================================================================
// Add a sibling node to this particle node
void OptiMass::ParticleNode::AddSibling( ParticleNode* ptlNode ){
    siblingPtls.push_back(ptlNode);
}
// Get final state of this particle state
std::vector< std::string > OptiMass::ParticleNode::getFinalStates(){
    std::vector< std::string > output;
    if(siblingPtls.empty()){
        output.push_back(label);
        return output;
    }else{
        for(auto it = siblingPtls.begin(); it != siblingPtls.end(); ++it){
            std::vector< std::string > vec_final_states = (*it)->getFinalStates();
            for(auto it_str = vec_final_states.begin(),
                     it_str_end = vec_final_states.end();
                    it_str != it_str_end;
                    ++it_str
            ){
                output.push_back(*it_str);
            }
        }
        return output;
    }
}
// Get particle label
std::string OptiMass::ParticleNode::GetLabel(){
    return label;
}
// Get string which represents tree structure from this node
std::string OptiMass::ParticleNode::toString(){
    std::string output;
    output += label + "\n";
    for(auto it = siblingPtls.begin(); it != siblingPtls.end(); ++it){
        std::string s((*it)->toString());
        std::vector< std::string > str_lines = split(s,'\n');
        for(auto it_line = str_lines.begin(),
                 it_line_end = str_lines.end();
                it_line != it_line_end;
                ++it_line
        ){
            output += "|" + *it_line + "\n";
        }
    }
    return output;
}



// ========================================================================
// ========================================================================

bool OptiMass::ptlIn(ParticleSubelements<int> vec_visible_index, ParticleSubelements<int> vec_invisible_index){
	bool res = true;
	for(auto it_invisible_index = vec_invisible_index.visible.begin(),
             it_invisible_index_end = vec_invisible_index.visible.end();
            it_invisible_index != it_invisible_index_end;
            ++it_invisible_index
    ){
		bool buf = false;
		for(auto it_visible_index = vec_visible_index.visible.begin(),
                 it_visible_index_end = vec_visible_index.visible.end();
            it_visible_index != it_visible_index_end;
            ++it_visible_index
        ){
    	    if(*it_visible_index == *it_invisible_index)
				buf = true;
		}
		res = res && buf;
	}
	for(auto it_invisible_index = vec_invisible_index.invisible.begin(),
             it_invisible_index_end = vec_invisible_index.invisible.end();
            it_invisible_index != it_invisible_index_end;
            ++it_invisible_index
    ){
		bool buf = false;
		for(auto it_visible_index = vec_visible_index.invisible.begin(),
                 it_visible_index_end = vec_visible_index.invisible.end();
            it_visible_index != it_visible_index_end;
            ++it_visible_index
        ){
      	    if(*it_visible_index == *it_invisible_index)
				buf = true;
		}
		res = res && buf;
	}
	return res;
}

OptiMass::ParticleSubelements<int> OptiMass::ptlRemove(ParticleSubelements<int> vec_visible_index, ParticleSubelements<int> vec_invisible_index,int ptl_invisible_index){
    ParticleSubelements<int> output = vec_visible_index;
	if(ptlIn(vec_visible_index,vec_invisible_index)){
		output.visible = {} ;
		output.invisible = {} ;
		//if(vec_visible_index.size()<=2){
		//	output[2] = {};
		//}else{
		output.optimize = vec_visible_index.optimize;
		//}
	    for(auto it_visible_index = vec_visible_index.visible.begin(),
                 it_visible_index_end = vec_visible_index.visible.end();
            it_visible_index != it_visible_index_end;
            ++it_visible_index
        ){
			bool buf = true;
	        for(auto it_invisible_index = vec_invisible_index.visible.begin(),
                     it_invisible_index_end = vec_invisible_index.visible.end();
                    it_invisible_index != it_invisible_index_end;
                    ++it_invisible_index
            ){
      	        if(*it_visible_index == *it_invisible_index)
					buf = false;
			}
			if(buf){
				output.visible.push_back(*it_visible_index);
			}
		}
	    for(auto it_visible_index = vec_visible_index.invisible.begin(),
                 it_visible_index_end = vec_visible_index.invisible.end();
            it_visible_index != it_visible_index_end;
            ++it_visible_index
        ){
			bool buf = true;
	        for(auto it_invisible_index = vec_invisible_index.invisible.begin(),
                     it_invisible_index_end = vec_invisible_index.invisible.end();
                    it_invisible_index != it_invisible_index_end;
                    ++it_invisible_index
            ){
      	        if(*it_visible_index == *it_invisible_index)
					buf = false;
			}
		    if(buf){
				output.invisible.push_back(*it_visible_index);
			}
		}
		output.optimize.push_back(ptl_invisible_index);
	}else{
	}
	return output;
}

void OptiMass::printPtlDict(PtlIndexDict dict){
    for (auto it = dict.begin(),
              it_end = dict.end();
            it != it_end;
            ++it
    ){
    	std::cout << it->first << " is | ";
       	for (auto it2 = it->second.visible.begin(),
                  it2_end = it->second.visible.end();
                it2 != it2_end;
                ++it2
        ){
    		std::cout << *it2 << " ";
    	}
    	std::cout << "| ";
       	for (auto it2 = it->second.invisible.begin(),
                  it2_end = it->second.invisible.end();
                it2 != it2_end;
                ++it2
        ){
    		std::cout << *it2 << " ";
    	}
    	std::cout << "| ";
//    	if(it->second.size() > 2){
           	for (auto it2 = it->second.optimize.begin(),
                      it2_end = it->second.optimize.end();
                    it2 != it2_end;
                    ++it2
            ){
        		std::cout << *it2 << " ";
        	}
        	std::cout << "| ";
//    	}
    	std::cout << std::endl;
    }
}
