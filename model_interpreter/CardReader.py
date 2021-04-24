#!/usr/bin/env python2
import xml.etree.ElementTree as ET
from xml.etree.ElementTree import Element, Comment, SubElement, ElementTree, dump
import re
import os, sys

def indent(elem, level=0):
    """ Cleanup indentation in XML """
    i = "\n" + level*(" "*4)
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + (" "*4)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i



class CardReader():
    """ Class for handling input card """
    
    """ Shortcut to element group """
    elem_card = None
    elem_decay_chain = None
    elem_particle_target = None
    elem_particle_label = None
    elem_particle_property = None
    elem_particle_invisible_subsystem = None
    elem_constraints = None
    elem_alm = None
    elem_minuit = None
    elem_calc = None
    elem_tree = None
    elem_tree_effective = None

    postfix_rules = [
            [".M()","GetSubsystemMass","double"],
            [".M2()","GetSubsystemMassSquare","double"],
            [".Px()","GetSubsystemPx","double"],
            [".Py()","GetSubsystemPy","double"],
            [".Pz()","GetSubsystemPz","double"],
            [".E()","GetSubsystemE","double"],
            [".Pt()","GetSubsystemPt","double"],
            [".V()","GetSubsystemMomentum","TLorentzVector"]
        ]
    mass_function_rules = {
            "M" : "CalcMassFromIndice",
            "M2" : "CalcMassSquareFromIndice",
            "MT2" : "CalcMTSquareFromIndice"
        }
    group_function_rules = {
            "max" : "CalcMaximum" ,
            "mean" : "CalcMean" 
        }
    
    def __init__(self,filename):
        """ Constructor for CardReader """
        """ Creating XML document model from preexisting file """
        
        doc = ET.parse(filename)
        self.elem_card = doc.getroot()
        
        # class name as the file name of the input process card #rev
        self.__class_name = filename.split("/")[-1].strip(".xml")

        self.elem_decay_chain = self.elem_card.find("DecayChains")
        self.elem_particle_target = self.elem_card.find("ParticleMassFunction")
        self.elem_particle_label = self.elem_card.find("ParticleLabels")
        self.elem_particle_property = self.elem_card.find("ParticleProperties")
        self.elem_particle_invisible_subsystem = self.elem_card.find("ParticleInvisibleSubsystem")
        self.elem_constraints = self.elem_card.find("Constraints")
        self.elem_alm = self.elem_card.find("ALM")
        self.elem_minuit = self.elem_card.find("Minuit")
        self.elem_calc = self.elem_card.find("Calc")
        self.elem_tree = Element("Process")
        self.elem_tree_effective = Element("ProcessEffective")

        self.__interpret_particle_name()
        self.__interpret_process_syntax() #rev:0-1
        self.__interpret_particle_label()
        self.__generate_effective_tree()

        #parental map #rev:done
        self.parents_map = {c: p for p in self.elem_tree.getiterator() for c in p}
        self.parents_map_effective = {c: p for p in self.elem_tree_effective.getiterator() for c in p}

        #label list and set from elem_tree 
        self.labels_elem_tree = [ i.get('label') for i in self.elem_tree.getiterator() ]
        self.labels_elem_tree = [ x for x in self.labels_elem_tree if x is not None and x.strip() ]
        
        #label list and set from elem_particle_label
        self.labels_elem_particle_label = [ i.get('label') for i in self.elem_particle_label.getiterator() ]
        self.labels_elem_particle_label = [ x for x in self.labels_elem_particle_label if x is not None and x.strip() ]

        #Check consistency of Particle elem in 'elem_tree' and 'elem_particle_label'  #rev
        # 1) discriminative particle labels in elem_tree #syntax_error:1 #rev:check 
        if len(set(self.labels_elem_tree)) != len(self.labels_elem_tree):

            str_decay = [ x.text.strip() for x in self.elem_decay_chain.findall("DecayChain") ]
            message = '\n (!) Syntax Error (!)\n' \
                    + ' ==================== \n' \
                    + ' The DecayChains - {} of the process - [{}], has duplicated particle labels in the tree graph syntax.'.format(str_decay, self.get_class_name())
            sys.exit(message)

        # 2) discriminative particle labels in elem_particle_label #syntax_error:2 #rev:check
        if len(set(self.labels_elem_particle_label)) != len(self.labels_elem_particle_label):

            message = '\n (!) Syntax Error (!)\n' \
                    + ' ==================== \n' \
                    + ' The ParticleLabels - {} of the process - [{}], has duplicated particle labels in the tree'.format(self.labels_elem_particle_label, self.get_class_name())
            sys.exit(message)

        # 3) injective particle labels from 'labels_elem_tree' 
        #    to 'labels_elem_particle_label' #syntax_error:3 #rev:check
        if not set(self.labels_elem_tree).issubset(self.labels_elem_particle_label):

            str_decay = [ x.text.strip() for x in self.elem_decay_chain.findall("DecayChain") ]
            message = '\n (!) Syntax Error (!)\n' \
                    + ' ==================== \n' \
                    + ' The Particles in the ParticleLabels branch - {} of the process - [{}], has missing particle labels, not enough to cover the particles appeared in the DecayChains - {}.'.format(self.labels_elem_particle_label, self.get_class_name(), str_decay)
            sys.exit(message)

        # 4) Check (name, ...) between ParticleLabels and ParticleProperties
        #rev:TBD ...

        # 5) Check the registered particles in the ParticleInvisibleSubsystem are the root particles of the DecayChains #syntax_error:8        
        subsystems = []
        str_decay = [ x.text.strip() for x in self.elem_decay_chain.findall("DecayChain") ]
        for id_subsystem, elem_subsystem in enumerate(self.elem_particle_invisible_subsystem.findall("Subsystem")):

            ptl_subsystem = []

            for elem in elem_subsystem.findall("Particle"):

                label = elem.get("label")
                if self.get_parent(label)[0] is not None:
                    message = '\n (!) Syntax Error (!)\n' \
                            + ' ==================== \n' \
                            + ' The particle - (currently, {}) for {}-th invisible subsystem of the process - [{}], should be chosen as the one of the root particles representing each decay chain appeared in the DecayChains - {}.'.format(label, id_subsystem, self.get_class_name(), str_decay)
                    sys.exit(message)

                ptl_subsystem.append(elem.get("label"))

            subsystems.append(ptl_subsystem)

        # 6) Check if the total number of registered root particles in the ParticleInvisibleSubsystem = number total 'DecayChain'(s) in the 'DecayChains' #syntax_error:9
        n_root_ptls = self.get_nelem_of_nested_list(subsystems)
        if n_root_ptls != len(self.elem_decay_chain):
            message = '\n (!) Syntax Error (!)\n' \
                    + ' ==================== \n' \
                    + ' The root particles - {} consisting of the invisible subsystems of the process - [{}], should be matched one-to-one with the roots of the whole DecayChains defined as - {}.'.format(subsystems, self.get_class_name(), str_decay)
            sys.exit(message)


        pass # end of __init__


    def get_final_state_visibles_by_invisible_subsystem(self): #rev
        # rev) (invisible_subsystem) VS (indep_semi_inv_system) VS (indep_PT_conserv_system)

        vis_list_tot_system = []
        for root_list_a_subsystem in self.get_invisible_subsystem():
            vis_list_a_subsystem = []
            for root in root_list_a_subsystem:
                for ptl in self.get_siblings_list(root):
                    ptl_ref = self.get_ptl_ref_elem_tree(ptl)
                    if ptl_ref.find('Particle') == None and ptl_ref.get('invisible') != 'True':
                        vis_list_a_subsystem.append(ptl)
            vis_list_tot_system.append(vis_list_a_subsystem)

        return vis_list_tot_system


    def get_invisibles_by_invisible_subsystem(self): #rev
        # rev) (invisible_subsystem) VS (indep_semi_inv_system) VS (indep_PT_conserv_system)

        inv_list_tot_system = []
        for root_list_a_subsystem in self.get_invisible_subsystem():
            inv_list_a_subsystem = []
            for root in root_list_a_subsystem:
                for ptl in self.get_siblings_list(root):
                    ptl_ref = self.get_ptl_ref_elem_tree(ptl)
                    if ptl_ref.find('Particle') == None and ptl_ref.get('invisible') == 'True':
                        inv_list_a_subsystem.append(ptl)
            inv_list_tot_system.append(inv_list_a_subsystem)

        return inv_list_tot_system


    def get_ptls_by_invisible_subsystem(self): #rev:done

        ptl_list_tot_system = []
        for root_list_a_subsystem in self.get_invisible_subsystem():
            ptl_list_a_subsystem = []
            for root in root_list_a_subsystem:
                for ptl in self.get_siblings_list(root):
                    ptl_list_a_subsystem.append(ptl)
            ptl_list_tot_system.append(ptl_list_a_subsystem)
        
        return ptl_list_tot_system


    def get_nelem_of_nested_list(self, elem): #rev:new,done
        count = 0
        if isinstance(elem, list):
            for elem_sub in elem:
                count += self.get_nelem_of_nested_list(elem_sub)
        else:
            count += 1

        return count


    def get_ptl_ref_elem_tree(self, ptl_label):

        for elem in self.elem_tree.getiterator():
            if elem.get('label') == ptl_label:
                return elem
        else:
            raise ValueError(' => There is no particle element with the label {}. '.format(ptl_label))


    def get_parent(self, ptl_label): #rev

        ptl_ref = self.get_ptl_ref_elem_tree(ptl_label)
        parent_ref = self.parents_map[ptl_ref]
        parent_label = parent_ref.get('label')

        return parent_label, parent_ref 


    def trace_add_parents(self, parents_list=None): #rev
        """
        * Args
        1) parents_list: a list with at least one offstring particle label(s) from the self.elem_tree
        
        * Return: a list of parent particles of the oldest offstring particle (rightmost)   
        """

        oldest_label = parents_list[-1]
        oldest_ref = self.get_ptl_ref_elem_tree(oldest_label)
        #oldest_ref = self.elem_particle_label.findall('./Particle[@label="{}"]'.format(oldest_label))
        
        parent_ref = self.parents_map[oldest_ref]
        parent_label = parent_ref.get('label')

        if parent_label is not None:
            parents_list.append(parent_label)
            self.trace_add_parents(parents_list=parents_list)

        return parents_list


    def get_parents_list(self, ptl_label=None): #rev
        """
        * Args
        1) ptl_label: a particle label from which its parental particles are to be traced and added.

        * Return: a list of parental particle labels, such as [ptl_label, ..., highest ancester ptl label]
        """

        if ptl_label is None or ptl_label not in self.labels_elem_tree:
            message = ' => Not proper particle label input for the arg - "ptl_label", within the candidate set - {}.'.format(self.labels_elem_tree)
            sys.exit(message)

        parents_list = [ptl_label]
        self.trace_add_parents(parents_list=parents_list)

        return parents_list        


    def get_siblings_list(self, ptl_label=None): #rev
        """
        * Args
        1) ptl_label: a particle label from which its sibling particles are to be traced and added.

        * Return: a list of all siblings from the ptl_label, ~ [ptl_label, ... ]
        """
        
        if ptl_label not in self.labels_elem_tree:
            message = ' => The input particle label - {} is not in the elements of decay chains - {}.'.format(ptl_label, self.labels_elem_tree)
            sys.exit(message)

        ptl_ref = self.get_ptl_ref_elem_tree(ptl_label)
        ptl_list = []
        for elem in ptl_ref.getiterator():
            ptl_list.append(elem.get('label'))

        return ptl_list


    #def str_split_keeping_parenthesis(arg): #rev:previous one
    def str_split_into_effective_branches(self, arg): #rev:done
        """
        Given the arg in the syntax for a 'DecayChain' of .xml process file,
        this function returns the list of decaying branch syntax string, from the highest 2 levels
            Lv1) initial mother ptl decay to N daughter ptls 
            Lv2) N daughter ptl decay processes, each of which includes its own sub-decay process.
        In total, len(output list) = 1+N. 
        """
        # Stripping out the highest decay branch
        arg.strip('() ')

        # Check the consistency of DecayBranch syntax: 
        # 1) seperators : N(',')=N('(')=N(')') #syntax_error:4 #rev:check 
        if not ( (arg.count(',') == arg.count('(')) \
                and (arg.count('(') == arg.count(')')) ):
            
            message =  '\n  (!) Syntax Error (!)\n' \
              + '  ==================== \n' \
              + ' : The DecayChain [{}] of the process [{}] is not properly using the separaters -[",", "(", ")"] for describing a tree graph.\n'.format(arg, self.get_class_name()) \
              + '\n   In particular,\n' \
              + '   1) All of the (sibling/subsequent) decay branches should be separated by [",":comma].\n' \
              + '   2) Except for the highest one, all of subsequent decay branches should be started with "(", and closed with adding ")" right after the scripts for their potential sub-sub-...-sequantial decay branches inside, recursively.\n' \
              + '   3) All in all, N[","] = N["("] = N[")"] \n\n'
            sys.exit(message)
       
        # Init. 
        level = 0
        output = []
        buf = ""

        for char in arg:
            if level < 0:
                print(' (!) Negative Parenthesis Level Error (!) ')
                return None        
            elif (char != '(' and char != ')' and char != ',' ) or (char == ',' and level > 0):
                buf += char
            
            elif char == ',' and level == 0:
            #    print("===============")
                output.append(buf.strip())
                buf = ""

            elif char == '(' and level == 0:
                level += 1
            elif char == '(' and level > 0:
                buf += char
                level += 1
            elif char == ')' and level == 1:
                level -= 1
            elif char == ')' and level > 1:
                buf += char
                level -= 1
        output.append(buf.strip())

        # Check 2)
        #   A child decaying subsequently, should be appeared once in advance, as a child of parent. #syntax_error:5
        if len(output)>1:
            branch_mom = output[0]
            cs = branch_mom.split('-')[1].split()
            branch_cs = output[1:]
            cs_branch = [ str_b.split()[0] for str_b in branch_cs ]
            if not set(cs_branch).issubset(set(cs)):
                message = '\n (!) Syntax Error (!)\n' \
                        + ' ==================== \n' \
                        + ' : The DecayChain ("{}") of the process ("{}") includes a sub-branch whose origin cannot be specified in its higher level decay. '.format(arg, self.get_class_name())
                sys.exit(message)

        return output
    
    def __interpret_particle_name(self):
        for ptl_label in self.elem_particle_label.findall("Particle"):
            for ptl_name in self.elem_particle_property.findall("Particle"):
                if ptl_label.get("name") != None and ptl_label.get("name") == ptl_name.get("name"):
                    ptl_label.attrib.update(ptl_name.attrib)

    def __interpret_process_syntax(self):
        for syntax in self.get_process_syntax():
            self.elem_tree.append(self.__process_parser(syntax)) #rev:0

    def __interpret_particle_label(self):
        for elem_ptl in self.elem_tree.iter("Particle"):
            for ptl_label in self.elem_particle_label.findall("Particle"):
                if ptl_label.get("label") != None and ptl_label.get("label") == elem_ptl.get("label"):
                    elem_ptl.attrib.update(ptl_label.attrib)

    def __single_process_parser(self, arg):
        """ Return Particle element, parsing arg into decay chain tree representation in XML
            arg should be single process syntax."""

        procPtl = arg.split('-')
        initPtl = procPtl.pop(0).strip(' ')
    
        # Check process syntax: 
        #   N of initial decaying particle should be one for each single process #syntax_error:6 #rev:check
        if len(initPtl.split()) > 1:
            message = '\n (!) Syntax Error (!)\n' \
                    + ' ==================== \n' \
                    + ' : The decay branch ({})) in the process ("{}") has multiple initial particles (should be one). '.format(arg, self.get_class_name())
            sys.exit(message)
                

        motherPtlElem = Element('Particle')
        motherPtlElem.set("label",initPtl)

        #siblingPtls = procPtl.pop(0).split(' ') #rev:done
        siblingPtls = procPtl[0].split()

        # Check process syntax: 
        #   N of final state particle should be >=2 for each single process #syntax_error:7 #rev:check
        if len(siblingPtls) < 2:
            message = '\n (!) Syntax Error (!)\n' \
                    + ' ==================== \n' \
                    + ' : The decay branch ({})) in the process ("{}") has only one final state particle (should be >=2). '.format(arg, self.get_class_name())
            sys.exit(message)
                

        for procPtlBuf in siblingPtls:
            SubElement(motherPtlElem,'Particle').attrib["label"] = procPtlBuf
            #a = SubElement(motherPtlElem,'Particle')
            #a.attrib["label"] = procPtlBuf
            #if(procPtlBuf != ''):
            #    SubElement(motherPtlElem,'Particle').attrib["label"] = procPtlBuf
        
        return motherPtlElem

    def __process_parser(self, arg):
        """ Return Particle element, parsing arg into decay chain tree representation in XML"""
        # translate arg into 
        buf = self.str_split_into_effective_branches(arg)
        
        # if arg is single process, pass arg to singleProcessParser
        if(buf[0].strip('() ') == arg.strip('() ')):
            return self.__single_process_parser(arg.strip('() '))
        # if 
        else:
            motherPtlElem = self.__single_process_parser(buf.pop(0))

            subPtlElems = map(lambda st:self.__process_parser(st),buf)

            for elem_siblings in motherPtlElem.findall("Particle"):
                for elem_siblings_detail in subPtlElems:
                    if elem_siblings.get("label") != None and elem_siblings.get("label") == elem_siblings_detail.get("label"):
                        elem_siblings.extend(elem_siblings_detail.findall("Particle"))
                        sib_ext = [ sib.get('label') for sib in elem_siblings_detail.findall("Particle")]
            return motherPtlElem

    def __generate_effective_tree(self):        
        for syntax in self.get_process_syntax():
            self.elem_tree_effective.append(self.__process_parser(syntax))
        for elem_ptl in self.elem_tree_effective.iter("Particle"):
            for ptl_label in self.elem_particle_label.findall("Particle"):
                if ptl_label.get("label") != None and ptl_label.get("label") == elem_ptl.get("label"):
                    elem_ptl.attrib.update(ptl_label.attrib)
        for elem in self.elem_tree_effective.iter("Particle"):
            if elem.get("optimize_target") == "True":
                invisible_found = False
                for elem_ptl_sub in elem.iter("Particle"):
                    if elem_ptl_sub.find("Particle") == None and elem_ptl_sub.get("invisible") == "True" and not invisible_found:
                        invisible_found = True
                        print ("On effective node: " + elem.get("label") + ", invisible: " + elem_ptl_sub.get("label") + " detected. Removing subelements of "  + elem.get("label") + " on effective tree")
                    elif elem_ptl_sub.find("Particle") == None and elem_ptl_sub.get("invisible") == "True" and invisible_found:
                        print ("Error! " + elem_ptl.get("label") + " optimize target with more than two invisible particle selected.")
                if invisible_found:
                    elem.attrib["invisible"] = "True"
                    for elem_sub in elem.findall("Particle"):
                        elem.remove(elem_sub)

    def __interpret_mass_elem(self,elem_input):
        output = ""
        if elem_input.tag == "Particle":        
            output = "new MassFunctionParticle( dynamic_cast<ProcessTree*>(&process_tree_), \"" + elem_input.get("label") +"\", &" + self.mass_function_rules[elem_input.get("mass_function")] + " )"
        elif elem_input.tag == "ParticleGroup":
            elem_iterator = elem_input.findall("Particle")
            elem_iterator.extend(elem_input.findall("ParticleGroup"))
            output = "new MassFunctionGroup( dynamic_cast<ProcessTree*>(&process_tree_) , &" + self.group_function_rules[elem_input.get("group_function")] + " , {\n"
            for elem in elem_iterator:
                if elem.get("mass_function") == None:
                    elem.attrib["mass_function"] = elem_input.get("mass_function")
                output += "    dynamic_cast<MassFunctionInterface*>( " + self.__interpret_mass_elem(elem).replace("\n","\n    ") + " ),\n"
            output = output.rstrip(",\n");
            output += "\n} )"
        return output

    def get_class_name(self):
        #return self.elem_card.attrib["classname"] #rev
        return self.__class_name

    def get_description(self):
        return self.elm_card.attrib["description"]

    def get_process_syntax(self):
        output = []
        for elem in self.elem_decay_chain.findall("DecayChain"):
            output.append(elem.text.strip())
        return output

    def get_invisibles(self):
        output = []
        for ptl in self.elem_particle_label.findall("Particle"):
            if(ptl.get('invisible') == 'True'):
                output.append(ptl.attrib["label"])
        return output

    def get_target_particles(self):
        output = []
#        for ptl in self.elem_particle_target.findall("Particle"):
        for ptl in self.elem_particle_target.iter("Particle"):
            output.append(ptl.attrib["label"])
        return output

    def get_optimize_particles(self):
        output = []
        for elem_ptl in self.elem_tree.iter("Particle"):
            if(elem_ptl.find("Particle") != None and elem_ptl.get("optimize_target") == "True"):
                invisible_found = False
                for elem_ptl_sub in elem_ptl.iter("Particle"):
                    if elem_ptl_sub.find("Particle") == None and elem_ptl_sub.get("invisible") == "True" and not invisible_found:
                        output.append(elem_ptl.get("label"))
                        elem_ptl_sub.attrib["disabled"] = "True"
                        invisible_found = True
                    elif elem_ptl_sub.find("Particle") == None and elem_ptl_sub.get("invisible") == "True" and invisible_found:
                        print ("Error! " + elem_ptl.get("label") + " optimize target with more than two invisible particle selected.")
            elif elem_ptl.find("Particle") == None and elem_ptl.get("disabled") != "True":
                for elem_label_info in self.elem_particle_label.findall("Particle"):
                    if elem_ptl.get("label") ==  elem_label_info.get("label") and elem_label_info.get('invisible')=="True":
                        output.append(elem_ptl.get("label"))
        return output

    def get_invisible_subsystem(self):
        output = []
        for elem_subsystem in self.elem_particle_invisible_subsystem.findall("Subsystem"):
            output_subsystem = []
            for elem_input in elem_subsystem.findall("Particle"):
                if elem_input.get("label") != None:
                    output_subsystem.append(elem_input.get("label"))
        #        for elem_on_tree in self.elem_tree_effective.findall("Particle"):
        #            if elem_on_tree.get("label") == elem_input.get("label"):
        #                for elem_final in elem_on_tree.iter("Particle"):
        #                    if elem_final.find("Particle") == None and elem_final.get("invisible") == "True":
        #                        print (elem_final.get("label"))
        #                        if not (elem_final.get("label") in output_subsystem):
        #                            output_subsystem.append(elem_final.get("label"))
        #    print (str(output_subsystem))
            output.append(output_subsystem)
        return output
    
    def get_mass_dict(self):
        output = {}        
        for elem in self.elem_particle_label.findall("Particle"):
            if elem.get("mass") != None:
                output[elem.attrib["label"]] = elem.attrib["mass"]
            else:
                output[elem.attrib["label"]] = "0"
        return output

    def get_alm_dict(self):
        output = {}        
        for elem in self.elem_alm.findall("param"):
            output[elem.attrib["name"]] = elem.attrib["value"]
        return output

    def get_minuit_dict(self):
        output = {}        
        for elem in self.elem_minuit.findall("param"):
            output[elem.attrib["name"]] = elem.attrib["value"]
        return output

    def get_mass_function(self):
        output = "mass_interface_ = dynamic_cast<MassFunctionInterface*>( " + self.__interpret_mass_elem(self.elem_particle_target.find("ParticleGroup")) + "\n);\n"
        return output

    def get_constraints_syntax(self):
        output = []
        for elem in self.elem_constraints.findall("Constraint"):
            output.append(elem.text.strip())
        return output

    def get_constraints_params(self):
        output = set([])
        for elem in self.elem_constraints.findall("Constraint"):
            for rule in self.postfix_rules:
                for cmd in re.findall('([a-zA-Z0-9]+' + rule[0].replace('.','\.').replace('(','\(').replace(')','\)') + ')',elem.text.strip()):
                    output.add(rule[2] + ' ' + cmd.replace('.','_').rstrip('()') + ' = process_tree_.' + rule[1] + '(\"' + cmd.replace(rule[0],'') + '\")' )
        return list(output)

    def get_constraints(self):
        output = []
        for elem in self.elem_constraints.findall("Constraint"):
            res = elem.text.strip()
            for rule in self.postfix_rules:
                res = re.sub('([a-zA-Z0-9]+)' + rule[0].replace('.','\.').replace('(','\(').replace(')','\)'), '\\1' + rule[0].replace('.','_').rstrip('()') , res)
                #res = re.sub('([a-zA-Z0-9]+)' + rule[0].replace('.','\.').replace('(','\(').replace(')','\)'), rule[1] + '(\"\\1\")', res)
            output.append(res)
        return output

    def print_code(self):
        for elem in self.elem_card.iter("code"):
            print (elem.text.strip())

    def print_tree(self):
        indent(self.elem_tree)
        dump(self.elem_tree)

    def print_tree_effective(self):
        indent(self.elem_tree_effective)
        dump(self.elem_tree_effective)
    
    def dump(self):
        dump(self.elem_card)


# rev
    def get_final_state_visibles(self):
        output = []
        for elem_ptl in self.elem_tree.iter("Particle"):
            if elem_ptl.find("Particle") == None and (elem_ptl.attrib["label"] not in self.get_invisibles()) :
#            if elem_ptl.find("Particle") == None :
                output.append(elem_ptl.attrib["label"])
        return output

# rev
    def get_node_particles(self):
        output = []
        for elem_ptl in self.elem_tree.iter("Particle"):
            if elem_ptl.find("Particle") != None :
                output.append(elem_ptl.attrib["label"])
        return output

# rev
    def get_reconstructed_particles(self):
        output = []
        for elem_ptl in self.elem_tree.iter("Particle"):
            if not (elem_ptl.find("Particle") == None and (elem_ptl.attrib["label"] not in self.get_invisibles()) ) :
                output.append(elem_ptl.attrib["label"])
        return output
        
    def get_constraints_formula(self):
        output = []
        for elem in self.elem_constraints.findall("Constraint"):
            res = elem.text.strip()
            output.append(res)
        return output


