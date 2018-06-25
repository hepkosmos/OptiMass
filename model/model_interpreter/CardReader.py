import xml.etree.ElementTree as ET
from xml.etree.ElementTree import Element, Comment, SubElement, ElementTree, dump
import re

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

def str_split_keeping_parenthesis(arg):
    """ Returns string list, splitted by comma while maintaing lowest parenthesis block. """
    if( arg.count('(') != arg.count(')') ) :
        print('Parenthesis does not match')
        return None
    #stripping needed
    arg.strip('() ')
    level = 0
    output = []
    buf = ""
    for char in arg:
        #print(char + " received at level " + str(level) + " "+ str(char == ',') + str(char == ',' and level == 0))
        if level < 0:
            print('negative parenthesis level error')
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
    return output

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
        
#        print("reading card from input file")
        doc = ET.parse(filename)
        self.elem_card = doc.getroot()
        
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
        self.__interpret_process_syntax()
        self.__interpret_particle_label()
        self.__generate_effective_tree()
    
    def __interpret_particle_name(self):
        for ptl_label in self.elem_particle_label.findall("Particle"):
            for ptl_name in self.elem_particle_property.findall("Particle"):
                if ptl_label.get("name") != None and ptl_label.get("name") == ptl_name.get("name"):
                    ptl_label.attrib.update(ptl_name.attrib)

    def __interpret_process_syntax(self):
        for syntax in self.get_process_syntax():
            self.elem_tree.append(self.__process_parser(syntax))

    def __interpret_particle_label(self):
        for elem_ptl in self.elem_tree.iter("Particle"):
            for ptl_label in self.elem_particle_label.findall("Particle"):
                if ptl_label.get("label") != None and ptl_label.get("label") == elem_ptl.get("label"):
                    elem_ptl.attrib.update(ptl_label.attrib)

    def __single_process_parser(self, arg):
        """ Return Particle element, parsing arg into decay chain tree representation in XML
            arg should be single process syntax."""
        procPtl = arg.split('-')
        initPtl = procPtl.pop(0).split(' ').pop(0)
        siblingPtls = procPtl.pop(0).split(' ')
        motherPtlElem = Element('Particle')
        motherPtlElem.set("label",initPtl)
        for procPtlBuf in siblingPtls:
            if(procPtlBuf != ''):
                SubElement(motherPtlElem,'Particle').attrib["label"] = procPtlBuf
        return motherPtlElem

    def __process_parser(self, arg):
        """ Return Particle element, parsing arg into decay chain tree representation in XML"""
        # translate arg into 
        buf = str_split_keeping_parenthesis(arg)
        
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
                        print "On effective node: " + elem.get("label") + ", invisible: " + elem_ptl_sub.get("label") + " detected. Removing subelements of "  + elem.get("label") + " on effective tree"
                    elif elem_ptl_sub.find("Particle") == None and elem_ptl_sub.get("invisible") == "True" and invisible_found:
                        print "Error! " + elem_ptl.get("label") + " optimize target with more than two invisible particle selected."
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
        return self.elem_card.attrib["classname"]

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
            #print(elem_ptl.get("optimize_target"))
            if(elem_ptl.find("Particle") != None and elem_ptl.get("optimize_target") == "True"):
                invisible_found = False
                for elem_ptl_sub in elem_ptl.iter("Particle"):
                    if elem_ptl_sub.find("Particle") == None and elem_ptl_sub.get("invisible") == "True" and not invisible_found:
                        output.append(elem_ptl.get("label"))
                        elem_ptl_sub.attrib["disabled"] = "True"
                        invisible_found = True
                    elif elem_ptl_sub.find("Particle") == None and elem_ptl_sub.get("invisible") == "True" and invisible_found:
                        print "Error! " + elem_ptl.get("label") + " optimize target with more than two invisible particle selected."
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
        #            print "Afsda : " + elem_input.get("label")
                    output_subsystem.append(elem_input.get("label"))
        #        for elem_on_tree in self.elem_tree_effective.findall("Particle"):
        #            if elem_on_tree.get("label") == elem_input.get("label"):
        #                for elem_final in elem_on_tree.iter("Particle"):
        #                    if elem_final.find("Particle") == None and elem_final.get("invisible") == "True":
        #                        print elem_final.get("label")
        #                        if not (elem_final.get("label") in output_subsystem):
        #                            output_subsystem.append(elem_final.get("label"))
        #    print str(output_subsystem)
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
            print elem.text.strip()
    def print_tree(self):
        indent(self.elem_tree)
        dump(self.elem_tree)
    def print_tree_effective(self):
        indent(self.elem_tree_effective)
        dump(self.elem_tree_effective)
    
    def dump(self):
        dump(self.elem_card)


# (!Revision-NEW!)
    def get_final_state_visibles(self):
        output = []
        for elem_ptl in self.elem_tree.iter("Particle"):
            if elem_ptl.find("Particle") == None and (elem_ptl.attrib["label"] not in self.get_invisibles()) :
#            if elem_ptl.find("Particle") == None :
                output.append(elem_ptl.attrib["label"])
        return output

# (!Revision-NEW!)
    def get_node_particles(self):
        output = []
        for elem_ptl in self.elem_tree.iter("Particle"):
            if elem_ptl.find("Particle") != None :
                output.append(elem_ptl.attrib["label"])
        return output

# (!Revision-NEW!)
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
