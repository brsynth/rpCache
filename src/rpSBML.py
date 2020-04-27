import libsbml
from hashlib import md5
import os
import logging
import copy

## @package RetroPath SBML writer
# Documentation for SBML representation of the different model
#
# To exchange between the different workflow nodes, the SBML (XML) format is used. This
# implies using the libSBML library to create the standard definitions of species, reactions, etc...
# Here we also define our own annotations that are used internally in that we call BRSYNTH nodes.
# The object holds an SBML object and a series of methods to write and access BRSYNTH related annotations


##################################################################
############################### rpSBML ###########################
##################################################################


## libSBML reader for RetroPath
# Converts an SBML object (or file) into the internal format
#
class rpSBML:
    ## Constructor
    #
    # @param model libSBML model object
    # @param docModel libSBML Document object
    # @param nameSpaceModel libSBML name space (not required)
    def __init__(self, modelName, document=None, path=None):
        self.logger = logging.getLogger(__name__)
        self.logger.info('Started instance of rpSBML')
        self.modelName = modelName
        self.document = document
        if self.document==None:
            self.model = None
        else:
            self.model = self.document.getModel()
        self.path = path
        self.miriam_header = {'compartment': {'go': 'go/GO:', 'mnx': 'metanetx.compartment/', 'bigg': 'bigg.compartment/', 'seed': 'seed/', 'name': 'name/'}, 'reaction': {'mnx': 'metanetx.reaction/', 'rhea': 'rhea/', 'reactome': 'reactome/', 'bigg': 'bigg.reaction/', 'sabiork': 'sabiork.reaction/', 'ec': 'ec-code/', 'biocyc': 'biocyc/', 'lipidmaps': 'lipidmaps/', 'uniprot': 'uniprot/'}, 'species': {'mnx': 'metanetx.chemical/', 'chebi': 'chebi/CHEBI:', 'bigg': 'bigg.metabolite/', 'hmdb': 'hmdb/', 'kegg_c': 'kegg.compound/', 'kegg_d': 'kegg.drug/', 'biocyc': 'biocyc/META:', 'seed': 'seed.compound/', 'metacyc': 'metacyc.compound/', 'sabiork': 'sabiork.compound/', 'reactome': 'reactome/R-ALL-'}}
        self.header_miriam = {'compartment': {'go': 'go', 'metanetx.compartment': 'mnx', 'bigg.compartment': 'bigg', 'seed': 'seed', 'name': 'name'}, 'reaction': {'metanetx.reaction': 'mnx', 'rhea': 'rhea', 'reactome': 'reactome', 'bigg.reaction': 'bigg', 'sabiork.reaction': 'sabiork', 'ec-code': 'ec', 'biocyc': 'biocyc', 'lipidmaps': 'lipidmaps', 'uniprot': 'uniprot'}, 'species': {'metanetx.chemical': 'mnx', 'chebi': 'chebi', 'bigg.metabolite': 'bigg', 'hmdb': 'hmdb', 'kegg.compound': 'kegg_c', 'kegg.drug': 'kegg_d', 'biocyc': 'biocyc', 'seed.compound': 'seed', 'metacyc.compound': 'metacyc', 'sabiork.compound': 'sabiork', 'reactome': 'reactome'}}

    #######################################################################
    ############################# PRIVATE FUNCTIONS ####################### 
    #######################################################################

    '''
    # define Python user-defined exceptions
    class Error(Exception):
        """Base class for other exceptions"""
        pass

    #We do this so that it is not stopping
    class libSBMLError(Error):
        pass
    '''

    ## Check the libSBML calls
    #
    # Check that the libSBML python calls do not return error INT and if so, display the error. Taken from: http://sbml.org/Software/libSBML/docs/python-api/create_simple_model_8py-example.html
    #
    # @param value The SBML call
    # @param message The string that describes the call
    def _checklibSBML(self, value, message):
        if value is None:
            self.logger.error('LibSBML returned a null value trying to ' + message + '.')
            #raise SystemExit('LibSBML returned a null value trying to ' + message + '.')
            #raise libSBMLError('LibSBML returned a null value trying to ' + message + '.')
            #raise self.libSBMLError
            raise AttributeError
        elif type(value) is int:
            if value==libsbml.LIBSBML_OPERATION_SUCCESS:
                return
            else:
                err_msg = 'Error encountered trying to ' + message + '.' \
                        + 'LibSBML returned error code ' + str(value) + ': "' \
                        + libsbml.OperationReturnValue_toString(value).strip() + '"'
                self.logger.error(err_msg)
                #raise self.libSBMLError
                #raise libSBMLError(err_msg)
                #raise SystemExit(err_msg)
                raise AttributeError
        else:
            #self.logger.info(message)
            return None


    ## String to SBML ID
    #
    # Convert any String to one that is compatible with the SBML meta_id formatting requirements
    #
    # @param name The input string
    def _nameToSbmlId(self, name):
        IdStream = []
        count = 0
        end = len(name)
        if '0' <= name[count] and name[count] <= '9':
            IdStream.append('_')
        for count in range(0, end):
            if (('0' <= name[count] and name[count] <= '9') or
                    ('a' <= name[count] and name[count] <= 'z') or
                    ('A' <= name[count] and name[count] <= 'Z')):
                IdStream.append(name[count])
            else:
                IdStream.append('_')
        Id = ''.join(IdStream)
        if Id[len(Id) - 1] != '_':
            return Id
        return Id[:-1]


    ## String to hashed ID
    #
    # Hash an input string and then pass it to _nameToSbmlId()
    #
    # @param input string
    def _genMetaID(self, name):
        return self._nameToSbmlId(md5(str(name).encode('utf-8')).hexdigest())


    ## compare two dictionarry of lists and return the 
    #
    def _compareXref(self, current, toadd):
        toadd = copy.deepcopy(toadd)
        for database_id in current:
            try:
                list_diff = [i for i in toadd[database_id] if i not in current[database_id]]
                if not list_diff:
                    toadd.pop(database_id)
                else:
                    toadd[database_id] = list_diff
            except KeyError:
                pass
        return toadd


    ######################################################################
    ####################### Annotations ##################################
    ######################################################################


    ## Returns a default annotation string
    #
    # @param meta_id String or None Default meta ID
    #
    def _defaultBothAnnot(self, meta_id):
        return '''<annotation>
  <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/">
    <rdf:Description rdf:about="#'''+str(meta_id or '')+'''">
      <bqbiol:is>
        <rdf:Bag>
        </rdf:Bag>
      </bqbiol:is>
    </rdf:Description>
    <rdf:BRSynth rdf:about="#'''+str(meta_id or '')+'''">
      <brsynth:brsynth xmlns:brsynth="http://brsynth.eu">
      </brsynth:brsynth>
    </rdf:BRSynth>
  </rdf:RDF>
</annotation>'''


    ## Returns a default annotation string
    #
    # @param meta_id String or None Default meta ID
    #
    def _defaultBRSynthAnnot(self, meta_id):
        return '''<annotation>
  <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/">
    <rdf:BRSynth rdf:about="#'''+str(meta_id or '')+'''">
      <brsynth:brsynth xmlns:brsynth="http://brsynth.eu">
      </brsynth:brsynth>
    </rdf:BRSynth>
  </rdf:RDF>
</annotation>'''


    ## Returns a default annotation string
    #
    # @param meta_id String or None Default meta ID
    #
    def _defaultMIRIAMAnnot(self, meta_id):
        return '''<annotation>
  <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/">
    <rdf:Description rdf:about="#'''+str(meta_id or '')+'''">
      <bqbiol:is>
        <rdf:Bag>
          <rdf:li rdf:resource="http://identifiers.org/uniprot/test2i33"/>
        </rdf:Bag>
      </bqbiol:is>
    </rdf:Description>
  </rdf:RDF>
</annotation>'''



    ## Either add or update the value of a BRSynth annotation
    #
    # @sbase_obj libSBML object that may be compartment, reaction or species
    #
    def addUpdateBRSynth(self, sbase_obj, annot_header, value, units=None, isAlone=False, isList=False, isSort=True, meta_id=None):
        if isList:
            annotation = '''<annotation>
      <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/" xmlns:bqmodel="http://biomodels.net/model-qualifiers/">
        <rdf:BRSynth rdf:about="#adding">
          <brsynth:brsynth xmlns:brsynth="http://brsynth.eu">
            <brsynth:'''+str(annot_header)+'''>'''
            if isSort:
                for name in sorted(value, key=value.get, reverse=True):
                    if isAlone:
                        annotation += '<brsynth:'+str(name)+'>'+str(value[name])+'</brsynth:'+str(name)+'>'
                    else:
                        if units:
                            annotation += '<brsynth:'+str(name)+' units="'+str(units)+'" value="'+str(value[name])+'" />'
                        else:
                            annotation += '<brsynth:'+str(name)+' value="'+str(value[name])+'" />'
            else:
                for name in value:
                    if isAlone:
                        annotation += '<brsynth:'+str(name)+'>'+str(value[name])+'</brsynth:'+str(name)+'>'
                    else:
                        if units:
                            annotation += '<brsynth:'+str(name)+' units="'+str(units)+'" value="'+str(value[name])+'" />'
                        else:
                            annotation += '<brsynth:'+str(name)+' value="'+str(value[name])+'" />'
            annotation += '''
            </brsynth:'''+str(annot_header)+'''>
          </brsynth:brsynth>
        </rdf:BRSynth>
      </rdf:RDF>
    </annotation>'''
        else:
            #### create the string
            annotation = '''<annotation>
      <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/" xmlns:bqmodel="http://biomodels.net/model-qualifiers/">
        <rdf:BRSynth rdf:about="#adding">
          <brsynth:brsynth xmlns:brsynth="http://brsynth.eu">'''
            if isAlone:
                annotation += '<brsynth:'+str(annot_header)+'>'+str(value)+'</brsynth:'+str(annot_header)+'>'
            else:
                if units:
                    annotation += '<brsynth:'+str(annot_header)+' units="'+str(units)+'" value="'+str(value)+'" />'
                else:
                    annotation += '<brsynth:'+str(annot_header)+' value="'+str(value)+'" />'
            annotation += '''
          </brsynth:brsynth>
        </rdf:BRSynth>
      </rdf:RDF>
    </annotation>'''
        annot_obj = libsbml.XMLNode.convertStringToXMLNode(annotation)
        if annot_obj==None:
            self.logger.error('Cannot conver this string to annotation object: '+str(annotation))
            return False
        #### retreive the annotation object
        brsynth_annot = None
        obj_annot = sbase_obj.getAnnotation()
        if not obj_annot:
            #sbase_obj.setAnnotation(self._defaultBRSynthAnnot(meta_id))
            sbase_obj.setAnnotation(libsbml.XMLNode.convertStringToXMLNode(self._defaultBRSynthAnnot(meta_id)))
            obj_annot = sbase_obj.getAnnotation()
            if not obj_annot:
                self.logger.error('Cannot update BRSynth annotation')
                return False
        brsynth_annot = obj_annot.getChild('RDF').getChild('BRSynth').getChild('brsynth')
        if not brsynth_annot:
             self.logger.error('Cannot find the BRSynth annotation')
             return False
        #add the annotation and replace if it exists
        if brsynth_annot.getChild(annot_header).toXMLString()=='':
            toWrite_annot = annot_obj.getChild('RDF').getChild('BRSynth').getChild('brsynth').getChild(annot_header)
            self._checklibSBML(brsynth_annot.addChild(toWrite_annot), 'Adding annotation to the brsynth annotation')
        else:
            self._checklibSBML(brsynth_annot.removeChild(brsynth_annot.getIndex(annot_header)),
                'Removing annotation '+str(annot_header))
            toWrite_annot = annot_obj.getChild('RDF').getChild('BRSynth').getChild('brsynth').getChild(annot_header)
            self._checklibSBML(brsynth_annot.addChild(toWrite_annot), 'Adding annotation to the brsynth annotation')
        return True


    ## Function to update or create a MIRIAM annotation for a compartment, reaction or species
    #
    # @sbase_obj libSBML object that may be compartment, reaction or species
    # @type_param String Type of MIRIAM to add. May be 'compartment', 'reaction', 'species'
    # @xref type of anno
    #
    def addUpdateMIRIAM(self, sbase_obj, type_param, xref, meta_id=None):
        if not type_param in ['compartment', 'reaction', 'species']:
            self.logger.error('type_param must be '+str(['compartment', 'reaction', 'species'])+' not '+str(type_param))
            return False
        miriam_annot = None
        isReplace = False
        try:
            miriam_annot = sbase_obj.getAnnotation().getChild('RDF').getChild('Description').getChild('is').getChild('Bag')
            miriam_elements = self.readMIRIAMAnnotation(sbase_obj.getAnnotation())
            if not miriam_elements:
                #tmp_annot = libsbml.XMLNode.convertStringToXMLNode(self._defaultMIRIAMAnnot(meta_id))
                isReplace = True
                if not meta_id:
                    meta_id = self._genMetaID('tmp_addUpdateMIRIAM')
                #miriam_annot = libsbml.XMLNode.convertStringToXMLNode(self._defaultMIRIAMAnnot(meta_id))
                miriam_annot_1 = libsbml.XMLNode.convertStringToXMLNode(self._defaultBothAnnot(meta_id))
                miriam_annot = miriam_annot_1.getChild('RDF').getChild('Description').getChild('is').getChild('Bag')
            else:
                miriam_elements = None
        except AttributeError:
            try:
                #Cannot find MIRIAM annotation, create it
                isReplace = True
                if not meta_id:
                    meta_id = self._genMetaID('tmp_addUpdateMIRIAM')
                miriam_annot = libsbml.XMLNode.convertStringToXMLNode(self._defaultMIRIAMAnnot(meta_id))
                miriam_annot = miriam_annot.getChild('RDF').getChild('Description').getChild('is').getChild('Bag')
            except AttributeError:
                self.logger.error('Fatal error fetching the annotation')
                return False
        #compile the list of current species
        inside = {}
        for i in range(miriam_annot.getNumChildren()):
            single_miriam = miriam_annot.getChild(i)
            if single_miriam.getAttributes().getLength()>1:
                self.logger.error('MIRIAM annotations should never have more than 1: '+str(single_miriam.toXMLString()))
                continue
            single_miriam_attr = single_miriam.getAttributes()
            if not single_miriam_attr.isEmpty():
                try:
                    db = single_miriam_attr.getValue(0).split('/')[-2]
                    v = single_miriam_attr.getValue(0).split('/')[-1]
                    inside[self.header_miriam[type_param][db]].append(v)
                except KeyError:
                    try:
                        db = single_miriam_attr.getValue(0).split('/')[-2]
                        v = single_miriam_attr.getValue(0).split('/')[-1]
                        inside[self.header_miriam[type_param][db]] = [v]
                    except KeyError:
                        self.logger.error('Cannot find the self.header_miriram entry '+str(db))
                        return False
            else:
                self.logger.warning('Cannot return MIRIAM attribute')
                pass
        #add or ignore 
        toadd = self._compareXref(inside, xref)
        for database_id in toadd:
            for species_id in toadd[database_id]:
                #not sure how to avoid having it that way
                try:
                    #determine if the dictionnaries 
                    annotation = '''<annotation>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/" xmlns:bqmodel="http://biomodels.net/model-qualifiers/">
<rdf:Description rdf:about="#tmp">
  <bqbiol:is>
    <rdf:Bag>'''
                    if type_param=='species':
                        if database_id=='kegg' and species_id[0]=='C':
                            annotation += '''
          <rdf:li rdf:resource="http://identifiers.org/'''+self.miriam_header[type_param]['kegg_c']+str(species_id)+'''"/>'''
                        elif database_id=='kegg' and species_id[0]=='D':
                            annotation += '''
          <rdf:li rdf:resource="http://identifiers.org/'''+self.miriam_header[type_param]['kegg_d']+str(species_id)+'''"/>'''
                        else:
                            annotation += '''
          <rdf:li rdf:resource="http://identifiers.org/'''+self.miriam_header[type_param][database_id]+str(species_id)+'''"/>'''
                    else:
                        annotation += '''
          <rdf:li rdf:resource="http://identifiers.org/'''+self.miriam_header[type_param][database_id]+str(species_id)+'''"/>'''
                    annotation += '''
    </rdf:Bag>
  </bqbiol:is>
</rdf:Description>
</rdf:RDF>
</annotation>'''
                    toPass_annot = libsbml.XMLNode.convertStringToXMLNode(annotation)
                    toWrite_annot = toPass_annot.getChild('RDF').getChild('Description').getChild('is').getChild('Bag').getChild(0)
                    miriam_annot.insertChild(0, toWrite_annot)
                except KeyError:
                    #WARNING need to check this
                    self.logger.warning('Cannot find '+str(database_id)+' in self.miriam_header for '+str(type_param))
                    continue
        if isReplace:
            ori_miriam_annot = sbase_obj.getAnnotation()
            self._checklibSBML(ori_miriam_annot.getChild('RDF').getChild('Description').getChild('is').removeChild(0), 'Removing annotation "is"')
            self._checklibSBML(ori_miriam_annot.getChild('RDF').getChild('Description').getChild('is').addChild(miriam_annot), 'Adding annotation to the brsynth annotation')
        return True


    #####################################################################
    ########################## INPUT/OUTPUT #############################
    #####################################################################


    ## Open an SBML using libSBML 
    #
    # Situation where an SBML is passed to add the heterologous pathway
    #
    # @param inFile String Path to the input SBML file
    def readSBML(self, inFile):
        if not os.path.isfile(inFile):
            self.logger.error('Invalid input file')
            raise FileNotFoundError
        document = libsbml.readSBMLFromFile(inFile)
        self._checklibSBML(document, 'reading input file')
        errors = document.getNumErrors()
        #display the errors in the log accordning to the severity
        for err in [document.getError(i) for i in range(document.getNumErrors())]:
            #TODO if the error is related to packages not enabled (like groups or fbc) activate them
            if err.isFatal:
                self.logger.error('libSBML reading error: '+str(err.getShortMessage()))
                raise FileNotFoundError
            else:
                self.logger.warning('libSBML reading warning: '+str(err.getShortMessage()))
        model = document.getModel()
        if not model:
            loging.error('Either the file was not read correctly or the SBML is empty')
            raise FileNotFoundError
        self.document = document
        self.model = model
        #enabling the extra packages if they do not exists when reading a model
        if not self.model.isPackageEnabled('groups'):
            self._checklibSBML(self.model.enablePackage(
                'http://www.sbml.org/sbml/level3/version1/groups/version1',
                'groups',
                True),
                    'Enabling the GROUPS package')
            self._checklibSBML(self.document.setPackageRequired('groups', False), 'enabling groups package')
        if not self.model.isPackageEnabled('fbc'):
            self._checklibSBML(self.model.enablePackage(
                'http://www.sbml.org/sbml/level3/version1/fbc/version2',
                'fbc',
                True),
                    'Enabling the FBC package')
            self._checklibSBML(self.document.setPackageRequired('fbc', False), 'enabling FBC package')


    ## Export a libSBML model to file
    #
    # Export the libSBML model to an SBML file
    # 
    # @param model libSBML model to be saved to file
    # @param model_id model id, note that the name of the file will be that
    # @param path Non required parameter that will define the path where the model will be saved
    def writeSBML(self, path):
        ####### check the path #########
        #need to determine where are the path id's coming from
        p = None
        if path:
            if path[-1:]=='/':
                path = path[:-1]
            if not os.path.isdir(path):
                if self.path:
                    p = self.path
                else:
                    self.logger.error('The output path is not a directory: '+str(path))
                    return False
            else:
                p = path
        else:
            p = self.path
        ########## check and create folder #####
        if not os.path.exists(p):
            os.makedirs(p)
        libsbml.writeSBMLToFile(self.document, p+'/'+str(self.modelName)+'.sbml')
        return True


    #####################################################################
    ########################## FindCreate ###############################
    #####################################################################


    ## Find the objective (with only one reaction associated) based on the reaction ID and if not found create it
    #
    #
    #
    def findCreateObjective(self, reactions, coefficients, isMax=True, objective_id=None):
        fbc_plugin = self.model.getPlugin('fbc')
        self._checklibSBML(fbc_plugin, 'Getting FBC package')
        for objective in fbc_plugin.getListOfObjectives():
            if objective.getId()==objective_id:
                return objective_id
            if not set([i.getReaction() for i in objective.getListOfFluxObjectives()])-set(reactions):
                return objective.getId()
        #Cannot find a valid objective create it
        if not objective_id:
            self.logger.warning('Creating a new objective')
            objective_id = 'obj_rpFBA'+'_'.join(reactions)
        self.createMultiFluxObj(objective_id,
                                reactions,
                                coefficients,
                                isMax)
        return objective_id


    #####################################################################
    ########################## READ #####################################
    #####################################################################


    ## Return the reaction ID's and the pathway annotation
    #
    # TODO: replace the name of this function with readRPpathwayIDs
    def readRPpathwayIDs(self, pathway_id='rp_pathway'):
        groups = self.model.getPlugin('groups')
        rp_pathway = groups.getGroup(pathway_id)
        self._checklibSBML(rp_pathway, 'retreiving groups rp_pathway')
        toRet = []
        for member in rp_pathway.getListOfMembers():
            toRet.append(member.getIdRef())
        return toRet


    ## Read the reaction rules from the BRSYNTH annotation
    #
    #@param pathway_id default='rp_pathway' unique ID (per SBML) where the heterologous pathways are stored
    #@return toRet dictionnary with the reaction rule and rule_id as key
    def readRPrules(self, pathway_id='rp_pathway'):
        toRet = {}
        for reacId in self.readRPpathwayIDs(pathway_id):
            reac = self.model.getReaction(reacId)
            brsynth_annot = self.readBRSYNTHAnnotation(reac.getAnnotation())
            if not brsynth_annot['rule_id']=='' and not brsynth_annot['smiles']=='':
                toRet[brsynth_annot['rule_id']] = brsynth_annot['smiles'].replace('&gt;', '>')
        return toRet


    ## Return the species annitations 
    #
    #
    def readRPspecies(self, pathway_id='rp_pathway'):
        reacMembers = {}
        for reacId in self.readRPpathwayIDs(pathway_id):
            reacMembers[reacId] = {}
            reacMembers[reacId]['products'] = {}
            reacMembers[reacId]['reactants'] = {}
            reac = self.model.getReaction(reacId)
            for pro in reac.getListOfProducts():
                reacMembers[reacId]['products'][pro.getSpecies()] = pro.getStoichiometry()
            for rea in reac.getListOfReactants():
                reacMembers[reacId]['reactants'][rea.getSpecies()] = rea.getStoichiometry()
        return reacMembers


    ## Return the species
    #
    #
    def readUniqueRPspecies(self, pathway_id='rp_pathway'):
        rpSpecies = self.readRPspecies()
        toRet = []
        for i in rpSpecies:
            for y in rpSpecies[i]:
                for z in rpSpecies[i][y]:
                    if not z in toRet:
                        toRet.append(z)
        return toRet
        #reacMembers = self.readRPspecies(pathway_id)
        #return set(set(ori_rp_path['products'].keys())|set(ori_rp_path['reactants'].keys()))


    ## Return the MIRIAM annotations of species
    #
    #
    def readTaxonAnnotation(self, annot):
        try:
            toRet = {}
            bag = annot.getChild('RDF').getChild('Description').getChild('hasTaxon').getChild('Bag')
            for i in range(bag.getNumChildren()):
                str_annot = bag.getChild(i).getAttrValue(0)
                if str_annot=='':
                    self.logger.warning('This contains no attributes: '+str(bag.getChild(i).toXMLString()))
                    continue
                dbid = str_annot.split('/')[-2].split('.')[0]
                if len(str_annot.split('/')[-1].split(':'))==2:
                    cid = str_annot.split('/')[-1].split(':')[1]
                else:
                    cid = str_annot.split('/')[-1]
                if not dbid in toRet:
                    toRet[dbid] = []
                toRet[dbid].append(cid)
            return toRet
        except AttributeError:
            return {}


    ## Return the MIRIAM annotations of species
    #
    #
    def readMIRIAMAnnotation(self, annot):
        try:
            toRet = {}
            bag = annot.getChild('RDF').getChild('Description').getChild('is').getChild('Bag')
            for i in range(bag.getNumChildren()):
                str_annot = bag.getChild(i).getAttrValue(0)
                if str_annot=='':
                    self.logger.warning('This contains no attributes: '+str(bag.getChild(i).toXMLString()))
                    continue
                dbid = str_annot.split('/')[-2].split('.')[0]
                if len(str_annot.split('/')[-1].split(':'))==2:
                    cid = str_annot.split('/')[-1].split(':')[1]
                else:
                    cid = str_annot.split('/')[-1]
                if not dbid in toRet:
                    toRet[dbid] = []
                toRet[dbid].append(cid)
            return toRet
        except AttributeError:
            return {}


    ## Takes for input a libSBML annotatio object and returns a dictionnary of the annotations
    #
    def readBRSYNTHAnnotation(self, annot):
        toRet = {'dfG_prime_m': {},
                 'dfG_uncert': {},
                 'dfG_prime_o': {},
                 'path_id': None,
                 'step_id': None,
                 'sub_step_id': None,
                 'rule_score': None,
                 'smiles': None,
                 'inchi': None,
                 'inchikey': None,
                 'selenzyme': None,
                 'rule_id': None,
                 'rule_ori_reac': None,
                 'rule_score': None,
                 'global_score': None}
        bag = annot.getChild('RDF').getChild('BRSynth').getChild('brsynth')
        for i in range(bag.getNumChildren()):
            ann = bag.getChild(i)
            if ann=='':
                self.logger.warning('This contains no attributes: '+str(ann.toXMLString()))
                continue
            if ann.getName()=='dfG_prime_m' or ann.getName()=='dfG_uncert' or ann.getName()=='dfG_prime_o' or ann.getName()[0:4]=='fba_' or ann.getName()=='flux_value':
                toRet[ann.getName()] = {
                        'units': ann.getAttrValue('units'),
                        'value': float(ann.getAttrValue('value'))}
            elif ann.getName()=='path_id' or ann.getName()=='step_id' or ann.getName()=='sub_step_id':
                try:
                    toRet[ann.getName()] = int(ann.getAttrValue('value'))
                except ValueError:
                    toRet[ann.getName()] = None
            elif ann.getName()=='rule_score' or ann.getName()=='global_score' or ann.getName()[:5]=='norm_':
                try:
                    toRet[ann.getName()] = float(ann.getAttrValue('value'))
                except ValueError:
                    toRet[ann.getName()] = None
            elif ann.getName()=='smiles':
                toRet[ann.getName()] = ann.getChild(0).toXMLString().replace('&gt;', '>')
            #lists in the annotation
            elif ann.getName()=='selenzyme' or ann.getName()=='rule_ori_reac':
                toRet[ann.getName()] = {}
                for y in range(ann.getNumChildren()):
                    selAnn = ann.getChild(y)
                    try:
                        toRet[ann.getName()][selAnn.getName()] = float(selAnn.getAttrValue('value'))
                    except ValueError:
                        toRet[ann.getName()][selAnn.getName()] = selAnn.getAttrValue('value')
            else:
                toRet[ann.getName()] = ann.getChild(0).toXMLString()
        return toRet


    ## Function to return the products and the species associated with a reaction
    #
    # @return Dictionnary with right==product and left==reactants
    def readReactionSpecies_old(self, reaction, isID=False):
        #TODO: check that reaction is either an sbml species; if not check that its a string and that
        # it exists in the rpsbml model
        toRet = {'left': {}, 'right': {}}
        #reactants
        for i in range(reaction.getNumReactants()):
            reactant_ref = reaction.getReactant(i)
            reactant = self.model.getSpecies(reactant_ref.getSpecies())
            if isID:
                toRet['left'][reactant.getId()] = int(reactant_ref.getStoichiometry())
            else:
                toRet['left'][reactant.getName()] = int(reactant_ref.getStoichiometry())
        #products
        for i in range(reaction.getNumProducts()):
            product_ref = reaction.getProduct(i)
            product = self.model.getSpecies(product_ref.getSpecies())
            if isID:
                toRet['right'][product.getId()] = int(product_ref.getStoichiometry())
            else:
                toRet['right'][product.getName()] = int(product_ref.getStoichiometry())
            toRet['reversible'] = reaction.getReversible()
        return toRet


    ## Function to return the products and the species associated with a reaction
    #
    # @return Dictionnary with right==product and left==reactants
    def readReactionSpecies(self, reaction):
        #TODO: check that reaction is either an sbml species; if not check that its a string and that
        # it exists in the rpsbml model
        toRet = {'left': {}, 'right': {}}
        #reactants
        for i in range(reaction.getNumReactants()):
            reactant_ref = reaction.getReactant(i)
            toRet['left'][reactant_ref.getSpecies()] = int(reactant_ref.getStoichiometry())
        #products
        for i in range(reaction.getNumProducts()):
            product_ref = reaction.getProduct(i)
            toRet['right'][product_ref.getSpecies()] = int(product_ref.getStoichiometry())
        return toRet


    #####################################################################
    ######################### INQUIRE ###################################
    #####################################################################


    ## Function to find out if the model already contains a species according to its name
    #
    #
    def speciesExists(self, speciesName, compartment_id='MNXC3'):
        if speciesName in [i.getName() for i in self.model.getListOfSpecies()] or speciesName+'__64__'+compartment_id in [i.getId() for i in self.model.getListOfSpecies()]:
            return True
        return False


    ## Function to determine if a species CAN be a product of any reaction.
    #
    # Note that this is only determines if a species can possibly be produced, but does not
    # guarantee it
    #
    # @param species_id String ID of the species
    # @param ignoreReactions List Default is empty, ignore specific reactions
    def isSpeciesProduct(self, species_id, ignoreReactions=[]):
        #return all the parameters values
        param_dict = {i.getId(): i.getValue() for i in self.model.parameters}
        for reaction in self.model.getListOfReactions():
            if reaction.getId() not in ignoreReactions:
                #check that the function is reversible by reversibility and FBC bounds
                if reaction.reversible:
                    reaction_fbc = reaction.getPlugin('fbc')
                    #strict left to right
                    if param_dict[reaction_fbc.getLowerFluxBound()]>=0 and param_dict[reaction_fbc.getUpperFluxBound()]>0:
                        if species_id in [i.getSpecies() for i in reaction.getListOfProducts()]:
                            return True
                    #can go both ways
                    elif param_dict[reaction_fbc.getLowerFluxBound()]<0 and param_dict[reaction_fbc.getUpperFluxBound()]>0:
                        if species_id in [i.getSpecies() for i in reaction.getListOfProducts()]:
                            return True
                        elif species_id in [i.getSpecies() for i in reaction.getListOfReactants()]:
                            return True
                    #strict right to left
                    elif param_dict[reaction_fbc.getLowerFluxBound()]<0 and param_dict[reaction_fbc.getUpperFluxBound()]<=0 and param_dict[reaction_fbc.getLowerFluxBound()]<param_dict[reaction_fbc.getUpperFluxBound()]:
                        if species_id in [i.getSpecies() for i in reaction.getListOfReactants()]:
                            return True
                    else:
                        self.logger.warning('isSpeciesProduct does not find the directionailty of the reaction for reaction: '+str(species_id))
                        return True
                else:
                    #if the reaction is not reversible then product are the only way to create it
                    if species_id in [i.getSpecies() for i in reaction.getListOfProducts()]:
                        return True
        return False


    #########################################################################
    ################### CONVERT BETWEEEN FORMATS ############################
    #########################################################################


    ## Really used to complete the monocomponent reactions   
    #{'rule_id': 'RR-01-503dbb54cf91-49-F', 'right': {'TARGET_0000000001': 1}, 'left': {'MNXM2': 1, 'MNXM376': 1}, 'pathway_id': 1, 'step': 1, 'sub_step': 1, 'transformation_id': 'TRS_0_0_17'}
    #
    def outPathsDict(self, pathway_id='rp_pathway'):
        pathway = {}
        for member in self.readRPpathwayIDs(pathway_id):
            #TODO: need to find a better way
            reaction = self.model.getReaction(member)
            brsynthAnnot = self.readBRSYNTHAnnotation(reaction.getAnnotation())
            speciesReac = self.readReactionSpecies(reaction)
            step = {'reaction_id': member,
                    'reaction_rule': brsynthAnnot['smiles'],
                    'rule_score': brsynthAnnot['rule_score'],
                    'rule_id': brsynthAnnot['rule_id'],
                    'rule_ori_reac': brsynthAnnot['rule_ori_reac'],
                    'right': speciesReac['right'],
                    'left': speciesReac['left'],
                    'path_id': brsynthAnnot['path_id'],
                    'step': brsynthAnnot['step_id'],
                    'sub_step': brsynthAnnot['sub_step_id']}
            pathway[brsynthAnnot['step_id']] = step
        return pathway


    #########################################################################
    ############################# COMPARE MODELS ############################
    #########################################################################


    ## Find out if two libSBML Species or Reactions come from the same species
    #
    # Compare two dictionnaries and if any of the values of any of the same keys are the same then the 
    # function return True, and if none are found then return False
    #
    # @param libSBML Annotation object for one of the 
    # @return Boolean to determine if they are the same
    def compareBRSYNTHAnnotations(self, source_annot, target_annot):
        source_dict = self.readBRSYNTHAnnotation(source_annot)
        target_dict = self.readBRSYNTHAnnotation(target_annot)
        #list the common keys between the two
        for same_key in list(set(list(source_dict.keys())).intersection(list(target_dict.keys()))):
            if source_dict[same_key]==target_dict[same_key]:
                return True
        return False


    ## Find out if two libSBML Species or Reactions come from the same species
    #
    # Compare two dictionnaries and if any of the values of any of the same keys are the same then the 
    # function return True, and if none are found then return False
    #
    # @param libSBML Annotation object for one of the 
    # @return Boolean to determine if they are the same
    def compareMIRIAMAnnotations(self, source_annot, target_annot):
        source_dict = self.readMIRIAMAnnotation(source_annot)
        target_dict = self.readMIRIAMAnnotation(target_annot)
        #list the common keys between the two
        for com_key in set(list(source_dict.keys()))-(set(list(source_dict.keys()))-set(list(target_dict.keys()))):
            #compare the keys and if same is non-empty means that there 
            #are at least one instance of the key that is the same
            if bool(set(source_dict[com_key]) & set(target_dict[com_key])):
                return True
        return False


    ## Compare an annotation and a dictionnary structured
    #
    #
    def compareAnnotations_annot_dict(self, source_annot, target_dict):
        source_dict = self.readMIRIAMAnnotation(source_annot)
        #list the common keys between the two
        for com_key in set(list(source_dict.keys()))-(set(list(source_dict.keys()))-set(list(target_dict.keys()))):
            #compare the keys and if same is non-empty means that there 
            #are at least one instance of the key that is the same
            if bool(set(source_dict[com_key]) & set(target_dict[com_key])):
                return True
        return False


    ## Compare two dictionnaries sutructured as dict
    #
    #
    def compareAnnotations_dict_dict(self, source_dict, target_dict):
        #list the common keys between the two
        for com_key in set(list(source_dict.keys()))-(set(list(source_dict.keys()))-set(list(target_dict.keys()))):
            #compare the keys and if same is non-empty means that there 
            #are at least one instance of the key that is the same
            if bool(set(source_dict[com_key]) & set(target_dict[com_key])):
                return True
        return False


    ## Function to compare two SBML's RP pathways
    #
    # Function that compares the annotations of reactions and if not found, the annotations of all
    # species in that reaction to try to recover the correct ones. Because we are working with
    # intermediate cofactors for the RP generated pathways, the annotation crossreference will
    # not work. Best is to use the cross-reference to the original reaction
    #
    def compareRPpathways(self, measured_sbml):
        #return all the species annotations of the RP pathways
        try:
            meas_rp_species = measured_sbml.readRPspecies()
            found_meas_rp_species = measured_sbml.readRPspecies()
            for meas_step_id in meas_rp_species:
                meas_rp_species[meas_step_id]['annotation'] = measured_sbml.model.getReaction(meas_step_id).getAnnotation()
                found_meas_rp_species[meas_step_id]['found'] = False
                for spe_name in meas_rp_species[meas_step_id]['reactants']:
                    meas_rp_species[meas_step_id]['reactants'][spe_name] = measured_sbml.model.getSpecies(spe_name).getAnnotation()
                    found_meas_rp_species[meas_step_id]['reactants'][spe_name] = False
                for spe_name in meas_rp_species[meas_step_id]['products']:
                    meas_rp_species[meas_step_id]['products'][spe_name] = measured_sbml.model.getSpecies(spe_name).getAnnotation()
                    found_meas_rp_species[meas_step_id]['products'][spe_name] = False
            rp_rp_species = self.readRPspecies()
            for rp_step_id in rp_rp_species:
                rp_rp_species[rp_step_id]['annotation'] = self.model.getReaction(rp_step_id).getAnnotation()
                for spe_name in rp_rp_species[rp_step_id]['reactants']:
                    rp_rp_species[rp_step_id]['reactants'][spe_name] = self.model.getSpecies(spe_name).getAnnotation()
                for spe_name in rp_rp_species[rp_step_id]['products']:
                    rp_rp_species[rp_step_id]['products'][spe_name] = self.model.getSpecies(spe_name).getAnnotation()
        except AttributeError:
            self.logger.error('TODO: debug, for some reason some are passed as None here')
            return False, {}
        #compare the number of steps in the pathway
        if not len(meas_rp_species)==len(rp_rp_species):
            self.logger.warning('The pathways are not of the same length')
            return False, {}
        ############## compare using the reactions ###################
        for meas_step_id in measured_sbml.readRPpathwayIDs():
            for rp_step_id in rp_rp_species:
                if self.compareMIRIAMAnnotations(rp_rp_species[rp_step_id]['annotation'], meas_rp_species[meas_step_id]['annotation']):
                    found_meas_rp_species[meas_step_id]['found'] = True
                    found_meas_rp_species[meas_step_id]['rp_step_id'] = rp_step_id
                    break
        ############## compare using the species ###################
        for meas_step_id in measured_sbml.readRPpathwayIDs():
            #if not found_meas_rp_species[meas_step_id]['found']:
            for rp_step_id in rp_rp_species:
                # We test to see if the meas reaction elements all exist in rp reaction and not the opposite
                #because the measured pathways may not contain all the elements
                ########## reactants ##########
                for meas_spe_id in meas_rp_species[meas_step_id]['reactants']:
                    for rp_spe_id in rp_rp_species[rp_step_id]['reactants']:
                        if self.compareMIRIAMAnnotations(meas_rp_species[meas_step_id]['reactants'][meas_spe_id], rp_rp_species[rp_step_id]['reactants'][rp_spe_id]):
                            found_meas_rp_species[meas_step_id]['reactants'][meas_spe_id] = True
                            break
                        else:
                            if self.compareBRSYNTHAnnotations(meas_rp_species[meas_step_id]['reactants'][meas_spe_id], rp_rp_species[rp_step_id]['reactants'][rp_spe_id]):
                                found_meas_rp_species[meas_step_id]['reactants'][meas_spe_id] = True
                                break
                ########### products ###########
                for meas_spe_id in meas_rp_species[meas_step_id]['products']:
                    for rp_spe_id in rp_rp_species[rp_step_id]['products']:
                        if self.compareMIRIAMAnnotations(meas_rp_species[meas_step_id]['products'][meas_spe_id], rp_rp_species[rp_step_id]['products'][rp_spe_id]):
                            found_meas_rp_species[meas_step_id]['products'][meas_spe_id] = True
                            break
                        else:
                            if self.compareBRSYNTHAnnotations(meas_rp_species[meas_step_id]['products'][meas_spe_id], rp_rp_species[rp_step_id]['products'][rp_spe_id]):
                                found_meas_rp_species[meas_step_id]['products'][meas_spe_id] = True
                                break
                ######### test to see the difference
                pro_found = [found_meas_rp_species[meas_step_id]['products'][i] for i in found_meas_rp_species[meas_step_id]['products']]
                rea_found = [found_meas_rp_species[meas_step_id]['reactants'][i] for i in found_meas_rp_species[meas_step_id]['reactants']]
                if pro_found and rea_found:
                    if all(pro_found) and all(rea_found):
                        found_meas_rp_species[meas_step_id]['found'] = True
                        found_meas_rp_species[meas_step_id]['rp_step_id'] = rp_step_id
                        break
        ################# Now see if all steps have been found ############
        if all(found_meas_rp_species[i]['found'] for i in found_meas_rp_species):
            found_meas_rp_species['measured_model_id'] = measured_sbml.model.getId()
            found_meas_rp_species['rp_model_id'] = self.model.getId()
            return True, found_meas_rp_species
        else:
            return False, {}


    #########################################################################
    ############################# MODEL APPEND ##############################
    #########################################################################


    ## Set a given reaction's upper and lower bounds
    #
    # Sets the upper and lower bounds of a reaction. Note that if the numerical values passed
    # are not recognised, new parameters are created for each of them
    #
    def setReactionConstraints(self,
                               reaction_id,
                               upper_bound,
                               lower_bound,
                               unit='mmol_per_gDW_per_hr',
                               is_constant=True):
        reaction = self.model.getReaction(reaction_id)
        if not reaction:
            self.logger.error('Cannot find the reaction: '+str(reaction_id))
            return False
        reac_fbc = reaction.getPlugin('fbc')
        self._checklibSBML(reac_fbc, 'extending reaction for FBC')
        ########## upper bound #############
        old_upper_value = self.model.getParameter(reac_fbc.getUpperFluxBound()).value
        upper_param = self.createReturnFluxParameter(upper_bound, unit, is_constant)
        self._checklibSBML(reac_fbc.setUpperFluxBound(upper_param.getId()),
            'setting '+str(reaction_id)+' upper flux bound')
        ######### lower bound #############
        old_lower_value = self.model.getParameter(reac_fbc.getLowerFluxBound()).value
        lower_param = self.createReturnFluxParameter(lower_bound, unit, is_constant)
        self._checklibSBML(reac_fbc.setLowerFluxBound(lower_param.getId()),
            'setting '+str(reaction_id)+' lower flux bound')
        return old_upper_value, old_lower_value


    ##### ADD SOURCE FROM ORPHAN #####
    #if the heterologous pathway from the self.model contains a sink molecule that is not included in the 
    # original model (we call orhpan species) then add another reaction that creates it
    #TODO: that transports the reactions that creates the species in the
    # extracellular matrix and another reaction that transports it from the extracellular matrix to the cytoplasm
    #TODO: does not work
    def fillOrphan(self, 
            rpsbml=None, 
            pathway_id='rp_pathway', 
            compartment_id='MNXC3',
            upper_flux_bound=999999,
            lower_flux_bound=10):
        if rpsbml==None:
            model = self.model
        else:
            model = rpsbml.model
        self.logger.info('Adding the orphan species to the GEM model')
        #only for rp species
        groups = model.getPlugin('groups')
        rp_pathway = groups.getGroup(pathway_id)
        reaction_id = sorted([(int(''.join(x for x in i.id_ref if x.isdigit())), i.id_ref) for i in rp_pathway.getListOfMembers()], key=lambda tup: tup[0], reverse=True)[0][1]
        #for reaction_id in [i.getId() for i in self.model.getListOfReactions()]:
        for species_id in set([i.getSpecies() for i in model.getReaction(reaction_id).getListOfReactants()]+[i.getSpecies() for i in model.getReaction(reaction_id).getListOfProducts()]):
            if rpsbml==None:
                isSpePro = self.isSpeciesProduct(species_id, [reaction_id])
            else:
                isSpePro = rpsbml.isSpeciesProduct(species_id, [reaction_id])
            if not isSpePro:
                #create the step
                createStep = {'rule_id': None,
                              'left': {species_id.split('__')[0]: 1},
                              'right': {},
                              'step': None,
                              'sub_step': None,
                              'path_id': None,
                              'transformation_id': None,
                              'rule_score': None,
                              'mnxr': None}
                #create the model in the 
                if rpsbml==None:
                    self.createReaction('create_'+species_id,
                                        upper_flux_bound,
                                        lower_flux_bound,
                                        createStep,
                                        compartment_id)
                else:
                    rpsbml.createReaction('create_'+species_id,
                                        upper_flux_bound,
                                        lower_flux_bound,
                                        createStep,
                                        compartment_id)



    ## Merge two models species and reactions using the annotations to recognise the same species and reactions
    #
    # The source mode has to have both the GROUPS and FBC packages enabled in its SBML. The course must have a groups
    #called rp_pathway. If not use the readSBML() function to create a model
    # We add the reactions and species from the rpsbml to the target_model
    # 
    # @param target_model input libsbml model object where we will add the reactions and species from self.model
    # @param pathway_id String default is rp_pathway, name of the pathway id of the groups object
    # @param addOrphanSpecies Boolean Default False
    # @param bilevel_obj Tuple of size 2 with the weights associated with the targetSink and GEM objective function
    #
    #def mergeModels(self, target_rpsbml, pathway_id='rp_pathway', fillOrphanSpecies=False, compartment_id='MNXC3'):
    def mergeModels(self, target_rpsbml):#, fillOrphanSpecies=False, compartment_id='MNXC3'):
        #target_rpsbml.model = target_document.getModel()
        #Find the ID's of the similar target_rpsbml.model species
        ################ MODEL FBC ########################
        if not target_rpsbml.model.isPackageEnabled('fbc'):
            self._checklibSBML(target_rpsbml.model.enablePackage(
                'http://www.sbml.org/sbml/level3/version1/fbc/version2',
                'fbc',
                True),
                    'Enabling the FBC package')
        target_fbc = target_rpsbml.model.getPlugin('fbc')
        source_fbc = self.model.getPlugin('fbc')
        #note sure why one needs to set this as False
        self._checklibSBML(self.document.setPackageRequired('fbc', False), 'enabling FBC package')
        ################ UNITDEFINITIONS ######
        #return the list of unit definitions id's for the target to avoid overwritting
        #WARNING: this means that the original unit definitions will be prefered over the new one
        target_unitDefID = [i.getId() for i in target_rpsbml.model.getListOfUnitDefinitions()]
        for source_unitDef in self.model.getListOfUnitDefinitions():
            if not source_unitDef.getId() in target_unitDefID: #have to compare by ID since no annotation
                #create a new unitDef in the target
                target_unitDef = target_rpsbml.model.createUnitDefinition()
                self._checklibSBML(target_unitDef, 'fetching target unit definition')
                #copy unitDef info to the target
                self._checklibSBML(target_unitDef.setId(source_unitDef.getId()),
                    'setting target unit definition ID')
                self._checklibSBML(target_unitDef.setAnnotation(source_unitDef.getAnnotation()),
                    'setting target unit definition Annotation')
                for source_unit in source_unitDef.getListOfUnits():
                    #copy unit info to the target unitDef
                    target_unit = target_unitDef.createUnit()
                    self._checklibSBML(target_unit, 'creating target unit')
                    self._checklibSBML(target_unit.setKind(source_unit.getKind()),
                        'setting target unit kind')
                    self._checklibSBML(target_unit.setExponent(source_unit.getExponent()),
                        'setting target unit exponent')
                    self._checklibSBML(target_unit.setScale(source_unit.getScale()),
                        'setting target unit scale')
                    self._checklibSBML(target_unit.setMultiplier(source_unit.getMultiplier()),
                        'setting target unit multiplier')
                target_unitDefID.append(source_unitDef.getId()) #add to the list to make sure its not added twice
        ################ COMPARTMENTS ###############
        # WARNING: Compare by MIRIAM annotations
        for source_compartment in self.model.getListOfCompartments():
            found = False
            source_annotation = source_compartment.getAnnotation()
            if not source_annotation:
                self.logger.warning('No annotation for the source of compartment '+str(source_compartment.getId()))
                continue
            if source_compartment.getId() in [i.getId() for i in target_rpsbml.model.getListOfCompartments()]:
                found = True
            else:
                for target_compartment in target_rpsbml.model.getListOfCompartments():
                    target_annotation = target_compartment.getAnnotation()
                    if not target_annotation:
                        self.logger.warning('No annotation for the target of compartment: '+str(target_compartment.getId()))
                        continue
                    if self.compareMIRIAMAnnotations(source_annotation, target_annotation):
                        found = True
                        break
            if not found:
                target_compartment = target_rpsbml.model.createCompartment()
                self._checklibSBML(target_compartment, 'Creating target compartment')
                self._checklibSBML(target_compartment.setMetaId(source_compartment.getMetaId()),
                        'setting target metaId')
                #make sure that the ID is different
                if source_compartment.getId()==target_compartment.getId():
                    self._checklibSBML(target_compartment.setId(source_compartment.getId()+'_sourceModel'),
                            'setting target id')
                else:
                    self._checklibSBML(target_compartment.setId(source_compartment.getId()),
                            'setting target id')
                self._checklibSBML(target_compartment.setName(source_compartment.getName()),
                        'setting target name')
                self._checklibSBML(target_compartment.setConstant(source_compartment.getConstant()),
                        'setting target constant')
                self._checklibSBML(target_compartment.setAnnotation(source_compartment.getAnnotation()),
                        'setting target annotation')
                self._checklibSBML(target_compartment.setSBOTerm(source_compartment.getSBOTerm()),
                        'setting target annotation')
        ################ PARAMETERS ###########
        #WARNING: here we compare by ID
        targetParametersID = [i.getId() for i in target_rpsbml.model.getListOfParameters()]
        for source_parameter in self.model.getListOfParameters():
            if not source_parameter.getId() in targetParametersID:
                target_parameter = target_rpsbml.model.createParameter()
                self._checklibSBML(target_parameter, 'creating target parameter')
                self._checklibSBML(target_parameter.setId(source_parameter.getId()), 'setting target parameter ID')
                self._checklibSBML(target_parameter.setSBOTerm(source_parameter.getSBOTerm()),
                    'setting target parameter SBO')
                self._checklibSBML(target_parameter.setUnits(source_parameter.getUnits()),
                    'setting target parameter Units')
                self._checklibSBML(target_parameter.setValue(source_parameter.getValue()),
                    'setting target parameter Value')
                self._checklibSBML(target_parameter.setConstant(source_parameter.getConstant()),
                    'setting target parameter ID')
        ################ FBC GENE PRODUCTS ########################
        #WARNING: here we compare by ID
        targetGenProductID = [i.getId() for i in target_fbc.getListOfGeneProducts()]
        for source_geneProduct in source_fbc.getListOfGeneProducts():
            if not source_geneProduct.getId() in targetGenProductID:
                target_geneProduct = target_fbc.createGeneProduct()
                self._checklibSBML(target_geneProduct, 'creating target gene product')
                self._checklibSBML(target_geneProduct.setId(source_geneProduct.getId()),
                    'setting target gene product id')
                self._checklibSBML(target_geneProduct.setLabel(source_geneProduct.getLabel()),
                    'setting target gene product label')
                self._checklibSBML(target_geneProduct.setName(source_geneProduct.getName()),
                    'setting target gene product name')
                self._checklibSBML(target_geneProduct.setMetaId(source_geneProduct.getMetaId()),
                    'setting target gene product meta_id')
        ############### FBC OBJECTIVES ############
        #WARNING: here we compare by ID
        targetObjectiveID = [i.getId() for i in target_fbc.getListOfObjectives()]
        sourceObjectiveID = [i.getId() for i in source_fbc.getListOfObjectives()]
        for source_objective in source_fbc.getListOfObjectives():
            if not source_objective.getId() in targetObjectiveID:
                target_objective = target_fbc.createObjective()
                self._checklibSBML(target_objective, 'creating target objective')
                self._checklibSBML(target_objective.setId(source_objective.getId()), 'setting target objective')
                self._checklibSBML(target_objective.setName(source_objective.getName()), 'setting target objective')
                self._checklibSBML(target_objective.setType(source_objective.getType()),
                        'setting target objective type')
                for source_fluxObjective in source_objective.getListOfFluxObjectives():
                    target_fluxObjective = target_objective.createFluxObjective()
                    self._checklibSBML(target_fluxObjective, 'creating target flux objective')
                    self._checklibSBML(target_fluxObjective.setName(source_fluxObjective.getName()),
                        'setting target flux objective name')
                    self._checklibSBML(target_fluxObjective.setCoefficient(source_fluxObjective.getCoefficient()),
                        'setting target flux objective coefficient')
                    self._checklibSBML(target_fluxObjective.setReaction(source_fluxObjective.getReaction()),
                        'setting target flux objective reaction')
                    self._checklibSBML(target_fluxObjective.setAnnotation(source_fluxObjective.getAnnotation()),
                        'setting target flux obj annotation from source flux obj')
                self._checklibSBML(target_objective.setAnnotation(source_objective.getAnnotation()),
                        'setting target obj annotation from source obj')
        ################ SPECIES ####################
        #Try to detect the current model species
        model_species_convert = {} #this is passed to be used for reaction ID conversion
        toAdd_model_species_ids = [i.getId() for i in self.model.getListOfSpecies()]
        model_species_ids = [i.getId() for i in self.model.getListOfSpecies()]
        targetModel_species_ids = [i.getId() for i in target_rpsbml.model.getListOfSpecies()]
        for model_species in self.model.getListOfSpecies():
            if model_species.getId() in targetModel_species_ids:
                toAdd_model_species_ids.remove(model_species.getId())
                continue
            model_species_annotation = model_species.getAnnotation()
            if not model_species_annotation:
                self.logger.error('Cannot find annotations for species: '+str(model_species.getId()))
                continue
            for targetModel_species in target_rpsbml.model.getListOfSpecies():
                targetModel_species_annotation = targetModel_species.getAnnotation()
                if not targetModel_species_annotation:
                    self.logger.warning('Cannot find annotations for species: '+str(targetModel_species.getId()))
                    continue
                if self.compareMIRIAMAnnotations(model_species_annotation, targetModel_species_annotation):
                    #if we find a match from MIRIAM and not from ID we assume that it supersedes the ID and need to replace it in reactions
                    model_species_convert[model_species.getId()] = targetModel_species.getId()
                    toAdd_model_species_ids.remove(model_species.getId())
                    break
        #add the new species
        for model_species_id in toAdd_model_species_ids:
            model_species = self.model.getSpecies(model_species_id)
            if not model_species:
                self.logger.error('Cannot retreive model species: '+str(model_species_id))
            self._checklibSBML(model_species, 'fetching source species')
            targetModel_species = target_rpsbml.model.createSpecies()
            self._checklibSBML(targetModel_species, 'creating species')
            self._checklibSBML(targetModel_species.setMetaId(model_species.getMetaId()),
                    'setting target metaId')
            self._checklibSBML(targetModel_species.setId(model_species.getId()),
                    'setting target id')
            self._checklibSBML(targetModel_species.setCompartment(model_species.getCompartment()),
                    'setting target compartment')
            self._checklibSBML(targetModel_species.setInitialConcentration(
                model_species.getInitialConcentration()),
                    'setting target initial concentration')
            self._checklibSBML(targetModel_species.setBoundaryCondition(
                model_species.getBoundaryCondition()),
                    'setting target boundary concentration')
            self._checklibSBML(targetModel_species.setHasOnlySubstanceUnits(
                model_species.getHasOnlySubstanceUnits()),
                    'setting target has only substance units')
            self._checklibSBML(targetModel_species.setBoundaryCondition(
                model_species.getBoundaryCondition()),
                    'setting target boundary condition')
            self._checklibSBML(targetModel_species.setConstant(model_species.getConstant()),
                'setting target constant')
            self._checklibSBML(targetModel_species.setAnnotation(model_species.getAnnotation()),
                'setting target annotation')
        ################ REACTIONS ###################
        # Here we go by ID and check if it exists. If it does we need to check the species are the same
        # and if not, add it changing the ID #model_species_convert
        toAdd_model_reaction_ids = [i.getId() for i in self.model.getListOfReactions()]
        model_reaction_ids = [i.getId() for i in self.model.getListOfReactions()]
        targetModel_reaction_ids = [i.getId() for i in target_rpsbml.model.getListOfReactions()]
        for model_reaction in self.model.getListOfReactions():
            if model_reaction.getId() in targetModel_reaction_ids:
                toAdd_model_reaction_ids.remove(model_reaction.getId())
                continue
            model_reaction_speciesID = []
            for reactant in model_reaction.getListOfReactants():
                if reactant.species in model_species_convert:
                    model_reaction_speciesID.append(model_species_convert[reactant.species])
                else:
                    model_reaction_speciesID.append(reactant.species)
            model_reaction_productsID = []
            for product in model_reaction.getListOfProducts():
                if product.species in model_species_convert:
                    model_reaction_productsID.append(model_species_convert[product.species])
                else:
                    model_reaction_productsID.append(product.species)
            for targetModel_reaction in target_rpsbml.model.getListOfReactions():
                #loop through target model reactions
                targetModel_reaction_speciesID = [i.species for i in targetModel_reaction.getListOfReactants()]
                targetModel_reaction_productsID = [i.species for i in targetModel_reaction.getListOfProducts()]
                if not set(model_reaction_speciesID)-set(targetModel_reaction_speciesID) and not set(model_reaction_productsID)-set(targetModel_reaction_productsID):
                    self.logging.info('The reactions species and products are the same')
                    toAdd_model_reaction_ids.remove(model_reaction.getId())
                    continue
        #add the new reactions
        for model_reaction_id in toAdd_model_reaction_ids:
            model_reaction = self.model.getReaction(model_reaction_id)
            if not model_reaction:
                self.logger.error('Cannot retreive model reaction: '+str(model_reaction_id))
            self._checklibSBML(model_reaction, 'fetching source reaction')
            target_reaction = target_rpsbml.model.createReaction()
            self._checklibSBML(target_reaction, 'create reaction')
            target_fbc = target_reaction.getPlugin('fbc')
            self._checklibSBML(target_fbc, 'fetching target FBC package')
            source_fbc = model_reaction.getPlugin('fbc')
            self._checklibSBML(source_fbc, 'fetching source FBC package')
            source_upperFluxBound = source_fbc.getUpperFluxBound()
            self._checklibSBML(source_upperFluxBound, 'fetching upper flux bound')
            self._checklibSBML(target_fbc.setUpperFluxBound(source_upperFluxBound),
                    'setting upper flux bound')
            source_lowerFluxBound = source_fbc.getLowerFluxBound()
            self._checklibSBML(source_lowerFluxBound, 'fetching lower flux bound')
            self._checklibSBML(target_fbc.setLowerFluxBound(source_lowerFluxBound),
                    'setting lower flux bound')
            self._checklibSBML(target_reaction.setId(model_reaction.getId()), 'set reaction id')
            self._checklibSBML(target_reaction.setName(model_reaction.getName()), 'set name')
            self._checklibSBML(target_reaction.setSBOTerm(model_reaction.getSBOTerm()),
                    'setting the reaction system biology ontology (SBO)') #set as process
            #TODO: consider having the two parameters as input to the function
            self._checklibSBML(target_reaction.setReversible(model_reaction.getReversible()),
                    'set reaction reversibility flag')
            self._checklibSBML(target_reaction.setFast(model_reaction.getFast()),
                    'set reaction "fast" attribute')
            self._checklibSBML(target_reaction.setMetaId(model_reaction.getMetaId()), 'setting species meta_id')
            self._checklibSBML(target_reaction.setAnnotation(model_reaction.getAnnotation()),
                    'setting annotation for source reaction')
            #Reactants
            for model_reaction_speciesID in [i.species for i in model_reaction.getListOfReactants()]:
                target_reactant = target_reaction.createReactant()
                self._checklibSBML(target_reactant, 'create target reactant')
                if model_reaction_speciesID in model_species_convert:
                    self._checklibSBML(target_reactant.setSpecies(
                        model_species_convert[model_reaction_speciesID]), 'assign reactant species')
                else:
                    self._checklibSBML(target_reactant.setSpecies(model_reaction_speciesID),
                        'assign reactant species')
                source_reactant = model_reaction.getReactant(model_reaction_speciesID)
                self._checklibSBML(source_reactant, 'fetch source reactant')
                self._checklibSBML(target_reactant.setConstant(source_reactant.getConstant()),
                        'set "constant" on species '+str(source_reactant.getConstant()))
                self._checklibSBML(target_reactant.setStoichiometry(source_reactant.getStoichiometry()),
                        'set stoichiometry ('+str(source_reactant.getStoichiometry)+')')
            #Products
            for model_reaction_productID in [i.species for i in model_reaction.getListOfProducts()]:
                target_reactant = target_reaction.createProduct()
                self._checklibSBML(target_reactant, 'create target reactant')
                if model_reaction_productID in model_species_convert:
                    self._checklibSBML(target_reactant.setSpecies(
                        model_species_convert[model_reaction_productID]), 'assign reactant product')
                else:
                    self._checklibSBML(target_reactant.setSpecies(model_reaction_productID),
                        'assign reactant product')
                source_reactant = model_reaction.getProduct(model_reaction_productID)
                self._checklibSBML(source_reactant, 'fetch source reactant')
                self._checklibSBML(target_reactant.setConstant(source_reactant.getConstant()),
                        'set "constant" on product '+str(source_reactant.getConstant()))
                self._checklibSBML(target_reactant.setStoichiometry(source_reactant.getStoichiometry()),
                        'set stoichiometry ('+str(source_reactant.getStoichiometry)+')')
        #### GROUPS #####
        #TODO loop through the groups to add them
        if not target_rpsbml.model.isPackageEnabled('groups'):
            self._checklibSBML(target_rpsbml.model.enablePackage(
                'http://www.sbml.org/sbml/level3/version1/groups/version1',
                'groups',
                True),
                    'Enabling the GROUPS package')
        #!!!! must be set to false for no apparent reason
        self._checklibSBML(self.document.setPackageRequired('groups', False), 'enabling groups package')
        source_groups = self.model.getPlugin('groups')
        self._checklibSBML(source_groups, 'fetching the source model groups')
        target_groups = target_rpsbml.model.getPlugin('groups')
        self._checklibSBML(target_groups, 'fetching the target model groups')
        '''
        if not pathway_id in [i.getId() for i in target_groups.getListOfGroups()]:
            self._checklibSBML(target_groups.addGroup(source_groups.getGroup(pathway_id)),
                    'copying the source groups "rp_pathway" to the target groups')
        '''
        for group in target_groups.getListOfGroups():
            self._checklibSBML(target_groups.addGroup(group),
                    'copying the source groups "rp_pathway" to the target groups')
        ###### TITLES #####
        target_rpsbml.model.setId(target_rpsbml.model.getId()+'__'+self.model.getId())
        target_rpsbml.model.setName(target_rpsbml.model.getName()+' merged with '+self.model.getId())
        '''
        if fillOrphanSpecies==True:
            self.fillOrphan(target_rpsbml, pathway_id, compartment_id)
        '''


    #########################################################################
    ############################# MODEL CREATION FUNCTIONS ##################
    #########################################################################


    ## Create libSBML model instance
    #
    # Function that creates a new libSBML model instance and initiates it with the appropriate packages. Creates a cytosol compartment
    #
    # @param name The name of the model
    # @param model_id The id of the mode
    # @param meta_id meta_id of the model. Default None means that we will generate a hash from the model_id
    def createModel(self, name, model_id, meta_id=None):
        ## sbmldoc
        self.sbmlns = libsbml.SBMLNamespaces(3,1)
        self._checklibSBML(self.sbmlns, 'generating model namespace')
        self._checklibSBML(self.sbmlns.addPkgNamespace('groups',1), 'Add groups package')
        self._checklibSBML(self.sbmlns.addPkgNamespace('fbc',2), 'Add FBC package')
        #sbmlns = libsbml.SBMLNamespaces(3,1,'groups',1)
        self.document = libsbml.SBMLDocument(self.sbmlns)
        self._checklibSBML(self.document, 'generating model doc')
        #!!!! must be set to false for no apparent reason
        self._checklibSBML(self.document.setPackageRequired('fbc', False), 'enabling FBC package')
        #!!!! must be set to false for no apparent reason
        self._checklibSBML(self.document.setPackageRequired('groups', False), 'enabling groups package')
        ## sbml model
        self.model = self.document.createModel()
        self._checklibSBML(self.model, 'generating the model')
        self._checklibSBML(self.model.setId(model_id), 'setting the model ID')
        model_fbc = self.model.getPlugin('fbc')
        model_fbc.setStrict(True)
        if meta_id==None:
            meta_id = self._genMetaID(model_id)
        self._checklibSBML(self.model.setMetaId(meta_id), 'setting model meta_id')
        self._checklibSBML(self.model.setName(name), 'setting model name')
        self._checklibSBML(self.model.setTimeUnits('second'), 'setting model time unit')
        self._checklibSBML(self.model.setExtentUnits('mole'), 'setting model compartment unit')
        self._checklibSBML(self.model.setSubstanceUnits('mole'), 'setting model substance unit')


    ## Create libSBML compartment 
    #
    # cytoplasm compartment TODO: consider seperating it in another function if another compartment is to be created
    #
    # @param model libSBML model object to add the compartment
    # @param size Set the compartement size
    # @return boolean Execution success
    #TODO: set the compName as None by default. To do that you need to regenerate the compXref to 
    #use MNX ids as keys instead of the string names
    def createCompartment(self, size, compId, compName, compXref, meta_id=None):
        comp = self.model.createCompartment()
        self._checklibSBML(comp, 'create compartment')
        self._checklibSBML(comp.setId(compId), 'set compartment id')
        if compName:
            self._checklibSBML(comp.setName(compName), 'set the name for the cytoplam')
        self._checklibSBML(comp.setConstant(True), 'set compartment "constant"')
        self._checklibSBML(comp.setSize(size), 'set compartment "size"')
        self._checklibSBML(comp.setSBOTerm(290), 'set SBO term for the cytoplasm compartment')
        if meta_id==None:
            meta_id = self._genMetaID(compId)
        self._checklibSBML(comp.setMetaId(meta_id), 'set the meta_id for the compartment')
        ############################ MIRIAM ############################
        self.addUpdateMIRIAM(comp, 'compartment', compXref, meta_id)


    ## Create libSBML unit definition
    #
    # Function that creates a unit definition (composed of one or more units)
    #
    # @param model libSBML model to add the unit definition
    # @param unit_id ID for the unit definition
    # @param meta_id meta_id for the unit definition. If None creates a hash from unit_id
    # @return Unit definition
    def createUnitDefinition(self, unit_id, meta_id=None):
        unitDef = self.model.createUnitDefinition()
        self._checklibSBML(unitDef, 'creating unit definition')
        self._checklibSBML(unitDef.setId(unit_id), 'setting id')
        if meta_id==None:
            meta_id = self._genMetaID(unit_id)
        self._checklibSBML(unitDef.setMetaId(meta_id), 'setting meta_id')
        #self.unitDefinitions.append(unit_id)
        return unitDef


    ## Create libSBML unit
    #
    # Function that created a unit
    #
    # @param unitDef libSBML unit definition
    # @param libsmlunit libSBML unit parameter
    # @param exponent Value for the exponent (ex 10^5 mol/sec)
    # @param scale Value for the scale 
    # @param multiplier Value for the multiplie
    # @return Unit
    def createUnit(self, unitDef, libsbmlunit, exponent, scale, multiplier):
        unit = unitDef.createUnit()
        self._checklibSBML(unit, 'creating unit')
        self._checklibSBML(unit.setKind(libsbmlunit), 'setting the kind of unit')
        self._checklibSBML(unit.setExponent(exponent), 'setting the exponenent of the unit')
        self._checklibSBML(unit.setScale(scale), 'setting the scale of the unit')
        self._checklibSBML(unit.setMultiplier(multiplier), 'setting the multiplier of the unit')


    ## Create libSBML parameters
    #
    # Parameters, in our case, are used for the bounds for FBA analysis. Unit parameter must be an instance of unitDefinition
    #
    # @param parameter_id SBML id
    # @param value Float value for this parameter
    # @param unit libSBML unit parameter
    # @param meta_id String Optional parameter for SBML meta_id
    # @return libSBML parameter object
    def createReturnFluxParameter(self, 
            value, 
            unit='mmol_per_gDW_per_hr', 
            is_constant=True,
            parameter_id=None,
            meta_id=None):
        if parameter_id:
            param_id = parameter_id
        else:
            if value>=0:
                param_id = 'B_'+str(round(abs(value), 4)).replace('.', '_')
            else:
                param_id = 'B__'+str(round(abs(value), 4)).replace('.', '_')
        if param_id in [i.getId() for i in self.model.getListOfParameters()]:
            return self.model.getParameter(param_id)
        else:
            newParam = self.model.createParameter()
            self._checklibSBML(newParam, 'Creating a new parameter object')
            self._checklibSBML(newParam.setConstant(is_constant), 'setting as constant')
            self._checklibSBML(newParam.setId(param_id), 'setting ID')
            self._checklibSBML(newParam.setValue(value), 'setting value')
            self._checklibSBML(newParam.setUnits(unit), 'setting units')
            self._checklibSBML(newParam.setSBOTerm(625), 'setting SBO term')
            if meta_id==None:
                meta_id = self._genMetaID(parameter_id)
            self._checklibSBML(newParam.setMetaId(meta_id), 'setting meta ID')
            #self.parameters.append(parameter_id)
            return newParam


    ## Create libSBML reaction
    #
    # Create a reaction. fluxBounds is a list of libSBML.UnitDefinition, length of exactly 2 with the first position that is the upper bound and the second is the lower bound. reactants_dict and reactants_dict are dictionnaries that hold the following parameters: name, compartment, stoichiometry
    #
    # @param name Name for the reaction
    # @param reaction_id Reaction ID
    # @param fluxUpperBounds FBC id for the upper flux bound for this reaction
    # @param fluxLowerBounds FBC id for the lower flux bound for this reaction
    # @param step 2D dictionnary with the following structure {'left': {'name': stoichiometry, ...}, 'right': {}}
    # @param reaction_smiles String smiles description of this reaction (added in BRSYNTH annotation)
    # @param compartment_id String Optinal parameter compartment ID
    # @param isTarget Boolean Flag to suppress the warning that the passed step is missing information. Used in this case for the target compound
    # @param hetero_group Groups Optional parameter object that holds all the heterologous pathways
    # @param meta_id String Optional parameter reaction meta_id
    # @return meta_id meta ID for this reaction
    #TODO as of now not generic, works when creating a new SBML file, but no checks if modifying existing SBML file
    def createReaction(self,
            reac_id,
            fluxUpperBound,
            fluxLowerBound,
            step,
            compartment_id,
            reaction_smiles=None,
            reacXref={},
            pathway_id=None,
            meta_id=None):
        reac = self.model.createReaction()
        self._checklibSBML(reac, 'create reaction')
        ################ FBC ####################
        reac_fbc = reac.getPlugin('fbc')
        self._checklibSBML(reac_fbc, 'extending reaction for FBC')
        #bounds
        upper_bound = self.createReturnFluxParameter(fluxUpperBound)
        lower_bound = self.createReturnFluxParameter(fluxLowerBound)
        self._checklibSBML(reac_fbc.setUpperFluxBound(upper_bound.getId()), 'setting '+str(reac_id)+' upper flux bound')
        self._checklibSBML(reac_fbc.setLowerFluxBound(lower_bound.getId()), 'setting '+str(reac_id)+' lower flux bound')
        #########################################
        #reactions
        self._checklibSBML(reac.setId(reac_id), 'set reaction id') #same convention as cobrapy
        #self._checklibSBML(reac.setName(str(reac_id)+), 'set name') #same convention as cobrapy
        self._checklibSBML(reac.setSBOTerm(185), 'setting the system biology ontology (SBO)') #set as process
        #TODO: consider having the two parameters as input to the function
        self._checklibSBML(reac.setReversible(True), 'set reaction reversibility flag')
        self._checklibSBML(reac.setFast(False), 'set reaction "fast" attribute')
        if meta_id==None:
            meta_id = self._genMetaID(reac_id)
        self._checklibSBML(reac.setMetaId(meta_id), 'setting species meta_id')
        #TODO: check that the species exist
        #reactants_dict
        for reactant in step['left']:
            spe = reac.createReactant()
            self._checklibSBML(spe, 'create reactant')
            #use the same writing convention as CobraPy
            self._checklibSBML(spe.setSpecies(str(reactant)+'__64__'+str(compartment_id)), 'assign reactant species')
            #TODO: check to see the consequences of heterologous parameters not being constant
            self._checklibSBML(spe.setConstant(True), 'set "constant" on species '+str(reactant))
            self._checklibSBML(spe.setStoichiometry(float(step['left'][reactant])),
                'set stoichiometry ('+str(float(step['left'][reactant]))+')')
        #TODO: check that the species exist
        #products_dict
        for product in step['right']:
            pro = reac.createProduct()
            self._checklibSBML(pro, 'create product')
            self._checklibSBML(pro.setSpecies(str(product)+'__64__'+str(compartment_id)), 'assign product species')
            #TODO: check to see the consequences of heterologous parameters not being constant
            self._checklibSBML(pro.setConstant(True), 'set "constant" on species '+str(product))
            self._checklibSBML(pro.setStoichiometry(float(step['right'][product])),
                'set the stoichiometry ('+str(float(step['right'][product]))+')')
        ############################ MIRIAM ############################
        self.addUpdateMIRIAM(reac, 'reaction', reacXref, meta_id)
        ###### BRSYNTH additional information ########
        if reaction_smiles:
            self.addUpdateBRSynth(reac, 'smiles', reaction_smiles, None, True, False, False, meta_id)
        if step['rule_id']:
            self.addUpdateBRSynth(reac, 'rule_id', step['rule_id'], None, True, False, False, meta_id)
        #TODO: need to change the name and content (to dict) upstream
        if step['rule_ori_reac']:
            #self.addUpdateBRSynthList(reac, 'rule_ori_reac', step['rule_ori_reac'], True, False, meta_id)
            self.addUpdateBRSynth(reac, 'rule_ori_reac', step['rule_ori_reac'], None, False, True, False, meta_id)
            #sbase_obj, annot_header, value, units=None, isAlone=False, isList=False, isSort=True, meta_id=None)
        if step['rule_score']:
            self.addUpdateBRSynth(reac, 'rule_score', step['rule_score'], None, False, False, False, meta_id)
        if step['path_id']:
            self.addUpdateBRSynth(reac, 'path_id', step['path_id'], None, False, False, False, meta_id)
        if step['step']:
            self.addUpdateBRSynth(reac, 'step_id', step['step'], None, False, False, False, meta_id)
        if step['sub_step']:
            self.addUpdateBRSynth(reac, 'sub_step_id', step['sub_step'], None, False, False, False, meta_id)
        #### GROUPS #####
        if not pathway_id==None:
            groups_plugin = self.model.getPlugin('groups')
            hetero_group = groups_plugin.getGroup(pathway_id)
            if not hetero_group:
                self.logger.warning('The pathway_id '+str(pathway_id)+' does not exist in the model')
            else:
                newM = hetero_group.createMember()
                self._checklibSBML(newM, 'Creating a new groups member')
                self._checklibSBML(newM.setIdRef(reac_id), 'Setting name to the groups member')


    ## Create libSBML reaction
    #
    # Create a reaction. fluxBounds is a list of libSBML.UnitDefinition, length of exactly 2 with the first position that is the upper bound and the second is the lower bound. reactants_dict and reactants_dict are dictionnaries that hold the following parameters: name, compartment_id, stoichiometry
    #
    # @param chemIdDictionnary containing all the cross references that we know of, can be empty)
    # @param chemXref Dictionnary containing all the cross references that we know of, can be empty
    # @param meta_id Name for the reaction
    # @param inchi String Inchi associated with this species
    # @param smiles String SMILES associated with this species
    # @param compartment_id String Set this species to belong to another compartment_id than the one globally set by self.compartment_id
    # @param charge Optional parameter describing the charge of the molecule of interest
    # @param chemForm Optional chemical formulae of the substrate (not SMILES or InChI)
    # @param dG Optinal Thermodynamics constant for this species
    # @param dG_uncert Optional Uncertainty associated with the thermodynamics of the reaction 
    def createSpecies(self,
            species_id,
            compartment_id,
            species_name=None,
            chemXref={},
            inchi=None,
            inchikey=None,
            smiles=None,
            species_group_id=None,
            meta_id=None):
            #TODO: add these at some point -- not very important
            #charge=0,
            #chemForm=''):
        spe = self.model.createSpecies()
        self._checklibSBML(spe, 'create species')
        ##### FBC #####
        spe_fbc = spe.getPlugin('fbc')
        self._checklibSBML(spe_fbc, 'creating this species as an instance of FBC')
        #spe_fbc.setCharge(charge) #### These are not required for FBA 
        #spe_fbc.setChemicalFormula(chemForm) #### These are not required for FBA
        #if compartment_id:
        self._checklibSBML(spe.setCompartment(compartment_id), 'set species spe compartment')
        #else:
        #    #removing this could lead to errors with xref
        #    self._checklibSBML(spe.setCompartment(self.compartment_id), 'set species spe compartment')
        #ID same structure as cobrapy
        #TODO: determine if this is always the case or it will change
        self._checklibSBML(spe.setHasOnlySubstanceUnits(False), 'set substance units')
        self._checklibSBML(spe.setBoundaryCondition(False), 'set boundary conditions')
        self._checklibSBML(spe.setConstant(False), 'set constant')
        #useless for FBA (usefull for ODE) but makes Copasi stop complaining
        self._checklibSBML(spe.setInitialConcentration(1.0), 'set an initial concentration')
        #same writting convention as COBRApy
        self._checklibSBML(spe.setId(str(species_id)+'__64__'+str(compartment_id)), 'set species id')
        if meta_id==None:
            meta_id = self._genMetaID(species_id)
        self._checklibSBML(spe.setMetaId(meta_id), 'setting reaction meta_id')
        if species_name==None:
            self._checklibSBML(spe.setName(species_id), 'setting name for the metabolite '+str(species_id))
        else:
            self._checklibSBML(spe.setName(species_name), 'setting name for the metabolite '+str(species_name))
        #this is setting MNX id as the name
        #this is setting the name as the input name
        ###### annotation ###
        self.addUpdateMIRIAM(spe, 'species', chemXref, meta_id)
        ###### BRSYNTH additional information ########
        if smiles:
            self.addUpdateBRSynth(spe, 'smiles', smiles, None, True, False, False, meta_id)
            #                   sbase_obj, annot_header, value, units=None, isAlone=False, isList=False, isSort=True, meta_id=None)
        if inchi:
            self.addUpdateBRSynth(spe, 'inchi', inchi, None, True, False, False, meta_id)
        if inchikey:
            self.addUpdateBRSynth(spe, 'inchikey', inchikey, None, True, False, False, meta_id)
        #### GROUPS #####
        #TODO: check that it actually exists
        if not species_group_id==None:
            groups_plugin = self.model.getPlugin('groups')
            hetero_group = groups_plugin.getGroup(species_group_id)
            if not hetero_group:
                self.logger.warning('The species_group_id '+str(species_group_id)+' does not exist in the model')
                #TODO: consider creating it if
            else:
                newM = hetero_group.createMember()
                self._checklibSBML(newM, 'Creating a new groups member')
                self._checklibSBML(newM.setIdRef(str(species_id)+'__64__'+str(compartment_id)), 'Setting name to the groups member')


    ## Create libSBML pathway
    #
    # Create the collection of reactions that constitute the pathway using the Groups package and create the custom IBIBSA annotations
    #
    # @param model libSBML model to add the unit definition
    # @param reaction_id Reaction ID
    # @param name Name for the reaction
    # @param fluxBounds list of size 2 that describe the FBC upper and lower bounds for this reactions flux
    # @param reactants list of species that are the reactants of this reaction
    # @param products list of species that are the products of this reaction
    # @param reaction_smiles String smiles description of this reaction (added in BRSYNTH annotation)
    # @return hetero_group The number libSBML groups object to pass to createReaction to categorise the new reactions
    def createPathway(self, pathway_id, meta_id=None):
        groups_plugin = self.model.getPlugin('groups')
        new_group = groups_plugin.createGroup()
        new_group.setId(pathway_id)
        if meta_id==None:
            meta_id = self._genMetaID(pathway_id)
        new_group.setMetaId(meta_id)
        new_group.setKind(libsbml.GROUP_KIND_COLLECTION)
        new_group.setAnnotation(self._defaultBRSynthAnnot(meta_id))


    ## Create libSBML gene
    #
    # Create the list of genes in the model including its custom BRSYNTH annotatons
    #
    # @param model libSBML model to add the unit definition
    # @param reac libSBML reaction object
    # @param step_id The step for the number of 
    # @return libSBML gene object
    def createGene(self, reac, step_id, meta_id=None):
        #TODO: pass this function to Pablo for him to fill with parameters that are appropriate for his needs
        geneName = 'RP'+str(step_id)+'_gene'
        fbc_plugin = self.model.getPlugin('fbc')
        #fbc_plugin = reac.getPlugin("fbc")
        gp = fbc_plugin.createGeneProduct()
        gp.setId(geneName)
        if meta_id==None:
            meta_id = self._genMetaID(str(geneName))
        gp.setMetaId(meta_id)
        gp.setLabel('gene_'+str(step_id))
        gp.setAssociatedSpecies('RP'+str(step_id))
        ##### NOTE: The parameters here require the input from Pablo to determine what he needs
        #gp.setAnnotation(self._defaultBothAnnot(meta_id))


    ## Create libSBML flux objective 
    # WARNING DEPRECATED -- use the createMultiFluxObj() with lists of size one to define an objective function
    # with a single reaction
    # Using the FBC package one can add the FBA flux objective directly to the model. This function sets a particular reaction as objective with maximization or minimization objectives
    #
    # @param model libSBML model to add the unit definition
    # @param fluxobj_id The id given to this particular objective
    # @param reactionName The name or id of the reaction that we are setting a flux objective
    # @param coefficient FBA coefficient 
    # @param isMax Boolean to determine if we are maximizing or minimizing the objective
    # @param meta_id Set the meta_id
    # @return Boolean exit code
    def createFluxObj(self, fluxobj_id, reactionName, coefficient, isMax=True, meta_id=None):
        fbc_plugin = self.model.getPlugin('fbc')
        target_obj = fbc_plugin.createObjective()
        #TODO: need to define inpiut metaID
        target_obj.setAnnotation(self._defaultBRSynthAnnot(meta_id))
        target_obj.setId(fluxobj_id)
        if isMax:
            target_obj.setType('maximize')
        else:
            target_obj.setType('minimize')
        fbc_plugin.setActiveObjectiveId(fluxobj_id) # this ensures that we are using this objective when multiple
        target_flux_obj = target_obj.createFluxObjective()
        target_flux_obj.setReaction(reactionName)
        target_flux_obj.setCoefficient(coefficient)
        if meta_id==None:
            meta_id = self._genMetaID(str(fluxobj_id))
        target_flux_obj.setMetaId(meta_id)
        target_flux_obj.setAnnotation(self._defaultBRSynthAnnot(meta_id))


    ## Create libSBML flux objective 
    #
    # Using the FBC package one can add the FBA flux objective directly to the model. This function sets a particular reaction as objective with maximization or minimization objectives
    #
    # @param model libSBML model to add the unit definition
    # @param fluxobj_id The id given to this particular objective
    # @param reactionName The name or id of the reaction that we are setting a flux objective
    # @param coefficient FBA coefficient 
    # @param isMax Boolean to determine if we are maximizing or minimizing the objective
    # @param meta_id Set the meta_id
    # @return Boolean exit code
    def createMultiFluxObj(self, fluxobj_id, reactionNames, coefficients, isMax=True, meta_id=None):
        if not len(reactionNames)==len(coefficients):
            self.logger.error('The size of reactionNames is not the same as coefficients')
            return False
        fbc_plugin = self.model.getPlugin('fbc')
        target_obj = fbc_plugin.createObjective()
        target_obj.setAnnotation(self._defaultBRSynthAnnot(meta_id))
        target_obj.setId(fluxobj_id)
        if isMax:
            target_obj.setType('maximize')
        else:
            target_obj.setType('minimize')
        fbc_plugin.setActiveObjectiveId(fluxobj_id) # this ensures that we are using this objective when multiple
        for reac, coef in zip(reactionNames, coefficients):
            target_flux_obj = target_obj.createFluxObjective()
            target_flux_obj.setReaction(reac)
            target_flux_obj.setCoefficient(coef)
            if meta_id==None:
                meta_id = self._genMetaID(str(fluxobj_id))
            target_flux_obj.setMetaId(meta_id)
            target_flux_obj.setAnnotation(self._defaultBRSynthAnnot(meta_id))


    ##############################################################################################
    ############################### Generic Model ################################################
    ##############################################################################################


    ## Generate a generic model 
    #
    # Since we will be using the same type of parameters for the RetroPath model, this function
    # generates a libSBML model with parameters that will be mostly used
    #
    #
    #
    def genericModel(self, 
            modelName, 
            model_id, 
            compXref, 
            compartment_id,
            upper_flux_bound=999999,
            lower_flux_bound=0):
        self.createModel(modelName, model_id)
        # mmol_per_gDW_per_hr -- flux
        unitDef = self.createUnitDefinition('mmol_per_gDW_per_hr')
        self.createUnit(unitDef, libsbml.UNIT_KIND_MOLE, 1, -3, 1)
        self.createUnit(unitDef, libsbml.UNIT_KIND_GRAM, 1, 0, 1)
        self.createUnit(unitDef, libsbml.UNIT_KIND_SECOND, 1, 0, 3600)
        # kj_per_mol -- thermodynamics
        gibbsDef = self.createUnitDefinition('kj_per_mol')
        self.createUnit(gibbsDef, libsbml.UNIT_KIND_JOULE, 1, 3, 1)
        self.createUnit(gibbsDef, libsbml.UNIT_KIND_MOLE, -1, 1, 1)
        ### set the bounds
        upBound = self.createReturnFluxParameter(upper_flux_bound)
        lowBound = self.createReturnFluxParameter(lower_flux_bound)
        #compartment
        #TODO: create a new compartment 
        #self.createCompartment(1, 'MNXC3', 'cytoplasm', compXref)
        #try to recover the name from the Xref
        try:
            name = compXref['name'][0]
        except KeyError:
            name = compartment_id+'_name'
        self.createCompartment(1, compartment_id, name, compXref)

