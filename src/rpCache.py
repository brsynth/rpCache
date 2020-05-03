import os
from rdkit.Chem import MolFromSmiles, MolFromInchi, MolToSmiles, MolToInchi, MolToInchiKey, AddHs
import csv
import logging
from pickle import load as pickle_load
from pickle import loads as pickle_loads
from pickle import dumps as pickle_dumps
import gzip
import urllib.request
import re
import tarfile
import shutil
import redis
from RedisDict import RedisDict
from argparse import ArgumentParser as argparse_ArgumentParser
import sys
import time
from itertools import chain as itertools_chain
sys.path.insert(0, '/home/rpCache')

#######################################################
################### rpCache  ##########################
#######################################################

from collections import deque
try:
    from reprlib import repr
except ImportError:
    pass

def total_size(o, handlers={}, verbose=False):
    """ Returns the approximate memory footprint an object and all of its contents.

    Automatically finds the contents of the following builtin containers and
    their subclasses:  tuple, list, deque, dict, set and frozenset.
    To search other containers, add handlers to iterate over their contents:

        handlers = {SomeContainerClass: iter,
                    OtherContainerClass: OtherContainerClass.get_elements}

    """
    dict_handler = lambda d: chain.from_iterable(d.items())
    all_handlers = {tuple: iter,
                    list: iter,
                    deque: iter,
                    dict: dict_handler,
                    set: iter,
                    frozenset: iter,
                   }
    all_handlers.update(handlers)     # user handlers take precedence
    seen = set()                      # track which object id's have already been seen
    default_size = getsizeof(0)       # estimate sizeof object without __sizeof__

    def sizeof(o):
        if id(o) in seen:       # do not double count the same object
            return 0
        seen.add(id(o))
        s = getsizeof(o, default_size)

        if verbose:
            print(s, type(o), repr(o), file=stderr)

        for typ, handler in all_handlers.items():
            if isinstance(o, typ):
                s += sum(map(sizeof, handler(o)))
                break
        return s

    return sizeof(o)


## Class to generate the cache
#
# Contains all the functions that parse different files, used to calculate the thermodynamics and the FBA of the
#the other steps. These should be called only when the files have changes
class rpCache:



    ## Cache constructor
    #
    # @param self The object pointer
    # @param inputPath The path to the folder that contains all the input/output files required
    # @param db Mode of storing objects ('file' or 'redis')
    def __init__(self, db='file', print_infos=False):

        self.store_mode = db

        if self.store_mode!='file':
            self.redis = redis.StrictRedis(host=self.store_mode, port=6379, db=0, decode_responses=True)
            self.deprecatedMNXM_mnxm = RedisDict('deprecatedMNXM_mnxm', self.redis)
            self.deprecatedMNXR_mnxr = RedisDict('deprecatedMNXR_mnxr', self.redis)
            self.mnxm_strc = RedisDict('mnxm_strc', self.redis)
            self.chemXref = RedisDict('chemXref', self.redis)
            self.rr_reactions = RedisDict('rr_reactions', self.redis)
            self.chebi_mnxm = RedisDict('chebi_mnxm', self.redis)
            # rpReader attributes
            self.inchikey_mnxm = RedisDict('inchikey_mnxm', self.redis)
            self.compXref = RedisDict('compXref', self.redis)
            self.name_compXref = RedisDict('name_compXref', self.redis)
        else:
            self.deprecatedMNXM_mnxm = None
            self.deprecatedMNXR_mnxr = None
            self.mnxm_strc = None
            self.chemXref = None
            self.rr_reactions = None
            self.chebi_mnxm = None
            # rpReader attributes
            self.inchikey_mnxm = None
            self.compXref = None
            self.name_compXref = None

        #given by Thomas
        self.logger = logging.getLogger(__name__)
        self.logger.info('Started instance of rpCache')

        self.print = print_infos

        # Common attribues
        self.convertMNXM = {'MNXM162231': 'MNXM6',
                            'MNXM84': 'MNXM15',
                            'MNXM96410': 'MNXM14',
                            'MNXM114062': 'MNXM3',
                            'MNXM145523': 'MNXM57',
                            'MNXM57425': 'MNXM9',
                            'MNXM137': 'MNXM588022'}

        self.dirname = os.path.dirname(os.path.abspath( __file__ ))+"/.."


        if not self.loadCache():
            raise ValueError


    #####################################################
    ################# ERROR functions ###################
    #####################################################

    ## Error function for the convertion of structures
    #
    class Error(Exception):
        pass


    ## Error function for the convertion of structures
    #
    class DepictionError(Error):
        def __init__(self, message):
            #self.expression = expression
            self.message = message


    ##########################################################
    ################## Private Functions #####################
    ##########################################################

    ## Private function to fetch the required data, parse them and generate the pickle
    #
    #  Opens the previously generated cache to the object memory
    #
    # @param The oject pointer
    # @return Boolean detemining the success of the function or not
    def loadCache(self, fetchInputFiles=False):

        attr_names = {
        #   KEY                  : [attribute(s) name(s) list, args list to the function]
            'deprecatedMNXM_mnxm': [['deprecatedMNXM_mnxm'], ['chem_xref.tsv']],
            'deprecatedMNXR_mnxr': [['deprecatedMNXR_mnxr'], ['reac_xref.tsv']],
            'mnxm_strc': [['mnxm_strc'], ['rr_compounds.tsv', 'chem_prop.tsv']],
            'chemXref': [['chemXref'], ['chem_xref.tsv']],
            'chebi_mnxm': [['chebi_mnxm'], []],
            'rr_reactions': [['rr_reactions'], ['rules_rall.tsv']],
            'inchikey_mnxm': [['inchikey_mnxm'], []],
            'compXref': [['compXref', 'name_compXref'], ['comp_xref.tsv']]
        }

        start_time = end_time = 0
        # For each attribute name
        for attr_name in attr_names:
            # For each real attribute
            for attr in attr_names[attr_name][0]:
                # If cache is not already loaded
                if not self.cache_loaded(attr):
                    # Then, for each corresponding input file
                    ####################### Fetch the files if necessary ######################
                    for input in attr_names[attr_name][1]:
                        filename = self.dirname+'/input_cache/'+input
                        if not os.path.isfile(filename):
                            print("Downloading "+input+"...", end = '', flush=True)
                            start_time = time.time()
                            self.fetch_input_file(input, self.dirname+'/input_cache/')
                            end_time = time.time()
#                            print(" (%.2fs)" % (end_time - start_time))
                        else:
                            print(filename+" already downloaded ", end = '', flush=True)
                        self.print_OK(end_time-start_time)
                    ###################### Populate the cache #################################
                    start_time = time.time()
                    self.populate_cache(attr_names[attr_name], self.dirname+'/cache')
                    end_time = time.time()
                    self.print_OK(end_time-start_time)

                elif not self.store_mode=='file':

                    print(" ".join(attr_names[attr_name][0])+" already in db ", end = '', flush=True)
                    self.print_OK()

            # Load cache from file
            if self.store_mode=='file':
                for i in range(len(attr_names[attr_name][0])):
                    _attr_name = attr_names[attr_name][0][i]
                    filename = self.dirname+'/cache/'+_attr_name+'.pickle'
                    print("Loading "+_attr_name+" from "+filename+"...", end = '', flush=True)
                    data = self.load_cache_from_file(filename)
                    setattr(self, _attr_name, data)
                    self.print_OK()


        return True



    def cache_loaded(self, attr):
        if self.store_mode=='file':
            return os.path.isfile(self.dirname+'/cache/'+attr+'.pickle')
        else:
            return getattr(self, attr).exists()


    def fetch_input_file(self, file, dir):

        #################### make the local folders ############################
        # input_cache
        if not os.path.isdir(dir):
            os.mkdir(dir)

        url = 'https://www.metanetx.org/cgi-bin/mnxget/mnxref/'

        # 3xCommon + rpReader
        if file in ['reac_xref.tsv', 'chem_xref.tsv', 'chem_prop.tsv', 'comp_xref.tsv']:
            urllib.request.urlretrieve(url+file, dir+'/'+file)

        #TODO: need to add this file to the git or another location
        if file in ['rr_compounds.tsv', 'rxn_recipes.tsv']:
            urllib.request.urlretrieve('https://retrorules.org/dl/this/is/not/a/secret/path/rr02',
                                       dir+'/rr02_more_data.tar.gz')
            tar = tarfile.open(dir+'/rr02_more_data.tar.gz', 'r:gz')
            tar.extractall(dir)
            tar.close()
            shutil.move(dir+'/rr02_more_data/compounds.tsv',
                        dir+'/rr_compounds.tsv')
            shutil.move(dir+'/rr02_more_data/rxn_recipes.tsv',
                        dir)
            os.remove(dir+'rr02_more_data.tar.gz')
            shutil.rmtree(dir+'rr02_more_data')

        if file=='rules_rall.tsv':
            urllib.request.urlretrieve('https://retrorules.org/dl/preparsed/rr02/rp3/hs',
                                       dir+'/retrorules_rr02_rp3_hs.tar.gz')
            tar = tarfile.open(dir+'/retrorules_rr02_rp3_hs.tar.gz', 'r:gz')
            tar.extractall(dir)
            tar.close()
            shutil.move(dir+'/retrorules_rr02_rp3_hs/retrorules_rr02_flat_all.tsv', dir+'/rules_rall.tsv')
            os.remove(dir+'/retrorules_rr02_rp3_hs.tar.gz')
            shutil.rmtree(dir+'/retrorules_rr02_rp3_hs')




    def populate_cache(self, attributes, dir):

        # cache
        if not os.path.isdir(dir):
            os.mkdir(dir)

        data = self.gen_cache(attributes[0], [self.dirname+'/input_cache/'+input_file for input_file in attributes[1]])

        for i in range(len(data)):
            _attr_name = attributes[0][i]
            print("Storing "+_attr_name+" to "+self.store_mode+"...", end = '', flush=True)
            method = getattr(self, 'store_cache_to_'+self.store_mode)
            method(dir+'/'+_attr_name+'.pickle', data[i])



    def print_OK(self, time=-1):
        sys.stdout.write("\033[0;32m") # Green
        print(" OK", end = '', flush=True)
        sys.stdout.write("\033[0;0m") # Reset
        if time!=-1: print(" (%.2fs)" % time, end = '', flush=True)
        print()

    def print_FAILED(self):
        sys.stdout.write("\033[1;31m") # Red
        print(" Failed")
        sys.stdout.write("\033[0;0m") # Reset
        print()

    def gen_cache(self, attr_names, args):
        try:
            results = []
            print("Generating "+" ".join(attr_names)+"...", end = '', flush=True)
            # Choose method according to attribute name
            method = getattr(self, '_m_'+attr_names[0])
            # Apply method and expand 'args' list as arguments
            # Put results in a list
            results = [method(*args)]
            if type(results[0]) is tuple:
                results = list(itertools_chain(results[0]))
            self.print_OK()
            return results
        except:
            self.print_FAILED()
            raise


    def load_cache_from_file(self, filename, gzip=False):
        if gzip:
            return pickle_load(gzip.open(filename, 'rb'))
        else:
            return pickle_load(open(filename, 'rb'))

    def store_cache_to_file(self, filename, data, gzip=False):
        pickle_obj = pickle_dumps(data)
        if gzip:
            filename += '.gz'
            with gzip.open(filename, "wb") as f:
            	f.write(pickle_obj)
        else:
            with open(filename, "wb") as f:
             	f.write(pickle_obj)

    def store_cache_to_db(self, attr_name, data):
        setattr(self, attr_name, RedisDict(attr_name, self.redis, data))



    ## Convert chemical depiction to others type of depictions
    #
    # Usage example:
    # - convert_depiction(idepic='CCO', otype={'inchi', 'smiles', 'inchikey'})
    # - convert_depiction(idepic='InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3', itype='inchi', otype={'inchi', 'smiles', 'inchikey'})
    #  @param self The object pointer
    #  @param idepic String depiction to be converted, str
    #  @param itype type of depiction provided as input, str
    #  @param otype types of depiction to be generated, {"", "", ..}
    #  @return odepic generated depictions, {"otype1": "odepic1", ..}
    def _convert_depiction(self, idepic, itype='smiles', otype={'inchikey'}):
        # Import (if needed)
        if itype == 'smiles':
            rdmol = MolFromSmiles(idepic, sanitize=True)
        elif itype == 'inchi':
            rdmol = MolFromInchi(idepic, sanitize=True)
        else:
            raise NotImplementedError('"{}" is not a valid input type'.format(itype))
        if rdmol is None:  # Check imprt
            raise self.DepictionError('Import error from depiction "{}" of type "{}"'.format(idepic, itype))
        # Export
        odepic = dict()
        for item in otype:
            if item == 'smiles':
                odepic[item] = MolToSmiles(rdmol)  # MolToSmiles is tricky, one mays want to check the possible options..
            elif item == 'inchi':
                odepic[item] = MolToInchi(rdmol)
            elif item == 'inchikey':
                odepic[item] = MolToInchiKey(rdmol)
            else:
                raise NotImplementedError('"{}" is not a valid output type'.format(otype))
        return odepic


    ## Function to create a dictionnary of old to new chemical id's
    #
    #  Generate a one-to-one dictionnary of old id's to new ones. Private function
    #
    # TODO: check other things about the mnxm emtry like if it has the right structure etc...
    def _checkMNXMdeprecated(self, mnxm):
        try:
            return self.deprecatedMNXM_mnxm[mnxm]
        except (KeyError, TypeError):
            return mnxm


    ## Function to create a dictionnary of old to new reaction id's
    #
    # TODO: check other things about the mnxm emtry like if it has the right structure etc...
    def _checkMNXRdeprecated(self, mnxr):
        try:
            return self.deprecatedMNXR_mnxr[mnxr]
        except (KeyError, TypeError):
            return mnxr


    #################################################################
    ################## Public functions #############################
    #################################################################


    #[TODO] merge the two functions
    ## Function to parse the chem_xref.tsv file of MetanetX
    #
    #  Generate a dictionnary of old to new MetanetX identifiers to make sure that we always use the freshest id's.
    # This can include more than one old id per new one and thus returns a dictionnary. Private function
    #
    #  @param self Object pointer
    #  @param chem_xref_path Input file path
    #  @return Dictionnary of identifiers
    #TODO: save the self.deprecatedMNXM_mnxm to be used in case there rp_paths uses an old version of MNX
    # def deprecatedMNXM(self, chem_xref_path):
    #     deprecatedMNXM_mnxm = {}
    #     with open(chem_xref_path) as f:
    #         c = csv.reader(f, delimiter='\t')
    #         for row in c:
    #             if not row[0][0]=='#':
    #                 mnx = row[0].split(':')
    #                 if mnx[0]=='deprecated':
    #                     deprecatedMNXM_mnxm[mnx[1]] = row[1]
    #         deprecatedMNXM_mnxm.update(self.convertMNXM)
    #         deprecatedMNXM_mnxm['MNXM01'] = 'MNXM1'
    #     return deprecatedMNXM_mnxm

    def _deprecatedMNX(self, xref_path):
        deprecatedMNX_mnx = {}
        with open(xref_path) as f:
            c = csv.reader(f, delimiter='\t')
            for row in c:
                if not row[0][0]=='#':
                    mnx = row[0].split(':')
                    if mnx[0]=='deprecated':
                        deprecatedMNX_mnx[mnx[1]] = row[1]
        return deprecatedMNX_mnx

    def _m_deprecatedMNXM_mnxm(self, chem_xref_path):
        deprecatedMNXM_mnxm = {}
        deprecatedMNXM_mnxm = self._deprecatedMNX(chem_xref_path)
        deprecatedMNXM_mnxm.update(self.convertMNXM)
        deprecatedMNXM_mnxm['MNXM01'] = 'MNXM1'
        return deprecatedMNXM_mnxm

    ## Function to parse the reac_xref.tsv file of MetanetX
    #
    #  Generate a dictionnary of old to new MetanetX identifiers to make sure that we always use the freshest id's.
    # This can include more than one old id per new one and thus returns a dictionnary. Private function
    #
    #  @param self Object pointer
    #  @param reac_xref_path Input file path
    #  @return Dictionnary of identifiers
    def _m_deprecatedMNXR_mnxr(self, reac_xref_path):
        return self._deprecatedMNX(reac_xref_path)


    ## Function to parse the chemp_prop.tsv file from MetanetX and compounds.tsv from RetroRules. Uses the InchIkey as key to the dictionnary
    #
    #  Generate a dictionnary gaving the formula, smiles, inchi and inchikey for the components
    #
    #  @param self Object pointer
    #  @param chem_prop_path Input file path
    #  @return mnxm_strc Dictionnary of formula, smiles, inchi and inchikey
    def _m_mnxm_strc(self, rr_compounds_path, chem_prop_path):
        mnxm_strc = {}
        for row in csv.DictReader(open(rr_compounds_path), delimiter='\t'):
            tmp = {'formula':  None,
                    'smiles': None,
                    'inchi': row['inchi'],
                    'inchikey': None,
                    'mnxm': self._checkMNXMdeprecated(row['cid']),
                    'name': None}
            try:
                resConv = self._convert_depiction(idepic=tmp['inchi'], itype='inchi', otype={'smiles','inchikey'})
                for i in resConv:
                    tmp[i] = resConv[i]
            except self.DepictionError as e:
                self.logger.warning('Could not convert some of the structures: '+str(tmp))
                self.logger.warning(e)
            mnxm_strc[tmp['mnxm']] = tmp
        with open(chem_prop_path) as f:
            c = csv.reader(f, delimiter='\t')
            for row in c:
                if not row[0][0]=='#':
                    mnxm = self._checkMNXMdeprecated(row[0])
                    tmp = {'formula':  row[2],
                            'smiles': row[6],
                            'inchi': row[5],
                            'inchikey': row[8],
                            'mnxm': mnxm,
                            'name': row[1]}
                    for i in tmp:
                        if tmp[i]=='' or tmp[i]=='NA':
                            tmp[i] = None
                    if mnxm in mnxm_strc:
                        mnxm_strc[mnxm]['formula'] = row[2]
                        mnxm_strc[mnxm]['name'] = row[1]
                        if not mnxm_strc[mnxm]['smiles'] and tmp['smiles']:
                            mnxm_strc[mnxm]['smiles'] = tmp['smiles']
                        if not mnxm_strc[mnxm]['inchikey'] and tmp['inchikey']:
                            mnxm_strc[mnxm]['inchikey'] = tmp['inchikey']
                    else:
                        #check to see if the inchikey is valid or not
                        otype = set({})
                        if not tmp['inchikey']:
                            otype.add('inchikey')
                        if not tmp['smiles']:
                            otype.add('smiles')
                        if not tmp['inchi']:
                            otype.add('inchi')
                        itype = ''
                        if tmp['inchi']:
                            itype = 'inchi'
                        elif tmp['smiles']:
                            itype = 'smiles'
                        else:
                            self.logger.warning('No valid entry for the convert_depiction function')
                            continue
                        try:
                            resConv = self._convert_depiction(idepic=tmp[itype], itype=itype, otype=otype)
                            for i in resConv:
                                tmp[i] = resConv[i]
                        except self.DepictionError as e:
                            self.logger.warning('Could not convert some of the structures: '+str(tmp))
                            self.logger.warning(e)
                        mnxm_strc[tmp['mnxm']] = tmp
        return mnxm_strc


    ## Function to parse the chem_xref.tsv file of MetanetX
    #
    #  Generate a dictionnary of all cross references for a given chemical id (MNX) to other database id's
    #
    #  @param self Object pointer
    #  @param chem_xref_path Input file path
    #  @return a The dictionnary of identifiers
    #TODO: save the self.deprecatedMNXM_mnxm to be used in case there rp_paths uses an old version of MNX
    def _m_chemXref(self, chem_xref_path):
        chemXref = {}
        with open(chem_xref_path) as f:
            c = csv.reader(f, delimiter='\t')
            for row in c:
                if not row[0][0]=='#':
                    mnx = self._checkMNXMdeprecated(row[1])
                    if len(row[0].split(':'))==1:
                        dbName = 'mnx'
                        dbId = row[0]
                    else:
                        dbName = row[0].split(':')[0]
                        dbId = ''.join(row[0].split(':')[1:])
                        if dbName=='deprecated':
                            dbName = 'mnx'
                    #mnx
                    if not mnx in chemXref:
                        chemXref[mnx] = {}
                    if not dbName in chemXref[mnx]:
                        chemXref[mnx][dbName] = []
                    if not dbId in chemXref[mnx][dbName]:
                        chemXref[mnx][dbName].append(dbId)
                    ### DB ###
                    if not dbName in chemXref:
                        chemXref[dbName] = {}
                    if not dbId in chemXref[dbName]:
                        chemXref[dbName][dbId] = mnx
        return chemXref


    ## Function to parse the chem_xref.tsv file of MetanetX
    #
    #  Generate a dictionnary of all cross references for a given chemical id (MNX) to other database id's
    #
    #  @param self Object pointer
    #  @param chem_xref_path Input file path
    #  @return a The dictionnary of identifiers
    #TODO: save the self.deprecatedMNXM_mnxm to be used in case there rp_paths uses an old version of MNX
#    def _m_chebi_mnxm(self, chemXref):
    def _m_chebi_mnxm(self):
        chebi_mnxm = {}
        for mnxm in self.chemXref:
            if 'chebi' in self.chemXref[mnxm]:
                for c in self.chemXref[mnxm]['chebi']:
                    chebi_mnxm[c] = mnxm
        return chebi_mnxm


    ## Function to parse the rules_rall.tsv from RetroRules
    #
    #  Extract from the reactions rules the ruleID, the reactionID, the direction of the rule directed to the origin reaction
    #
    #  @param self The object pointer.
    #  @param path The input file path.
    #  @return rule Dictionnary describing each reaction rule
    def _m_rr_reactions(self, rules_rall_path):
        rr_reactions = {}
        try:
            #with open(rules_rall_path, 'r') as f:
            #    reader = csv.reader(f, delimiter = '\t')
            #    next(reader)
            #    rule = {}
            #    for row in reader:
            for row in csv.DictReader(open(rules_rall_path), delimiter='\t'):
                #NOTE: as of now all the rules are generated using MNX
                #but it may be that other db are used, we are handling this case
                #WARNING: can have multiple products so need to seperate them
                products = {}
                for i in row['Product_IDs'].split('.'):
                    mnxm = self._checkMNXMdeprecated(i)
                    if not mnxm in products:
                        products[mnxm] = 1
                    else:
                        products[mnxm] += 1
                try:
                    #WARNING: one reaction rule can have multiple reactions associated with them
                    #To change when you can set subpaths from the mutliple numbers of
                    #we assume that the reaction rule has multiple unique reactions associated
                    if row['# Rule_ID'] not in rr_reactions:
                        rr_reactions[row['# Rule_ID']] = {}
                    if row['# Rule_ID'] in rr_reactions[row['# Rule_ID']]:
                        self.logger.warning('There is already reaction '+str(row['# Rule_ID'])+' in reaction rule '+str(row['# Rule_ID']))
                    rr_reactions[row['# Rule_ID']][row['Reaction_ID']] = {'rule_id': row['# Rule_ID'], 'rule_score': float(row['Score_normalized']), 'reac_id': self._checkMNXRdeprecated(row['Reaction_ID']), 'subs_id': self._checkMNXMdeprecated(row['Substrate_ID']), 'rel_direction': int(row['Rule_relative_direction']), 'left': {self._checkMNXMdeprecated(row['Substrate_ID']): 1}, 'right': products}
                except ValueError:
                    self.logger.error('Problem converting rel_direction: '+str(row['Rule_relative_direction']))
                    self.logger.error('Problem converting rule_score: '+str(row['Score_normalized']))
            return rr_reactions
        except FileNotFoundError as e:
                self.logger.error('Could not read the rules_rall file ('+str(rules_rall_path)+')')
                return {}


    def _m_inchikey_mnxm(self):
        inchikey_mnxm = {}
        for mnxm in self.mnxm_strc:
            inchikey = self.mnxm_strc[mnxm]['inchikey']
            if not inchikey: inchikey = 'NO_INCHIKEY'
            if not inchikey in inchikey_mnxm:
                inchikey_mnxm[inchikey] = []
            inchikey_mnxm[inchikey].append(mnxm)
        return inchikey_mnxm

    # rpReader
    ## Function to parse the compXref.tsv file of MetanetX
    #
    #  Generate a dictionnary of compartments id's (MNX) to other database id's
    #
    #  @param self Object pointer
    #  @param chem_xref_path Input file path
    #  @return a The dictionnary of identifiers
    #TODO: save the self.deprecatedMNXM_mnxm to be used in case there rp_paths uses an old version of MNX
    def _m_compXref(self, compXref_path):
        compXref = {}
        name_compXref = {}
        try:
            with open(compXref_path) as f:
                c = csv.reader(f, delimiter='\t')
                #not_recognised = []
                for row in c:
                    #cid = row[0].split(':')
                    if not row[0][0]=='#':
                        #collect the info
                        mnxc = row[1]
                        if len(row[0].split(':'))==1:
                            dbName = 'mnx'
                            dbCompId = row[0]
                        else:
                            dbName = row[0].split(':')[0]
                            dbCompId = ''.join(row[0].split(':')[1:])
                            dbCompId = dbCompId.lower()
                        if dbName=='deprecated':
                            dbName = 'mnx'
                        #create the dicts
                        if not mnxc in compXref:
                            compXref[mnxc] = {}
                        if not dbName in compXref[mnxc]:
                            compXref[mnxc][dbName] = []
                        if not dbCompId in compXref[mnxc][dbName]:
                            compXref[mnxc][dbName].append(dbCompId)
                        #create the reverse dict
                        if not dbCompId in name_compXref:
                            name_compXref[dbCompId] = mnxc
        except FileNotFoundError:
            self.logger.error('compXref file not found')
            return {}
        return compXref,name_compXref


def add_arguments(parser):
    parser.add_argument('-sm', '--store_mode', type=str, default='file',
                        help='data storage mode: file or db')
    parser.add_argument('-p', '--print', type=bool, default=False,
                        help='print additional informations')
    return parser

def build_parser():
    return add_arguments(argparse_ArgumentParser('Python script to pre-compute data'))

def entrypoint(params=sys.argv[1:]):
    parser = build_parser()

    args = parser.parse_args(params)

    rpcache = rpCache(args.store_mode, args.print)

##
#
#
if __name__ == "__main__":

    entrypoint()
