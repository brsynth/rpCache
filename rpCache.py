import os
from rdkit.Chem import MolFromSmiles, MolFromInchi, MolToSmiles, MolToInchi, MolToInchiKey, AddHs
import csv
import logging
import pickle
import gzip
import urllib.request
import re
import tarfile
import shutil
import redis
from sys import stdout as sys_stdout
from sys import getsizeof as sys_getsizeof
import time

#######################################################
################### rpCache  ##########################
#######################################################

from sys import getsizeof, stderr
from itertools import chain
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
    def __init__(self, db='file', print=False):
        self.store_mode = db
        if self.store_mode!='file':
            self.redis = redis.StrictRedis(host=self.store_mode, port=6379, db=0)
        #given by Thomas
        self.logger = logging.getLogger(__name__)
        self.logger.info('Started instance of rpCache')

        self.print = print_size

        # Common attribues
        self.convertMNXM = {'MNXM162231': 'MNXM6',
                            'MNXM84': 'MNXM15',
                            'MNXM96410': 'MNXM14',
                            'MNXM114062': 'MNXM3',
                            'MNXM145523': 'MNXM57',
                            'MNXM57425': 'MNXM9',
                            'MNXM137': 'MNXM588022'}
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

        dirname = os.path.dirname(os.path.abspath( __file__ ))

        #################### make the local folders ############################
        # input_cache
        if not os.path.isdir(dirname+'/input_cache'):
            os.mkdir(dirname+'/input_cache')
        # cache
        if not os.path.isdir(dirname+'/cache'):
            os.mkdir(dirname+'/cache')

        # ###################### Fetch the files if necessary ######################
        url = 'https://www.metanetx.org/cgi-bin/mnxget/mnxref/'

        # 3xCommon + rpReader
        for file in ['reac_xref.tsv', 'chem_xref.tsv', 'chem_prop.tsv', 'comp_xref.tsv']:
            if not os.path.isfile(dirname+'/input_cache/'+file) or fetchInputFiles:
                print("Downloading "+file+"...", end = '', flush=True)
                start = time.time()
                urllib.request.urlretrieve(url+file, dirname+'/input_cache/'+file)
                end = time.time()
                print(" (%.2fs)" % (end - start))


        #TODO: need to add this file to the git or another location
        for file in ['rr_compounds.tsv', 'rxn_recipes.tsv']:
            if not os.path.isfile(dirname+'/input_cache/'+file) or fetchInputFiles:
                print("Downloading "+file+"...", end = '', flush=True)
                start = time.time()
                urllib.request.urlretrieve('https://retrorules.org/dl/this/is/not/a/secret/path/rr02',
                                           dirname+'/input_cache/rr02_more_data.tar.gz')
                end = time.time()
                print(" (%.2fs)" % (end - start))
                tar = tarfile.open(dirname+'/input_cache/rr02_more_data.tar.gz', 'r:gz')
                tar.extractall(dirname+'/input_cache/')
                tar.close()
                shutil.move(dirname+'/input_cache/rr02_more_data/compounds.tsv',
                            dirname+'/input_cache/rr_compounds.tsv')
                shutil.move(dirname+'/input_cache/rr02_more_data/rxn_recipes.tsv',
                            dirname+'/input_cache/')
                os.remove(dirname+'/input_cache/rr02_more_data.tar.gz')
                shutil.rmtree(dirname+'/input_cache/rr02_more_data')

        if not os.path.isfile(dirname+'/input_cache/rules_rall.tsv') or fetchInputFiles:
            print("Downloading rules_rall.tsv...", end = '', flush=True)
            start = time.time()
            urllib.request.urlretrieve('https://retrorules.org/dl/preparsed/rr02/rp3/hs',
                                       dirname+'/input_cache/retrorules_rr02_rp3_hs.tar.gz')
            end = time.time()
            print(" (%.2fs)" % (end - start))
            tar = tarfile.open(dirname+'/input_cache/retrorules_rr02_rp3_hs.tar.gz', 'r:gz')
            tar.extractall(dirname+'/input_cache/')
            tar.close()
            shutil.move(dirname+'/input_cache/retrorules_rr02_rp3_hs/retrorules_rr02_flat_all.tsv', dirname+'/input_cache/rules_rall.tsv')
            os.remove(dirname+'/input_cache/retrorules_rr02_rp3_hs.tar.gz')
            shutil.rmtree(dirname+'/input_cache/retrorules_rr02_rp3_hs')

        ###################### Populate the cache #################################
        input_cache = dirname+'/input_cache'
        attributes = {
            'deprecatedMNXM_mnxm': [input_cache+'/chem_xref.tsv'],
            'deprecatedMNXR_mnxr': [input_cache+'/reac_xref.tsv'],
            'mnxm_strc': [input_cache+'/rr_compounds.tsv', input_cache+'/chem_prop.tsv'],
            'chemXref': [input_cache+'/chem_xref.tsv'],
            'chebi_mnxm': [],
            'rr_reactions': [input_cache+'/rules_rall.tsv'],
            'inchikey_mnxm': []
        }
        for attribute in attributes:
            start = time.time()
            if not getattr(self, attribute):
                self._processAttribute(attribute, dirname, attributes[attribute])
                if self.print:
                    print(" ("+str(total_size(getattr(self,attribute)))+" bytes)", end = '', flush=True)
            end = time.time()
            if self.print:
                print(" (%.2fs)" % (end - start))

        attributes = {
            'compXref': [input_cache+'/comp_xref.tsv']
        }
        for attribute in attributes:
            start = time.time()
            if not getattr(self, attribute):
                self._processAttribute2([attribute, 'name_'+attribute], dirname, attributes[attribute])
                if self.print:
                    print(" ("+str(total_size(getattr(self,attribute)))+" "+str(total_size(getattr(self,'name_'+attribute)))+" bytes)", end = '', flush=True)
            end = time.time()
            if self.print:
                print(" (%.2fs)" % (end - start))

        return True


    def _processAttribute(self, attribute, dirname, args):
        try:
            # Check if attribute is set
            if self.checkAttribute(attribute, dirname):
                print("Loading "+attribute+" from "+self.store_mode+"...", end = '', flush=True)
                result = self.loadAttribute(attribute, dirname)
            else:
                print("Generating "+attribute+" to "+self.store_mode+"...", end = '', flush=True)
                # Choose method according to attribute name
                method = getattr(self, '_m_'+attribute)
                # Apply method and expand 'args' list as arguments
                result = method(*args)
                # Store pickle
                self.storeAttribute(attribute, result, dirname)
            sys_stdout.write("\033[0;32m") # Green
            print(" OK", end = '', flush=True)
            sys_stdout.write("\033[0;0m") # Reset
        except:
            sys_stdout.write("\033[1;31m") # Red
            print(" Failed")
            sys_stdout.write("\033[0;0m") # Reset
            raise
        # Set attribute to value
        setattr(self, attribute, result)


    # Process with two outputs method
    def _processAttribute2(self, attributes, dirname, args):
        try:
            results = []
            # Check if attributes are set
            check = True
            i=0
            while check and i<len(attributes):
                check = check and self.checkAttribute(attributes[i], dirname)
                i+=1
            if check:
                print("Loading "+" ".join(attributes)+" from "+self.store_mode+"...", end = '', flush=True)
                for i in range(len(attributes)):
                    results += [self.loadAttribute(attributes[i], dirname)]
            else:
                print("Generating "+" ".join(attributes)+" to "+self.store_mode+"...", end = '', flush=True)
                # Choose method according to attribute name
                method = getattr(self, '_m_'+attributes[0])
                # Apply method and expand 'args' list as arguments
                # Put results in a list
                results = method(*args)
                for i in range(len(results)):
                    # Store pickle
                    self.storeAttribute(attributes[i], results[i], dirname)
            sys_stdout.write("\033[0;32m") # Green
            print(" OK", end = '', flush=True)
            sys_stdout.write("\033[0;0m") # Reset
        except:
            sys_stdout.write("\033[1;31m") # Red
            print(" Failed")
            sys_stdout.write("\033[0;0m") # Reset
            raise
        # Set attribute to value
        for i in range(len(results)):
            setattr(self, attributes[i], results[i])


    def storePickle(self, pickle_key, pickle_obj, dirname='./', gzip=False):
        if self.store_mode=='redis':
            self._storePickleToDB(pickle_key, pickle_obj)
        else:
            self._storePickleToFile(pickle_key, pickle_obj, dirname, gzip)

    def storeAttribute(self, attribute, data, dirname='./', gzip=False):
        self.storePickle(attribute+'.pickle', pickle.dumps(data), dirname, gzip)

    def _storePickleToFile(self, pickle_key, pickle_obj, dirname, gzip):
        filename = dirname+'/cache/'+pickle_key
        if gzip:
            filename += '.gz'
            with gzip.open(filename, "wb") as f:
            	f.write(pickle_obj)
        else:
            with open(filename, "wb") as f:
            	f.write(pickle_obj)

    def _storePickleToDB(self, pickle_key, pickle_obj):
        self.redis.set(pickle_key, pickle_obj)
# redisClient.hset(hashName, 1, "Cheesecake")


    def loadPickle(self, pickle_key, dirname='./', gz=False):
        if self.store_mode=='redis':
            return self._loadPickleFromDB(pickle_key, gz)
        else:
            return self._loadPickleFromFile(pickle_key, dirname, gz)

    def loadAttribute(self, attribute, dirname='./', gz=False):
        return self.loadPickle(attribute+'.pickle', dirname, gz)

    def _loadPickleFromFile(self, pickle_key, dirname, gz=False):
        filename = dirname+'/cache/'+pickle_key
        if gz:
            filename += '.gz'
            return pickle.load(gzip.open(filename, 'rb'))
        else:
            return pickle.load(open(filename, 'rb'))

    # def _loadPickleFromFile(self, filename, gz=False):
    #     if gzip:
    #         return pickle.load(gzip.open(filename, 'rb'))
    #     else:
    #         return pickle.load(open(filename, 'rb'))

    def _loadPickleFromDB(self, pickle_key, gz=False):
        return pickle.loads(self.redis.get(pickle_key))

    def setPickle(self, attribute, data):
        attribute = data

    def getPickle(self, attribute):
        return attribute

    def checkAttribute(self, attribute, dirname):
        return self.checkPickle(attribute+'.pickle', dirname)

    def checkPickle(self, pickle_key, dirname):
        if self.store_mode=='redis':
            return self.redis.exists(pickle_key)
        else:
            return os.path.isfile(dirname+'/cache/'+pickle_key)


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
        deprecatedMNX_mnx = self._deprecatedMNX(chem_xref_path)
        deprecatedMNX_mnx.update(self.convertMNXM)
        deprecatedMNX_mnx['MNXM01'] = 'MNXM1'
        return deprecatedMNX_mnx

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
        try:
            #with open(rules_rall_path, 'r') as f:
            #    reader = csv.reader(f, delimiter = '\t')
            #    next(reader)
            #    rule = {}
            #    for row in reader:
            rule = {}
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
                    if row['# Rule_ID'] not in rule:
                        rule[row['# Rule_ID']] = {}
                    if row['# Rule_ID'] in rule[row['# Rule_ID']]:
                        self.logger.warning('There is already reaction '+str(row['# Rule_ID'])+' in reaction rule '+str(row['# Rule_ID']))
                    rule[row['# Rule_ID']][row['Reaction_ID']] = {'rule_id': row['# Rule_ID'], 'rule_score': float(row['Score_normalized']), 'reac_id': self._checkMNXRdeprecated(row['Reaction_ID']), 'subs_id': self._checkMNXMdeprecated(row['Substrate_ID']), 'rel_direction': int(row['Rule_relative_direction']), 'left': {self._checkMNXMdeprecated(row['Substrate_ID']): 1}, 'right': products}
                except ValueError:
                    self.logger.error('Problem converting rel_direction: '+str(row['Rule_relative_direction']))
                    self.logger.error('Problem converting rule_score: '+str(row['Score_normalized']))
        except FileNotFoundError as e:
                self.logger.error('Could not read the rules_rall file ('+str(rules_rall_path)+')')
                return {}
        return rule


    def _m_inchikey_mnxm(self):
        inchikey_mnxm = {}
        for mnxm in self.mnxm_strc:
            if not self.mnxm_strc[mnxm]['inchikey'] in inchikey_mnxm:
                inchikey_mnxm[self.mnxm_strc[mnxm]['inchikey']] = []
            inchikey_mnxm[self.mnxm_strc[mnxm]['inchikey']].append(mnxm)
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
        name_pubDB_xref = {}
        compName_mnxc = {}
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
                        if not mnxc in name_pubDB_xref:
                            name_pubDB_xref[mnxc] = {}
                        if not dbName in name_pubDB_xref[mnxc]:
                            name_pubDB_xref[mnxc][dbName] = []
                        if not dbCompId in name_pubDB_xref[mnxc][dbName]:
                            name_pubDB_xref[mnxc][dbName].append(dbCompId)
                        #create the reverse dict
                        if not dbCompId in compName_mnxc:
                            compName_mnxc[dbCompId] = mnxc
        except FileNotFoundError:
            self.logger.error('compXref file not found')
            return {}
        return name_pubDB_xref, compName_mnxc


if __name__ == "__main__":
    rpcache = rpCache('file')
