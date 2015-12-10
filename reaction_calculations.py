from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from itertools import product
from itertools import permutations
from copy import copy
import pybel
import numpy as np

def spam(n):
    out = []
    for perm in getPerms(n):
        elem = [int(i) for i in list(perm)]
        out.append(elem)
    return out

def getPerms(n):
    for i in getCandidates(n):
        for perm in set(permutations(i)):
            yield ''.join(perm)

def getCandidates(n):
    for i in range(0, n+1):
        res = "1"*i+"0" *(n-i)
        yield res

# Adapted code from
# https://github.com/rdkit/rdkit/issues/626
def GetStereoIsomers(mol, maxNum=12):
    out = []

    chiralCenters = Chem.FindMolChiralCenters(mol, includeUnassigned=True)

    # keep only unassigned chiral centers
    chiralCenters = [c for c in chiralCenters if c[1] == "?"]

    # return the molecule object if no unassigned centers were found
    if chiralCenters == []:
        return [mol]

    # All bit permutations with number of bits equals numbers of chiralCenters
    elements = spam(len(chiralCenters))

    for isoId, element in enumerate(elements):
        for centerId, i in enumerate(element):
            atomId = chiralCenters[centerId][0]
            if i == 0:
                mol.GetAtomWithIdx(atomId).SetChiralTag(Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW)
            elif i == 1:
                mol.GetAtomWithIdx(atomId).SetChiralTag(Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW)
        outmol = copy(mol)
        out.append(outmol)
        if len(out) >= maxNum:
            break
    return out


# Iterates over all values in list of lists
# From http://stackoverflow.com/questions/6340351/python-iterating-through-list-of-list
def traverse(o, tree_types=(list, tuple)):
    if isinstance(o, tree_types):
        for value in o:
            for subvalue in traverse(value, tree_types):
                yield subvalue
    else:
        yield o

def convert_sugar_forms(molecule):
    rxn1 = '[O:1]1[C:2]([C:8])[C:3][C:4][C:5][C:6]1[O:7][H]>>[H][O:1][C:2]([C:8])[C:3][C:4][C:5][C:6]=[O:7]'
    rxn2 = '[O:6]1[C:2]([O:7][H])([C:1])[C:3][C:4][C:5]1>>[C:1][C:2](=[O:7])[C:3][C:4][C:5][O:6][H]'
    rxn3 = '[O:1]1[C:2]([H:8])[C:3][C:4][C:5][C:6]1[O:7][H]>>[H][O:1][C:2]([H:8])[C:3][C:4][C:5][C:6]=[O:7]'
    rxn4 = '[O:6]1[C:2]([O:7][H])([H:1])[C:3][C:4][C:5]1>>[H:1][C:2](=[O:7])[C:3][C:4][C:5][O:6][H]'
    sugarrxns = [AllChem.ReactionFromSmarts(rxn1), 
                 AllChem.ReactionFromSmarts(rxn2),
                 AllChem.ReactionFromSmarts(rxn3),
                 AllChem.ReactionFromSmarts(rxn4)]
    rxnproducts = [molecule]
    seen = set()
    seen.add(Chem.MolToSmiles(Chem.RemoveHs(molecule), isomericSmiles=True))    
    for sugarrxn in sugarrxns:
        prods = sugarrxn.RunReactants((molecule,))
        for p in traverse(prods):
            
            smilesprod = Chem.MolToSmiles(p, isomericSmiles=True)
            if smilesprod not in seen:
                pmol = Chem.MolFromSmiles(smilesprod)
                seen.add(smilesprod)
                rxnproducts.append(pmol)
    return rxnproducts

def getMaxTC(fps1, fps2):
    if fps1 == [] or fps2 == []:
        return 0.
    tcs = []
    for pairfps in product(fps1, fps2):
        tc = DataStructs.TanimotoSimilarity(pairfps[0], pairfps[1])
        if tc == 1.0:
            return tc
        tcs.append(tc)
    return np.max(tcs)

class ligandData(object):
    
    def __init__(self, smilesfile=None, smirksfile=None, reactions={}, smiles={},
                 molecules={}, delimiter=None, linearforms=True):
        self.smiles = smiles
        self.molecules = molecules
        self.reactions = reactions
        self.fingerprints = None
        self.linearforms = linearforms
        if smilesfile is not None:
            self.read_in_molecules(smilesfile=smilesfile,
                                   delimiter=delimiter,
                                   linearforms=linearforms)
            print 'Number of molecules: %d' % len(self.molecules)
            self.make_fingerprints()
        if smirksfile is not None:
            self.read_in_smirks(smirksfile)

        
    def read_in_smiles(self, smilesfilename, delimiter=None):
        with open(smilesfilename, 'r') as handle:
            smiles = {}
            lines = handle.readlines()
            for line in lines:
                if delimiter is not None:
                    fields = line.strip().split(delimiter)
                else:
                    fields = line.strip().split()
                fields = [f.strip() for f in fields]
                smiles_string = pybel.readstring('smi', fields[0]).write('can')
                smiles[fields[1]] = smiles_string
            self.smiles = smiles


    def read_in_molecules(self, smilesfile, delimiter=None, linearforms=True):
        with open(smilesfile, 'r') as handle:
            molecules = {}
            lines = handle.readlines()
            smilesdict = {}
            for line in lines:
                if delimiter is not None:
                    fields = line.strip().split(delimiter)
                else:
                    fields = line.strip().split()
                smiles = fields[0].strip()
                mol = Chem.MolFromSmiles(smiles)

                #if linearforms:
                #    # Convert to linear forms
                #    identifier = fields[1].strip()
                #    linearmols = convert_sugar_conversions(mol)
                #else:
                #    # Preserve form from database
                #    identifier = fields[1].split('_')[0]
                #    linearmols = [mol]

                #if identifier not in molecules.keys():
                #    if len(linearmols) > 1:
                #        molecules[identifier] = linearmols[1]
                #        smiles = Chem.MolToSmiles(linearmols[1], isomericSmiles=True)
                #    else:
                #        molecules[identifier] = mol
                #    smiles = pybel.readstring('smi', smiles).write('can')
                #    smilesdict[identifier] = smiles
                identifier = fields[1].strip()
                if identifier not in molecules.keys():
                    molecules[identifier] = mol
                    smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
                    smilesdict[identifier] = smiles
            self.molecules = molecules
            self.smiles = smilesdict

    def read_in_smirks(self, smirksfilename):
        with open(smirksfilename, 'r') as handle:
            reactions = {}
            lines = handle.readlines()
            for line in lines:
                fields = lines.strip().split()
                reactions[fields[1]] = fields[2]
            self.reactions = reactions

    def make_fingerprints(self):
        fingerprints = {}
        for k in self.smiles.keys():
            s = self.smiles[k]
            rdkmol = Chem.MolFromSmiles(s)
            fp = AllChem.GetMorganFingerprintAsBitVect(rdkmol, 2, useChirality=True)
            fingerprints[k] = fp
        self.fingerprints = fingerprints

    def compute_match(self):
        idxs = []
        dictlist = []
        for f1 in self.fingerprints.keys():
            ds = {}
            for f2 in self.fingerprints.keys():
                fpA = self.fingerprints[f1]
                fpB = self.fingerprints[f2]
                sim = DataStructs.TanimotoSimilarity(fpA, fpB)
                ds[f2] = sim
            idxs.append(f1)
            dictlist.append(ds)
        return idxs, dictlist

    def reaction_sets(self, reactionkey):
        idxs = []
        count = 0
        edges = []

        isofpdict = {}
        for j in self.molecules.keys():
            m =  self.molecules[j]
            m = Chem.AddHs(m)
            if self.linearforms:
                ms = convert_sugar_forms(m)
            else:
                ms = [m]
            refs = []
            for refmol in ms:
                refs.extend(GetStereoIsomers(refmol))
            refs = [Chem.RemoveHs(rdkmol) for rdkmol in refs]
            refs = [Chem.MolFromSmiles(Chem.MolToSmiles(ref, isomericSmiles=True)) for ref in refs]
            fps = [AllChem.GetMorganFingerprintAsBitVect(rdkmol, 2, useChirality=True) for rdkmol in refs]
            isofpdict[j] = fps

        for k in self.molecules.keys():
            count += 1
            idxs.append(k)
            fps = []
            m = self.molecules[k]
            reactantmol = Chem.AddHs(m)
            if self.linearforms:
                rs = convert_sugar_forms(reactantmol)
            else:
                rs = [reactantmol]

            for reaction in self.reactions[reactionkey]:
                rdk_rxn = AllChem.ReactionFromSmarts(reaction)
                
                for r in rs:
                    rxn_products = rdk_rxn.RunReactants((r,))
                    
                    for rxn_p in traverse(rxn_products):
                        rxn_p_noH = Chem.RemoveHs(rxn_p)
                        pstereos = GetStereoIsomers(rxn_p_noH)
                        smis = [Chem.MolToSmiles(p, isomericSmiles=True) for p in pstereos]
                        pmols = [Chem.MolFromSmiles(s) for s in smis]
                        cfps = [AllChem.GetMorganFingerprintAsBitVect(pm, 2, useChirality=True) for pm in pmols]
                        fps.extend(cfps)
                
            for j in isofpdict.keys():
                score = getMaxTC(isofpdict[j], fps)
                if score == 1.0:
                    edges.append([k, j])

        return edges

    def write_possible_reaction_sets(self, textfile):
        with open(textfile, 'w') as handle:
            for key in self.reactions.keys():
                print key
                edges = self.reaction_sets(key)
                for edge in edges:
                    handle.write('%s, %s, %s\n' % (key, edge[0], edge[1]))


