from rdkit import Chem

import os

class Filter(object):

    def __init__(self, *args, **kwargs):
        self.poised_filters = {'Amides(1)': ['[#7:1][C;x0:2]=[O:3]>>[#7:1].Cl[C:2]=[O:3]'],
 'Benzimidazole(11)': ['[#7:1]1[#6:9][#7:2][#6:4]2[#6:3]1[#6:8][#6:7][#6:6][#6:5]2>>[#7:1][c:3]1[c:8][c:7][c:6][c:5][c:4]1[#7:2].Cl[#6:9]=O'],
 'Benzoxazole(13)': ['[#7:1]1[#6:9][#8:2][#6:4]2[#6:3]1[#6:8][#6:7][#6:6][#6:5]2>>[#7:1][c:3]1[c:8][c:7][c:6][c:5][c:4]1[#8:2].Cl[#6:9]=O'],
 'Ester_Coupling(10)': ['[#8:1][C;x0:2]=[O:3]>>[#8:1].Cl[C:2]=[O:3]'],
 'Ether_Coupling(9)1': ['[CH0:2]-[#8R0:1]>>[#8:1].[#6:2]Br',
  '[CH1:2]-[#8R0:1]>>[#8:1].[#6:2]Br',
  '[CH2R0:2]-[#8:1]>>[#8:1].[#6:2]Br',
  '[CH1R0:2]-[#8:1]>>[#8:1].[#6:2]Br',
  '[CH2:2]-[#8R0:1]>>[#8:1].[#6:2]Br',
  '[CH3:2]-[#8R0:1]>>[#8:1].[#6:2]Br',
  '[CH0R0:2]-[#8:1]>>[#8:1].[#6:2]Br',
  '[n:3][c:2]-[#8:1]>>[#8:1].[n:3][#6:2]Br'],
 'Indole(14)': ['[c:10]1[c:9][nH:1][c:3]2[c:8]1[c:7][c:6][c:5][c:4]2>>N[N:1][c:3]1[c:8]([2H])[c:7][c:6][c:5][c:4]1.[C:10][C:9]=O'],
 'N-Alkylation(8)1': ['[CH0:2]-[#7R0:1]>>[#7:1].[#6:2]Br',
  '[CH1:2]-[#7R0:1]>>[#7:1].[#6:2]Br',
  '[CH2R0:2]-[#7:1]>>[#7:1].[#6:2]Br',
  '[CH1R0:2]-[#7:1]>>[#7:1].[#6:2]Br',
  '[CH2:2]-[#7R0:1]>>[#7:1].[#6:2]Br',
  '[CH3:2]-[#7R0:1]>>[#7:1].[#6:2]Br',
  '[CH0R0:2]-[#7:1]>>[#7:1].[#6:2]Br'],
 'Oxadiazole(15)': ['[#6:6][c:4]1[n:5][o:3][c:1][n:2]1>>[O:2]-[C:1]=[O:3].[#6:6][C:4]#[N:5]'],
 'Reductive_Amination(7)1': ['[CH1:2]-[#7R0:1]>>[#7:1].[#6:2]=O',
  '[CH2R0:2]-[#7:1]>>[#7:1].[#6:2]=O',
  '[CH1R0:2]-[#7:1]>>[#7:1].[#6:2]=O',
  '[CH2:2]-[#7R0:1]>>[#7:1].[#6:2]=O',
  '[CH3:2]-[#7R0:1]>>[#7:1].[#6:2]=O'],
 'SNAr(2)': ['[c:2][N:1][#6:3]>>[#6:3]-[#7:1].[c:2]Br'],
 'Sonogashira(5)': ['[#6;a:1][C:2]#[C:3]>>[#6;a:1]Br.[C:2]#[C:3]'],
 'Sulfonamide(6)': ['[#7:1][S:2](=[O:3])=[O:4]>>[#7:1].Cl[S:2](=[O:3])=[O:4]'],
 'Suzuki_Coupling(4)': ['[#6;a:1]-[#6;a:2]>>[#6;a:2]Br.[#6;a:1]-[#5](-[#8])-[#8]'],
 'Triazole(12)': ['[#6:6][n:1]1[c:3][c:4][n:2][n:5]1>>[#6:6][#7:1][#7:5]=[#7:2].[C:3]#[C:4]'],
 'Urea(3)': ['[#7:1][C;x0]([#7:2])=O>>[#7:1].[#7:2]']}
        self.starts = {}
        for key in self.poised_filters:
            self.starts[key] = [Chem.MolFromSmarts(x.split(">>")[0]) for x in self.poised_filters[key]]

    def pass_filter(self, mol):
        """
        Filter a given molecule given a series of SMARTS patterns.
        :param mol: an input RDKit Molecule
        :param patterns: a dict of SMARTS patterns
        :return: a list of the reactions each mol can do.
        """

        out_list = []
        for key in self.starts:
            for patt in self.starts[key]:
                if mol.HasSubstructMatch(patt):
                    out_list.append(key)
        return out_list

    def get_writers(self, dir_base):
        """
        Get all the writers of the SD files
        :param output_path:
        :return:
        """
        out_d = {}
        for x in self.starts:
            out_d[x] = Chem.SDWriter(os.path.join(dir_base, x))
        return out_d
