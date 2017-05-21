# $Id: ShowFeats.py 537 2007-08-20 14:54:35Z landrgr1 $
#
# Created by Greg Landrum Aug 2006
#
#
from __future__ import print_function
# # --- ----  --- ----  --- ----  --- ----  --- ----  --- ----
import sys, os, getopt
from rdkit import RDConfig
from rdkit import Chem
from rdkit.Chem import AllChem

_version = "0.3.2"

_usage = """
   ShowFeats [optional args] <filenames>
  if "-" is provided as a filename, data will be read from stdin (the console)
"""

_welcomeMessage = "This is ShowFeats version %s" % (_version)

import math
#set up the logger:
from rdkit import RDLogger as logging
logger = logging.logger()
logger.setLevel(logging.INFO)

from rdkit import Geometry
from rdkit.Chem.Features import FeatDirUtilsRD as FeatDirUtils

_featColors = {
  'Donor': (0, 1, 1),
  'Acceptor': (1, 0, 1),
  'NegIonizable': (1, 0, 0),
  'PosIonizable': (0, 0, 1),
  'ZnBinder': (1, .5, .5),
  'Aromatic': (1, .8, .2),
  'LumpedHydrophobe': (.5, .25, 0),
  'Hydrophobe': (.5, .25, 0),
}


def _getVectNormal(v, tol=1e-4):
  if math.fabs(v.x) > tol:
    res = Geometry.Point3D(v.y, -v.x, 0)
  elif math.fabs(v.y) > tol:
    res = Geometry.Point3D(-v.y, v.x, 0)
  elif math.fabs(v.z) > tol:
    res = Geometry.Point3D(1, 0, 0)
  else:
    raise ValueError('cannot find normal to the null vector')
  res.Normalize()
  return res


_canonArrowhead = None
viewer = None

def _buildCanonArrowhead(headFrac, nSteps, aspect):
  global _canonArrowhead
  startP = RDGeometry.Point3D(0, 0, headFrac)
  _canonArrowhead = [startP]

  scale = headFrac * aspect
  baseV = RDGeometry.Point3D(scale, 0, 0)
  _canonArrowhead.append(baseV)

  twopi = 2 * math.pi
  for i in range(1, nSteps):
    v = RDGeometry.Point3D(scale * math.cos(i * twopi), scale * math.sin(i * twopi), 0)
    _canonArrowhead.append(v)


_globalArrowCGO = []
_globalSphereCGO = []
# taken from pymol's cgo.py
BEGIN = 2
END = 3
TRIANGLE_FAN = 6
COLOR = 6
VERTEX = 4
NORMAL = 5
SPHERE = 7
CYLINDER = 9
ALPHA = 25



def _cgoArrowhead(viewer, tail, head, radius, color, label, headFrac=0.3, nSteps=10, aspect=.5):
  global _globalArrowCGO
  delta = head - tail
  normal = _getVectNormal(delta)
  delta.Normalize()

  dv = head - tail
  dv.Normalize()
  dv *= headFrac
  startP = head

  normal *= headFrac * aspect

  cgo = [BEGIN, TRIANGLE_FAN, COLOR, color[0], color[1], color[2], NORMAL, dv.x, dv.y, dv.z, VERTEX,
         head.x + dv.x, head.y + dv.y, head.z + dv.z]
  base = [BEGIN, TRIANGLE_FAN, COLOR, color[0], color[1], color[2], NORMAL, -dv.x, -dv.y, -dv.z,
          VERTEX, head.x, head.y, head.z]
  v = startP + normal
  cgo.extend([NORMAL, normal.x, normal.y, normal.z])
  cgo.extend([VERTEX, v.x, v.y, v.z])
  base.extend([VERTEX, v.x, v.y, v.z])
  for i in range(1, nSteps):
    v = FeatDirUtils.ArbAxisRotation(360. / nSteps * i, delta, normal)
    cgo.extend([NORMAL, v.x, v.y, v.z])
    v += startP
    cgo.extend([VERTEX, v.x, v.y, v.z])
    base.extend([VERTEX, v.x, v.y, v.z])

  cgo.extend([NORMAL, normal.x, normal.y, normal.z])
  cgo.extend([VERTEX, startP.x + normal.x, startP.y + normal.y, startP.z + normal.z])
  base.extend([VERTEX, startP.x + normal.x, startP.y + normal.y, startP.z + normal.z])
  cgo.append(END)
  base.append(END)
  cgo.extend(base)

  #viewer.server.renderCGO(cgo,label)
  _globalArrowCGO.extend(cgo)


def ShowArrow(viewer, tail, head, radius, color, label, transparency=0, includeArrowhead=True):
  global _globalArrowCGO
  if transparency:
    _globalArrowCGO.extend([ALPHA, 1 - transparency])
  else:
    _globalArrowCGO.extend([ALPHA, 1])
  _globalArrowCGO.extend([CYLINDER,
                          tail.x,
                          tail.y,
                          tail.z,
                          head.x,
                          head.y,
                          head.z,
                          radius * .10,
                          color[0],
                          color[1],
                          color[2],
                          color[0],
                          color[1],
                          color[2], ])

  if includeArrowhead:
    _cgoArrowhead(viewer, tail, head, radius, color, label)


def ShowMolFeats(mol, factory, radius=0.5, confId=-1, showOnly=True, name='',
                 transparency=0.0, colors=None, excludeTypes=[], useFeatDirs=True, featLabel=None,
                 dirLabel=None, includeArrowheads=True, writeFeats=False, showMol=True,
                 featMapFile=False):
  global _globalSphereCGO
  out_list = []
  import json
  if not name:
    if mol.HasProp('_Name'):
      name = mol.GetProp('_Name')
    else:
      name = 'molecule'
  if not colors:
    colors = _featColors

  molFeats = factory.GetFeaturesForMol(mol)
  if not featLabel:
    featLabel = '%s-feats' % name
    # viewer.server.resetCGO(featLabel)
  if not dirLabel:
    dirLabel = featLabel + "-dirs"
    # viewer.server.resetCGO(dirLabel)

  for i, feat in enumerate(molFeats):
    family = feat.GetFamily()
    if family in excludeTypes:
      continue
    pos = feat.GetPos(confId)
    color = colors.get(family, (.5, .5, .5))
    nm = '%s(%d)' % (family, i + 1)

    if transparency:
      _globalSphereCGO.extend([ALPHA, 1 - transparency])
    else:
      _globalSphereCGO.extend([ALPHA, 1])
    _globalSphereCGO.extend([COLOR, color[0], color[1], color[2], SPHERE, pos.x, pos.y, pos.z,
                             radius])

    aidText = ''
    if writeFeats:
      aidText = ' '.join([str(x + 1) for x in feat.GetAtomIds()])
      print('%s\t%.3f\t%.3f\t%.3f\t1.0\t# %s' % (family, pos.x, pos.y, pos.z, aidText))


    if featMapFile:
      print("  family=%s pos=(%.3f,%.3f,%.3f) weight=1.0" % (family, pos.x, pos.y, pos.z), end='',
            file=featMapFile)

    if useFeatDirs:
      ps = []
      if family == 'Aromatic':
        ps, fType = FeatDirUtils.GetAromaticFeatVects(
          mol.GetConformer(confId), feat.GetAtomIds(), pos, scale=1.0)
      elif family == 'Donor':
        aids = feat.GetAtomIds()
        if len(aids) == 1:
          featAtom = mol.GetAtomWithIdx(aids[0])
          hvyNbrs = [x for x in featAtom.GetNeighbors() if x.GetAtomicNum() != 1]
          if len(hvyNbrs) == 1:
            ps, fType = FeatDirUtils.GetDonor1FeatVects(mol.GetConformer(confId), aids, scale=1.0)
          elif len(hvyNbrs) == 2:
            ps, fType = FeatDirUtils.GetDonor2FeatVects(mol.GetConformer(confId), aids, scale=1.0)
          elif len(hvyNbrs) == 3:
            ps, fType = FeatDirUtils.GetDonor3FeatVects(mol.GetConformer(confId), aids, scale=1.0)
      elif family == 'Acceptor':
        aids = feat.GetAtomIds()
        if len(aids) == 1:
          featAtom = mol.GetAtomWithIdx(aids[0])
          hvyNbrs = [x for x in featAtom.GetNeighbors() if x.GetAtomicNum() != 1]
          if len(hvyNbrs) == 1:
            ps, fType = FeatDirUtils.GetAcceptor1FeatVects(
              mol.GetConformer(confId), aids, scale=1.0)
          elif len(hvyNbrs) == 2:
            ps, fType = FeatDirUtils.GetAcceptor2FeatVects(
              mol.GetConformer(confId), aids, scale=1.0)
          elif len(hvyNbrs) == 3:
            ps, fType = FeatDirUtils.GetAcceptor3FeatVects(
              mol.GetConformer(confId), aids, scale=1.0)
      this_l = [family, pos.x, pos.y, pos.z]
      for tail, head in ps:
        ShowArrow(viewer, tail, head, radius, color, dirLabel, transparency=transparency,
                  includeArrowhead=includeArrowheads)
        vect = head - tail
        print('dir=(%.3f,%.3f,%.3f)' % (vect.x, vect.y, vect.z),)

        if vect:
          this_l.extend([vect.x, vect.y, vect.z])
    out_list.append(this_l)
    if featMapFile:
      aidText = ' '.join([str(x + 1) for x in feat.GetAtomIds()])
      print('# %s' % (aidText), file=featMapFile)
  return out_list


def get_feats(molecule):
  fdef = open(os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef'), 'r').read()
  factory = AllChem.BuildFeatureFactoryFromString(fdef)
  return ShowMolFeats(molecule, factory, transparency=0.25,
                   showOnly=False, includeArrowheads=True,
                   writeFeats=True, showMol=True)



if __name__ == '__main__':

  mol = Chem.MolFromSmiles('Cc1sc2c(c1C)C(c1ccc(Cl)cc1)=NC(CC(=O)OC(C)(C)C)c1[nH]nc(C)[n+]1-2')
  AllChem.EmbedMolecule(mol)
  print(get_feats(mol))
