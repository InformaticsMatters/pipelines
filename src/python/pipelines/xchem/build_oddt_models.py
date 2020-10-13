#!/usr/bin/env python

# Copyright 2020 Informatics Matters Ltd.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Run this to generate the RFScore and NNScore models.
The files RFScore_v1_pdbbind2016.pickle, RFScore_v2_pdbbind2016.pickle, RFScore_v3_pdbbind2016.pickle and
NNScore_pdbbind2016.pickle are generated.
If you want them to be re-generated they must first be deleted.
"""


from oddt.virtualscreening import virtualscreening as vs

ligands = '../../data/mpro/hits-17.sdf.gz'
protein = '../../data/mpro/Mpro-x0387_0.pdb'
pipeline=vs()
print('Loading')
pipeline.load_ligands('sdf', ligands)
print('Scoring with rfscore')
pipeline.score(function='rfscore_v1', protein=protein)
pipeline.score(function='rfscore_v2', protein=protein)
pipeline.score(function='rfscore_v3', protein=protein)
print('Scoring with nnscore')
pipeline.score(function='nnscore', protein=protein)
# print('Scoring with plecscore')
# pipeline.score(function='pleclinear', protein=protein)
# pipeline.score(function='plecnn', protein=protein)
# pipeline.score(function='plecrf', protein=protein)
print('Writing')
pipeline.write('sdf', 'scored.sdf')
print('Done')