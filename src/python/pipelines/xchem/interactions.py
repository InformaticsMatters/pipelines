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
Classes and utility functions for calculating interactions using ODDT.
"""

import sys, json, math

I_TYPE_HBOND = 'HydrogenBond'
I_TYPE_HALOGEN = 'HalogenBond'
I_TYPE_HYDROPHOBIC = 'Hydrophobic'
I_TYPE_SALT_BRIDGE = 'SaltBridge'
I_TYPE_PI_STACKING = 'PiStacking'
I_TYPE_PI_CATION = 'PiCation'

I_SUFFIX = 'Interaction'
I_NAME_HBOND = I_TYPE_HBOND + I_SUFFIX
I_NAME_HALOGEN = I_TYPE_HALOGEN + I_SUFFIX
I_NAME_HYDROPHOBIC = I_TYPE_HYDROPHOBIC + I_SUFFIX
I_NAME_SALT_BRIDGE = I_TYPE_SALT_BRIDGE + I_SUFFIX
I_NAME_PI_STACKING = I_TYPE_PI_STACKING + I_SUFFIX
I_NAME_PI_CATION = I_TYPE_PI_CATION + I_SUFFIX


class Interaction:
    def __init__(self, canon_residue, protein_pos, ligand_pos, distance, ligand_atom):
        self.canon_residue = canon_residue
        self.protein_pos = protein_pos
        self.ligand_pos = ligand_pos
        self.distance = distance
        self.ligand_atom = ligand_atom
        self.score_value = None
        self.score_name = None

    def compare(self, other):
        p_dist = distance_between_points(self.protein_pos, other.protein_pos)
        l_dist = distance_between_points(self.ligand_pos, other.ligand_pos)
        score = 0
        # score += max(0.0, (2.0 - p_dist) / 4.0)
        score += max(0.0, (2.0 - l_dist) / 2.0)
        # if score > 0:
        #     print('       ', score, l_dist, self.ligand_pos, other.ligand_pos)

        return score

    def compare_and_store(self, other, name):
        sc = self.compare(other)
        if not self.score_value or sc > self.score_value:
            self.score_value = sc
            self.score_name = name
        return sc

    def __eq__(self, other):
        if not isinstance(other, Interaction):
            return False
        else:
            return self.canon_residue == other.canon_residue and self.protein_pos == other.protein_pos and self.ligand_pos == other.ligand_pos

    def __hash__(self):
        return hash(self.canon_residue, self.protein_pos, self.ligand_pos)


class InteractionType:
    def __init__(self, interaction_type, interactions):
        self.interaction_type = interaction_type
        self.interactions = interactions

    def addInteraction(self, interaction):
        self.interactions.append(interaction)

    def compare(self, other):
        if self.interaction_type != other.interaction_type:
            print('WARNING: comparing interactions of different types:', self.interaction_type, other.interaction_type)
        scores = []
        data = []
        for inter1 in self.interactions:
            for inter2 in other.interactions:
                if inter1.canon_residue == inter2.canon_residue:
                    score = inter1.compare(inter2)
                    scores.append(score)
                    if score > 0:
                        # print('  ', self.interaction_type, inter1.canon_residue, score)
                        data.append((inter1.canon_residue, score))
        return sum(scores), len(scores), data

    def asText(self):
        content = []
        for i in self.interactions:
            if i.ligand_atom is not None:
                s = "%s, [%s, %s, %s], [%s, %s, %s], [%.3f, %s]" % (
                i.canon_residue, i.protein_pos[0], i.protein_pos[1], i.protein_pos[2],
                i.ligand_pos[0], i.ligand_pos[1], i.ligand_pos[2],
                i.distance, i.ligand_atom)
            else:
                s = "%s, [%s, %s, %s], [%s, %s, %s], [%.3f]" % (
                i.canon_residue, i.protein_pos[0], i.protein_pos[1], i.protein_pos[2],
                i.ligand_pos[0], i.ligand_pos[1], i.ligand_pos[2], i.distance)
            if i.score_value:
                s += ", [%s, %s]" % (i.score_value, i.score_name)
            content.append(s)
        return "\n".join(content)


class InteractionSet:
    def __init__(self, id, ligand_name, interaction_types):
        self.id = id
        self.ligand_name = ligand_name
        self.interaction_types = interaction_types

    def add(self, interaction_type):
        self.interaction_types.append(interaction_type)

    def compare(self, other):
        scores = []
        matches = []
        for itype1 in self.interaction_types:
            for itype2 in other.interaction_types:
                if itype1.interaction_type == itype2.interaction_type:
                    score, count, data = itype1.compare(itype2)
                    if count:
                        # print('  ', itype1.interaction_type, score, count)
                        matches.append((itype1.interaction_type, data))
                    scores.append(score)
        return sum(scores), len(scores), matches


class InteractionEncoder(json.JSONEncoder):
    def default(self, i):
        if isinstance(i, Interaction):
            return i.__dict__
        elif isinstance(i, InteractionType):
            return {i.interaction_type: i.interactions}
        elif isinstance(i, InteractionSet):
            d = {} #collections.OrderedDict()
            d['id'] = i.id
            d['ligand_name'] = i.ligand_name
            for t in i.interaction_types:
                d[t.interaction_type] = t.interactions
            return d
        else:
            return super().default(i)


def distance_between_points(pos1, pos2):
    return math.sqrt((pos1[0] - pos2[0]) ** 2 + (pos1[1] - pos2[1]) ** 2 + (pos1[2] - pos2[2]) ** 2)


def from_json(text):
    results = []
    records = json.loads(text)
    for record in records:
        itypes = []
        name = None
        for key in record:
            if key == 'id':
                id = record['id']
            elif key == 'ligand_name':
                name = record['ligand_name']
            elif key.endswith('Interaction'):
                itype = InteractionType(key, [])
                itypes.append(itype)
                values = record[key]
                for value in values:
                    inter = Interaction(value['canon_residue'], value['protein_pos'], value['ligand_pos'],
                                        value['distance'], value.get('ligand_atom', None))
                    itype.addInteraction(inter)
            else:
                print('WARNING: unexpected field %s' % (key))

        iset = InteractionSet(id, name, itypes)
        results.append(iset)
    return results


def compare_interactions(reference_file, test_file):
    with open(reference_file, "r") as ref:
        txt = ref.read()
        ref_data = from_json(txt)
        print('Found', len(ref_data), 'reference items')
    with open(test_file, "r") as test:
        txt = test.read()
        test_data = from_json(txt)
        print('Found', len(test_data), 'test items')
    for i, iset1 in enumerate(test_data, start=1):
        canonical_sites = {}
        for j, iset2 in enumerate(ref_data, start=1):
            # print('Comparing', i, j, iset1.ligand_name, iset2.ligand_name)
            score, count, matches = iset1.compare(iset2)
            if score:
                for match in matches:
                    if match[1]:
                        # print('Matched', i, j, iset1.ligand_name, iset2.ligand_name, score)
                        for inter in match[1]:
                            # print('  ', match[0], inter[0], inter[1])
                            type_canon = (match[0], inter[0])
                            if type_canon in canonical_sites:
                                current_best = canonical_sites[type_canon]
                                if inter[1] > current_best[0][1]:
                                    canonical_sites[type_canon] = (inter, j, iset2.ligand_name)
                            else:
                                canonical_sites[type_canon] = (inter, j, iset2.ligand_name)
        total_score = 0
        for type_canon in canonical_sites:
            inter_j_name = canonical_sites[type_canon]
            total_score += inter_j_name[0][1]

        print('Comparing', i, iset1.ligand_name, total_score)
        for type_canon in canonical_sites:
            inter_j_name = canonical_sites[type_canon]
            print('  ', type_canon[0], type_canon[1], inter_j_name[0][1], inter_j_name[1], inter_j_name[2])


def test_read_write_json():
    with open("report.json", "r") as f:
        txt = f.read()
        data = from_json(txt)
        print('Found', len(data), 'items')
    with open("report2.json", 'w') as report:
        json.dump(data, report, cls=InteractionEncoder)
    for iset1 in data:
        for iset2 in data:
            print('Comparing', iset1.ligand_name, iset2.ligand_name)
            score, count = iset1.compare(iset2)


def main():
    compare_interactions(sys.argv[1], sys.argv[2])


if __name__ == "__main__":
    main()
