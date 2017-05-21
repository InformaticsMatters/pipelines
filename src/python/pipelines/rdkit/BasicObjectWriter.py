#!/usr/bin/env python

# Copyright 2017 Informatics Matters Ltd.
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

import uuid, json

class BasicObjectWriter():

    def __init__(self, file):
        if type(file) == str:
            self.file = open(file, 'w')
        else:
            self.file = file
        self.count = 0

    def writeHeader(self):
        self.file.write('[')

    def writeFooter(self):
        self.file.write(']')

    def write(self, dictOfValues, objectUUID=None):

        d = {}

        if objectUUID:
            d['uuid'] = objectUUID
        else:
            d['uuid'] = str(uuid.uuid4())

        d['values'] = dictOfValues
        json_str = json.dumps(d)
        if self.count > 0:
            self.file.write(',\n')
        self.file.write(json_str)
        self.count += 1

    def close(self):
        if self.file:
            self.file.close()


