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

class TsvWriter():

    def __init__(self, file, headersAsOrderedDict):
        if type(file) == str:
            self.file = open(file, 'w')
        else:
            self.file = file
        self.headersAsOrderedDict = headersAsOrderedDict

    def write(self, dictOfValues):
        h = self.headersAsOrderedDict
        count = 0
        for k in h:
            if count > 0:
                self.file.write('\t')
            if k in dictOfValues:
                self.file.write(str(dictOfValues[k]))
            count += 1
        self.file.write('\n')

    def writeHeader(self):
        d = self.headersAsOrderedDict
        count = 0
        for k in d:
            if count > 0:
                self.file.write('\t')
            self.file.write(str(d[k]))
            count += 1
        self.file.write('\n')

    def writeFooter(self):
        pass

    def close(self):
        if self.file:
            self.file.close()




