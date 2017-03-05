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

import sys, gzip, json

class StreamJsonListLoader():
    """
    When you have a big JSON file containing a list, such as

    [{ ... },{ ... },{ ... }, ... ]

    And it's too big to be practically loaded into memory and parsed by json.load,
    This class comes to the rescue. It lets you lazy-load the large json list.

    This code is borrowed from this stackoverlow question:
    http://stackoverflow.com/questions/6886283/how-i-can-i-lazily-read-multiple-json-objects-from-a-file-stream-in-python
    """

    def __init__(self, filename_or_stream):
        if type(filename_or_stream) == str:
            self.stream = open(filename_or_stream)
        else:
            self.stream = filename_or_stream

        if not self.stream.read(1) == '[':
            raise NotImplementedError('Only JSON-streams of lists (that start with a [) are supported.')

    def __iter__(self):
        return self

    def next(self):
        read_buffer = self.stream.read(1)
        while True:
            try:
                json_obj = json.loads(read_buffer)

                if not self.stream.read(1) in [',',']']:
                    raise Exception('JSON seems to be malformed: object is not followed by comma (,) or end of list (]).')
                return json_obj
            except ValueError:
                next_char = self.stream.read(1)
                read_buffer += next_char
                while next_char != '}':
                    next_char = self.stream.read(1)
                    if next_char == '':
                        raise StopIteration
                    read_buffer += next_char

    def close(self):
        if self.stream:
            self.stream.close()


def main():
    loader = StreamJsonListLoader(sys.stdin)
    count = 0
    for j in loader:
        #print "got some json", j
        count += 1
    loader.close()
    print "Found",count,"records"


if __name__ == "__main__":
    main()
