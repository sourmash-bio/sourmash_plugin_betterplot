"""xyz plugin description"""

usage="""
   sourmash scripts xyz
"""

epilog="""
See https://github.com/xyz for more examples.

Need help? Have questions? Ask at http://github.com/sourmash-bio/sourmash/issues!
"""

import argparse
import sourmash

from sourmash.index import LinearIndex
from sourmash.logging import debug_literal
from sourmash.plugins import CommandLinePlugin

from sourmash.save_load import (Base_SaveSignaturesToLocation,
                                _get_signatures_from_rust)


###

#
# load_from plugin
#

def load_sketches(location, *args, **kwargs):
    if location and location.endswith('.xyz'):
        # ... add your code here ...
        # return LinearIndex(siglist)
        pass

load_sketches.priority = 5


#
# save_to plugin - supports output to .xyz
#

class SaveSignatures_XYZ(Base_SaveSignaturesToLocation):
    "Save signatures to a location."
    def __init__(self, location):
        super().__init__(location)
        self.keep = []

    @classmethod
    def matches(self, location):
        # match anything that is not None or ""
        if location:
            return location.endswith('.xyz')

    def __repr__(self):
        return f"SaveSignatures_XYZ('{self.location}')"

    def open(self):
        pass

    def close(self):
        # actually do the writing...
        pass

    def add(self, ss):
        super().add(ss)
        self.keep.append(ss)

#
# CLI plugin - supports 'sourmash scripts xyz'
#

class Command_XYZ(CommandLinePlugin):
    command = 'xyz'             # 'scripts <command>'
    description = __doc__       # output with -h
    usage = usage               # output with no args/bad args as well as -h
    epilog = epilog             # output with -h
    formatter_class = argparse.RawTextHelpFormatter # do not reformat multiline

    def __init__(self, subparser):
        super().__init__(subparser)
        # add argparse arguments here.
        debug_literal('RUNNING cmd_xyz.__init__')

    def main(self, args):
        # code that we actually run.
        super().main(args)
        print('RUNNING cmd', self, args)
