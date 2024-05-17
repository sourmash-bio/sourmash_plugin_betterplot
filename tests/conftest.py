import os
import sys

import pytest

from sourmash_tst_utils import TempDirectory, RunnerContext
sys.stdout = sys.stderr

@pytest.fixture
def runtmp():
    with TempDirectory() as location:
        yield RunnerContext(location)
