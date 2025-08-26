import pytest
import io

import sys
import os

# Add the parent rep to the sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))


from bin import pdb as pdb

def test_write_a_pdb_line():
    # Create fake file
    fake_file = io.StringIO()
    # Write
    pdb.write_a_pdb_line(fake_file, 1, 'EDG', 1, [-8.0, -9.0, 7.310], -1.0)

    # Read the fake file
    fake_file.seek(0)  # Seek the start of the file
    content = fake_file.read()

    wanted_line = f"ATOM      1   H1 EDG     1      -8.000  -9.000   7.310  1.00 -1.00\n"
    assert content == wanted_line
