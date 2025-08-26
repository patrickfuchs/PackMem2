import pytest

import sys
import os

# Add the parent rep to the sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))


from  bin import dico as d

def test_get_radius():
    radius = {"DMPG C13": 2.01, "DMPG H13A": 1.34, "DMPG H13B": 1.34, "DMPG OC3": 1.77, "DMPG C12": 2.28, "DMPG P": 2.15, "DMPG HX": 1.34}
    tested_ouput = d.get_radius(radius, "DMPG", "C13")
    wanted_output = 2.01
    assert(tested_ouput == wanted_output)

def test_get_aliphatic():
    aliph = {"DMPG C13": 'n', "DMPG H13A": 'n', "DMPG H13B": 'n', "DMPG OC3": 'n', "DMPG C12": 'n', "DMPG P": 'n', "DMPG C22": 'a'}
    tested_ouput = d.get_aliphatic(aliph, "DMPG", "C13")
    wanted_output = 'n'
    assert(tested_ouput == wanted_output)