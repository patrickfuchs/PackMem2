import pytest
import numpy as np

import sys
import os

# Add the parent rep to the sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))


from bin import dico as d


def test_get_value():
    radius = {"DMPG C13": 2.01, "DMPG H13A": 1.34, "DMPG H13B": 1.34, "DMPG OC3": 1.77, "DMPG C12": 2.28, "DMPG P": 2.15, "DMPG HX": 1.34}
    tested_ouput = d.get_value(radius, "DMPG", "C13")
    wanted_output = 2.01
    assert(tested_ouput == wanted_output)

    aliph = {"DMPG C13": 'n', "DMPG H13A": 'n', "DMPG H13B": 'n', "DMPG OC3": 'n', "DMPG C12": 'n', "DMPG P": 'n', "DMPG C22": 'a'}
    tested_ouput = d.get_value(aliph, "DMPG", "C13")
    wanted_output = 'n'
    assert tested_ouput == wanted_output

def test_del_key_dict():
    dictionary = {1 : "okay",
            2: "nope",
            3: "okay",
            4: "okay",
            5: "nope"}
    list_del = [2, 5]
    tested_ouput = d.del_key_dict(dictionary, list_del)
    wanted_output = {1 : "okay",
                    3: "okay",
                    4: "okay"}
    assert tested_ouput == wanted_output

def test_max_value_dict():
    dictionary = {1: np.array([2.34, 3.34, 4.34, 5.34, 6.34, 7.34, 8.34]),
                  2: np.array([2.7, 3.7, 4.7, 5.7, 6.7, 7.7, 8.7, 9.7,
                               10.7, 11.7])}
    tested_out = d.max_value_dict(dictionary)
    wanted_out = 11.7
    assert tested_out == wanted_out
