import os
import sys
import string
import copy
from math import *
import warnings
import CifFile
import subprocess
from utils import *
from uctools import *
from ESPInterfaces import *
from elementdata import *

cf = CifFile.ReadCif("cifs/alpha-Mn.cif")
print cf
