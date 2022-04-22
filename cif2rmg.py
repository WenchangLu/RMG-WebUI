import os
import sys
import string
import copy
from math import *
from datetime import datetime
from optparse import OptionParser, OptionGroup
import warnings
import CifFile
import subprocess
from utils import *
from uctools import *
from ESPInterfaces import *
from elementdata import *

#################################################################
# Open and read CIF file
cif_file = sys.argv[1]
if cif_file:
    if not os.path.exists(cif_file):
        sys.stderr.write("***Error: The file "+cif_file+" could not be found.\n")
        sys.exit(2)
    cif_file_name = cif_file.split("/")[-1]
    cf = CifFile.ReadCif(cif_file)

##############################################
# Get blocks
cfkeys = list(cf.keys())
cb = cf.get(cfkeys[0])
# Get reference data
ref = ReferenceData()
ref.getFromCIF(cb)
# Get cell data
cd = CellData()
cd.force = True
cd.getFromCIF(cb)


##############################################
# Generate cell
cd.conventional()

inputcell = copy.copy(cd)

############################################################################################
# Output file mode (overwrite or append?)
outputfile = "rmg_input"

################################################################################################
docstring = ''
rmginput = RMGFile(cd, docstring)
print(rmginput) 
