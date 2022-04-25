import os
import sys
import string
import copy
import math
from datetime import datetime
from optparse import OptionParser, OptionGroup
import warnings
import CifFile
import subprocess
from utils import *
from uctools import *

class rmg_interface():
    def cif2cell(self, cif_file=None): 
        #################################################################
        # Open and read CIF file
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
        if self.reducetoprim:
            cd.primitive()
        else:
            cd.conventional()

        inputcell = copy.copy(cd)

        self.cell = cd


        self.ibrav = 0
        t = LatticeMatrix(self.cell.latticevectors)
        for i in range(3):
            for j in range(3):
                t[i][j] = self.cell.latticevectors[i][j]*self.cell.lengthscale
        ortho = abs(self.cell.a - t[0][0]) < 1.0e-5 
        ortho = ortho and abs(self.cell.b - t[1][1]) < 1.0e-5
        ortho = ortho and abs(self.cell.c - t[2][2]) < 1.0e-5

        if ortho: self.ibrav = 8
        system = self.cell.crystal_system()
        setting = self.cell.spacegroupsetting
        if system == 'cubic':
            if self.cell.primcell:
                if setting == 'P':
                    self.ibrav = 1
                elif setting == 'F':
                    self.ibrav = 2
                elif setting == 'I':
                    self.ibrav = 3
            else:
                self.ibrav = 1
        if system == 'hexagonal':
            if self.cell.primcell:
                if setting == 'P':
                    self.ibrav = 4
            #elif setting == 'R':
            #    return 5

    def cell2rmg(self):
        self.cell.newunit("bohr")
        # set up species list
        tmp = set([])
        for a in self.cell.atomdata:
            for b in a:
                tmp.add(b.spcstring())
        self.species = list(tmp)
        filestring = ""
        #
        # some default input options
        filestring += """
atomic_coordinate_type = "Absolute"  
write_eigvals_period = "10"  
input_wave_function_file = "Wave/wave"  
output_wave_function_file = "Wave/wave"  
"""
        brav_type = {
            0:"None",
            1:"Cubic Primitive",
            2:"Cubic Face Centered",
            3:"Cubic Body Centered",
            4:"Hexagonal Primitive",
            8:"Orthorhombic Primitive"
        }
        filestring += 'bravais_lattice_type="%s"  \n'%brav_type[self.ibrav]
        if self.ibrav !=0:
            filestring += 'a_length="%16.8f"  \n'%self.cell.a
            filestring += 'b_length="%16.8f"  \n'%self.cell.b
            filestring += 'c_length="%16.8f"  \n'%self.cell.c
        else:
            filestring += 'lattice_vector="  \n'
            filestring += str(t)
            filestring += '"  \n'

        filestring += 'atoms="  \n'
        atom_format = "%s  %.12e %.12e %.12e"
        for a in self.cell.atomdata:
            for b in a:
                t = Vector(mvmult3(self.cell.latticevectors,b.position))
                for i in range(3):
                    t[i] = self.cell.lengthscale*t[i]
                filestring += atom_format%(b.spcstring(),t[0], t[1], t[2])
                filestring += "  1 1 1 0.0 0.0 0.0  \n"
        filestring += '"  \n'

        return filestring

    def __init__(self, filename, filetype):
        self.reducetoprim = True
        if filetype == ".cif":
            self.cif2cell(filename)

        self.rmginput = self.cell2rmg()
