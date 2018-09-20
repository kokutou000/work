import numpy as np
import sys
import subprocess
import shutil

"""
  2018.08.28
  this script file should be on 
      data/sh***/output/auto_output.py
"""

#from code_output import funcs as fnc
import code_output.funcs as fnc

### make dirs
fnc.code_mkdirs()
####### memo dirname --- 2018.07.27
### dirname[0] ---> binfile
#fnc.code_binfile()
### dirname[1] ---> 2d
### dirname[2] ---> 2d/xy
### dirname[10] ---> 2d/yz
#fnc.code_2dxy()
#fnc.code_2dyz()
### dirname[18] ---> 3dbin
#fnc.code_3dbin()
### dirname[19] ---> r3dbline
#fnc.code_r3dbline_calc()
#fnc.code_r3dbline_plot()
### dirname[30] ---> calc_kappa
#fnc.code_calckappa()
### dirname[31] ---> vapor
#fnc.code_vapor_makefile()
### dirname[33] ---> kappa2
fnc.code_kappa2()
####################################
print("sucess")

sys.exit()

