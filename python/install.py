#!/usr/bin/env python

# copy LAMMPS src/liblammps.so and lammps.py to system dirs

from __future__ import print_function

instructions = """
Syntax: python install.py [-h] [pydir]
        pydir = target dir for lammps.py and liblammps.so
                default = Python site-packages dir
"""

import sys,os,shutil

if (len(sys.argv) > 1 and sys.argv[1] == "-h") or len(sys.argv) > 2:
  sys.exit(instructions)

if len(sys.argv) == 2: pydir = sys.argv[1]
else: pydir = ""

# copy lammps.py to pydir if it exists
# if pydir not specified, install in site-packages via distutils setup()

if pydir:
  if not os.path.isdir(pydir):
    sys.exit("ERROR: pydir %s does not exist" % pydir)

  print("Copying python/lammps.py to %s\n" % pydir)
  try:
    shutil.copyfile(os.path.join("..", "python", "lammps.py"), os.path.join(pydir, 'lammps.py'))
  except shutil.Error:
    pass # source and destination are identical

  shlib = os.path.join("..", "src", "liblammps.dylib")
  if not os.path.isfile(shlib):
    shlib = os.path.join("..", "src", "liblammps.dll")
  if not os.path.isfile(shlib):
    shlib = os.path.join("..", "src", "liblammps.so")
  if not os.path.isfile(shlib):
    sys.exit("Cannot find LAMMPS shared library file %s" % shlib)

  print("Copying %s to %s\n" % (shlib, pydir))
  try:
     shutil.copyfile(shlib, os.path.join(pydir,"liblammps.so") )
  except shutil.Error:
    pass # source and destination are identical
  sys.exit()

print("Installing lammps.py in Python site-packages dir")

os.chdir(os.path.join('..', 'python')       # in case invoked via make in src dir

# extract version string from header
fp = open(os.path.join('..', 'src', 'version.h'), 'r')
txt=fp.read().split('"')[1].split()
verstr=txt[0]+txt[1]+txt[2]
fp.close()

from distutils.core import setup
from distutils.sysconfig import get_python_lib
import site
tryuser=False

shlib = os.path.join("..", "src", "liblammps.dylib")
if not os.path.isfile(shlib):
  shlib = os.path.join("..", "src", "liblammps.dll")
if not os.path.isfile(shlib):
  shlib = os.path.join("..", "src", "liblammps.so")
if not os.path.isfile(shlib):
  sys.exit("Cannot find LAMMPS shared library file %s" % shlib)

try:
  sys.argv = ["setup.py","install"]    # as if had run "python setup.py install"
  setup(name = "lammps",
        version = verstr,
        author = "Steve Plimpton",
        author_email = "sjplimp@sandia.gov",
        url = "http://lammps.sandia.gov",
        description = "LAMMPS molecular dynamics library",
        py_modules = ["lammps"],
        data_files = [(get_python_lib(), [shlib])])
except:
  tryuser=True
  print ("Installation into global site-packages dir failed.\nTrying user site dir %s now." % site.USER_SITE)


if tryuser:
  try:
    sys.argv = ["setup.py", "install", "--user"]    # as if had run "python setup.py install --user"
    setup(name = "lammps",
    version = verstr,
    author = "Steve Plimpton",
    author_email = "sjplimp@sandia.gov",
    url = "http://lammps.sandia.gov",
    description = "LAMMPS molecular dynamics library",
    py_modules = ["lammps"],
    data_files = [(site.USER_SITE, [shlib])])
  except:
    print("Installation into user site package dir failed.\nGo to ../python and install manually.")
