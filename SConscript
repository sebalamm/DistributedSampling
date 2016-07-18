#/******************************************************************************
# * SConscript
# *
# * Source of the sampling routine
# ******************************************************************************
# * Copyright (C) 2016 Sebastian Lamm <lamm@ira.uka.de>
# *
# * This program is free software: you can redistribute it and/or modify it
# * under the terms of the GNU General Public License as published by the Free
# * Software Foundation, either version 3 of the License, or (at your option)
# * any later version.
# *
# * This program is distributed in the hope that it will be useful, but WITHOUT
# * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
# * more details.
# *
# * You should have received a copy of the GNU General Public License along with
# * this program.  If not, see <http://www.gnu.org/licenses/>.
# *****************************************************************************/

# The main SConscript file for the code.
import platform
import sys
import os

# Get the current platform.
SYSTEM = platform.uname()[0]

Import('env')

# SYSTEMNAME=platform.platform(aliased=0, terse=0)
# if "supermuc" in os.environ['MYSYSTEMNAME']:
#     env['CXX'] = 'mpiCC'
#     env['CC']  = 'mpicc'
#     env.Append(LINKFLAGS = '-openmp')
#     env.Append(CXXFLAGS = '-msse4.2 -openmp -w')
#     env.Append(CCFLAGS = '-msse4.2 -openmp -w')
# else: 

env['CXX'] = 'mpicxx'
env['CC']  = 'mpicc'
env.Append(CXXFLAGS = '-msse4.2')
env.Append(CCFLAGS = '-msse4.2')

if SYSTEM == 'Darwin':
    env['OMPI_CXX'] = 'g++-4.8'
    env['OMPI_CC'] = 'gcc-4.8'

if env['program'] == 'library':
    stocc = env.Library('stocc', [
        'extern/mersenne64/userintf.cpp',
        'extern/mersenne64/mersenne.cpp',
        'extern/stocc64/stoc1.cpp',
        'extern/stocc64/stoc2.cpp',
        'extern/stocc64/stoc3.cpp',
        'extern/stocc64/wnchyppr.cpp',
        'extern/stocc64/fnchyppr.cpp'
    ])

if env['program'] == 'run_methodD':
    env.Program('run_methodD', ['app/run_methodD.cpp'], LIBS=['libargtable2','gomp','stocc'])

if env['program'] == 'run_methodH':
    env.Program('run_methodH', ['app/run_methodH.cpp'], LIBS=['libargtable2','gomp','stocc'])

if env['program'] == 'run_methodR':
    env.Program('run_methodR', ['app/run_methodR.cpp'], LIBS=['libargtable2','gomp','stocc'])

if env['program'] == 'run_methodP':
    env.Program('run_methodP', ['app/run_methodP.cpp'], LIBS=['libargtable2','gomp','stocc'])

if env['program'] == 'run_experiments':
    env.Program('run_experiments', ['app/run_experiments.cpp'], LIBS=['libargtable2','gomp','stocc'])
