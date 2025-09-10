from init_settings import init_settings
from Include.probeData import probeData
from Include.postProcessing import postProcessing
#--------------------------------------------------------------------------
#                         CU Multi-GNSS SDR  
# (C) Updated by Yafeng Li, Nagaraj C. Shivaramaiah and Dennis M. Akos
# Based on the original work by Darius Plausinaitis, Peter Rinder, 
# Nicolaj Bertelsen and Dennis M. Akos
#--------------------------------------------------------------------------

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
# USA.
#--------------------------------------------------------------------------
#
# Script initializes settings and environment of the software receiver.
# Then the processing is started.

#--------------------------------------------------------------------------
# CVS record:
# $Id: init.m,v 1.14.2.21 2006/08/22 13:46:00 dpl Exp $

# Clean up the environment first =========================================
# clear all; close all; clc;

#--- Include folders with functions ---------------------------------------
# addpath ('../IF_Data_Set')    # IF data sets for each SDR receivers

# Print startup ==========================================================
# print(['\n',
#     'Welcome to:  softGNSS\n\n',
#     'An open source GNSS SDR software project initiated by:\n\n',
#     '              Danish GPS Center/Aalborg University\n\n',
#     'The code was improved by GNSS Laboratory/University of Colorado.\n\n',
#     'The software receiver softGNSS comes with ABSOLUTELY NO WARRANTY;\n',
#     'for details please read license details in the file license.txt. This\n',
#     'is free software, and  you  are  welcome  to  redistribute  it under\n',
#     'the terms described in the license.\n\n'])
# print('                   -------------------------------\n\n')

# Initialize constants, settings =========================================
settings = init_settings()

try:
    print(f'Probing data ({settings.fileName})...')
    probeData(settings)
except Exception as e:
    print(str(e))
    print('  (change settings in "init_settings.py" to reconfigure)')
    exit()

print('  Raw IF data plotted ')
print('  (change settings in "init_settings.py" to reconfigure)')
print(' ')
gnssStart = input('Enter "1" to initiate GNSS processing or "0" to exit : ')

if gnssStart.strip() == "1":
    postProcessing(settings)

