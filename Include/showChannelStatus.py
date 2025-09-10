def show_channel_status(channel, settings):
    """
    Prints the status of all channels in a table.

    show_channel_status(channel, settings)

    Inputs:
        channel     - data for each channel. It is used to initialize and
                      at the processing of the signal (tracking part).
        settings    - receiver settings
    --------------------------------------------------------------------------
    #                           SoftGNSS v3.0
    # 
    # Copyright (C) Peter Rinder and Nicolaj Bertelsen
    # Written by Peter Rinder Nicolaj Bertelsen and Darius Plausinaitis
    # Based on Peter Rinder and Nicolaj Bertelsen
    #--------------------------------------------------------------------------
    #This program is free software; you can redistribute it and/or
    #modify it under the terms of the GNU General Public License
    #as published by the Free Software Foundation; either version 2
    #of the License, or (at your option) any later version.
    #
    #This program is distributed in the hope that it will be useful,
    #but WITHOUT ANY WARRANTY; without even the implied warranty of
    #MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    #GNU General Public License for more details.
    #
    #You should have received a copy of the GNU General Public License
    #along with this program; if not, write to the Free Software
    #Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
    #USA.
    #--------------------------------------------------------------------------
    """
    print('\n*=========*=====*===============*===========*=============*========*')
    print('| Channel | PRN |   Frequency   |  Doppler  | Code Offset | Status |')
    print('*=========*=====*===============*===========*=============*========*')

    for channel_nr in range(settings.numberOfChannels):
        ch = channel[channel_nr]
        if ch['status'] != '-':
            doppler = ch['acquiredFreq'] - settings.IF
            print("|      {:2d} | {:3d} |  {:>11.5e} | {:>9.0f} | {:>11d} |   {:>3s}  |".format(
                channel_nr + 1,
                ch['PRN'],
                ch['acquiredFreq'],
                doppler,
                ch['codePhase'],
                ch['status']
            ))
        else:
            print("|      {:2d} | --- |  ------------ |   -----   |    ------   |  Off   |".format(
                channel_nr + 1
            ))

    print('*=========*=====*===============*===========*=============*========*')
    print()
