# -*- coding: utf-8 -*-
# Copyright 2013 Pierre de Buyl
#
# This file is part of md_tools
#
# md_tools is free software and is licensed under the modified BSD license (see
# LICENSE file).

# This file contains code from the nMOLDYN project. The complete nMOLDYN license
# "CeCILL" is available at licenses/LICENSE_nMOLDYN-3.0.10.txt, in the
# distribution root directory.
# nMOLDYN is Copyright by E. Pellegrini, K. Hinsen and G.R. Kneller
# http://dirac.cnrs-orleans.fr/nMOLDYN/

# The ScientificPython modules
from Scientific import N
from Scientific.FFT import fft, inverse_fft

def correlation(inputSeries1, inputSeries2 = None):
    """Returns the numerical correlation between two signals.

    @param inputSeries1: the first signal.
    @type inputSeries1: NumPy array   
    
    @param inputSeries2: if not None, the second signal otherwise the correlation will be an autocorrelation.
    @type inputSeries2: NumPy array or None
    
    @return: the result of the numerical correlation.
    @rtype: NumPy array

    @note: if |inputSeries1| is a multidimensional array the correlation calculation is performed on
    the first dimension.

    @note: The correlation is computed using the FCA algorithm.
    """

    # The signal must not be empty.
    if len(inputSeries1) <= 0:
        raise Error('One or both time series are empty.')

    # The length of inputSeries1 is stored in inputSeries1Length
    inputSeries1Length = len(inputSeries1)

    # extendedLength = 2*len(inputSeries1)
    extendedLength = 2*inputSeries1Length

    # The FCA algorithm:

    # 1) computation of the FFT of inputSeries1 zero-padded until extendedLength
    # The computation is done along the 0-axis
    FFTSeries1 = fft(inputSeries1,extendedLength,0)        

    if inputSeries2 is None:
            # Autocorrelation case
        FFTSeries2 = FFTSeries1
    else:
        # 2) computation of the FFT of inputSeries2 zero-padded until extendedLength
        # The computation is  done along the 0-axis
        FFTSeries2 = fft(inputSeries2,extendedLength,0)

    # 3) Product between FFT(inputSeries1)* and FFT(inputSeries2)
    FFTSeries1 = N.conjugate(FFTSeries1)*FFTSeries2

    # 4) inverse FFT of the product
    # The computation is done along the 0-axis
    FFTSeries1 = inverse_fft(FFTSeries1,len(FFTSeries1),0)

    # This refers to (1/(N-m))*Sab in the published algorithm.
    # This is the correlation function defined for positive indexes only.
    if len(FFTSeries1.shape) == 1:
        corr = FFTSeries1.real[:inputSeries1Length] / (inputSeries1Length-N.arange(inputSeries1Length))
    else:
        corr = N.add.reduce(FFTSeries1.real[:inputSeries1Length],1) / (inputSeries1Length-N.arange(inputSeries1Length))

    return corr

def calc_msd(series):
        """Calculates the atomic term.
        
        @param atomIndex: the index of the selected atom.
        @type atomIndex: integer
        
        @param trajectory: the trajectory.
        @type trajectory: MMTK.Trajectory.Trajectory object
        """
                
        # series = 2D Numeric array. The positions of the selected atom |at| from the first step to the
        # last step with the selected step increment.
        # series = trajectory.readParticleTrajectory(atomIndex, first = self.first, last = self.last, skip = self.skip).array
                            
        # dsq is the squared norm of the position for each time step
        # dsq refers to DSQ(k) in the published algorithm
        dsq = N.add.reduce(series * series,1)

        # sum_dsq1 is the cumulative sum of dsq
        sum_dsq1 = N.add.accumulate(dsq)

        # sum_dsq1 is the reversed cumulative sum of dsq
        sum_dsq2 = N.add.accumulate(dsq[::-1])

        # sumsq refers to SUMSQ in the published algorithm
        sumsq = 2.*sum_dsq1[-1]

        # this line refers to the instruction SUMSQ <-- SUMSQ - DSQ(m-1) - DSQ(N - m) of the published algorithm
        # In this case, msd is an array because the instruction is computed for each m ranging from 0 to len(traj) - 1
        # So, this single instruction is performing the loop in the published algorithm
        Saabb  = sumsq - N.concatenate(([0.], sum_dsq1[:-1])) - N.concatenate(([0.], sum_dsq2[:-1]))

        # Saabb refers to SAA+BB/(N-m) in the published algorithm
        # Sab refers to SAB(m)/(N-m) in the published algorithm
        Saabb = Saabb / (len(dsq) - N.arange(len(dsq)))
        Sab   = 2.*correlation(series)

        atomicMSD = Saabb - Sab
                                
        return atomicMSD
