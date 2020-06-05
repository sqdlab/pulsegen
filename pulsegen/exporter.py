'''
Created on 21/05/2013

@author: Matthias Baur
'''

# TODO: Add wfm_count, wfm_adv, marker_masks to sequence files

import numpy as np
import os
import shutil
import logging
import sys
from time import sleep

def to_bytes(s):
    '''Convert s to bytes, assuming latin-1 encoding'''
    if sys.version_info[0] >= 3:
        return bytes(s, encoding='latin-1')
    else:
        return s

def num_pad(number):
    return "{:=05d}".format(number)

def rescale(x, from_lowhigh, to_lowhigh):
    '''
    Rescale `x` to be in the range `to_lowhigh`
    when its values are in the range `from_lowhigh`
    
    Paramteres
    ----------
    x : float or numpy array of floats
        Input
    from_lowhigh: tuple
        The range of the numbers `x` (min, max)
    to_lowhigh: tuple
        The range of the output numbers (min, max)
        
    Returns
    -------
    r : float or numpy array of floats
        Rescaled value of `x`
    
    Example:
    --------
    >>> rescale(np.array([-1, 0, 1]), (-1, 1), (0, 10))
    array([  0.,   5.,  10.])
    
    '''    
    if np.any(x < from_lowhigh[0]) or np.any(x > from_lowhigh[1]):
        x = np.clip(x, *from_lowhigh)
        logging.warning('Input samples clipped.')
    scaling = 1.0*((to_lowhigh[1] - to_lowhigh[0])
                   /(from_lowhigh[1] - from_lowhigh[0]))
    
    return scaling*(x - from_lowhigh[0]) + to_lowhigh[0]

def pat_name_format(fileprefix, chlabel=None, number=None):
    if chlabel is not None and number is not None:
        return '{0:s}Ch{1:d}{2:s}'.format(fileprefix, chlabel, num_pad(number))
    elif chlabel is not None:
        return '{0:s}Ch{1:d}'.format(fileprefix, chlabel)
    elif number is not None:
        return '{0:s}{2:s}'.format(fileprefix, num_pad(number))
    else:
        return fileprefix
    
# def export(awg_type, directory, fileprefix, wfm_sequence, marker_sequence,
#            lowhigh):
#     if awg_type == "Tektronix_5014":
#         exporter = Tektronix5014Export(directory, fileprefix)
#     elif awg_type == "Agilent_N824x":
#         exporter =AgilentN8241Export(directory, fileprefix)
#     else:
#         raise ValueError('Unknown AWG model `{0}`'.format(awg_type))
#         
#     exporter.export_sequence(wfm_sequence, marker_sequence, lowhigh)

def export(awg_sequence, directory, fileprefix, rmdir=True):
    '''
    Export an AWGSequence
    
    Parameters:
    -----------
    awg_sequence: AWGSequence
        The awg sequence object
    directory: String
        Directory where the pattern files should be exported to
    fileprefix: String
        Prefix of the pattern files to be saved. Does not include the `.pat`
        or `.seq` respectively.
    rmdir: Bool
        If True, empty existing sequence directory before exporting
    '''
    awg_model = awg_sequence.awg.model
    if awg_model == "Tektronix_5014":
        exporter = Tektronix5014Export(directory, fileprefix)
    elif awg_model == "Agilent_N824x":
        exporter =AgilentN8241Export(directory, fileprefix)
    elif awg_model == "Agilent_81180A":
        exporter =Agilent81180AExport(directory, fileprefix)
    else:
        raise ValueError('Unknown AWG model `{0}`'.format(awg_model))
    
    waveforms = [seq.waveforms for seq in awg_sequence.sampled_sequences]
    markers = [seq.markers for seq in awg_sequence.sampled_sequences]
    lowhighs = [ch.lowhigh for ch in awg_sequence.channels]
    
    exporter.export_sequence(waveforms, markers, lowhighs, rmdir)

def deduplicate(arrs):
    '''
    Find duplicate arrays in arrs
    
    Implements a naive O(N**2) comparison of the current element 
    with previous elements.
    
    Input:
        arrs - iterable of arrays
    Output:
        unique - unique elements
        indices - indices into unique that will restore the original array
    '''
    unique = []
    indices = []
    inverse = []
    for arr in arrs:
        found = False
        for idx, seen in enumerate(unique):
            if np.all(np.array(seen) == np.array(arr)):
                found = True
                break
        if found:
            inverse.append(idx)
        else:
            indices.append(len(inverse))
            inverse.append(len(unique))
            unique.append(arr)
    return unique, indices, inverse
    
class Tektronix5014Export(object):
    '''
    Tektronix 5014 Exporter
    
    '''
    ## Even though Tek5014 only has 14 bits for the data,
    ## the pattern file format requires the data to be saved as 16 bit integers.
    ## When uploaded to the AWG, the two least significant bits are dropped.
    wfm_mask = 2**16 - 1
    marker1_mask = 2**0
    marker2_mask = 2**1
    lowhigh = (0, wfm_mask)
    uint24 = np.dtype('<u2, u1')
    
    def __init__(self, directory, fileprefix):
        self.dir = directory
        self.fileprefix = fileprefix
        self._suffix = ".pat"
    
    def _pack_pattern(self, waveform, marker1, marker2):
        pattern_words = np.bitwise_and(waveform, Tektronix5014Export.wfm_mask)
        marker1 = np.array(np.round(marker1), dtype=bool)
        marker2 = np.array(np.round(marker2), dtype=bool)
        marker_bytes = marker1 * Tektronix5014Export.marker1_mask
        marker_bytes = np.bitwise_or(
            marker_bytes, marker2 * Tektronix5014Export.marker2_mask)
        pattern = np.array(
            list(zip(pattern_words, marker_bytes)), dtype=Tektronix5014Export.uint24)
        return pattern
    
    def export_sequence_file(self, ch_labels):
        '''
        Save the sequence file. Details about the format of this sequence file
        can be found in the programmer manual of the Tektronix AWG400.
        
        '''
        def string_format(ch1_label, ch2_label, ch3_label, ch4_label):
            string = ('"' + pat_name_format(self.fileprefix, 1, ch1_label)
                      + self._suffix + '",'
                      '"' + pat_name_format(self.fileprefix, 2, ch2_label)
                      + self._suffix + '",' 
                      '"' + pat_name_format(self.fileprefix, 3, ch3_label)
                      + self._suffix + '",' 
                      '"' + pat_name_format(self.fileprefix, 4, ch4_label)
                      + self._suffix + '"')
            return to_bytes(string + ',1,1,0,0\r\n')
        
        if not os.path.isdir(self.dir):
            os.mkdir(self.dir)
        
        file_name = os.path.join(self.dir, self.fileprefix + '.seq')
        with open(file_name, 'wb', 1) as seq_file:
            seq_file.write(b"MAGIC 3004\r\n")
            seq_file.write(b"LINES %d\r\n"%(len(ch_labels[0])))
            for labels in zip(*ch_labels):
                seq_file.write(string_format(*labels))
    
    def export_pattern(self, waveform, marker1, marker2, lowhigh,
                       ch_idx=None, number=None):
        '''
        Save the waveform and markers in a pattern file. 
        Details about the pattern file format can be found in the
        programmer manual of the Tektronix AWG400.
        
        '''
        waveform_rescaled = rescale(
            waveform, lowhigh, Tektronix5014Export.lowhigh)
        waveform_rescaled_rounded = np.array(
            np.round(waveform_rescaled), dtype=int)
        pattern = self._pack_pattern(
            waveform_rescaled_rounded, marker1, marker2)
        
        fileprefix = pat_name_format(self.fileprefix, ch_idx + 1, number)
        
        filename = os.path.join(self.dir, fileprefix + self._suffix)
        
        num_bytes = 3 * len(pattern)
        num_digits = len(str(num_bytes))
        header = b"MAGIC 2003\r\n#%d%d" %(num_digits, num_bytes)
        
        if not os.path.isdir(self.dir):
            os.makedirs(self.dir)
        with open(filename, 'wb') as pat_obj:
            pat_obj.write(header)
            pat_obj.write(pattern.data)
    
    def export_sequence(self, wfm_sequences, marker_sequences, lowhighs, rmdir=True):
        '''
        Save all waveforms and markers in pattern files and create a 
        sequence file.
        
        Parameters:
        -----------
        wfm_sequences: List
            List of sequences, where each entry with index `idx` defines
            the sequence to be played on the AWG channel `idx`.
        marker_sequences: List
            List of marker sequences, where each entry with index `idx` defines
            the marker sequence to be played on the AWG markers corresponding to
            the AWG channel `idx`.
        lowhighs: List
            Each entry defines the low and high levels of the corresponding
            AWG channel, used to rescale the waveform. See also the help of
            the rescale function.
        rmdir: Bool
            If True, empty existing sequence directory before exporting
        '''
        num_channels = len(wfm_sequences)
        if num_channels != 4 or num_channels != len(marker_sequences):
            raise ValueError('`wfm_sequence` and `marker_sequence` should have'
                             ' exactly four channels')

        num_patterns = [len(wfm_sequence) for wfm_sequence in wfm_sequences]
        if min(num_patterns) != max(num_patterns):
            raise ValueError('sequences must have equal length on all channels' 
                             'of the device.')
        
        if not os.path.isdir(self.dir):
            os.makedirs(self.dir)
        else:
            if rmdir:
                shutil.rmtree(self.dir)
                # TODO: Nicer way to prevent permission denied error than sleep?
                sleep(0.1)
                os.makedirs(self.dir)

        # file names for all channels
        seq_labelss = []
        for ch_idx in range(num_channels):
            wfm_sequence = wfm_sequences[ch_idx]
            marker_sequence = marker_sequences[ch_idx]

            if len(wfm_sequence) != len(marker_sequence):
                raise ValueError('`wfm_sequence` and `marker_sequence` must'
                                 ' have equal length')

            wfm_unique, _, wfm_inverse = deduplicate(wfm_sequence)
            marker_unique, _, marker_inverse = deduplicate(marker_sequence)
            # the sequence file does not allow separate specification of the marker
            # and waveform files. we have to write both out when either is different
            seq_unique, seq_labels, seq_inverse = deduplicate(zip(wfm_inverse, 
                                                                  marker_inverse))

            for seq_idxs, seq_label in zip(seq_unique, seq_labels):
                wfm_idx, marker_idx = seq_idxs
                self.export_pattern(wfm_unique[wfm_idx], 
                                    marker_unique[marker_idx][0],
                                    marker_unique[marker_idx][1], 
                                    lowhighs[ch_idx], ch_idx, seq_label)
            seq_labelss.append([seq_labels[idx] for idx in seq_inverse])

        self.export_sequence_file(seq_labelss)

class AgilentN8241Export(object):
    '''
    Agilent N8241/2 Exporter
    
    '''
    lowhigh = (-1, 1)
    marker1_mask = 2**6
    marker2_mask = 2**7
    
    ## TODO: Do we have to add a blank waveform and blank markers at the
    ## beginning as a workaround for the problem that the N8242A AWG starts
    ## the first pattern file when triggered with the start software trigger.
    ## Adding the blank waveform makes sure that the first real waveform
    ## is triggered with the hardware trigger instead of the software trigger.
    ## This of course depends on how the AWG is addressed in QTLab.
    ## Maybe one should take care of that in QTLab rather than here.
    
    def __init__(self, directory, fileprefix):
        self.dir = directory
        self.fileprefix = fileprefix
        self._suffix = ".bin"
        
    def _pack_markers(self, marker1, marker2):
        marker1 = np.array(np.round(marker1), dtype=bool)
        marker2 = np.array(np.round(marker2), dtype=bool)
        markers = marker1 * AgilentN8241Export.marker1_mask
        markers = np.bitwise_or(
            markers, marker2 * AgilentN8241Export.marker2_mask)
        return np.array(markers, dtype=np.int8)
        
    def export_pattern(self, waveform, marker1, marker2, lowhigh,
                       ch_idx=None, number=None):
        '''
        Save the waveform and markers in a pattern file.
        
        Waveform and markers use two separate files.
        
        The waveform is saved in a binary file, where each data point is
        represented by a 32bit wide floating point number between -1 and 1.
        
        The markers are saved in another binary file, where each data point is
        a 8 bit integer, where the last two bits define the markers. Bit 6 is
        marker one, and bit 7 is marker 2 (counting the bits from 0 - 7).
        
        '''
        waveform = np.array(waveform, np.float32)
        
        if len(marker2) != len(marker1):
            raise ValueError('markers must have same length')
        
        if len(waveform) < 128:
            raise ValueError('waveform must have at least 128 samples')
        
        remainder = len(waveform) % 8
        if remainder != 0:
            logging.warning(('Waveform length is not a multiple of 8.'
                            ' Waveform is automatically padded at the end with'
                            ' the last element to fulfill this requirement.'
                            ))
            waveform = np.pad(waveform, (0, 8 - remainder), 'edge')
            marker1 = np.pad(marker1, (0, 1), 'edge')
            marker2 = np.pad(marker2, (0, 1), 'edge')
            
        if 8 * len(marker1) != len(waveform):
            raise ValueError('len(waveform) must be equal to 8 * len(marker)')
            
        fileprefix = pat_name_format(self.fileprefix, ch_idx + 1, number)
        filename = os.path.join(self.dir, fileprefix + self._suffix)
        
        if not os.path.isdir(self.dir):
            os.makedirs(self.dir)
        with open(filename, 'wb') as wfm_file:
            wfm_file.write(waveform)
        
        markers = self._pack_markers(marker1, marker2)
        marker_filename = os.path.join(
            self.dir, 'marker_' + fileprefix + self._suffix)
        with open(marker_filename, 'wb') as marker_file:
            marker_file.write(markers.data)
            
    def export_sequence_file(self, ch_labels):
        '''
        Save the sequence file. This sequence file is intended to look similar
        to the one used for the Tektronix 5014 arbitrary waveform generator.
        
        '''
        def string_format(ch1_label, ch2_label):
            string = ('"' + pat_name_format(self.fileprefix, 1, ch1_label)
                      + self._suffix + '",'
                      '"' + pat_name_format(self.fileprefix, 2, ch2_label)
                      + self._suffix + '"'
                      +'\r\n')
            return to_bytes(string)
        
        seq_filename = os.path.join(self.dir, self.fileprefix + '.seq')
        if not os.path.isdir(self.dir):
            os.mkdir(self.dir)
        with open(seq_filename, 'wb') as seq_file:
            for ch1_label, ch2_label in zip(*ch_labels):
                seq_file.write(string_format(ch1_label, ch2_label))
                
    def export_sequence(self, wfm_sequences, marker_sequences, lowhighs, rmdir=True):
        '''
        Save all waveforms and markers in pattern files and create a 
        sequence file.
        
        Parameters:
        -----------
        wfm_sequences: List
            List of sequences, where each entry with index `idx` defines
            the sequence to be played on the AWG channel `idx`.
        marker_sequences: List
            List of marker sequences, where each entry with index `idx` defines
            the marker sequence to be played on the AWG markers corresponding to
            the AWG channel `idx`.
        lowhighs: List
            Each entry defines the low and high levels of the corresponding
            AWG channel, used to rescale the waveform. See also the help of
            the rescale function.
        rmdir: Bool
            If True, empty existing sequence directory before exporting
        '''
        num_channels = len(wfm_sequences)
        if num_channels != 2 or num_channels != len(marker_sequences):
            raise ValueError('`wfm_sequence` and `marker_sequence` should have'
                             ' exactly two channels')
        
        num_patterns = [len(wfm_sequence) for wfm_sequence in wfm_sequences]
        if min(num_patterns) != max(num_patterns):
            raise ValueError('sequences must have equal length on all channels' 
                             'of the device.')
        
        if not os.path.isdir(self.dir):
            os.makedirs(self.dir)
        else:
            if rmdir:
                shutil.rmtree(self.dir)
                # TODO: Get rid of this sleep
                sleep(0.1)
                os.makedirs(self.dir)
        
        # file names for all channels
        seq_labelss = []
        for ch_idx in range(num_channels):
            wfm_sequence = wfm_sequences[ch_idx]
            marker_sequence = marker_sequences[ch_idx]
            
            if len(wfm_sequence) != len(marker_sequence):
                raise ValueError('`wfm_sequence` and `marker_sequence` must'
                                 ' have equal length')
                
            wfm_unique, _, wfm_inverse = deduplicate(wfm_sequence)
            marker_unique, _, marker_inverse = deduplicate(marker_sequence)
            # the sequence file does not allow separate specification of the marker
            # and waveform files. we have to write both out when either is different
            seq_unique, seq_labels, seq_inverse = deduplicate(zip(wfm_inverse, 
                                                                  marker_inverse))

            for seq_idxs, seq_label in zip(seq_unique, seq_labels):
                wfm_idx, marker_idx = seq_idxs
                self.export_pattern(wfm_unique[wfm_idx], 
                                    marker_unique[marker_idx][0],
                                    marker_unique[marker_idx][1], 
                                    lowhighs[ch_idx], ch_idx, seq_label)
            seq_labelss.append([seq_labels[idx] for idx in seq_inverse])

        self.export_sequence_file(seq_labelss)
        
        

class Agilent81180AExport(object):
    '''
    Agilent 81180A Exporter
    
    '''
    lowhigh = (-1, 1)
    marker1_mask = 2**0
    marker2_mask = 2**1
    
    def __init__(self, directory, fileprefix):
        self.dir = directory
        self.fileprefix = fileprefix
        self._suffix = ".bin"
        
    def _pack_markers(self, marker1, marker2):
        marker1 = np.array(np.round(marker1), dtype=bool)
        marker2 = np.array(np.round(marker2), dtype=bool)
        markers = marker1 * Agilent81180AExport.marker1_mask
        markers = np.bitwise_or(
            markers, marker2 * Agilent81180AExport.marker2_mask)
        return np.array(markers, dtype=np.uint8)
        
    def export_pattern(self, waveform, marker1, marker2, lowhigh,
                       ch_idx=None, number=None):
        '''
        Save the waveform and markers in a pattern file.
        
        Waveform and markers use two separate files.
        
        The waveform is saved in a binary file, where each data point is
        represented by a 32bit wide floating point number between -1 and 1.
        
        The markers are saved in another binary file, where each data point is
        a 8 bit integer, where the last two bits define the markers. Bit 6 is
        marker one, and bit 7 is marker 2 (counting the bits from 0 - 7).
        
        '''
        waveform = np.array(waveform, np.float32)
        
        if len(marker2) != len(marker1):
            raise ValueError('markers must have same length')
        
        if len(waveform) < 320:
            raise ValueError('waveform must have at least 320 samples')
        
        remainder = len(waveform) % 32
        if remainder != 0:
            logging.warning(('Waveform length is not a multiple of 32.'
                            ' Waveform is automatically padded at the end with'
                            ' the last element to fulfill this requirement.'
                            ))
            waveform = np.pad(waveform, (0, 32 - remainder), 'edge')
            marker1 = np.pad(marker1, (0, 8 - len(marker1)%8), 'edge')
            marker2 = np.pad(marker2, (0, 8 - len(marker2)%8), 'edge')
            
        if 4 * len(marker1) != len(waveform):
            raise ValueError('len(waveform) must be equal to 4 * len(marker)')
            
        fileprefix = pat_name_format(self.fileprefix, ch_idx + 1, number)
        filename = os.path.join(self.dir, fileprefix + self._suffix)
        
        if not os.path.isdir(self.dir):
            os.makedirs(self.dir)
        with open(filename, 'wb') as wfm_file:
            wfm_file.write(waveform)
        
        markers = self._pack_markers(marker1, marker2)
        marker_filename = os.path.join(
            self.dir, 'marker_' + fileprefix + self._suffix)
        with open(marker_filename, 'wb') as marker_file:
            marker_file.write(markers.data)
            
    def export_sequence_file(self, ch_labels):
        '''
        Save the sequence file. This sequence file is intended to look similar
        to the one used for the Tektronix 5014 arbitrary waveform generator.
        
        '''
        def string_format(ch1_label, ch2_label):
            string = ('"' + pat_name_format(self.fileprefix, 1, ch1_label)
                      + self._suffix + '",'
                      '"' + pat_name_format(self.fileprefix, 2, ch2_label)
                      + self._suffix + '"'
                      +'\r\n')
            return string
        
        seq_filename = os.path.join(self.dir, self.fileprefix + '.seq')
        if not os.path.isdir(self.dir):
            os.mkdir(self.dir)
        with open(seq_filename, 'wb') as seq_file:
            for ch1_label, ch2_label in zip(*ch_labels):
                seq_file.write(string_format(ch1_label, ch2_label))
                
    def export_sequence(self, wfm_sequences, marker_sequences, lowhighs, rmdir=True):
        '''
        Save all waveforms and markers in pattern files and create a 
        sequence file.
        
        Parameters:
        -----------
        wfm_sequences: List
            List of sequences, where each entry with index `idx` defines
            the sequence to be played on the AWG channel `idx`.
        marker_sequences: List
            List of marker sequences, where each entry with index `idx` defines
            the marker sequence to be played on the AWG markers corresponding to
            the AWG channel `idx`.
        lowhighs: List
            Each entry defines the low and high levels of the corresponding
            AWG channel, used to rescale the waveform. See also the help of
            the rescale function.
        rmdir: Bool
            If True, empty existing sequence directory before exporting
        '''
        num_channels = len(wfm_sequences)
        if num_channels != 2 or num_channels != len(marker_sequences):
            raise ValueError('`wfm_sequence` and `marker_sequence` should have'
                             ' exactly two channels')
        
        num_patterns = [len(wfm_sequence) for wfm_sequence in wfm_sequences]
        if min(num_patterns) != max(num_patterns):
            raise ValueError('sequences must have equal length on all channels' 
                             'of the device.')
        
        if not os.path.isdir(self.dir):
            os.makedirs(self.dir)
        else:
            if rmdir:
                shutil.rmtree(self.dir)
                # TODO: Get rid of this sleep
                sleep(0.1)
                os.makedirs(self.dir)
        
        # file names for all channels
        seq_labelss = []
        for ch_idx in range(num_channels):
            wfm_sequence = wfm_sequences[ch_idx]
            marker_sequence = marker_sequences[ch_idx]
            
            if len(wfm_sequence) != len(marker_sequence):
                raise ValueError('`wfm_sequence` and `marker_sequence` must'
                                 ' have equal length')
                
            wfm_unique, _, wfm_inverse = deduplicate(wfm_sequence)
            marker_unique, _, marker_inverse = deduplicate(marker_sequence)
            # the sequence file does not allow separate specification of the marker
            # and waveform files. we have to write both out when either is different
            seq_unique, seq_labels, seq_inverse = deduplicate(zip(wfm_inverse, 
                                                                  marker_inverse))

            for seq_idxs, seq_label in zip(seq_unique, seq_labels):
                wfm_idx, marker_idx = seq_idxs
                self.export_pattern(wfm_unique[wfm_idx], 
                                    marker_unique[marker_idx][0],
                                    marker_unique[marker_idx][1], 
                                    lowhighs[ch_idx], ch_idx, seq_label)
            seq_labelss.append([seq_labels[idx] for idx in seq_inverse])

        self.export_sequence_file(seq_labelss)