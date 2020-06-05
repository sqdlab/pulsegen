'''
Created on 04/07/2013

@author: Matthias Baur
'''
import numpy as np
from . import pulse
import os
import logging
from scipy import signal
from .config import cfg
from .optimalcontrol import OptimalControl
from .exporter import export
from copy import deepcopy


class AWG(object):
    def __init__(self, model, ch_count, marker_count, granularity):
        self.model = model
        self.ch_count = ch_count
        self.marker_count = marker_count # per channel
        self.granularity = granularity

AWG_MAP = {
    'Tektronix_5014': AWG(model='Tektronix_5014',
                          ch_count=4, marker_count=2, granularity=1),
    'Tektronix_520': AWG(model='Tektronix_520',
                         ch_count=2, marker_count=2, granularity=1),
    'Agilent_N824x': AWG(model='Agilent_N824x',
                         ch_count=2, marker_count=2, granularity=8),
    'Agilent_81180A': AWG(model='Agilent_81180A',
                         ch_count=2, marker_count=2, granularity=4),
    }


class Filter(object):
    def __init__(self, window, cutoff, sampling_rate=None):
        self.window = window
        self.cutoff = cutoff
        self.Fs = sampling_rate
        
    def get_fir_filter(self):
        
        if self.window == "gaussian":
            # std_time * std_freq = 1/(2*pi)
            std = self.Fs / self.cutoff # Fs in Hz, self.cutoff in rad/s
            length = int(6*std)
            b = signal.gaussian(length, std)
            # Normalize the fir_filter such that filtered output
            # has same amplitude
            return b/sum(b)
        else:
            raise NotImplementedError(
                'fir window `{0:s}` not implemented.'.format(self.window))
        
    def filter(self, input):
        b = self.get_fir_filter()
        a = [1]
        return signal.filtfilt(b, a, input)


class Segment(object):
    '''
    classdocs
    
    '''
    def __init__(self, waveform=None, markers=None):
        '''
        Constructor
        
        '''
        if waveform is None:
            self.waveform = []
        else:
            self.waveform = waveform
        
        if markers is None:
            self.markers = []
        else:
            self.markers = markers


class Sampler(object):
    # TODO: Make fixed point array
    # TODO: Remove sampler, and add sampling to Channel class
    def __init__(self, channel):
        
        self._sampling_freqs = [channel.sampling_freq]
        self._sampling_freqs.extend(
            [channel.sampling_freq / channel.awg.granularity
             for __ in range(channel.awg.marker_count)]
            )
        self._channel = channel
        self._cache = {}
    
    def sample(self, segment):
        sampled_wfm = self._sample_train(
            segment.waveform, self._sampling_freqs[0],
            self._channel.fixed_points[0] + self._channel.delay
            )
        sampled_markers = []
        
        for idx in range(self._channel.awg.marker_count):
            try:
                marker = segment.markers[idx]
            except IndexError:
                marker = []
            
            marker_fixed_pt = (self._channel.fixed_points[idx+1]
                               + self._channel.marker_delays[idx])
            sampled_markers.append(
                self._sample_train(marker, self._sampling_freqs[idx+1],
                                   marker_fixed_pt, marker=True)
                )
        
        return Segment(waveform=sampled_wfm, markers=sampled_markers)
        
    def _sample_train(self, train, sampling_freq, fixed_point, marker=False):
        sampled_train = np.zeros(
            int(self._channel.pattern_length*sampling_freq), dtype=np.float32)
        end_time = fixed_point
            
        for pulse in reversed(train):
            (sampled_train, end_time) = self._add_sampled_pulse(
                pulse, end_time, sampled_train, sampling_freq, marker=marker)
        
        return sampled_train
    
    def _add_sampled_pulse(self, mypulse, end_time, pattern, sampling_freq,
                           marker=False):
        '''Sample the pulse `mypulse` and add it to the pattern `pattern`'''
        cache_id = mypulse.cache_id()
        
        if marker:
            (sampled_pulse, start_time, start_idx, stop_idx) = \
                pulse.sample_pulse(mypulse, end_time, sampling_freq)
        
        elif cache_id in self._cache and mypulse.optimal_control:
            sampled_pulse = self._cache[cache_id] * mypulse.amp
            start_idx = int(np.ceil(
                (end_time - self._channel.optimal_control.get_prepend()
                 - mypulse.len) * sampling_freq))
            
            stop_idx = start_idx + len(sampled_pulse)
            start_time = end_time - mypulse.len - mypulse.sep
        
        else:
            if (mypulse.optimal_control
                and (self._channel.optimal_control is not None)):
                if cache_id is not None:
                    pulse_amp = mypulse.amp
                    mypulse.amp = 1
                    
                sampled_pulse = self._channel.optimal_control(mypulse)
                
                start_idx = int(np.ceil(
                    (end_time - self._channel.optimal_control.get_prepend()
                     - mypulse.len) * sampling_freq))
                
                stop_idx = start_idx + len(sampled_pulse)
                start_time = end_time - mypulse.len - mypulse.sep
                
                if cache_id is not None:
                    print('caching: ' + str(mypulse))
                    self._cache[cache_id] = sampled_pulse
                    mypulse.amp = pulse_amp
                    sampled_pulse = sampled_pulse * pulse_amp
                
            else:
                (sampled_pulse, start_time, start_idx, stop_idx) = \
                    pulse.sample_pulse(mypulse, end_time, sampling_freq)
                
#                 if self._channel.filter is not None:
#                     filter = self._channel.filter
#                     pad = filter.padding
#                     sampled_pulse = np.pad(
#                         sampled_pulse, (pad, pad), 'constant')
#                     sampled_pulse = filter.filter(sampled_pulse)
#                     start_idx -= pad
#                     stop_idx += pad

        pattern[start_idx : stop_idx] += sampled_pulse
        return (pattern, start_time)


class Channel(object):
    '''
    classdocs
    
    '''
    def __init__(self, awg, sampling_freq, fixed_points, pattern_length, 
                 filter=None, optimal_control=None, delay=0, marker_delays=0,
                 lowhigh=(-1, 1)):
        self.awg = awg
        self.sampling_freq = sampling_freq
        self.filter = filter
        self.optimal_control = optimal_control
        self.delay = delay
        self.marker_delays = marker_delays
        self.pattern_length = pattern_length
        self.lowhigh = lowhigh
        
        if not isinstance(fixed_points, (list, tuple)):
            self.fixed_points = [fixed_points
                                 for __ in range(self.awg.marker_count + 1)]
        else:
            self.fixed_points = fixed_points


class Sequence(object):
    ADV_MODE_MAP = {
        "auto": 0,
        "continuous": 1,
        "play_one_repetition": 2,
        "play_all_repetitions": 3}
    
    MARKER_MASK_MAP = {
        "off": 0,
        "start": 1,
        "repeat": 2,
        "gate": 3}
    
    COUNT = 1
    ADV_MODE = ADV_MODE_MAP['play_one_repetition']
    MARKER_MASK = MARKER_MASK_MAP['off']
    
    def __init__(self, segments=None, counts=None, adv_modes=None,
                 marker_masks=None):
        self.segments = []
        self.counts = []
        self.adv_modes = []
        self.marker_masks = []
        self._markers_counter = 0
        self._wfm_counter = 0
        
        if segments is not None:
            self.segments = segments
        if counts is not None:
            self.counts = counts
        if adv_modes is not None:
            self.adv_modes = adv_modes
        if marker_masks is not None:
            self.marker_masks = marker_masks
    
    def __len__(self):
        return len(self.segments)
    
    @property
    def waveforms(self):
        return [segment.waveform for segment in self.segments]
    
    @property
    def markers(self):
        return [segment.markers for segment in self.segments]
    
    def append_segment(self, segment, count=COUNT, adv_mode=ADV_MODE,
                       marker_mask=MARKER_MASK):
        self.segments.append(segment)
        self.counts.append(count)
        self.adv_modes.append(adv_mode)
        self.marker_masks.append(marker_mask)
    
    def set_segment(self, idx, segment, count=COUNT, adv_mode=ADV_MODE,
                    marker_mask=MARKER_MASK):
        self.segments[idx] = segment
        self.counts[idx] = count
        self.adv_modes[idx] = adv_mode
        self.marker_masks[idx] = marker_mask
    
    def set_pulses(self, idx, pulse_train, mode=None):
        if (mode is None) or (mode == 'replace'):
            self.segments[idx].waveform = pulse_train
        elif mode == 'append':
            self.segments[idx].waveform.extend(pulse_train)
        elif mode == 'prepend':
            self.segments[idx].waveform = (pulse_train + 
                                           self.segments[idx].waveform)
        else:
            raise ValueError('Unsupported mode {0}.'.format(mode))
    
    def set_markers(self, idx, marker_trains, mode=None):
        if (mode is None) or (mode == 'replace'):
            self.segments[idx].markers = marker_trains
        elif (mode == 'append') or (mode == 'prepend'):
            self.segments[idx].markers = [
                (mt1 + mt2) if (mode == 'append') else (mt2 + mt1)
                for mt1, mt2 in zip(self.segments[idx].markers, marker_trains)
            ]
        else:
            raise ValueError('Unsupported mode {0}.'.format(mode))
    
    def append_pulses(self, pulse_train, count=COUNT, adv_mode=ADV_MODE,
                      marker_mask=MARKER_MASK):
        if self._wfm_counter >= self._markers_counter:
            self.append_segment(
                Segment(waveform=pulse_train), count, adv_mode, marker_mask)
        else:
            # TODO: How to handle count, adv_mode and marker_mask here?
            self.segments[self._markers_counter].wfm = pulse_train
        
        self._wfm_counter += 1
    
    def append_markers(self, marker_trains, count=COUNT, adv_mode=ADV_MODE,
                       marker_mask=MARKER_MASK):
        if self._markers_counter >= self._wfm_counter:
            self.append_segment(
                Segment(markers=marker_trains), count, adv_mode, marker_mask)
        else:
            # TODO: How to handle count, adv_mode and marker_mask here?
            self.segments[self._markers_counter].markers = marker_trains
        
        self._markers_counter += 1
    
    def sample(self, sampler):
        sampled_seq = Sequence()
        for param in zip(self.segments, self.counts, self.adv_modes,
                         self.marker_masks):
            sampled_seq.append_segment(sampler.sample(param[0]), *param[1:])
        
        return sampled_seq


class AWGSequence(object):
    def __init__(self, awg, id=None, channels=None, sequences=None):
        self.id = id
        self.awg = awg
        self.sampled_sequences = []
        
        if channels is None:
            self.channels = []
        else:
            self.channels = channels
        
        if sequences is None:
            self.sequences = []
        else:
            self.sequences = sequences
        
        if (channels is not None) and (sequences is not None):
            if len(sequences) != len(channels):
                raise ValueError('`sequences` and `channels` must have'
                                 ' equal length.')
        
        for ch in self.channels:
            self.sampling_freqs = [ch.sampling_freq]
            self.sampling_freqs.extend(
                [ch.sampling_freq/self.awg.granularity
                 for __ in self.awg.marker_count]
                )
    
    def sample(self):
        max_segments = max([len(seq) for seq in self.sequences])
        
        self.sampled_sequences = []
        for ch_idx, ch in enumerate(self.channels):
            seq = self.sequences[ch_idx]
            if len(seq) != max_segments:
                
                # Just to prevent too many warning messages. Most channels
                # are zero in normal experiments
                if len(seq) != 0:
                    logging.warning(__name__ + 'Number of segments not equal'
                                  ' for all sequences of AWG sequence {0}.'
                                  ' Appending segments with no waveform and'
                                  ' no markers to make all sequences of equal'
                                  ' length.'.format(self.id))
                
                [seq.append_segment(Segment()) for __ in
                 range(max_segments - len(seq))]
            
            sampler = Sampler(ch)
            sampled_seq = seq.sample(sampler)
            self.sampled_sequences.append(sampled_seq)
    
    def append_sequence(self, sequence):
        self.sequences.append(sequence)
    
    def append_channel(self, channel):
        self.channels.append(channel)


class MultiAWGSequence(object):
    def __init__(self):
        
        #self._cfgs = cfg.get_config()
        self._channel_cfgs = cfg.get_channel_config()
        self._channel_pair_cfgs = cfg.get_channel_pair_config()

        self._sequences = []
        self._channels = []
        self._awg_sequences = []
        self._ch_awg_map = {}
        
        self._init_channels()
        
    def _init_channels(self):
        for ch_cfg in self._channel_cfgs:
            if ch_cfg['use_fir']:
                filter = Filter(ch_cfg['fir_window'], 
                                ch_cfg['fir_cutoff'],
                                ch_cfg['sampling_rate'])
            else:
                filter = None
            
            if (#(ch_cfg['channel_type'] == 'single') and
                ch_cfg['use_optimal_control']
                ):
                oc = OptimalControl(ch_cfg['sampling_rate'],
                                    ch_cfg['sampling_rate_multiplier'],
                                    ch_cfg['ir_tlim'],
                                    ch_cfg['prepend'],
                                    ch_cfg['append'],
                                    ch_cfg['use_pinv'], filter)
                oc.load(ch_cfg['impulse_response'])
            else:
                oc = None
              
            channel = Channel(AWG_MAP[ch_cfg['awg_models']],
                              ch_cfg['sampling_rate'],
                              ch_cfg['fixed_point'],
                              ch_cfg['pattern_length'],
                              filter=filter,
                              optimal_control=oc,
                              delay=ch_cfg['delay'],
                              marker_delays=ch_cfg['marker_delays'],
                              lowhigh=ch_cfg['lowhigh'])
            
            self._channels.append(channel)
            self._sequences.append(Sequence())
    
    @property
    def channels(self):
        return self._channels
    
    @property
    def sampled_sequences(self):
        sequences = []
        for awg_seq in self._awg_sequences:
            sequences.extend(awg_seq.sampled_sequences)
        
        return sequences
    
    @property
    def sequences(self):
        return self._sequences
    
    def _create_pulse_list(self, pulse_train, cfg, chpair=False):
        if chpair:
            pulse_list = []
            [pulse_list.append(p.get_pulse(cfg)) for p in pulse_train]
            pulse_list = list(map(list, zip(*pulse_list)))
        else:
            pulse_list = []
            for p in pulse_train:
                if isinstance(p, pulse.Pulse):
                    pulse_list.append(p.get_pulse(cfg))
                else:
                    pulse_list.append(p)
          
        return pulse_list
    
    def _get_marker_pulse_lists(self, marker_trains, ch):
        pulse_list = []
        for train in marker_trains:
            pulse_list.append(self._create_pulse_list(
                    train, self._channel_cfgs[ch], chpair=False))
          
        return pulse_list
    
    def _get_pulse_lists(self, pulse_train, ch=None, chpair=None):
        if ch is not None:
            pulse_list = self._create_pulse_list(
                pulse_train, self._channel_cfgs[ch], chpair=False)
          
        elif chpair is not None:
            pulse_list = self._create_pulse_list(
                pulse_train, self._channel_pair_cfgs[chpair], chpair=True)
          
        else:
            ValueError("Either of `chpair` and `ch` must be specified")
          
        return pulse_list
    
    def append_pulses(self, pulse_train, ch=None, chpair=None, **kwargs):
        pulse_lists = self._get_pulse_lists(pulse_train, ch=ch, chpair=chpair)
        
        if chpair is not None:
            for ii, pulse_list in enumerate(pulse_lists):
                self._sequences[2*chpair + ii].append_pulses(
                                                pulse_list, **kwargs)
        
        elif ch is not None:
            self._sequences[ch].append_pulses(pulse_lists, **kwargs)
    
    def append_markers(self, marker_trains, ch, **kwargs):
        pulse_list = self._get_marker_pulse_lists(marker_trains, ch)
        self._sequences[ch].append_markers(pulse_list, **kwargs)
    
    def set_pulses(self, pulse_train, segment_idx, ch=None, chpair=None, mode=None):
        """
        Update analog pulses of a certain segment and channel.
        
        Parameters
        ----------
        pulse_train : [`Pulse`]
            Analog pulses to set/prepend/append.
        segment_idx : `int`
            Segment modified.
        ch : `int`
            Channel modified.
        chpair : `int`
            Channel pair modified. Takes precedence over `ch`.
        mode : `str`, default 'replace'
            One of 'replace', 'prepend', append'.
        """
        pulse_lists = self._get_pulse_lists(pulse_train, ch=ch, chpair=chpair)
        
        if chpair is not None:
            for ii, pulse_list in enumerate(pulse_lists):
                self._sequences[2*chpair + ii].set_pulses(
                    segment_idx, pulse_list, mode
                )
        
        elif ch is not None:
            self._sequences[ch].set_pulses(segment_idx, pulse_lists, mode)
      
    def set_markers(self, marker_trains, segment_idx, ch, mode=None):
        """
        Update markers of a certain segment and channel.
        
        Parameters
        ----------
        marker_trains : [[`Pulse`], [`Pulse`]]
            Marker pulses to set/prepend/append.
        segment_idx : `int`
            Segment modified.
        ch : `int`
            Channel modified.
        mode : `str`, default 'replace'
            One of 'replace', 'prepend', append'.
        """
        pulse_list = self._get_marker_pulse_lists(marker_trains, ch)
        self._sequences[ch].set_markers(segment_idx, pulse_list, mode)
    
    def equalize(self, mode='prepend'):
        """
        Equalize the duration of all channel sequences by prepending 
        appropriate `Spacer` pulses.
        
        Parameters
        ----------
        mode : `str`, default 'prepend'
            'prepend' to pad left, 'append' to pad right.
        """
        # determine length of all segments with the same index
        seq_count = len(self._sequences)
        seg_count = len(self)
        pulse_durations = np.zeros((seq_count, seg_count))
        marker_durations = np.zeros((seq_count, seg_count, 2))
        for seq_idx, sequence in enumerate(self.sequences):
            for seg_idx, segment in enumerate(sequence.segments):
                duration = sum(pulse.duration for pulse in segment.waveform)
                pulse_durations[seq_idx, seg_idx] = duration
                for marker_idx, marker_train in enumerate(segment.markers):
                    duration = sum(pulse.duration for pulse in marker_train)
                    marker_durations[seq_idx, seg_idx, marker_idx] = duration
        max_durations = np.max((pulse_durations.max(0), 
                              marker_durations.max(2).max(0)), axis=0)
        # pad sequences
        if max_durations.max() == 0.:
            return
        pulse_paddings = max_durations - pulse_durations
        marker_paddings = max_durations[:, np.newaxis] - marker_durations
        for seg_idx in range(seg_count):
            for seq_idx, sequence in enumerate(self.sequences):
                # pad analog pulse train
                if len(self.sequences[seq_idx]) <= seg_idx:
                    self.sequences[seq_idx].append_segment(Segment())
                pulse_padding = pulse_paddings[seq_idx, seg_idx]
                if pulse_padding != 0.:
                    spacer = pulse.spacer(pulse_padding)
                    self.set_pulses([spacer], seg_idx, ch=seq_idx, mode=mode)
                # pad marker pulse trains
                marker_padding = marker_paddings[seq_idx, seg_idx, :]
                if np.any(marker_padding != 0.):
                    spacer = [[pulse.spacer(padding)] 
                              for padding in marker_padding]
                    self.set_markers(spacer, seg_idx, ch=seq_idx, mode=mode)
    
    def _distribute_sequences(self):
        # TODO: Check somewhere whether each sequence in an awg sequence has
        # correct amount of waveforms, probably that should be handled
        # by the awg sequence class
        awg_ch_counter = 1
        awg_counter = 0
        
        self._awg_sequences = []
        for ch_idx, ch in enumerate(self._channels):
            current_awg = AWG_MAP[self._channel_cfgs[ch_idx]['awg_models']]
            
            if ch_idx == 0:
                previous_awg = current_awg
            
            if ((awg_ch_counter > current_awg.ch_count)
                | (current_awg != previous_awg)):
                
                awg_ch_counter = 1
            
            if awg_ch_counter == 1:
                awg_seq = AWGSequence(current_awg, id=awg_counter)
                self._awg_sequences.append(awg_seq)
                awg_counter += 1
            
            awg_seq.append_sequence(self._sequences[ch_idx])
            awg_seq.append_channel(ch)
            self._ch_awg_map[ch] = (awg_counter, awg_ch_counter)
            awg_ch_counter += 1
            previous_awg = current_awg
    
    def sample(self):
        self._distribute_sequences()
        [awg_seq.sample() for awg_seq in self._awg_sequences]
    
    def export(self, directory, fileprefix, rmdir=True):
        for idx, awg_seq in enumerate(self._awg_sequences):
            awg_directory = os.path.join(directory, 'AWG_{0:0=2d}'.format(idx))
            export(awg_seq, awg_directory, fileprefix, rmdir)
            
    def __len__(self):
        ''' 
            return number of waveforms of the sequence.
        '''
        if not len(self.sequences):
            return 0
        return max([sum(s.counts) for s in self.sequences])

    def append_tomography_sequence(self, qubit_param_list):
        """
        Construct a quorum of tomography pulses from a list of qubit control parameters
        Expects [[q1_amp, q1_sigma, q1_omega, q1_chpair], [q2_amp, q2_sigma, q2_omega, q2_chpair], ...[]]
        Returns the "outer product" of the tomography manipulation pulse, where each pulse is determined
        by the parameters in the qubit_parameters_list. The first qubit in the list will be the one which
        has the tomography sequencing change the fastest. The last in the list will be the qubit that has
        the manipulation pulses repeat the most before moving on to the next pulse.
        There is no assumption that the qubit control parameters will define pulses with consistent duration.
        All pulses are padded to match the maximum duration using a Spacer pulse added at then end.
    
        Modify the provided sequence by appending the each element of the appendix sequence to the original sequence waveform.
    
        """
        
        num_qubits = len(qubit_param_list)
        quorum_length = 4**num_qubits
    
        # Create an empty sequence holder
        tomo_seq = MultiAWGSequence()
        
        for di in range(num_qubits):
            
            # extract the values from the supplied list
            # ### This is the part that will need to change when the qubit parameters are stored as "Instruments"
            
            amplitude = qubit_param_list[di][0]
            sigma = qubit_param_list[di][1]
            omega = qubit_param_list[di][2]
            chpair = qubit_param_list[di][3]
            # mixer_calibration ...
            
            ident = {'envelope':pulse.GaussianPulse, 
            'amplitude': 0,
            'sigma':sigma,
            'phase':0, 
            'omega':omega}
    #            'mixer_calibration':mixer_calibration(mixer_ratio, mixer_phase)}    
            
            pix = ident.copy()
            pix.update({'amplitude':amplitude})
            pihalfx = ident.copy()
            pihalfx.update({'amplitude':amplitude/2})
            piy = ident.copy()
            piy.update({'amplitude':amplitude,'phase':-np.pi/2})
            pihalfy = ident.copy()
            pihalfy.update({'amplitude':amplitude/2,'phase':-np.pi/2})      
            
            quorum_template = [ident, pihalfx, pihalfy, pix]
    
            for q in range(quorum_length):
                # work out which manipulation pulse (from the quorum) goes in the corresponding segment slot
                # this is 01230123... or 00001111... etc
                index = q/(4**di) % 4
                # append the manipulation pulse to the channel pair corresponding to that qubit            
                tomo_seq.append_pulses(pulse_train=[pulse.Pulse(pulse.GaussianPulse)(**quorum_template[index])], chpair=chpair)
                
        # After the whole sequence is built, we can scan the structure to find the longest tomography pulse used anywhere
        # We couldn't do this beforehand, because the Pulse is not evaluated until it is added the sequence
                
        max_duration = 0.0
        
        for sequence in tomo_seq.sequences:
            for segment in sequence.segments:
                duration = 0.0
                for waveform in segment.waveform:
                    duration += waveform.duration
                if duration > max_duration:
                    max_duration = duration
                            
        # The whole structure is now traversed again, appending a Spacer onto any segment with a duration shorter
        # than the maximum duration
                        
        for sequence in tomo_seq.sequences:
            if sequence.segments == []:
                # There is no manipulation required on this qubit, but we DO need to
                # create a Spacer to pad out each segment to the correct duration
                for segment in range(quorum_length):
                    sequence.append_pulses(pulse_train=[pulse.Spacer(length = max_duration)])
    
            else:
                # For all the segments, ensure that the durations are padded to the maximum duration
                for segment in sequence.segments:
                    duration = 0.0
                    for waveform in segment.waveform:
                        duration += waveform.duration
                    if duration < max_duration:
                        spacer = pulse.Spacer(length = max_duration-duration)
                        segment.waveform.append(spacer)
                
        # This is kind of like a tensor product for sequences
        # The situation is complicated in that the original waveform list is stored part of a segment which is part of a sequence,
        # and we have to drill down to get to the appropriate lists.
        
        # When the pulse is "sampled" (ie converted from a descriptive list to a series of samples that are programmed into the AWG)
        # the list that describes the pulses is traversed backwards. This means that the tomography pulses that have been appended
        # push the original pulses further back in time relative to the start of the measurement. We don't need to make any further
        # adjustments.
            
        # The outer loop effectively refers to AWG channel numbers,
        # so the two inputs have to have this outer loop unrolled together
        
        for this_sequence, appendix_seq in zip(self.sequences, tomo_seq.sequences):
            
            new_segment_list = []
    
            new_counts = []
            new_wfm_counter = 0
            
            new_marker_masks = []
            new_markers_counter = 0
    
            new_adv_modes = []
    
            # Each sequence has seperate lists for the pulse control structure
            # We need to unroll these together to be able to build extened versions.
            # This is all because the sequence uses a structure of lists rather than a list of structure
            
            for this_segment, this_count, this_marker_mask, this_adv_mode in zip(
                this_sequence.segments,
                this_sequence.counts,
                this_sequence.marker_masks,
                this_sequence.adv_modes):
    
                # The reason for the following assertion is that the tomography pulse is appended
                # to the segment and that whole segment will be repeated.
                assert (this_count == 1), "We can't adapt a sequence that uses replicated waveforms"
    
                # For this segment, we are going to build several new segments in the new segment list
                for segment_id, this_appendix_segment in enumerate(appendix_seq.segments):
                    # copy the original segment non-destructively (shallow copy is not sufficient)
                    new_segment = deepcopy(this_segment)
                    # Append the appendix waveform(s) and update the pulse count summaries
                    new_segment.waveform.extend(this_appendix_segment.waveform)
                    
                    # Special case for the sequence start marker. We only want the pulse to be copied for
                    # segment zero. Otherwise force the first marker to blank (if it exists at all)
    
                    if (segment_id > 0):
                        if new_segment.markers != []:
                            new_segment.markers[0] = []
                            
                    # add this to the new segment list we're building
                    new_segment_list.append(new_segment)
                    new_counts.append(this_count)
                    new_marker_masks.append(this_marker_mask)
                    new_adv_modes.append(this_adv_mode)
                    
                    new_wfm_counter     += 1 # Not sure what these are supposed to count so it's difficult to make it right.
                    new_markers_counter += 1 # Also there seems to be a bug in the original
                    
            # now attach the new extended lists back onto the original sequence channel
            this_sequence.segments = new_segment_list
            this_sequence.counts = new_counts
            this_sequence.adv_modes = new_adv_modes
            this_sequence.marker_masks = new_marker_masks
            this_sequence._markers_counter = new_markers_counter
            this_sequence._wfm_counter = new_wfm_counter
    





# The following definitions are outside the class scope. They are left here because it might be preferable to use them
# as conventional function calls rather as methods.

def construct_tomography_sequence(qubit_param_list):
    """
    Construct a quorum of tomography pulses from a list of qubit control parameters
    Expects [[q1_amp, q1_sigma, q1_omega, q1_chpair], [q2_amp, q2_sigma, q2_omega, q2_chpair], ...[]]
    Returns the "outer product" of the tomography manipulation pulse, where each pulse is determined
    by the parameters in the qubit_parameters_list. The first qubit in the list will be the one which
    has the tomography sequencing change the fastest. The last in the list will be the qubit that has
    the manipulation pulses repeat the most before moving on to the next pulse.
    There is no assumption that the qubit control parameters will define pulses with consistent duration.
    All pulses are padded to match the maximum duration using a Spacer pulse added at then end and all
    channels that are not manipulated explicitly are padded using Spacers of the corresponding duration.
    This sequence is then ready to append to an arbitrary measurement sequence.
    """

    num_qubits = len(qubit_param_list)
    quorum_length = 4**num_qubits    

    # Create an empty sequence holder
    tomo_seq = MultiAWGSequence()
    
    for di in range(num_qubits):
        
        # extract the values from the supplied list
        # ### This is the part that will need to change when the qubit parameters are stored as "Instruments"
        
        amplitude = qubit_param_list[di][0]
        sigma = qubit_param_list[di][1]
        omega = qubit_param_list[di][2]
        chpair = qubit_param_list[di][3]
        # mixer_calibration ...
        
        ident = {'envelope':pulse.GaussianPulse, 
        'amplitude': 0,
        'sigma':sigma,
        'phase':0, 
        'omega':omega}
#            'mixer_calibration':mixer_calibration(mixer_ratio, mixer_phase)}    
        
        pix = ident.copy()
        pix.update({'amplitude':amplitude})
        pihalfx = ident.copy()
        pihalfx.update({'amplitude':amplitude/2})
        piy = ident.copy()
        piy.update({'amplitude':amplitude,'phase':-np.pi/2})
        pihalfy = ident.copy()
        pihalfy.update({'amplitude':amplitude/2,'phase':-np.pi/2})      
        
        quorum_template = [ident, pihalfx, pihalfy, pix]

        for q in range(quorum_length):
            # work out which manipulation pulse (from the quorum) goes in the corresponding segment slot
            # this is 01230123... or 00001111... etc
            index = q/(4**di) % 4
            # append the manipulation pulse to the channel pair corresponding to that qubit            
            tomo_seq.append_pulses(pulse_train=[pulse.Pulse(pulse.GaussianPulse)(**quorum_template[index])], chpair=chpair)
            
    # After the whole sequence is built, we can scan the structure to find the longest tomography pulse used anywhere
    # We couldn't do this beforehand, because the Pulse is not evaluated until it is added the sequence
            
    max_duration = 0.0
    
    for sequence in tomo_seq.sequences:
        for segment in sequence.segments:
            duration = 0.0
            for waveform in segment.waveform:
                duration += waveform.duration
            if duration > max_duration:
                max_duration = duration
                        
    # The whole structure is now traversed again, appending a Spacer onto any segment with a duration shorter
    # than the maximum duration
                    
    for sequence in tomo_seq.sequences:
        if sequence.segments == []:
            # There is no manipulation required on this qubit, but we DO need to
            # create a Spacer to pad out each segment to the correct duration
            for segment in range(quorum_length):
                sequence.append_pulses(pulse_train=[pulse.Spacer(length = max_duration)])

        else:
            # For all the segments, ensure that the durations are padded to the maximum duration
            for segment in sequence.segments:
                duration = 0.0
                for waveform in segment.waveform:
                    duration += waveform.duration
                if duration < max_duration:
                    spacer = pulse.Spacer(length = max_duration-duration)
                    segment.waveform.append(spacer)
            
    return tomo_seq

def append_tomography_sequence( orig_sequence, quorum ):
    """
    Modify the provided sequence by appending the each element of the appendix sequence to the original sequence waveform.
    For example [[wav1],[wav2],[wav3]] [[app1],[app2]] becomes [[wav1,app1],[wav1,app2],[wav2,app1],[wav2,app2],[wav3,app1],[wav3,app2]].
    Both parameters are of type MultiAWGSequence.
    """
    # This is kind of like a tensor product for sequences
    # The situation is complicated in that the original waveform list is stored part of a segment which is part of a sequence,
    # and we have to drill down to get to the appropriate lists.
    
    # When the pulse is "sampled" (ie converted from a descriptive list to a series of samples that are programmed into the AWG)
    # the list that describes the pulses is traversed backwards. This means that the tomography pulses that have been appended
    # push the original pulses further back in time relative to the start of the measurement. We don't need to make any further
    # adjustments.
        
    # The outer loop effectively refers to AWG channel numbers,
    # so the two inputs have to have this outer loop unrolled together
    
    for this_sequence, appendix_seq in zip(orig_sequence.sequences, quorum.sequences):
        
        new_segment_list = []

        new_counts = []
        new_wfm_counter = 0
        
        new_marker_masks = []
        new_markers_counter = 0

        new_adv_modes = []

        # Each sequence has separate lists for the pulse control structure
        # We need to unroll these together to be able to build extended versions.
        # This is all because the sequence uses a structure of lists rather than a list of structure
        
        for this_segment, this_count, this_marker_mask, this_adv_mode in zip(
            this_sequence.segments,
            this_sequence.counts,
            this_sequence.marker_masks,
            this_sequence.adv_modes):

            # The reason for the following assertion is that the tomography pulse is appended
            # to the segment and that whole segment will be repeated.
            assert (this_count == 1), "We can't adapt a sequence that uses replicated waveforms"

            # For this segment, we are going to build several new segments in the new segment list
            for segment_id, this_appendix_segment in enumerate(appendix_seq.segments):
                # copy the original segment non-destructively (shallow copy is not sufficient)
                new_segment = deepcopy(this_segment)
                # Append the appendix waveform(s) and update the pulse count summaries
                new_segment.waveform.extend(this_appendix_segment.waveform)
                
                # Special case for the sequence start marker. We only want the pulse to be copied for
                # segment zero. Otherwise force the first marker to blank (if it exists at all)

                if (segment_id > 0):
                    if new_segment.markers != []:
                        new_segment.markers[0] = []
                        
                # add this to the new segment list we're building
                new_segment_list.append(new_segment)
                new_counts.append(this_count)
                new_marker_masks.append(this_marker_mask)
                new_adv_modes.append(this_adv_mode)
                
                new_wfm_counter     += 1 # Not sure what these are supposed to count so it's difficult to make it right.
                new_markers_counter += 1 # Also there seems to be a bug in the original
                
        # now attach the new extended lists back onto the original sequence channel
        this_sequence.segments = new_segment_list
        this_sequence.counts = new_counts
        this_sequence.adv_modes = new_adv_modes
        this_sequence.marker_masks = new_marker_masks
        this_sequence._markers_counter = new_markers_counter
        this_sequence._wfm_counter = new_wfm_counter
