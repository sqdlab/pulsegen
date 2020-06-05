'''
Created on 02/05/2013

@author: Matthias Baur
'''

# Requires numpy 1.7 or more recent
import numpy as np
from scipy import interpolate
import collections
from copy import copy
import warnings
import six

def near_equal(a, b, eps):
    if abs(a-b) < eps:
        return True
    return False

class OptimalControl(object):
    '''
    Create an OptimalControl object used to find the optimal input `x`
    of a linear device, characterized by an impulse response,
    to receive a given pulse shape `y` at the output.
    
    Since the device is linear, this problem can be described by
    the linear equation y = A x, where A is a not necessarily rectangular
    matrix, determined by the impulse response.
    
    Since A is generally not rectangular, we can not calculate the inverse
    of A to find x for a given y. Instead, we determine the Moore-Penrose
    pseudo-inverse A^+, which solves the problem in a least-squares sense,
    i.e. |A^+ y - x|^2 is minimal.
    
    Attributes:
    -----------
    sampling_rate: float
        Sampling rate of the vector `x`
    sampling_rate_multiplier: float
        The scaling between the sampling rate of vector `x` and vector `y`.
        The wanted pulse shape `y` at the output of the linear device is 
        sampled with a higher sampling rate to receive more accurate fits.
    ir_tlim: tuple (tmin, tmax)
        The start and stop time of the impulse response used to build A
        in seconds.
    prepend: float
        The pulse is zero padded before the pulse for `prepend` seconds.
    append: float
        The pulse is zero padded after the pulse for `append` seconds.
    use_pinv: boolean
        If True, the function load expects the pseudo-inverse as a
        parameter.
        If False, the impulse response is expected.
    filter: Filter
        Object of the sequence.Filter class, to low pass filter the 
        pulse shape given by `y`, prior to fitting the pulse. 
            
    '''

    def __init__(self, sampling_rate, sampling_rate_multiplier, ir_tlim, 
                 prepend, append, use_pinv, filter):
        for param in ['sampling_rate', 'sampling_rate_multiplier', 'ir_tlim',
                      'prepend', 'append', 'use_pinv']:
            setattr(self, "_" + param, locals()[param])
        
        self._filter = copy(filter)
        
        # TODO: Remove 'use_pinv' and only determine it from the file given
        # for impulse_response
        self._Ts = 1 / self._sampling_rate
        self._Ta = self._Ts / self._sampling_rate_multiplier
        if self._filter is not None:
            self._filter.Fs = self._filter.Fs * self._sampling_rate_multiplier
        self._impulse_response = []
        self._map_matrix = []
        self._pinv_map_matrix = []
        
    def _load_impulse_response(self, impulse_response):
        # From file
        if isinstance(impulse_response, six.string_types):
            # Load the impulse response data from a space separated txt file.
            # First column are the time points, second column the amplitude.
            data = np.loadtxt(impulse_response, dtype=np.float)
            data = data.T
            self._impulse_response = self._resample_impulse_response(data)
            return True
        
        # From list
        elif isinstance(impulse_response, np.ndarray):
            if len(impulse_response) != 2:
                raise ValueError(
                    'Invalid list size in creation of OptimalControl. '
                    'The list should have exactly two elements')
            else:
                self._impulse_response = \
                    self._resample_impulse_response(impulse_response)
                return True
        
        else:
            return False
         
    def _resample_impulse_response(self, impulse_response):
        # Resample with the analog sampling rate,
        # rescale and shift the impulse response
        
        impulse_response[1] = impulse_response[1] / np.max(impulse_response[1])
        # Set the baseline offset of the impulse response to 0
        rightbaseline = np.mean(
            impulse_response[1, -10:len(impulse_response[1])])
        leftbaseline = np.mean(impulse_response[1, 1:10])
        if abs(rightbaseline - leftbaseline) > 1.0e-3:
            raise ValueError("The difference between the end point values" 
                             "of the impulse response is too large."
                             "Remeasure the impulse response with a larger"
                             "time window and/or more averaging.")
        impulse_response[1] = (impulse_response[1]
                               - (leftbaseline + rightbaseline)/2)
        
        # Check whether the impulse_response is long enough,
        # otherwise pad it with zeros
        Ts_ir = impulse_response[0,1] - impulse_response[0, 0]
        if ((self._ir_tlim[0] < impulse_response[0, 0]) or 
            (self._ir_tlim[1] > impulse_response[0, -1])):
            warnings.warn("The starting time of the impulse response"
                           "is larger than ir_tlim[0]. Zero padding the "
                           "impulse_response.")
            #raise UserWarning("The starting time of the impulse response" 
            #                  "is larger than ir_tlim[0]. Zero padding the"
            #                 "impulse_response.")
            x, y = impulse_response[0], impulse_response[1]

            while self._ir_tlim[0] < x[0]:
                x = np.insert(x,0,impulse_response[0, 0] - Ts_ir)
                y = np.insert(y,0,0)
                impulse_response = np.array([x, y])

            while self._ir_tlim[1] > x[-1]:
                x = np.append(x, impulse_response[0, -1] + Ts_ir)
                y = np.append(y,0)
                impulse_response = np.array([x, y])

#             while self._ir_tlim[0] < impulse_response[0, 0]:
#                 impulse_response[0].insert(0, impulse_response[0, 0] - Ts_ir)
#                 impulse_response[1].insert(0, 0)
#                 
#             while self._ir_tlim[1] > impulse_response[0, -1]:
#                 impulse_response[0].append(impulse_response[0, -1] + Ts_ir)
#                 impulse_response[1].append(0)
        
        # Select the impulse_response values which are in the time frame
        # self._ir_tlim
        impulse_response = impulse_response[..., np.all(
            (impulse_response[0] > self._ir_tlim[0], 
             impulse_response[0] < self._ir_tlim[1]
             ), axis = 0)]
        
        # Resample the impulse response with the sampling rate Ta
        # of the quasi analog desired pulse
        f = interpolate.UnivariateSpline(
            impulse_response[0], impulse_response[1], s=0, k=3)
        times = np.arange(
            impulse_response[0, 0], impulse_response[0, -1], self._Ta)
        ypoints = f(times)
        
        # Rescale the impulse response such that the output signal
        # is amplified by a factor 1.
        ypoints = ypoints/np.sum(ypoints) * self._sampling_rate_multiplier
        
        resampled_impulse_response = np.array([times, ypoints])
        return resampled_impulse_response
    
    def load(self, impulse_response):
        '''
        Load the impulse response or the pseudo-inverse map matrix.
        
        The attribute use_pinv determines whether the impulse response
        or the pseudo-inverse matrix is expected.
        
        Parameters:
        -----------
        impulse_response: 2d_array or string
            path and filename of the impulse_response
            or the pseudo-inverse map matrix. For the
            impulse_response, this can also be a 2d_array.
                
        Return:
        -------
            None
        
        '''
        if self._use_pinv:
            if not self._load_pinv_map_matrix(impulse_response):
                raise ValueError(
                    "pseudo inverse matrix {0} was generated with different "
                    "optimal control parameters than assigned to the current "
                    "instance of OptimalControl".format(impulse_response))
            else:
                return True
        else:
            if not self._load_impulse_response(impulse_response):
                raise Exception("can not load impulse response `{0}`"
                                .format(impulse_response))
            else:
                return True
    
    def set_map_matrix(self, fitting_window_length):
        rows = int(np.round(fitting_window_length/self._Ta 
                            + len(self._impulse_response[1]) 
                            - self._sampling_rate_multiplier
                            )
                   )
        columns = int(np.round(fitting_window_length/self._Ts))
        map_matrix = np.zeros((rows, columns), np.float)
        
        for i in range(columns):
            rowmin = i*int(self._sampling_rate_multiplier)
            rowmax = (i*int(self._sampling_rate_multiplier)
                      + len(self._impulse_response[1]))
            map_matrix[rowmin : rowmax, i] = self._impulse_response[1]
        
        self._map_matrix = map_matrix
        
    def get_prepend(self):
        return self._prepend
    
    def get_output(self, waveform):
        '''
        Calculate y = A waveform.
        
        '''
        return self._map_matrix.dot(waveform)
    
    def get_impulse_response(self):
        return self._impulse_response
    
    def set_pinv_map_matrix(self, fitting_window_length):
        self.set_map_matrix(fitting_window_length)
        self._pinv_map_matrix = np.linalg.pinv(self._map_matrix)
        
    def save_current_pinv_map_matrix(self, filename):
        tosave = {"ir_tlim": [self._ir_tlim],
                  "Ta": [self._Ta],
                  "Ts": [self._Ts],
                  "pinv_map_matrix": [self._pinv_map_matrix]
                  }
        np.savez(filename, **tosave)
        
    def save_pinv_map_matrix(self, filename, impulse_response,
                             fitting_window_length):
        '''
        Calculate and save the pseudo-inverse map matrix.
        
        Parameters:
        -----------
        filename: string
            Path and filename without the filename extension where the 
            pseudo-inverse map matrix should be saved.
        impulse_response: 2d-array or string
            path and filename of the impulse_response
            or the pseudo-inverse map matrix. For the
            impulse_response, this can also be a 2d_array.
        fitting_window_length: float
            The length of the vector `y` in seconds. This, together
            with the attributes prepend and append, determines the 
            maximal pulse length that can be fitted.
        
        '''
        # Save current data to restore afterwards
        current_impulse_response = self._impulse_response
        current_map_matrix = self._map_matrix
        current_pinv_map_matrix = self._pinv_map_matrix
        current_use_pinv = self._use_pinv
        
        self._use_pinv = False
        
        print("loading impulse response")
        self.load(impulse_response)
        print("calculating pseudo inverse map matrix")
        self.set_pinv_map_matrix(fitting_window_length)
        print("saving pseudo inverse map matrix")
        self.save_current_pinv_map_matrix(filename)
        print("done")
        
        # restore data
        self._impulse_response = current_impulse_response
        self._map_matrix = current_map_matrix
        self._pinv_map_matrix = current_pinv_map_matrix
        self._use_pinv = current_use_pinv
        
    def _load_pinv_map_matrix(self, filename):
        data = np.load(filename)
        keywords = ["ir_tlim", "Ta", "Ts"]
        is_equal = []
        eps = 1.0e-14
        # check whether the loaded pinv_map_matrix was generated
        # with the same optimal control options as currently set
        # in this instance of OptimalControl
        for kw in keywords:
            value = data[kw][0]
            self_value = getattr(self, "_" + kw)
            if (isinstance(value, collections.Iterable)
                and not isinstance(value, six.string_types)):
                is_equal.append(all(map(
                    lambda a, b: near_equal(a, b, eps), value, self_value)))
            else:
                is_equal.append(near_equal(value, self_value, eps))
        
        # test = np.hstack(
        #    [near_equal(data[kw][0], 
        #     getattr(self, "_" + kw), 1.0e-15) for kw in keywords]
        #    )
        if all(is_equal):
            self._pinv_map_matrix = data['pinv_map_matrix'][0]
            return True
        
        self._pinv_map_matrix = []
        return False
        
    def get_pinv_map_matrix(self):
        return self._pinv_map_matrix

    def __call__(self, mypulse):
        '''
        Calculate the optimal pulse `x`.
        
        Parameters:
        -----------
        mypulse: PulseBase
            The pulse shape `y` that we want to have at the output of the
            linear device.
                
        Returns:
        --------
        out: numpy.ndarray
            The waveform `x` that solves the equation A x = y
            in the least-squares sense.
                
        '''
        if (mypulse.len == 0) or (mypulse.amp == 0):
            return np.array([])
        
        if self._use_pinv:
            A = self._pinv_map_matrix
            times = np.arange(0, A.shape[1])*self._Ta
            
            if times[-1] < self._prepend + mypulse.len + self._append:
                raise ValueError(
                    'The window length of {0:.2e} of the pseudo '
                    'inverse map matrix is too small for the '
                    'pulse {1} with length {2:.2e} and the current '
                    'prepend and append settings'
                    .format(times[-1], mypulse, mypulse.len))
            # The - self._ir_tlim[0] is introduced in order to shift the pulse
            # by the offset which will be introduced by the impulse response.
            # If the maximum of the impulse response is at time 0,
            # this is roughly self._ir_tlim[0]
            y = mypulse.function(
                    times, self._prepend + mypulse.len - self._ir_tlim[0])
            if self._filter is not None:
                y = self._filter.filter(y)
            x = np.dot(A, y)
        else:
            fitting_window_length = mypulse.len + self._prepend + self._append
            self.set_map_matrix(fitting_window_length)
            A = self._map_matrix
            times = np.arange(0, A.shape[0])*self._Ta
            # The - self._ir_tlim[0] is introduced in order to shift the pulse
            # by the offset which will be introduced by the impulse response.
            # If the maximum of the impulse response is at time 0,
            # this is roughly self._ir_tlim[0]
            y = mypulse.function(
                    times, self._prepend + mypulse.len - self._ir_tlim[0])
            if self._filter is not None:
                y = self._filter.filter(y)
            x = np.linalg.lstsq(A, y)[0]
        
        return x
    