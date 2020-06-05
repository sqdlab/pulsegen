'''
Created on 16/04/2013

@author: Matthias Baur
'''

import numpy as np
import abc
import copy
from inspect import getargspec
import collections
import functools
import six

# TODO: Make stop_time a necessary attribute of a pulse? This requires 
# modification in the sequence module though.


class Param(object):
    '''
    This class is used to specify a configuration file parameter. It emulates 
    a numeric type such that arithmetic operations can be carried out on
    the string specifying a parameter, which is later replaced from the
    by the value of a dictionary entry (config file) by calling the instance
    with the dictionary as a parameter.
    
    '''
    def __init__(self, key):
        '''
        Create a Param instance.
        
        Parameters:
        -----------
        key: String
            Key of the configuration parameter.
        
        '''
        self._function = lambda cfg: cfg[key]
    
    def _copy(self):
        return copy.deepcopy(self)
    
#     def _check_types(self, function):
#         ''' 
#         decorate a function to check if the inputs are of a specific type
#         '''
#         types = (int, float, Param)
#         @functools.wraps(function)
#         def decorated_function(*args, **kwargs):
#             for (a, t) in zip(args, types):
#                 if not isinstance(a, t):
#                     raise TypeError('unsupported operand type(s) for {0}:'
#                                      ' `Param` and `{1}`'
#                                      .format(function.__name__, type(a)))
#             return function(*args, **kwargs)
#         return decorated_function
       
    def _make_function(self, other):
        if hasattr(other, '__class__') and other.__class__.__name__.endswith('Param'):
            return other
        else:
            return lambda cfg: other

    def __add__(self, other):
        other = self._make_function(other)
        new = self._copy()
        new._function = lambda cfg: self._function(cfg) + other(cfg)
        return new
    
    __radd__ = __add__
    
    def __sub__(self, other):
        other = self._make_function(other)
        new = self._copy()
        new._function = lambda cfg: self._function(cfg) - other(cfg)
        return new
    
    def __rsub__(self, other):
        other = self._make_function(other)
        new = self._copy()
        new._function = lambda cfg: other(cfg) - self._function(cfg)
        return new
    
    def __mul__(self, other):
        other = self._make_function(other)
        new = self._copy()
        new._function = lambda cfg: self._function(cfg) * other(cfg)
        return new
    
    __rmul__ = __mul__
    
    def __pow__(self, other):
        other = self._make_function(other)
        new = self._copy()
        new._function = lambda cfg: self._function(cfg) ** other(cfg)
        return new
    
    def __rpow__(self, other):
        other = self._make_function(other)
        new = self._copy()
        new._function = lambda cfg: other(cfg) ** self._function(cfg)
        return new
    
    def __div__(self, other):
        other = self._make_function(other)
        new = self._copy()
        new._function = lambda cfg: self._function(cfg) / other(cfg)
        return new
    
    def __rdiv__(self, other):
        other = self._make_function(other)
        new = self._copy()
        new._function = lambda cfg: other(cfg) / self._function(cfg)
        return new
    
    def __call__(self, cfg):
        '''
        Evaluate the arithmetic operations carried out on Param with the
        parameter value cfg[key], where cfg is a dictionary.
        
        '''
        return self._function(cfg)


def gauss(t, sigma, mu):
    '''
    Calculate the Gaussian, element-wise
    
    Parameters:
    -----------
        t: array_like or float
            Time.
        sigma: float
            Standard deviation.
        mu: float
            Expectation value.
            
    Returns:
    --------
        out : ndarray or float
            Output array, element-wise Gaussian of `t`.
    
    '''
    return np.exp(-(np.power((t - mu), 2))/(2*sigma*sigma))

def gausstruncated(t, sigma, truncate, mu):
    '''
    Calculate the truncated Gaussian that is forced 
    to start and end at zero, element-wise
    
    Parameters:
    -----------
        t: array_like of float
            Time.
        sigma: float
            Standard deviation.
        truncate:
            The Gaussian is truncated at times +- truncate * sigma.
        mu: float
            Expectation value.
            
    Returns:
    --------
        out: ndarray or float
            Output array, element-wise truncated Gaussian at `t`, 
            according to the function:
    
        function(t, tstop): 
            (gauss(t, sigma, tstop - len/2) 
            - gauss(tstop, sigma, tstop - len/2)
            ) / (1 - gauss(tstop, sigma, tstop - len/2))
    
    '''
    return (gauss(t, sigma, mu) - gauss(mu + sigma*truncate, sigma, mu)
            )/(1 - gauss(mu + sigma*truncate, sigma, mu))

def gausstruncated_derivative(t, sigma, truncate, mu):
    '''
    Calculate the first time derivative of the truncated Gaussian.
    
    Parameters:
    -----------
        t: array_like of float
            Time.
        sigma: float
            Standard deviation.
        truncate:
            The Gaussian is truncated at times +- truncate * sigma.
        mu: float
            Expectation value.
            
    Returns:
    --------
        out: ndarray or float
            Output array, element-wise first time derivative of the 
            truncated Gaussian at `t`.
    
    '''
    return (
            -(gauss(t, sigma, mu)*(t - mu))
            / ((1 - gauss(mu + sigma*truncate, sigma, mu))* np.power(sigma, 2))
            )

def sample_pulse(pulse, end_time, sampling_rate):
    '''
    Sample a single pulse and return its waveform.
    
    Paramters
    ---------
    pulse : Pulse
        The pulse object to be sampled and added to the waveformn.
    end_time : float
        End time of the pulse in seconds. The start time is automatically
        determined by the length of the pulse.
    sampling_rate : float
        Sampling rate in Hz.
    
    Returns 
    --------
    waveformn : ndarray
        waveformn including the sampled pulse
    start_time : int
        Starting time of the added pulse
        
    '''
    start_time = end_time - pulse.len
    start_index = int(round(start_time * sampling_rate))
        # +1 because np.arange does not include the stop index.
        # In [162]: np.arange(0, 4)
        # Out[162]: array([0, 1, 2, 3])
        # The same is true for indexing an error, the stop index is neglected.
        # In [163]: [0, 1, 2, 3][0:4]
        # Out[163]: [0, 1, 2, 3]
        # TODO: Do we want a +1 or not here? I guess not, otherwise we will 
        # have overlapping pulses if we concatenate them.
    end_index = int(round(end_time * sampling_rate))
    
    if pulse.len < 0:
        start_index, end_index = end_index, start_index
   
    if (0 < abs(pulse.len) < (0.5/sampling_rate)):
        sampled_pulse = np.array([])
    
    elif pulse.amp == 0:
        sampled_pulse = np.zeros(end_index - start_index)
    
    else:
        times = np.arange(start_index, end_index, 1) / sampling_rate
        sampled_pulse = pulse.function(times, end_time)
        
    start_time = start_time - pulse.sep
    return(sampled_pulse, start_time, start_index, end_index)

@six.add_metaclass(abc.ABCMeta)
class PulseBase(object):
    '''
    Abstract pulse base class.
    
    Attributes:
    -----------
        amp: float [-1, 1]
            Maximal amplitude of the pulse.
        len: float
            Length of the pulse.
        function: method
            A method f defining the shape of the pulse.
            Call the method with f(t, tstop).
        sep: float
            Time in seconds that separates this pulse 
            from a pulse earlier in time.
        duration: float
            Total time duration of the pulse given by
            len + sep.
    
    Optional:
        optimal_control: bool
            If True, then optimal control routine is applied 
            to this pulse prior to sampling, which deconvolves 
            the impulse response of the arbitrary waveform generator 
            and cables. Otherwise, sample the pulse as is (default).
    
    '''
    
    def __init__(self, amplitude, length, separation, optimal_control):
        if not (-1 <= amplitude) & (amplitude <= 1) :
            raise ValueError("amplitude {0} should be a number "
                             "between -1 and 1".format(amplitude))
        self._amp = amplitude
        self._length = length
        self._separation = separation
        self._optimal_control = optimal_control
        self._duration = length + separation
        self._EPSILON = 1.0e-14
        self.tstop = None
        
    @property
    def amp(self):
        '''
        Get or set the amplitude of the pulse.
        
        Type: float [-1, 1]
        
        '''
        return self._amp
    
    @amp.setter
    def amp(self, amp):
        if not (-1 <= amp) & (amp <= 1) :
            raise ValueError("amplitude {0} should be a number "
                             "between -1 and 1".format(amp))
        self._amp = amp
        
    @property
    def sep(self):
        '''
        Get or set the pulse separation of the pulse, which is defined
        as the time that separates this pulse from a pulse earlier in time.
        
        Type: float
         
        '''
        return self._separation
    
    @sep.setter
    def sep(self, separation):
        self._separation = separation
        self._duration = self._length + separation
    
    @property
    def optimal_control(self):
        '''
        Get or set whether to send the pulse through
        the optimal pulse control routine.
        
        Type: Boolean
        
        '''
        return self._optimal_control
    
    @optimal_control.setter
    def optimal_control(self, optimal_control):
        self._optimal_control = optimal_control

    @property
    def function(self):
        '''Get the function defining the pulse.'''
        # This makes it more consistent where the pulse starts
        # and where it stops, by eliminating rounding errors < self.EPSILON.
        # Pulse starting at T = 0 has index 0 different from 0.
        def f(t, tstop):
            t = t + self._EPSILON
            return self._function(t, tstop)
        
        if self.tstop is None:
            return f
        else:
            return lambda t: f(t, self.tstop)
    
    @property
    def len(self):
        '''
        Get or set the length of the pulse.
        
        Type: float
        
        '''
        return self._length
    
#     @len.setter
#     def len(self, length):
#         self._length = length
#         self._duration = length + self._separation
        
    @property
    def duration(self):
        '''Get the total duration of the pulse.'''
        return self._duration
    
    def whoami(self):
        return self.__class__.__name__
    
    def __call__(self, t, tstop):
        return self._function(t, tstop)
    
    @abc.abstractmethod
    def __str__(self):
        '''Define the string representation of the object.'''
        return
    
    @abc.abstractmethod
    def _function(self, t, stop):
        '''Define the function of the pulse.'''
        return
    
    # This function is not defined as __hash__, as __hash__ should only
    # be defined by immutable objects.
    @abc.abstractmethod
    def cache_id(self):
        '''
        Define the cache id of a pulse used to identify a cached pulse.
        
        This is used to cache sampled pulses to make the pattern generation
        faster, because optimal pulse control is generally pretty slow. If the
        pulse is the same up to the amplitude, we can reuse previously sampled
        pulses.
        
        Always include the output of whoami() function into the function which
        returns the id to get no dublicates.
         
        Return:
        -------
        None if not cacheable, int otherwise
        '''
        return


class SquarePulse(PulseBase):
    '''
    Create a square pulse.
    
    Attributes:
    -----------
        Same as PulseBase class
    
    '''
    def __init__(self, amplitude, length, separation, optimal_control=False):
        PulseBase.__init__(self, amplitude, length, separation, optimal_control)
        
    def _function(self, t, tstop):
        if self._length > 0:
            return self._amp * np.piecewise(
                t, [((tstop - self._length) < t) & (t < tstop)], [1, 0])
        if self._length < 0:
            return self._amp * np.piecewise(
                t, [((tstop - self._length) > t) & (t > tstop)], [1, 0])
        else:
            return 0
        
    def __str__(self):
        return "Square({0}, {1}, {2})".format(
            self._amp, self._length, self._separation)
    
    def cache_id(self):
        return hash((self.whoami(), self.len, self.sep, self.optimal_control))


class GaussianPulse(PulseBase):
    '''
    Create a truncated Gaussian pulse that is forced 
    to start and end at zero.
    
    function(t, tstop): 
        (gauss(t, sigma, tstop - len/2) 
        - gauss(tstop, sigma, tstop - len/2)
        ) / (1 - gauss(tstop, sigma, tstop - len/2))
    
    Attributes:
    -----------
        Same as PulseBase class, and following additional
        attributes.
        
        sigma: float
            The standard deviation of the Gaussian pulse.
        truncate: float [0, inf)
            The Gaussian is truncated at times +- truncate * sigma.

    '''
    def __init__(self, amplitude, sigma, truncate, separation,
                 optimal_control=False):
        self._sigma = sigma
        if truncate < 0:
            raise ValueError("truncate should be a positive number")
        self._trun = truncate
        PulseBase.__init__(
            self, amplitude, 2 * sigma * truncate, separation, optimal_control)
        
    def _function(self, t, tstop):
        f = lambda t: gausstruncated(
                t, self._sigma, self._trun, tstop - self._length/2)
        return self._amp * np.piecewise(
            t, [((tstop - self._length) < t) & (t < tstop)], [f, 0])
               
    def __str__(self):
        return "Gauss({0:.2f}, {1:.1f} ns, {2:.1f}, {3:.1f} ns)".format(
            self._amp, self._sigma*1.0e9, self._trun, self._separation * 1.0e9)
    
    def cache_id(self):
        return None


class ModulatedSinePulse(PulseBase):
    '''
    Create an amplitude, phase, and frequency modulated sine pulse.
    
    Attributes:
    -----------
        Same as PulseBase class, and following additional
        attributes.
        
        pulse: some pulse
            A pulse defining the envelope of the modulated sine pulse.
        omega: float or function reference
            The frequency of the sine in rad/s. Can be a float or a 
            function f(t,tstop) which defines the frequency modulation, 
            and takes the time t and the stop time tstop as input parameters.
        phase: float or function reference
            The phase of the sine in rad. Can be a float or a function 
            f(t,tstop) which defines the phase modulation, and takes 
            the time t and the stop time tstop as input parameters.

    '''
    
    def __init__(self, pulse, omega, phase, optimal_control=False):
        if not isinstance(pulse, PulseBase):
            raise TypeError("pulse must be an instance of a PulseBase subclass")
        self._pulse = pulse
        self._omega = omega
        self._phase = phase
#         self._fromfunction = False
        # amplitude set to 1 instead of self._pulse.amp -- amp is already part of the envelope function
        PulseBase.__init__(self, 1, self._pulse.len,
                           self._pulse.sep, optimal_control)
        
    def _function(self, t, tstop):        
        if type(self._omega).__name__ == 'function':
            omega = self._omega(t, tstop)
        else:
            omega = self._omega
        
        if type(self._phase).__name__ == 'function':
            phase = self._phase(t, tstop)
        else:
            phase = self._phase
        
        envelope = self._pulse.function(t, tstop)
        return self._amp * envelope * np.sin(omega*t + phase)
    
    def __str__(self):
        if (type(self._phase).__name__ == 'function'
            and type(self._omega).__name__ != 'function'):
            return "ModulatedSinePulse({0}, {1:.1f} MHz*2*pi, {2})".format(
                self._pulse, self._omega/(2*np.pi)*1.0e-6, self._phase)
        
        if (type(self._phase).__name__ != 'function'
            and type(self._omega).__name__ == 'function'):
            return "ModulatedSinePulse({0}, {0}, %.2f*pi)".format(
                self._pulse, self._omega, self._phase)
        
        if (type(self._phase).__name__ == 'function'
            and type(self._omega).__name__ == 'function'):
            return "ModulatedSinePulse({0}, {1}, {2})".format(
                self._pulse, self._omega, self._phase)
        
        else:
            return ("ModulatedSinePulse({0}, {1:.1f} MHz*2*pi, {2:.2f}*pi)"
                    .format(self._pulse, self._omega/(2*np.pi)*1.0e-6,
                            self._phase/(np.pi))
                    )

    def cache_id(self):
        return None

class ZeroPulse(PulseBase):
    '''
    Create a pulse with zero amplitude. Compared to the `Spacer`, the 
    `ZeroPulse` is intended to support easier line up with normal pulses
    as it includes a pulse separation attribute.
    
    Attributes:
    -----------
        amp: constant at zero
        len: float
            Length of the spacer pulse
        sep: float
            Time in seconds that separates this pulse 
            from a pulse earlier in time.
        optimal_control: constant False
    
    '''
    def __init__(self, length, separation):
        PulseBase.__init__(self, 0, length, separation, optimal_control=False)
        
    @PulseBase.amp.setter
    def amp(self, amp):
        raise AttributeError("`ZeroPulse` object does not support"
                             " the assignment of the `amp` attribute.")
    
    @PulseBase.optimal_control.setter
    def optimal_control(self, optimal_control):
        raise AttributeError("`ZeroPulse` object does not support"
                             " the assignment of the `optimal_control`"
                             " attribute.")
        
    def _function(self, t, tstop):
        return 0
        
    def __str__(self):
        return "Zero({0}, {1})".format(self._length, self._separation)

    def cache_id(self):
        return hash((self.whoami(), self._length, self._separation))

class Spacer(PulseBase):
    '''
    Create a pulse with zero amplitude and zero pulse separation
    
    Attributes:
    -----------
        amplitude: constant at zero
        length: float
            Length of the spacer pulse
        separation: constant at zero
        optimal_control: constant False
    
    '''
    def __init__(self, length):
        PulseBase.__init__(self, 0, length, 0, optimal_control=False)
        
    @PulseBase.amp.setter
    def amp(self, amp):
        raise AttributeError("`Spacer` object does not support the assignment"
                             " of the `amp` attribute.")
    
    @PulseBase.sep.setter
    def sep(self, separation):
        raise AttributeError("`Spacer` object does not support the assignment"
                             " of the `separation` attribute.")
    
    @PulseBase.optimal_control.setter
    def optimal_control(self, optimal_control):
        raise AttributeError("`Spacer` object does not support the assignment"
                             " of the `optimal_control` attribute.")
        
    def _function(self, t, tstop):
        return 0
        
    def __str__(self):
        return "Spacer({0})".format(self._length)
    
    def cache_id(self):
        return hash((self.whoami(), self._length))

class ArbitraryPulse(PulseBase):
    '''
    Create a pulse with an arbitrary function.
    
    Attributes:
    -----------
        Same as `PulseBase` class and the following additional attributes.
    
        function: function reference
            An arbitrary function f(t,tstop), which defines the envelope of 
            the pulse and takes the time t and stop time tstop as input
            parameters. This can be a lambda function or a normal function.
        
    Example:
    --------
    This code creates a square pulse:
    
    >>> def f(length): 
            return lambda t, tstop: np.piecewise(
                t, [(tstop - length < t) & (t < tstop)], [1, 0])
    >>> length = 5.0e-9
    >>> myPulse = ArbitraryPulse(0.5, length, 1.0e-8, f(length))
    >>> times = np.arange(4994, 5001, 1)*1.0e-9
    >>> myPulse.function(times, 5.0e-6)
    array([ 0. ,  0.5,  0.5,  0.5,  0.5,  0.5,  0. ])
    
    '''
    def __init__(self, amplitude, length, separation, function,
                 optimal_control=False):
        self._arbitrary_function = function
        PulseBase.__init__(self, amplitude, length, separation, optimal_control)

    def _function(self, t, tstop):
        return self._amp * self._arbitrary_function(t, tstop)
    
    @PulseBase.function.setter
    def function(self, function):
        self._function = function
        
    def __str__(self):
        return "ArbitraryPulse({0}, {1}, {2}, {3})".format(
            self._amp, self._length, self._separation, self._arbitrary_function)
    
    def cache_id(self):
        return hash((self.whoami(), self._length, self._separation, 
                     self._arbitrary_function, self._optimal_control))

class DRAGGaussPulse(PulseBase):
    '''
    Create a DRAGGaussPulse, referred to as the Y-only correction in 
    Phys. Rev. A 83, 012308 (2011)
    
    Attributes:
    -----------
        Same as `GaussianPulse` class, and following additional attributes.
        
        anharmonicity: float
            The anharmonicity of the qubit, defined as \omega_{21} - \omega_{10}
            in rad/s units.
    
        qscale: float
            Scaling factor of the q quadrature.
        
        ret: string
            Specify which functions to return. Possible values:
                AmpPhase (default), IQ
            
    '''
    def __init__(self, amplitude, sigma, truncate, anharmonicity,
                 qscale, separation, optimal_control=False, ret="AmpPhase"):
        self._sigma = sigma
        if truncate < 0:
            raise ValueError("truncate should be a positive number")
        self._trun = truncate
        self._anharmonicity = anharmonicity
        self._qscale = qscale
        self._ret = ret
        PulseBase.__init__(
            self, amplitude, 2 * sigma * truncate, separation, optimal_control)
    
    def _function(self, t, tstop):
        f_i = lambda t: gausstruncated(
            t, self._sigma, self._trun, tstop - self._length/2)
        
        f_q = lambda t: gausstruncated_derivative(
            t, self._sigma, self._trun, tstop - self._length/2
            ) / self._anharmonicity
        
        i_quadrature = np.piecewise(
            t, [((tstop - self._length) < t) & (t < tstop)], [f_i, 0])
        
        q_quadrature = - self._qscale * np.piecewise(
            t, [((tstop - self._length) < t) & (t < tstop)], [f_q, 0])
        
        if self._ret == "IQ":
            return self._amp * np.array([i_quadrature, q_quadrature])
        
        elif self._ret == "AmpPhase":
            signal = i_quadrature + 1j*q_quadrature
            return self._amp * np.array([np.abs(signal), np.angle(signal)])
        
        else:
            raise ValueError("Unknown `ret` value {0:s}".format(self._ret))
    
    @property
    def function(self):
        return lambda *args, **kwargs: PulseBase.function.fget(self)(*args, **kwargs)[0]
        
    @property
    def phase(self):
        return lambda *args, **kwargs: PulseBase.function.fget(self)(*args, **kwargs)[1]
    
    @property
    def qscale(self):
        return self._qscale
    
    def __str__(self):
        return ("DRAGGaussPulse({0:.3f}, {1:.1f} ns, {2:.2f}, {3:.1f} MHz*2*pi,"
                "{4:.2f}, {5:.1f} ns)"
                .format(self._amp,
                        self._sigma*1.0e9,
                        self._trun,
                        self._anharmonicity/(2*np.pi)*1.0e-6,
                        self._qscale,
                        self._separation*1.0e9)
                )
    
    def cache_id(self):
        return None


class Pulse(object):
    '''
    Pulse factory
    
    This class is used to create pulse instances to be used together with the
    pulse configuration and sequence module. It creates pulse objects with 
    pulse attributes given by the user as keyword arguments, or by the 
    pulse configuration file if the keyword argument for an attribute is
    missing.
    
    It is also used to create several dependent pulse objects, for example
    the two pulse objects for the I and Q channel of an upconversion mixer when
    used in the single sideband upconversion mode.
    
    '''
    def __init__(self, envelope=Param("pulse_shape"), ptype="SSB", **kwargs):
        '''
        Create a Pulse object
        
        Keyword arguments:
        -----------------
        envelope: a pulse object
            Specify the envelope of the pulse. Default value is the value of
            the configuration file entry `pulse_shape`.

        ptype: string (default "SSB")
            "single": A single pulse object given by `envelope` is created.
            
            "SSB": Two ModulatedSinePulse objects used for the I and Q channel
                   of a single sideband upconversion mixer are created.
                   The envelope of the pulse is given by `envelope`, and the
                   sine frequency given by the keyword `omega` (in rad/s) or 
                   the configuration file entry `if_freq` (in Hz).
                   
                   Additionaly to the attributes of ModulatedSinePulse, the
                   keyword argument `mixer_calibration` must be specified.
        
        kwargs:
            Keyword arguements defining the attributes of `envelope` and 
            `ptype` specific keywords.
            
            ptype "SSB":
                `mixer_calibration`: [(Iamp, Iphase), (Qamp, Qphase)]
                The amplitude of the I and Q quadratures are multiplied
                with Iamp and Qamp, respectively.
                Iphase and Qphase are added to the phase of the I and Q
                quadratures.
        
        '''
        self._envelope = envelope
        self._ptype = ptype
        self._kwargs = kwargs

    # TODO: Do we really need this copy here? Maybe only copy in get_pulse?
    def __call__(self, **kwargs):
        selfcopy = copy.deepcopy(self)
        selfcopy._kwargs.update(kwargs)
        return selfcopy
        
    def get_pulse(self, cfg):
        envelope = self._envelope
        if type(envelope).__name__.endswith('Param'):
            envelope = envelope(cfg)
        # All keyword arguments values set by the user that are Param instance
        # are replaced with the value from the config file
        kwargs = {}
        for key, value in self._kwargs.items():
            if type(value).__name__.endswith('Param'):
                kwargs[key] = value(cfg)
            else:
                kwargs[key] = value
        
        # Get all arguments of the constructor of the pulse class defined
        # by envelope
        pulse_argspec = getargspec(vars(envelope)["__init__"])
        pulse_arguments = pulse_argspec.args[1:]
        
        # Check whether any of the pulse arguments have been defined
        # by the user with **kwargs, and add them to pulse_kwargs
        pulse_kwords_user_set = set(kwargs).intersection(pulse_arguments)
        pulse_kwargs = dict(
            [(key, kwargs[key]) for key in pulse_kwords_user_set])
        
        # For all arguments not defined by the user, look up the value
        # in the cfg dictionary except the one that have default values
        arguments_user_set = set(kwargs)
        # remove the optional parameters with a default value
        if pulse_argspec.defaults is not None:
            pulse_arguments = pulse_arguments[:-len(pulse_argspec.defaults)]
            
        pulse_kwords_not_user_set = \
            set(pulse_arguments).difference(arguments_user_set)
        for key in pulse_kwords_not_user_set:
            pulse_kwargs[key] = cfg[key]
        
        # If the mixer is used in the single sideband mode,
        # we apply modulated sine waves to the I and Q port
        # of the upconversion mixer
        if self._ptype == "SSB":
            # Create the new envelope pulse
            i_envelope_pulse = envelope(**pulse_kwargs)
            q_envelope_pulse = envelope(**pulse_kwargs)
            
            if (envelope.__name__ == 'Spacer'
                or envelope.__name__ == 'ZeroPulse'):
                (i_pulse, q_pulse) = (i_envelope_pulse, q_envelope_pulse)
            
            else:
                pulse_kwargs.clear()
                # get the parameters for ModulatedSinePulse
                pulse_arguments = getargspec(
                    vars(ModulatedSinePulse)["__init__"]).args[1:]
                pulse_kwords_user_set = \
                    set(kwargs).intersection(pulse_arguments)
                pulse_kwargs = dict(
                    [(key, kwargs[key]) for key in pulse_kwords_user_set])
                
                if not "omega" in pulse_kwargs:
                    pulse_kwargs["omega"] = 2*np.pi*cfg["if_freq"]

                if "mixer_calibration" in kwargs:
                    mixer_calibration = kwargs["mixer_calibration"]
                else:
                    mixer_calibration = cfg["mixer_calibration"]
                    
                i_envelope_pulse.amp = (i_envelope_pulse.amp
                                        * mixer_calibration[0][0])
                q_envelope_pulse.amp = (q_envelope_pulse.amp
                                        * mixer_calibration[1][0])
                phase = pulse_kwargs.pop("phase")
                i_phase = mixer_calibration[0][1] + phase
                q_phase = mixer_calibration[1][1] + phase
                
                if envelope.__name__ == 'DRAGGaussPulse':
                    i_phase_copy = i_phase
                    q_phase_copy = q_phase
                    i_phase = lambda t, tstop: (
                        i_phase_copy + i_envelope_pulse.phase(t, tstop))
                    q_phase = lambda t, tstop: (
                        q_phase_copy + q_envelope_pulse.phase(t, tstop))

                i_pulse = ModulatedSinePulse(
                    i_envelope_pulse, phase=i_phase, **pulse_kwargs)
                q_pulse = ModulatedSinePulse(
                    q_envelope_pulse, phase=q_phase, **pulse_kwargs)
            
            return [i_pulse, q_pulse]

        elif self._ptype == "single":
            return envelope(**pulse_kwargs)

        else:
            raise ValueError('unknown pulse type %s' %self._ptype)
            return []
    
    def __repr__(self):
        return self.__str__()
        
    def __str__(self):
        return "Pulse({0}, {1}, {2})".format(
            self._envelope, self._ptype, self._kwargs)


mwpulse = Pulse(Param("pulse_shape"), ptype="SSB")
identity = Pulse(Param("pulse_shape"), ptype="SSB", amplitude=0., phase=0.)
pix = Pulse(Param("pulse_shape"), ptype="SSB", amplitude=Param("pi_amp"), phase=0.)
piminusx = Pulse(Param("pulse_shape"), ptype="SSB", amplitude=Param("pi_amp"), phase=np.pi)
piy = Pulse(Param("pulse_shape"), ptype="SSB", amplitude=Param("pi_amp"), phase=-np.pi/2)
piminusy = Pulse(Param("pulse_shape"), ptype="SSB", amplitude=Param("pi_amp"), phase=np.pi/2)
pihalfx = Pulse(Param("pulse_shape"), ptype="SSB", amplitude=Param("pi_half_amp"), phase=0.)
pihalfminusx = Pulse(Param("pulse_shape"), ptype="SSB", amplitude=Param("pi_half_amp"), phase=np.pi)
pihalfy = Pulse(Param("pulse_shape"), ptype="SSB", amplitude=Param("pi_half_amp"),
                phase=-np.pi/2)
pihalfminusy = Pulse(Param("pulse_shape"), ptype="SSB", amplitude=Param("pi_half_amp"),
                     phase=np.pi/2)

def mwspacer(length, **kwargs):
    return Pulse(Spacer, ptype="SSB")(length=length, **kwargs)
    
def mwdummy(length, **kwargs):
    return Pulse(ZeroPulse, ptype="SSB")(length=length, **kwargs)

def spacer(length, **kwargs):
    return Pulse(Spacer, ptype="single")(length=length, **kwargs)
    
def dummy(length, **kwargs):
    return Pulse(ZeroPulse, ptype="single")(length=length, **kwargs)

def rot(amp, phase, **kwargs):
    return Pulse(Param("pulse_shape"), ptype="SSB")(amplitude=amp, phase=phase,
                                             **kwargs)

def marker(length, **kwargs):
    return Pulse(SquarePulse, ptype="single")(amplitude=1, length=length,
                                              **kwargs)

def pattern_start_marker():
    return [marker(1.0e-7), spacer(Param("fixed_point") - 1.1e-7)]

def meas_marker():
    return [marker(Param("fixed_point") - 100e-9)]
