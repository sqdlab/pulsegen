########################################################
####												####
####		Pattern Configuration File				####
####												####
####	All parameters are given in units of		####
####	seconds, radians or no units.				####
####												####
####	If only one value is set to a parameter,	####
####	this value is distributed to all channel	####
####	pairs. If the value is a list, entry i is	####
####	assigned to channel pair i.					####
####												####
########################################################
from math import pi
from pulsegen.pulse import SquarePulse, GaussianPulse, DRAGGaussPulse
#from .pulse import SquarePulse, GaussianPulse, DRAGGaussPulse

NUMBER_OF_CHANNEL_PAIRS = 2

config = {
#########################
### PULSE CALIBRATION ###
#########################

#"pulse_shape": [SquarePulse],#, SquarePulse],
"pulse_shape": [GaussianPulse, GaussianPulse],
"pi_amp": [0.9, 0.9],
"pi_half_amp": [0.45, 0.45],
"mixer_calibration": [[(1, 0.5*pi), (1, 0*pi)], 
					  [(1, 0.5*pi), (1, 0*pi)]],
#"if_freq": 0.,
"if_freq": 2*pi*100e6,
				
##########################
###   PULSE CONFIG   #####
##########################

"sigma": 20.0e-9,
"truncate": 2,
"qscale": 0.5,
"separation": 5.0e-9,
"length": 20.0e-9,

##########################
### QUBIT PARAMETERS  ####
##########################
"anharmonicity": -2*pi*254e6,
				
##########################
### PATTERN PARAMETERS ###
##########################

"pattern_length": 5.0e-6,
"fixed_point": 4.0e-6,

#########################
### OPTIMAL CONTROL   ###
#########################

"use_optimal_control": False,

"impulse_response": [None, None
#					['C:\\Software\\QTLab\\custom_lib\\pulsegen\\test\\pinv_map_matrix_1us_sample.npz',
#					'C:\\Software\\QTLab\\custom_lib\\pulsegen\\test\\pinv_map_matrix_1us_sample.npz']
					],
"sampling_rate_multiplier": 10,
"ir_tlim": (-1e-8, 2e-7),
"prepend": 50e-9,
"append": 100e-9,
"use_pinv": True,
# Currently the filter only applies for pulses with optimal control
"use_fir": [False, False],
# window, cutoff frequency is given in units of the nyquist sampling rate
"fir_parameters": ("gaussian", 2*pi*300e6),
					 
#########################
### CHANNEL SETTINGS  ###
#########################

"channel_type": ["mixer", "mixer"],
#"awg_models": "Tektronix_5014",
"awg_models": "Agilent_N824x",
#"awg_models": "Agilent_81180A",
"lowhigh": (-1, 1),
"delay": [[0, 0], [0, 0]],
"marker_delays": [[(0, -82e-9 + 5e-9),(0,0)], [(0,0),(0,0)]],
"sampling_rate": 1.2e9
}