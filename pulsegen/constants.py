'''
Created on 16/04/2013

@author: Matthias Baur
'''
import os

dirname, filename = os.path.split(os.path.abspath(__file__))

HOME_DIR = dirname
CONFIG_FILE = os.path.join(HOME_DIR, 'pulse_config.py')
#PATTERN_DIR = os.path.join(HOME_DIR, "patterns")
#IR_DIR = os.path.join(HOME_DIR, "impulse_response")
