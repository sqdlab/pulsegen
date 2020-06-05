'''
Created on 13/05/2013

@author: uqmbaue1
'''
from .constants import CONFIG_FILE
#import pulse_config
import imp
import copy


class _Config(object):
    '''
    classdocs
    '''
    def __init__(self, filename):
        '''
        Constructor
        '''
        self._filename = filename
        self._config = {}
        self._channel_config = []
        self._channel_pair_config = []
        
    def load(self):
        '''Load settings from file'''
        
        with open(self._filename, 'r') as cfg_file:
            pulse_config = imp.load_source('pulse_config', self._filename, 
                                           cfg_file)
            config = pulse_config.config
            chpairs = pulse_config.NUMBER_OF_CHANNEL_PAIRS
            #del pulse_config
        self.parse(config, chpairs, update=False)
    
    def parse(self, config, chpairs=None, update=True):
        '''
        Load settings from a dictionary
        
        Input:
            config (dict) - configuration dictionary
            chpairs (int) - number of channel pairs, optional when 
                update is True
            update (bool) - if True, update the current configuration,
                if False, replace it.
        '''
        if update:
            if ((chpairs is not None) and 
                (hasattr(self, '_number_of_channel_pairs')) and 
                (self._number_of_channel_pairs != chpairs)):
                raise ValueError('number of channel pairs must be unchanged'+
                                 'when updating the configuration.')
            self._config.update(config)
        else:
            if chpairs is None:
                raise ValueError('number of channel pairs may only be omitted '+
                                 'when update is True.')
            self._number_of_channel_pairs = chpairs
            self._config = config
        self._channel_config = [
            {} for __ in range(0, 2 * self._number_of_channel_pairs)]
        self._channel_pair_config = [
            {} for __ in range(0, self._number_of_channel_pairs)]
        for key, value in self._config.items():
            for channel_pair in range(0, self._number_of_channel_pairs):
                if not isinstance(value, list):
                    self._channel_pair_config[channel_pair][key] = value
                    self._channel_config[2*channel_pair][key] = value
                    self._channel_config[2*channel_pair + 1][key] = value
                elif len(value) == self._number_of_channel_pairs:
                    channel_pair_value = value[channel_pair]
                    self._channel_pair_config[channel_pair][key] = \
                        channel_pair_value
                     
                    if not isinstance(channel_pair_value, list):
                        self._channel_config[2*channel_pair][key] = \
                            channel_pair_value
                        self._channel_config[2*channel_pair + 1][key] = \
                            channel_pair_value
                    elif (isinstance(channel_pair_value, list)
                          and len(channel_pair_value) == 2):
                        self._channel_config[2*channel_pair][key] = \
                            channel_pair_value[0]
                        self._channel_config[2*channel_pair + 1][key] = \
                            channel_pair_value[1]
                    else:
                        raise ValueError(
                            'for parameter {0}, wrong number of values {1}'
                            ' are assigned'.format(key, channel_pair_value))
                else:
                    raise ValueError(
                        'for parameter {0}, wrong number of values {1}'
                        ' are assigned'.format(key, value))
        
    def get_channel_config(self):
        return self._channel_config
    
    def get_channel_pair_config(self):
        return self._channel_pair_config
    
    def get_config(self):
        return self._config

cfg = _Config(CONFIG_FILE)
cfg.load()