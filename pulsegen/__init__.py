import logging
import importlib
import six

def reimport(module, package=__package__):
    '''import a module, reloading it if was already imported'''
    module = package + '.' + module
    if module in globals():
        logging.debug(__name__ + ': forcing reload of {0}'.format(module))
        six.moves.reload_module(globals()[module])
    else:
        globals()[module] = importlib.import_module(module, package)

#import ptplot_gui

reimport('pulse')
from .pulse import Param

from .pulse import (SquarePulse, ArbitraryPulse, Pulse, DRAGGaussPulse,
                    ModulatedSinePulse, ZeroPulse, Spacer, GaussianPulse)

from .pulse import (mwpulse, identity, pix, piy, piminusx, piminusy, pihalfx, 
                    pihalfy, pihalfminusx, pihalfminusy, mwspacer, mwdummy, 
                    marker, spacer, dummy, rot, pattern_start_marker, 
                    meas_marker)

reimport('config')

reimport('optimalcontrol')

reimport('exporter')

reimport('sequence')

from .sequence import MultiAWGSequence
