import logging

worstSeverity = logging.NOTSET

class TrackLogSeverity(logging.Handler) :
    '''Simplistic logging handler that simply stores the worst kind of log message seen'''

    def __init__(self) :
        logging.Handler.__init__(self)

    def emit(self, record) :
        global worstSeverity
        if record.levelno > worstSeverity :
            worstSeverity = record.levelno
