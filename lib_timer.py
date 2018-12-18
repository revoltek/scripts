import time, logging

class Timer(object):
    """
    context manager used to time the operations
    """

    def __init__(self, step = 'undefined', log = None):
        """
        log: is a logging istance to print the correct log format
        if nothing is passed, root is used
        """
        if log is None: self.log = logging
        else: self.log = log
        self.step = step

    def __enter__(self):
        self.log.debug("--> Starting \'" + self.step + "\'.")
        self.start = time.time()
        self.startcpu = time.clock()

    def __exit__(self, exit_type, value, tb):

        # if not an error
        if exit_type is None:
            self.log.debug("<-- Time for this step: %i s (cpu: %i s)." % ( ( time.time() - self.start), (time.clock() - self.startcpu) ))
