#!/usr/bin/env python3
# 24, BereshpolovaLab, University of Connecticut


# #########################################################
# CLogFileObject
# #########################################################

class CLogFileObject:

    def __init__(self, LogFileName):
        self.LogFile = None
        if LogFileName:
            try:
                self.LogFile = open(LogFileName, 'w')
            except IOError:
                print('Cannot open log-file! Something is wrong.')

    def __del__(self):
        if self.LogFile and not self.LogFile.closed:
            self.LogFile.close()
