#!/usr/bin/env python3
# 2024, BereshpolovaLab, University of Connecticut
# Developed by Victor Serdyukov, e-mails: Victor.Serdyukov@uconn.edu


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
