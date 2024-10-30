#!/usr/bin/env python3
# 24, BereshpolovaLab, University of Connecticut


###########################################################
# CPrintObject
###########################################################

class CPrintObject:

    def __init__(self, LogFile): self.LogFile = LogFile

#
# Get - functions
#

    def GetLogFile(self): return self.LogFile

#
# Print - functions
#

    def Print(self, List):
        String = ''.join(List)
        print(String)
        if self.LogFile and not self.LogFile.closed:
            self.LogFile.write(F'{String}\n')
