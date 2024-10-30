#!/usr/bin/env python3
# 24, BereshpolovaLab, University of Connecticut


import struct


#
# CPL_Base
#

class CPL_Base:

    def __init__(self, Size): self.Size = Size

#
# Get - functions
#

    def GetSize(self): return self.Size


#
# CPL_FileHeader
#

class CPL_FileHeader(CPL_Base):

    def __init__(self, Data):
        super().__init__(7504)
        ADFrequency, self.NumDSPCH, self.NumEventCH, self.NumSlowCH = \
            struct.unpack('iiii', Data[136:152])
        self.TimeStampPeriod = 1000.0 / ADFrequency

#
# Get - functions
#

    def GetNumDSPCH(self): return self.NumDSPCH

    def GetNumEventCH(self): return self.NumEventCH

    def GetNumSlowCH(self): return self.NumSlowCH

    def GetTimeStampPeriod(self): return self.TimeStampPeriod


#
# CPL_ChanHeader
#

class CPL_ChanHeader(CPL_Base):

    def __init__(self): super().__init__(1020)


#
# CPL_EventHeader
#

class CPL_EventHeader(CPL_Base):

    def __init__(self): super().__init__(296)


#
# CPL_SlowChannelHeader
#

class CPL_SlowChannelHeader(CPL_Base):

    def __init__(self): super().__init__(296)


#
# CPL_DataBlockHeader
#

class CPL_DataBlockHeader(CPL_Base):

    def __init__(self):
        super().__init__(16)
        self.Type, self.TimeStamp, self.Channel, self.Unit = \
            None, None, None, None

#
# Get - functions
#
    def GetType(self): return self.Type

    def GetTimeStamp(self): return self.TimeStamp

    def GetChannel(self): return self.Channel

    def GetUnit(self): return self.Unit

#
# Service - functions
#
    def ExtractData(self, Data, N1):
        N2, Length = N1 + self.Size, len(Data)
        if N2 > Length:
            N1 = 0
        else:
            self.Type, _, self.TimeStamp, self.Channel, self.Unit, WFN, \
                WWFN = struct.unpack('hHIhhhh', Data[N1:N2])
            N2 += 2 * WFN * WWFN
            N1 = 0 if N2 > Length else N2
        return N1
