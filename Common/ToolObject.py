#!/usr/bin/env python3
# 2024, BereshpolovaLab, University of Connecticut
# Developed by Victor Serdyukov, e-mails: Victor.Serdyukov@uconn.edu


import math
import numpy as np


###########################################################
# CToolObject
###########################################################

class CToolObject:

    def ReverseCorrelation0(self, SE, SpikeList, TW, Count):
        SN, StartCount, T1, T2 = 0, Count, \
            SE.GetTimeStamp() + TW[0], SE.GetTimeStamp() + TW[1]
        for SP in SpikeList[Count:]:
            if SP < T1:
                StartCount += 1
            if SP >= T1 and SP <= T2:
                SN += 1
            if SP >= T2:
                break
        if StartCount >= len(SpikeList):
            StartCount = None
        return SN, StartCount

    def Distance(self, P1, P2):
        return math.sqrt((P2[0] - P1[0]) ** 2 + (P2[1] - P1[1]) ** 2) \
            if P1[0] is not None and P1[1] is not None and \
            P2[0] is not None and P2[1] is not None \
            else None

    def GetDistanceString(self, P1, P2, PR=3):
        DistanceString = None
        if P1 is not None and P2 is not None:
            D = self.Distance(P1, P2)
            if D is not None:
                DistanceString = self.FMT(D, F'{{:>6.{PR}f}}')
        return DistanceString

    def FMT(self, v, PatterString='{:0.2f}'):
        return (v if isinstance(v, str) else PatterString.format(v))\
            if v else 'None'

    def Angle(self, P1, P2=None, K=1.0):
        ThetaDegrees = None
        if P1 is not None:
            if P1[0] is not None and P1[1] is not None:
                if P2 is not None:
                    if P2[0] is not None and P2[1] is not None:
                        P = np.asarray(P2) - np.asarray(P1)
                        ThetaDegrees = self.Angle(P, None, K)
                else:
                    ThetaRadians = np.arctan2(K * P1[0], P1[1])
                    ThetaDegrees = np.mod(np.degrees(ThetaRadians), 360)
                    if ThetaDegrees < 0:
                        ThetaDegrees += 180
        return ThetaDegrees

    def GetAngleString(self, P1, P2, K=1.0, PR=0):
        AngleString = None
        if P1 is not None and P2 is not None:
            Angle = self.Angle(P1, P2, K)
            if Angle is not None:
                AngleString = self.FMT(Angle, F'{{:>3.{PR}f}}')
        return AngleString
