#!/usr/bin/env python3
# 24, BereshpolovaLab, University of Connecticut


from scipy.ndimage.filters import gaussian_filter
import numpy as np

from Common import MADFPlotData as C_MADFPD


# #########################################################
# CDynamicsBasicPlotData
# #########################################################

class CDynamicsBasicPlotData(C_MADFPD.CMADFPlotData):

    def __init__(self, TimeWindow, PlotParameters, GridSize, LogFile):
        super().__init__(PlotParameters, GridSize, LogFile)
        self.TimeWindow, self.HotSpotMaximum, self.Strength = TimeWindow, 0, 0

    def CalculateStrength(self, StrengthMaskArray):
        for i, j in StrengthMaskArray:
            self.Strength += self.FilteredThresholdData2DNPArray[i, j]

#
# Get - functions
#

    def GetBaseLineInterval(self):
        return self.DictionaryParameter(
            self.PlotParameters, 'BaseLineInterval')

    def GetHotSpotMaximum(self): return self.HotSpotMaximum

    def GetStrength(self): return self.Strength


# #########################################################
# CDynamicsMainPlotData
# #########################################################

class CDynamicsMainPlotData(CDynamicsBasicPlotData):

    def __init__(
            self, Channel, Cluster, TimeWindow, PlotParameters,
            GridSize, LogFile):
        super().__init__(TimeWindow, PlotParameters, GridSize, LogFile)

        self.Channel, self.Cluster, self.StrengthMaskArray = \
            Channel, Cluster, None
        self.SetIJMaximumPosition(self.PlotParameters['IJMaximumPosition'])

#
# Get - functions
#

    def GetStrengthMaskArray(self): return self.StrengthMaskArray

#
# Analysis - functions
#

    def DetectMaximumPosition(self):
        self.HotSpotMaximum = self.FilteredThresholdData2DNPArray.max()
        if self.GetPlotParameters_MaximumAutoDetectFlag():
            Indexes = np.where(
                self.FilteredThresholdData2DNPArray == self.HotSpotMaximum)
            if len(Indexes[0]) == 1 and len(Indexes[1]) == 1:
                self.SetIJMaximumPosition(Indexes[0][0], Indexes[1][0])
            else:
                self.Print((
                    F'There are a few maximum positions in {self.Channel}.'
                    F'{self. Cluster}-', Parameters['Response'],
                    ' plot.\nPreset maximum position will be used.'))

    def CalculateStrengthMaskAndStrength(self):
        StrengthMaskArray = np.zeros(self.RawData2DNPArray.shape)
        I_MP, J_MP = self.GetIJMaximumPosition()
        R = self.GetPlotParameters_FieldRadius()
        for i in range(I_MP-R, I_MP+R+1):
            for j in range(J_MP-R, J_MP+R+1):
                if i >= 0 and i < StrengthMaskArray.shape[0] and\
                        j >= 0 and j < StrengthMaskArray.shape[1]:
                    StrengthMaskArray[i, j] = 1

        FieldTuple = ('FieldRemoveTuple', 'FieldAddTuple')
        for FTCount, KeyString in enumerate(FieldTuple):
            DictionaryParameter = self.DictionaryParameter(
                self.PlotParameters, KeyString)
            if DictionaryParameter:
                for i, j in DictionaryParameter:
                    i, j = i + I_MP, j + J_MP
                    if i >= 0 and i < StrengthMaskArray.shape[0] and\
                            j >= 0 and j < StrengthMaskArray.shape[1]:
                        StrengthMaskArray[i, j] = FTCount

        Indexes = np.where(StrengthMaskArray == 1)
        self.StrengthMaskArray = [
            (I_SMA, J_SMA) for I_SMA, J_SMA in zip(Indexes[0], Indexes[1])]
        self.CalculateStrength(self.StrengthMaskArray)


# #########################################################
# CDynamicsPlotData
# #########################################################

class CDynamicsPlotData(CDynamicsBasicPlotData):

    def __init__(self, TimeWindow, PlotParameters, GridSize, LogFile):
        super().__init__(TimeWindow, PlotParameters, GridSize, LogFile)

    def HotSpotMaximum_Strength(self, IJMaximumPosition, StrengthMaskArray):
        I, J = IJMaximumPosition
        self.HotSpotMaximum = self.FilteredThresholdData2DNPArray[I][J]
        self.CalculateStrength(StrengthMaskArray)
