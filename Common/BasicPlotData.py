#!/usr/bin/env python3
# 2024, BereshpolovaLab, University of Connecticut
# Developed by Victor Serdyukov, e-mails: Victor.Serdyukov@uconn.edu


import numpy as np

from scipy.ndimage.filters import gaussian_filter
from Common import ErrorObject as C_EO
from Common import ToolObject as C_TO


# #########################################################
# CBasicPlotData
# #########################################################

class CBasicPlotData(C_EO.CErrorObject, C_TO.CToolObject):

    def __init__(self, PlotParameters, GridSize, LogFile):
        super().__init__(LogFile)
        self.PlotParameters, self.GridSize, self.RawData2DNPArray, \
            self.FilteredData2DNPArray, self.FilteredThresholdData2DNPArray, \
            self.Maximum = PlotParameters, GridSize, np.zeros((
                GridSize[1], GridSize[0])), None, None, None
        self.TimeWindow = PlotParameters['TimeWindow']\
            if 'TimeWindow' in PlotParameters else None
        self.SetIJMaximumPosition(None)

#
# GetParameters - functions
#

    def GetPlotParameters(self): return self.PlotParameters

    def GetPlotParameters_Response(self):
        return self.DictionaryParameter(self.PlotParameters, 'Response')

#
# Get - functions
#

    def GetGridSize(self): return self.GridSize

    def GetMaximum(self): return self.Maximum

    def GetTimeWindow(self): return self.TimeWindow

    def GetTimeWindowPoint(self, Count): return self.TimeWindow[Count]

    def GetIJMaximumPosition(self): return self.IJMaximumPosition

    def GetRawData2DNPArray(self): return self.RawData2DNPArray

    def GetFilteredData2DNPArray(self): return self.FilteredData2DNPArray

    def GetRawData2DNPArray_Shape(self): return self.RawData2DNPArray.shape

    def GetIJMaximumPositionPoint(self, Count):
        return self.IJMaximumPosition[Count]

    def GetFilteredThresholdData2DNPArray(self):
        return self.FilteredThresholdData2DNPArray

    def GetFilteredData2DNPArray_Shape(self):
        return self.FilteredData2DNPArray.shape

    def GetFilteredThresholdData2DNPArray_Shape(self):
        return self.FilteredThresholdData2DNPArray.shape

#
# Set - functions
#

    def SetMaximum(self, Maximum): self.Maximum = Maximum

    def SetIJMaximumPosition(self, *args):
        if len(args) == 1:
            self.IJMaximumPosition = args[0]
        elif len(args) == 2:
            self.IJMaximumPosition = [args[0], args[1]]

#
# Analysis - functions
#

    def ReverseCorrelation(self, SE, SpikeList, Count):
        SN, StartCount = self.ReverseCorrelation0(
            SE, SpikeList, self.GetTimeWindow(), Count)
        self.RawData2DNPArray[SE.GetGridI(), SE.GetGridJ()] += SN
        return StartCount

    def FilterData(self, Filter, GaussianSigma):
        self.FilteredData2DNPArray = gaussian_filter(
            self.RawData2DNPArray, GaussianSigma)\
            if Filter == 'gaussian' else self.RawData2DNPArray
        self.FilteredThresholdData2DNPArray = self.FilteredData2DNPArray
        self.Maximum = self.FilteredThresholdData2DNPArray.max()

    def ThresholdData(self, Threshold, Maximum=None):
        if Maximum is None:
            Maximum = self.Maximum
        self.FilteredThresholdData2DNPArray[
            self.FilteredThresholdData2DNPArray < Maximum * Threshold] = 0
        self.Maximum = self.FilteredThresholdData2DNPArray.max()

    def FilterThresholdData(
            self, Filter, GaussianSigma, Threshold, Maximum=None):
        self.FilterData(Filter, GaussianSigma)
        self.ThresholdData(Threshold, Maximum)
