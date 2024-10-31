#!/usr/bin/env python3
# 2024, BereshpolovaLab, University of Connecticut
# Developed by Victor Serdyukov, e-mails: Victor.Serdyukov@uconn.edu


import tabulate

from DynamicsPackage import DynamicsChannelData as DP_DCD


# #########################################################
# CMSChannelData
# #########################################################

class CMSChannelData(DP_DCD.CDynamicsChannelData):

    def __init__(self, TimeWindowList, Parameters, CHP, GridSize, LogFile):
        super().__init__(TimeWindowList, Parameters, CHP, LogFile)

        self.MainPlotDataList, self.DynamicsPlotDataList = \
            self.CreateMPD_DPD_Slices(GridSize)

#
# Get slice - functions
#

    def GetMPD_Slice(self, Index=None): return self.GetMainPlotDataList()

    def GetDPD_Slice(self, Index=None): return self.GetDynamicsPlotDataList()

    def GetMPD_DPDL_Slices(self, Index=None):
        return self.MainPlotDataList, self.DynamicsPlotDataList

#
# Analysis - functions
#

    def FilterThresholdData(self): self.FilterThresholdData_Slice()

    def UpdateMaximums(self): self.UpdateMaximums_Slice()

    def MaximumPositions_HotSpotMaximum_Strength(self): self.OOMP_HSMS_Slice()

    def ReverseCorrelationOnOff(self, OnOffStimuliList):
        self.ReverseCorrelationOnOff_Slice(OnOffStimuliList)
