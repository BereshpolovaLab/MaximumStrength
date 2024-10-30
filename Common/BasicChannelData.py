#!/usr/bin/env python3
# 24, BereshpolovaLab, University of Connecticut


import tabulate

from Common import ErrorObject as C_EO


# ##########################################################
# CBasicChannelData
# ##########################################################

class CBasicChannelData(C_EO.CErrorObject):

    def __init__(self, CHP, LogFile):
        super().__init__(LogFile)
        self.CHP, self.SpikeList = CHP, []

#
# GetCHP - functions
#
    def GetCHP(self): return self.CHP

    def GetCHP_Channel(self):
        return self.DictionaryParameter(self.CHP, 'Channel')

    def GetCHP_Cluster(self):
        return self.DictionaryParameter(self.CHP, 'Cluster')

    def GetCHP_OOPP(self):
        return self.DictionaryParameter(self.CHP, 'OnOffPlotParameters')

    def GetCHP_Response(self):
        return self.DictionaryParameter(self.CHP, 'Response')

    def GetCHP_Threshold(self):
        return self.DictionaryParameter(self.CHP, 'Threshold')

#
# Get - functions
#
    def GetSpikeList(self): return self.SpikeList

    def GetSpikeListLength(self): return len(self.SpikeList)

    def GetStatus(self): return self.SpikeList != []

#
# Set - functions
#
    def SetSpikeList(self, List):
        self.SpikeList = List

#
# Analysis - functions
#

    def AddSpike(self, Channel, Cluster, TimeStamp):
        if self.GetCHP_Channel() == Channel and\
                self.GetCHP_Cluster() == Cluster:
            self.SpikeList.append(TimeStamp)

    def SortSpikes(self):
        self.SpikeList.sort()

    def PlotReverseCorrelationFilterData(
            self, PlotData, OnOffStimuliList, Filter, GaussianSigma):
        SL_Count, Count = 0, 0
        if PlotData.GetPlotParameters_Response() == 'Off':
            SL_Count = 1
        for SE in OnOffStimuliList[SL_Count]:
            Count = PlotData.ReverseCorrelation(SE, self.SpikeList, Count)
            if Count is None:
                break
        PlotData.FilterData(Filter, GaussianSigma)

    def ThresholdData(self, PlotData, Maximum=None):
        PlotData.ThresholdData(self.GetCHP_Threshold(), Maximum)

#
# Print-functions
#
    def PrintSpikeNumbers(self):
        HeaderList, Data2DList = ['CH.CL', 'Numeber of spikes'], [[
            F'{CD.GetCHP_Channel()}.{CD.GetCHP_Cluster()}',
            CD.GetSpikeListLength()] for CD in self.ChannelDataList]
        Tabulate = tabulate.tabulate(Data2DList, headers=HeaderList)
        self.Print((F'{Tabulate}\n\n'))
