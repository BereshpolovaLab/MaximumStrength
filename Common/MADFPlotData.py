#!/usr/bin/env python3
# 2024, BereshpolovaLab, University of Connecticut
# Developed by Victor Serdyukov, e-mails: Victor.Serdyukov@uconn.edu


from Common import BasicPlotData as C_BPD


# #########################################################
# CMADFBasicPlotData
# #########################################################

class CMADFPlotData(C_BPD.CBasicPlotData):

    def __init__(self, PlotParameters, GridSize, LogFile):
        super().__init__(PlotParameters, GridSize, LogFile)

#
# GetPlotParameters - functions
#

    def GetPlotParameters_MaximumAutoDetectFlag(self):
        return self.DictionaryParameter(
            self.PlotParameters, 'MaximumAutoDetectFlag')

    def GetPlotParameters_FieldRadius(self):
        return self.DictionaryParameter(
            self.PlotParameters, 'FieldRadius')
