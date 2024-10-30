#!/usr/bin/env python3se
# 24, BereshpolovaLab, University of Connecticut


from DynamicsPackage import DynamicsPlexonData as DP_DPD
from MaximumStrengthPackage import MSChannelData as MSP_MSCD


###########################################################
# CMSPlexonData
###########################################################

class CMSPlexonData(DP_DPD.CDynamicsPlexonData):

    def __init__(self, Parameters, LogFile):
        super().__init__(Parameters, MSP_MSCD.CMSChannelData, LogFile)

#
# Print - functions
#

    def PrintHotSpotMaximums(self):
        self.PrintCDParameters(
            MSP_MSCD.CMSChannelData.PrintHotSpotsMaximums)

    def PrintStrengths(self):
        self.PrintCDParameters(MSP_MSCD.CMSChannelData.PrintStrengths)

    def PrintSignIndexes(self):
        FuncList = [
            MSP_MSCD.CMSChannelData.GetMPD_Slice_HotSpotMaximum,
            MSP_MSCD.CMSChannelData.GetMPD_Slice_Strength,
            MSP_MSCD.CMSChannelData.GetDPD_Slice_HotSpotMaximum,
            MSP_MSCD.CMSChannelData.GetDPD_Slice_Strength]
        HeaderList = [
            'Main plot maximum', 'Main plot strength', 'Dynamics maximum',
            'Dynamics strength']
        self.PrintSignIndexes_Slice(FuncList, None, HeaderList)

#
# Show/Present - funstions
#

    def ShowPlots(self, PlotGraphicFileName):
        self.ShowPlots_Slice(PlotGraphicFileName)

    def ShowHotSpotMaximumFittingCurves(
            self, HotSpotMaximum_FittingCurve_GraphicFileName):
        self.MaximumMTSD = self.ShowCurves(
            'Maximum in the hot spot', 'Maximum',
            HotSpotMaximum_FittingCurve_GraphicFileName,
            self.GetParameters_HSMFC_Colors(),
            self.GetParameters_HSMFC_GFDPI(),
            MSP_MSCD.CMSChannelData.GetDPD_HotSpotMaximum)

    def ShowPowerFittingCurves(self, Power_FittingCurve_GraphicFileName):
        self.PowerMTSD = self.ShowCurves(
            'Strength of the hot spot field', 'Strength',
            Power_FittingCurve_GraphicFileName,
            self.GetParameters_PFC_Colors(),
            self.GetParameters_PFC_GFDPI(),
            MSP_MSCD.CMSChannelData.GetDPD_Strength)

    def ShowMaximumAndPowerMTSD(self):
        self.ShowPeakParameters(self.MaximumMTSD, 'Maximum peak parameters')
        self.ShowPeakParameters(
            self.PowerMTSD, 'Strength of the hot spot field peak parameters')
