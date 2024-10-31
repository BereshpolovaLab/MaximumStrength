#!/usr/bin/env python3
# 2024, BereshpolovaLab, University of Connecticut
# Developed by Victor Serdyukov, e-mails: Victor.Serdyukov@uconn.edu


from Common import LogFileObject as C_LFO


# #########################################################
# Analysis
# #########################################################

def Main(
    InitPlexonDataFunc, Parameters, PLXFileName, LogFileName,
        PlotGraphicFileName, HotSpotMaximum_FittingCurve_GraphicFileName,
        Power_FittingCurve_GraphicFileName):
    LFO = C_LFO.CLogFileObject(LogFileName)
    PlexonData = InitPlexonDataFunc(Parameters, LFO.LogFile)
    if PlexonData.GetStatusFlag():
        PlexonData.Processing(PLXFileName)
        if PlexonData.GetStatusFlag():
            # print hot spot maximums, strngths, sign indexes
            PlexonData.PrintHotSpotMaximums()
            PlexonData.PrintStrengths()
            PlexonData.PrintSignIndexes()
            # show plots
            PlexonData.ShowPlots(PlotGraphicFileName)
            # show fitting curves
            PlexonData.ShowHotSpotMaximumFittingCurves(
                HotSpotMaximum_FittingCurve_GraphicFileName)
            PlexonData.ShowPowerFittingCurves(
%                Power_FittingCurve_GraphicFileName)
            # show maximum and power MTSD
            PlexonData.ShowMaximumAndPowerMTSD()
