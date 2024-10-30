#!/usr/bin/env python3
# 24, BereshpolovaLab, University of Connecticut


import MaximumStrength_FirstChannel as MS_FTCH
import MaximumStrength_SecondChannel as MS_SCH

from DynamicsPackage import DynamicsMain as DP_DM
from MaximumStrengthPackage import MSPlexonData as MSP_MSPD


# #########################################################
# Parameters
# #########################################################

PLTS, HSMS, PWRS, FCS_SVGS = '_Plots.png', '_HotSpotMaximum', '_Power', \
    '_FittingCurves.svg'

InputDir = 'Data/'
OutputDir = 'Output/MaximumStrength/'
BasicFileName = '2022Jul26AFsn3-01'
PLXFileName = F'{InputDir}{BasicFileName}.plx'
LogFileName = F'{OutputDir}{BasicFileName}.txt'
PlotGraphicFileName = F'{OutputDir}{BasicFileName}{PLTS}'
HotSpotMaximum_FittingCurve_GraphicFileName = \
    F'{OutputDir}{BasicFileName}{HSMS}{FCS_SVGS}'
Power_FittingCurve_GraphicFileName = \
    F'{OutputDir}{BasicFileName}{PWRS}{FCS_SVGS}'

Parameters = {
    'ProtocolNumber': 11, 'VRR': 160, 'VRR_AutoDetectFlag': True,
    'PlexonDigitalCard': 0,
    'FileInterval': (None, None),
    'TimeWindows': {'T1': 0, 'T2': 100, 'Step': 1, 'Duration': 10},
    'PlotStep': 10,
    'Filter': 'gaussian',
    'GaussianSigma': 1,
    'Interpolation': 'bicubic',

    'Channels': (
        MS_FTCH.ChannelParameters, MS_SCH.ChannelParameters),

    # vertical - within channel, None - individual
    'MainPlotMode': 'vertical',
    # 'vertical', 'horizontal', 'group', None
    'DynamicsPlotMode': 'group',
    # print or not text data of maximums
    'OriginalMaximumShow': True,
    'ModifiedMaximumShow': False,

    # vertical - within dynamics window, channel, horizontal - ON of OFF,
    # None - individual, group - all channels
    'ShowPlotSubfield': {'I1': None, 'I2': None, 'J1': None, 'J2': None},
    'MaximumFieldShow': False,
    'StdThresholdCoefficient': 3,

    'PlotAdjustment': {
        'Top': 0.9, 'Left': 0.03, 'Right': 0.97, 'Bottom': 0.03,
        'HSpace': 0.0, 'WSpace': 0.35},
    'CurveAdjustment': {
        'Top': 0.9, 'Left': 0.03, 'Right': 0.97, 'Bottom': 0.03,
        'HSpace': 0.25, 'WSpace': 0.0},

    'TightLayoutFlag': False,
    'CurveWindowMaximized': False,
    'CurveWindowWidth': 5,
    'CurveWindowHeight': 5,
    'FittingCurveHBin': 20,
    'LinspacePointNumber': 200,
    'HotSpotMaximum_FittingCurve_Colors': ('green', 'yellow'),
    'Power_FittingCurve_Colors': ('red', 'blue'),
    'PlotGraphicFileDPI': 300,
    'HotSpotMaximum_FittingCurve_GraphicFileDPI': 300,
    'Power_FittingCurve_GraphicFileDPI': 300}


# #########################################################
# Analysis
# #########################################################

DP_DM.Main(
    MSP_MSPD.CMSPlexonData, Parameters, PLXFileName, LogFileName,
    PlotGraphicFileName, HotSpotMaximum_FittingCurve_GraphicFileName,
    Power_FittingCurve_GraphicFileName)
