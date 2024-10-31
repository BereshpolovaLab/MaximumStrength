#!/usr/bin/env python3
# 2024, BereshpolovaLab, University of Connecticut
# Developed by Victor Serdyukov, e-mails: Victor.Serdyukov@uconn.edu


# #########################################################
# Parameters
# #########################################################

ChannelParameters = {
    'Channel': 18,
    'Cluster': 1,
    'Threshold': 0.3,
    'TimeWindowPlotMode': 'OnOff',
    'MaximumIJPlotMode': 'OnOff',
    'FittingCurve': 'polyfit',
    'PolyfitDegree': 15,
    'MinMaxXFit': (None, None),
    'StdLatencyTimeDelay': 20,

    'OnOffPlotParameters': ({
        # ON
        'Response': 'On',
        'TimeWindow': (30, 60),
        'MaximumAutoDetectFlag': True,
        'IJMaximumPosition': (15, 16),
        'FieldRadius': 2,
        'FieldAddTuple': None,
        'FieldRemoveTuple': None,
        'FittingCurveEnable': True,
        'BaseLineInterval': (5, 24)}, {
        # OFF
        'Response': 'Off',
        'TimeWindow': (30, 60),
        'MaximumAutoDetectFlag': True,
        'IJMaximumPosition': (17, 18),
        'FieldRadius': 2,
        'FieldAddTuple': None,
        'FieldRemoveTuple': None,
        'FittingCurveEnable': True,
        'BaseLineInterval': (7, 28)})}
