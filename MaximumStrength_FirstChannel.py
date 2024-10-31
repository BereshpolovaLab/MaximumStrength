#!/usr/bin/env python3
# 2024, BereshpolovaLab, University of Connecticut
# Developed by Victor Serdyukov, e-mails: Victor.Serdyukov@uconn.edu


# #########################################################
# Parameters
# #########################################################

ChannelParameters = {
    'Channel': 19,
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
        'TimeWindow': (30, 50),
        'MaximumAutoDetectFlag': True,
        'IJMaximumPosition': (10, 11),
        'FieldRadius': 2,
        'FieldAddTuple': None,
        'FieldRemoveTuple': None,
        'FittingCurveEnable': True,
        'BaseLineInterval': (0, 20)}, {
        # OFF
        'Response': 'Off',
        'TimeWindow': (30, 50),
        'MaximumAutoDetectFlag': True,
        'IJMaximumPosition': (12, 13),
        'FieldRadius': 2,
        'FieldAddTuple': None,
        'FieldRemoveTuple': None,
        'FittingCurveEnable': True,
        'BaseLineInterval': (0, 20)})}
