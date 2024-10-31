#!/usr/bin/env python3
# 2024, BereshpolovaLab, University of Connecticut
# Developed by Victor Serdyukov, e-mails: Victor.Serdyukov@uconn.edu


import inspect
import numpy as np
import tabulate

from Common import BasicChannelData as C_BCD
from DynamicsPackage import DynamicsPlotData as DP_DPD


# #########################################################
# CDynamicsChannelData
# #########################################################

class CDynamicsChannelData(C_BCD.CBasicChannelData):

    def __init__(self, TimeWindowList, Parameters, CHP, LogFile):
        super().__init__(CHP, LogFile)

        self.TimeWindowList, self.Parameters, self.MainPlotDataList, \
            self.DynamicsPlotDataList = TimeWindowList, Parameters, [], []

#
# GetCHP - functions
#

    def GetCHP_MaximumIJPlotMode(self):
        return self.DictionaryParameter(self.CHP, 'MaximumIJPlotMode')

    def GetCHP_FittingCurve(self):
        return self.DictionaryParameter(self.CHP, 'FittingCurve')

    def GetCHP_PolyfitDegree(self):
        return self.DictionaryParameter(self.CHP, 'PolyfitDegree')

    def GetCHP_MinMaxXFit(self):
        return self.DictionaryParameter(self.CHP, 'MinMaxXFit')

    def GetCHP_StdLatencyTimeDelay(self):
        return self.DictionaryParameter(self.CHP, 'StdLatencyTimeDelay')

#
# GetParameters - functions
#

    def GetParameters_Filter(self):
        return self.DictionaryParameter(self.Parameters, 'Filter')

    def GetParameters_GaussianSigma(self):
        return self.DictionaryParameter(self.Parameters, 'GaussianSigma')

    def GetParameters_MainPlotMode(self):
        return self.DictionaryParameter(self.Parameters, 'MainPlotMode')

    def GetParameters_DynamicsPlotMode(self):
        return self.DictionaryParameter(self.Parameters, 'DynamicsPlotMode')

    def GetParameters_PlotStep(self):
        return self.DictionaryParameter(self.Parameters, 'PlotStep')

    def GetOriginalMaximumShow(self):
        return self.DictionaryParameter(self.Parameters, 'OriginalMaximumShow')

    def GetModifiedMaximumShow(self):
        return self.DictionaryParameter(self.Parameters, 'ModifiedMaximumShow')

#
# Get slice - functions
#

    # get main plot data slice
    def GetMPD_Slice(self, ScaleIndex=None):
        self.ShowError_DCFR()
        raise NotImplementedError()
        return None

    # get dyncamics plot data slice
    def GetDPD_Slice(self, ScaleIndex=None):
        self.ShowError_DCFR()
        raise NotImplementedError()
        return None

    # get main plot data and dyncamics plot data slices
    def GetMPD_DPDL_Slices(self, ScaleIndex=None):
        self.ShowError_DCFR()
        raise NotImplementedError()
        return None, None

#
# GetMPD_Slice - functions
#

    def GetMPD_Slice_HotSpotMaximum(self, ScaleIndex=None):
        return self.GetMPD_Slice_HotSpotMaximum_Strength(
            DP_DPD.CDynamicsBasicPlotData.GetHotSpotMaximum, ScaleIndex)

    def GetMPD_Slice_Strength(self, ScaleIndex=None):
        return self.GetMPD_Slice_HotSpotMaximum_Strength(
            DP_DPD.CDynamicsBasicPlotData.GetStrength, ScaleIndex)

    # get hot spot or strength for on-off responses
    def GetMPD_Slice_HotSpotMaximum_Strength(self, Func, ScaleIndex=None):
        MPDL, RV1, RV2 = self.GetMPD_Slice(ScaleIndex), None, None
        if MPDL[0] is not None and MPDL[1] is not None:
            RV1, RV2 = Func(MPDL[0]), Func(MPDL[1])
        return RV1, RV2

    def GetMPD_Slice_HotSpotMaximum_TitleString(self, ScaleIndex=None):
        MPDS = self.GetMPD_Slice(ScaleIndex)

        def IJMPPSF(RS, Count):
            String = ''
            if MPDS[Count] is not None:
                P0, P1 = MPDS[Count].GetIJMaximumPositionPoint(0), \
                    MPDS[Count].GetIJMaximumPositionPoint(1)
                String = F'{RS} - {P0} : {P1}'
            return String

        TitleString = 'Maximums in the hot spots [I:J]: '
        TitleString += IJMPPSF('On', 0)
        if MPDS[0] is not None and MPDS[1] is not None:
            TitleString += '; '
        TitleString += IJMPPSF('Off', 1)
        return TitleString

#
# GetDPD_Slice - functions
#

    # dynamics plot data hot spot maximum
    def GetDPD_HotSpotMaximum(self, ScaleIndex, OnOffCount, PlotCount):
        return self.GetDPD_PD(
            ScaleIndex, OnOffCount, PlotCount).GetHotSpotMaximum()

    # dynamics plot data strength
    def GetDPD_Strength(self, ScaleIndex, OnOffCount, PlotCount):
        return self.GetDPD_PD(
            ScaleIndex, OnOffCount, PlotCount).GetStrength()

    # get plot data of dynamics plot data slice
    def GetDPD_PD(self, ScaleIndex, OnOffCount, PlotCount):
        return self.GetDPD_Slice(ScaleIndex)[OnOffCount][PlotCount]

    # dynamics plot data hot spot maximum of dynamics plot data slice
    def GetDPD_Slice_HotSpotMaximum(self, ScaleIndex=None):
        return self.GetDPD_Slice_HotSpotMaximum_Strength(
            DP_DPD.CDynamicsBasicPlotData.GetHotSpotMaximum, ScaleIndex)

    # dynamics plot data strength of dynamics plot data slice
    def GetDPD_Slice_Strength(self, ScaleIndex=None):
        return self.GetDPD_Slice_HotSpotMaximum_Strength(
            DP_DPD.CDynamicsBasicPlotData.GetStrength, ScaleIndex)

    # dynamics plot data hot spot maximum and strength
    # of dynamics plot data slice
    def GetDPD_Slice_HotSpotMaximum_Strength(self, Func, ScaleIndex=None):

        def GetPeakHeight(DPDL, Count):
            RV = None
            Maximum = max([Func(DPD) for DPD in DPDL])
            T1, T2 = self.CHP['OnOffPlotParameters'][Count]['BaseLineInterval']
            List = [
                Func(DPD) for DPD in DPDL if DPD.GetTimeWindowPoint(0) >= T1
                and DPD.GetTimeWindowPoint(1) <= T2]
            Mean = np.mean(List) if List else 0
            RV = Maximum - Mean
            return RV

        OnDPDL, OffDPDL = self.GetDPD_Slice(ScaleIndex)
        RV1, RV2 = None, None

        if OnDPDL is not None and OffDPDL is not None:
            RV1, RV2 = GetPeakHeight(OnDPDL, 0), GetPeakHeight(OffDPDL, 1)
        return RV1, RV2

#
# Main plot data list functions
#

    def GetMainPlotDataList(self): return self.MainPlotDataList

#
# Dynamics plot data list functions
#

    def GetDynamicsPlotDataList(self): return self.DynamicsPlotDataList

#
# Get - functions
#

    def GetBaseLineIntervals(self, ScaleIndex=None):
        List = [None, None]
        for Count, PD in enumerate(self.GetMPD_Slice(ScaleIndex)):
            if PD is not None:
                List[Count] = PD.GetBaseLineInterval()
        return tuple(List)

#
# Creation - functions
#

    def CreateMPD_DPD_Slices(self, GridSize):
        MPDS, DPDS = [], []
        for Count, OOPP in enumerate(self.GetCHP_OOPP()):
            MPD = DP_DPD.CDynamicsMainPlotData(
                self.GetCHP_Channel(), self.GetCHP_Cluster(),
                OOPP['TimeWindow'], OOPP, GridSize, self.LogFile)\
                if OOPP is not None else None
            DPDL = [
                DP_DPD.CDynamicsPlotData(TW, OOPP, GridSize, self.LogFile)
                for TW in self.TimeWindowList] if OOPP is not None else None
            MPDS.append(MPD)
            DPDS.append(DPDL)
        return MPDS, DPDS

#
# Analysis functions
#

    def ReverseCorrelationOnOff(self, OnOffStimuliList):
        self.ShowError_DCFR()
        raise NotImplementedError()

    def FilterThresholdData(self):
        self.ShowError_DCFR()
        raise NotImplementedError()

    def UpdateMaximums(self):
        self.ShowError_DCFR()
        raise NotImplementedError()

    def MaximumPositions_HotSpotMaximum_Strength(self):
        self.ShowError_DCFR()
        raise NotImplementedError()

#
#  Analysis slice - functions
#

    def ReverseCorrelationOnOff_Slice(self, OnOffStimuliList, ScaleIndex=None):
        MPDS, DPDS = self.GetMPD_DPDL_Slices(ScaleIndex)
        for SL, MPD, DPDL in zip(OnOffStimuliList, MPDS, DPDS):
            if MPD is not None and DPDL is not None:
                self.RCLoop(self.MainRC, StimuliList=SL, PlotData=MPD)
                self.RCLoop(self.DynamicsRC, StimuliList=SL, PlotData=DPDL)

    def FilterThresholdData_Slice(self, ScaleIndex=None):

        def FTD(PD):
            PD.FilterThresholdData(
                self.GetParameters_Filter(),
                self.GetParameters_GaussianSigma(), self.GetCHP_Threshold())

        MPDS, DPDS = self.GetMPD_DPDL_Slices(ScaleIndex)
        for MPD, DPDL in zip(MPDS, DPDS):
            if MPD is not None and DPDL is not None:
                FTD(MPD)
                for DPD in DPDL:
                    FTD(DPD)

    def UpdateMaximums_Slice(self, ScaleIndex=None):

        def WarningFunc(PlotModeString, PlotModeValueString):
            self.Print((
                'Warning!!! \'', PlotModeString, '\'-parameter with \'',
                PlotModeValueString, '\'-value will be ignored for ',
                F'{self.GetCHP_Channel()}.{self.GetCHP_Cluster()}',
                ' channel because only one \'On\' or \'Off\' response data ',
                'is presnet.\n'))

        def SetMaximums(MPDL):
            Maximum = max([PD.GetMaximum() for PD in MPDL])
            for PD in MPDL:
                PD.SetMaximum(Maximum)

        if self.GetOriginalMaximumShow():
            self.PrintMaximums_Slice('Original', ScaleIndex)

        if self.GetParameters_MainPlotMode() == 'vertical':
            MPDS = self.GetMPD_Slice(ScaleIndex)
            if MPDS[0] is not None and MPDS[1] is not None:
                SetMaximums(MPDS)
            else:
                WarningFunc('MainPlotMode', 'vertical')

        if self.GetParameters_DynamicsPlotMode().lower() == 'horizontal':
            DPDS = self.GetDPD_Slice(ScaleIndex)
            for DPDL in DPDS:
                if DPDL is not None:
                    SetMaximums(DPDL)

        elif self.GetParameters_DynamicsPlotMode().lower() == 'vertical':
            DPDS = self.GetDPD_Slice(ScaleIndex)
            if DPDS[0] is not None and DPDS[1] is not None:
                for DPD1, DPD2 in zip(DPDS[0], DPDS[1]):
                    Maximum = max(DPD1.GetMaximum(), DPD2.GetMaximum())
                    for DPD in {DPD1, DPD2}:
                        DPD.SetMaximum(Maximum)
            else:
                WarningFunc('DynamicsPlotMode', 'vertical')

        elif self.GetParameters_DynamicsPlotMode().lower() == 'group':
            DPDS = self.GetDPD_Slice(ScaleIndex)
            if DPDS[0] is not None and DPDS[1] is not None:
                Maximum = max([DPD.GetMaximum() for DPD in DPDS[0] + DPDS[1]])
                for DPDL in DPDS:
                    for DPD in DPDL:
                        DPD.SetMaximum(Maximum)
            else:
                WarningFunc('DynamicsPlotMode', 'group')

        if self.GetModifiedMaximumShow():
            self.PrintMaximums_Slice('Modified', ScaleIndex)

    def OOMP_HSMS_Slice(self, ScaleIndex=None):

        def WarningFunc(RS):
            self.Print((
                'Warning!!! \'MaximumIJPlotMode\'-parameter with \'',
                RS, '\'-value will be ignored for ',
                F'{self.GetCHP_Channel()}.{self.GetCHP_Cluster()}',
                ' channel because \'', RS, '\' response data is None.\n'))

        MPDS, DPDS = self.GetMPD_DPDL_Slices(ScaleIndex)
        for MPD in MPDS:
            if MPD is not None:
                MPD.DetectMaximumPosition()
        if self.GetCHP_MaximumIJPlotMode().lower() == 'on':
            if MPDS[0] is not None and MPDS[1] is not None:
                MPDS[1].SetIJMaximumPosition(MPDS[0].GetIJMaximumPosition())
            elif MPDS[0] is None:
                self.WarningFunc2('On')
        elif self.GetCHP_MaximumIJPlotMode().lower() == 'off':
            if MPDS[0] is not None and MPDS[1] is not None:
                MPDS[0].SetIJMaximumPosition(MPDS[1].GetIJMaximumPosition())
            elif MPDS[0] is None:
                self.WarningFunc('Off')

        for MPD, DPDL in zip(MPDS, DPDS):
            if MPD is not None and DPDS is not None:
                MPD.CalculateStrengthMaskAndStrength()
                HotSpotPosition, StrengthMask = MPD.GetIJMaximumPosition(), \
                    MPD.GetStrengthMaskArray()
                for DP in DPDL:
                    DP.HotSpotMaximum_Strength(HotSpotPosition, StrengthMask)

#
#  Analysis service - functions
#

    def RCLoop(self, RCFunc, **KWArgs):
        Count = 0
        for SE in KWArgs['StimuliList']:
            KWArgs['StimuliEvent'], KWArgs['Count'] = SE, Count
            Count = RCFunc(**KWArgs)
            if Count is None:
                break

    def MainRC(self,  **KWArgs):
        return KWArgs['PlotData'].ReverseCorrelation(
            KWArgs['StimuliEvent'], self.SpikeList, KWArgs['Count'])

    def DynamicsRC(self, **KWArgs):
        Count = 0
        for DPD in KWArgs['PlotData']:
            KWArgs['PlotData'], KWArgs['Count'] = DPD, Count
            Count = self.MainRC(**KWArgs)
            if Count is None:
                break
        return Count

    def OOMP_HSMS(self, ScaleIndex):
        for MPD in MPDL:
            MPD.DetectMaximumPosition()
        if self.GetCHP_MaximumIJPlotMode().lower() == 'on':
            MPDL[1].SetIJMaximumPosition(MPDL[0].GetIJMaximumPosition())
        elif self.GetCHP_MaximumIJPlotMode().lower() == 'off':
            MPDL[0].SetIJMaximumPosition(MPDL[1].GetIJMaximumPosition())
        for MPD, DPD in zip(MPDL, DPDL):
            MPD.CalculateStrengthMaskAndStrength()
            for DP in DPD:
                DP.HotSpotMaximum_Strength(
                    MPD.GetIJMaximumPosition(), MPD.StrengthMaskArray)

#
# Print - function
#

    def PrintMaximums_Slice(self, TitleString, ScaleIndex=None):
        TitleString = 'Maximums - ' + TitleString
        self.PrintPlotData(
            lambda PD: PD.GetMaximum(), TitleString, ScaleIndex)

    def PrintHotSpotsMaximums(self, ScaleIndex=None):
        self.PrintPlotData(
            lambda PD: PD.GetHotSpotMaximum(),
            self.GetMPD_Slice_HotSpotMaximum_TitleString(ScaleIndex),
            ScaleIndex)

    def PrintStrengths(self, ScaleIndex):
        TitleString = 'Strength of the hot spot field'
        self.PrintPlotData(
            lambda PD: PD.GetStrength(), TitleString, ScaleIndex)

    def PrintPlotData(self, PlotDataFunc, TitleString, ScaleIndex):
        Count, String = None, ''
        MPDS, DPDS = self.GetMPD_DPDL_Slices(ScaleIndex)
        if MPDS[0] is not None and MPDS[1] is not None:

            if MPDS[0].GetTimeWindow() != MPDS[1].GetTimeWindow():
                t00, t01, t10, t11 = \
                    MPDS[0].GetTimeWindowPoint(0), \
                    MPDS[0].GetTimeWindowPoint(1), \
                    MPDS[1].GetTimeWindowPoint(0), \
                    MPDS[1].GetTimeWindowPoint(1)
                String = F'{t00}-{t01} / {t10}-{t11} ms'
            else:
                Count = 0
        elif MPDS[0] is not None and MPDS[1] is None:
            Count = 0
        elif MPDS[0] is None and MPDS[1] is not None:
            Count = 1
        if Count is not None:
            t0, t1 = MPDS[Count].GetTimeWindowPoint(0), \
                MPDS[Count].GetTimeWindowPoint(1)
            String = F'{t0} - {t1} ms'

        HeaderList = [F'{0} ms {String}'] + [
            F'{TW[0]} - {TW[1]} ms'
            for TWCount, TW in enumerate(self.TimeWindowList)
            if TWCount % self.GetParameters_PlotStep() == 0]

        RL, Data2DList = ('On', 'Off'), []
        for RS, MPD, DPDL in zip(RL, MPDS, DPDS):
            if MPD is not None and DPDS is not None:
                Data2DList.append([RS, PlotDataFunc(MPD)] + [
                    PlotDataFunc(DPD) for PlotCount, DPD in enumerate(DPDL)
                    if PlotCount % self.GetParameters_PlotStep() == 0])
        self.PrintData(TitleString, Data2DList, HeaderList)

    def PrintData(self, TitleString, Data2DList, HeaderList):
        TitleString, Tabulate = (
            F'{TitleString}\nChannel: {self.GetCHP_Channel()}, '
            F'Cluster: {self.GetCHP_Cluster()}'), \
            tabulate.tabulate(Data2DList, headers=HeaderList)
        self.Print((TitleString, '\n', Tabulate, '\n\n'))
