#!/usr/bin/env python3
# 2024, BereshpolovaLab, University of Connecticut
# Developed by Victor Serdyukov, e-mails: Victor.Serdyukov@uconn.edu


from matplotlib.path import Path
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
from numpy import asarray as ar
from scipy.ndimage.filters import gaussian_filter
from scipy.optimize import curve_fit
from scipy.optimize import leastsq
from scipy.optimize import minimize
from scipy.special import erf
from sklearn.metrics import r2_score
import matplotlib as mpl
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import tabulate


from Common import BasicPlexonData as C_BPD
from Common import ToolObject as C_TO


###########################################################
# CDynamicsPlexonData
###########################################################

class CDynamicsPlexonData(C_BPD.CBasicPlexonData, C_TO.CToolObject):

    def __init__(self, Parameters, ConstructorFunc, LogFile):
        super().__init__(Parameters, LogFile)
        self.MaximumMTSD, self.PowerMTSD = [], []

        self.CheckAndModifyParameters()
        if self.GetStatusFlag():
            self.TimeWindowList = self.GetTimeWindowList()

            self.ChannelDataList = [
                ConstructorFunc(
                    self.TimeWindowList, self.Parameters, CHP,
                    self.GridSize, LogFile)
                for CHP in self.GetParameters_Channels()]

#
# GetParameters - functions
#

    def GetParameters_TimeWindows(self):
        return self.DictionaryParameter(self.Parameters, 'TimeWindows')

    def GetParameters_PlotStep(self):
        return self.DictionaryParameter(self.Parameters, 'PlotStep')

    def GetParameters_ShowPlotSubfield(self):
        return self.DictionaryParameter(self.Parameters, 'ShowPlotSubfield')

    def GetParameters_MaximumFieldShow(self):
        return self.DictionaryParameter(self.Parameters, 'MaximumFieldShow')

    def GetParameters_TightLayoutFlag(self):
        return self.DictionaryParameter(self.Parameters, 'TightLayoutFlag')

# plot adjustment

    def GetParameters_PlotAdjustment_Top(self):
        return self.DictInDictParameter(
            self.Parameters, 'PlotAdjustment', 'Top')

    def GetParameters_PlotAdjustment_Left(self):
        return self.DictInDictParameter(
            self.Parameters, 'PlotAdjustment', 'Left')

    def GetParameters_PlotAdjustment_Right(self):
        return self.DictInDictParameter(
            self.Parameters, 'PlotAdjustment', 'Right')

    def GetParameters_PlotAdjustment_Bottom(self):
        return self.DictInDictParameter(
                self.Parameters, 'PlotAdjustment', 'Bottom')

    def GetParameters_PlotAdjustment_HSpace(self):
        return self.DictInDictParameter(
            self.Parameters, 'PlotAdjustment', 'HSpace')

    def GetParameters_PlotAdjustment_WSpace(self):
        return self.DictInDictParameter(
            self.Parameters, 'PlotAdjustment', 'WSpace')

# curve adjustment

    def GetParameters_CurveAdjustment_Top(self):
        return self.DictInDictParameter(
            self.Parameters, 'CurveAdjustment', 'Top')

    def GetParameters_CurveAdjustment_Left(self):
        return self.DictInDictParameter(
            self.Parameters, 'CurveAdjustment', 'Left')

    def GetParameters_CurveAdjustment_Right(self):
        return self.DictInDictParameter(
            self.Parameters, 'CurveAdjustment', 'Right')

    def GetParameters_CurveAdjustment_Bottom(self):
        return self.DictInDictParameter(
            self.Parameters, 'CurveAdjustment', 'Bottom')

    def GetParameters_CurveAdjustment_HSpace(self):
        return self.DictInDictParameter(
            self.Parameters, 'CurveAdjustment', 'HSpace')

    def GetParameters_CurveAdjustment_WSpace(self):
        return self.DictInDictParameter(
            self.Parameters, 'CurveAdjustment', 'WSpace')

# plot graphic file DPI

    def GetParameters_PlotGraphicFileDPI(self):
        return self.DictionaryParameter(self.Parameters, 'PlotGraphicFileDPI')

# hot spot maximum fitting curve

    def GetParameters_HSMFC_Colors(self):
        return self.DictionaryParameter(
            self.Parameters, 'HotSpotMaximum_FittingCurve_Colors')

    def GetParameters_HSMFC_GFDPI(self):
        return self.DictionaryParameter(
            self.Parameters, 'HotSpotMaximum_FittingCurve_GraphicFileDPI')

#  power fitting curve

    def GetParameters_PFC_Colors(self):
        return self.DictionaryParameter(
            self.Parameters, 'Power_FittingCurve_Colors')

    def GetParameters_PFC_GFDPI(self):
        return self.DictionaryParameter(
            self.Parameters, 'Power_FittingCurve_GraphicFileDPI')

# curve window

    def GetParameters_CW_Maximized(self):
        return self.DictionaryParameter(
            self.Parameters, 'CurveWindowMaximized')

    def GetParameters_CW_Width(self):
        return self.DictionaryParameter(self.Parameters, 'CurveWindowWidth')

    def GetParameters_CW_Height(self):
        return self.DictionaryParameter(self.Parameters, 'CurveWindowHeight')

# curve

    def GetParameters_LSPN(self):
        return self.DictionaryParameter(self.Parameters, 'LinspacePointNumber')

    def GetParameters_STC(self):
        return self.DictionaryParameter(
            self.Parameters, 'StdThresholdCoefficient')

    def GetParameters_FCHB(self):
        return self.DictionaryParameter(self.Parameters, 'FittingCurveHBin')

#
# Get - functions
#

    def GetTimeWindowList(self):
        TW = self.GetParameters_TimeWindows()
        return [
            (t1, t1 + self.DictionaryParameter(TW, 'Duration'))
            for t1 in range(
                self.DictionaryParameter(TW, 'T1'),
                self.DictionaryParameter(TW, 'T2') + 1,
                self.DictionaryParameter(TW, 'Step'))]

#
# Analysis - functions
#

    def CheckAndModifyOnOffParameters(self):

        def PrintError1():
            self.Print((
                'Error!!! THe length of \'OnOffPlotParameters\' ',
                'is bigger than 2.\n\n'))

        def PrintError2():
            self.Print((
                'Error!!! \'Response\'-parameters are not',
                '\'On\' and \'Off\'.'))

        def PrintError3(RS):
            self.Print((
                'Error!!! Both parameters in \'OnOffPlotParameters\' are ',
                '\'{0}\'.\n\n'.format(RS)))

        def PrintError4():
            self.Print((
                'Error!!! Please check \'OnOffPlotParameters\'-',
                'parameters in the configuration file. ',
                'Something is wrong.\n\n'))

        FlagList = []
        for CHP in self.GetParameters_Channels():
            Flag = False
            Length = len(self.DictionaryParameter(CHP, 'OnOffPlotParameters'))
            if Length > 2:
                PrintError()
            elif Length == 2:
                if 'Response' in CHP['OnOffPlotParameters'][0] and\
                        'Response' in CHP['OnOffPlotParameters'][1]:
                    RS1 = CHP['OnOffPlotParameters'][0]['Response'].lower()
                    RS2 = CHP['OnOffPlotParameters'][1]['Response'].lower()
                    if RS1 != RS2:
                        if (RS1 == 'on' and RS2 == 'off') or (
                                RS1 == 'off' and RS2 == 'on'):
                            Flag = True
                            if RS1 == 'off':
                                CHP['OnOffPlotParameters'] = (
                                    CHP['OnOffPlotParameters'][1],
                                    CHP['OnOffPlotParameters'][0])
                        else:
                            PrintError2()
                    else:
                        PrintError3(RS1)
            elif Length == 1:
                if 'Response' in CHP['OnOffPlotParameters'][0]:
                    RS = CHP['OnOffPlotParameters'][0]['Response'].lower()
                    if RS == 'on' or RS == 'off':
                        Flag = True
                        CHP['OnOffPlotParameters'] = (
                            CHP['OnOffPlotParameters'][0], None) if RS == 'on'\
                            else (None, CHP['OnOffPlotParameters'][0])
            FlagList.append(Flag)

        if False in FlagList:
            if not self.GetStatusFlag():
                PrintError4()
        else:
            self.SetStatusFlag(True)

    def CheckAndModifyOnOffTimeWindows(self):

        def GetTW(CHP, TWPM, Count, RS):

            def PrintError(CHP, RS):
                CH, CL = self.DictionaryParameter(CHP, 'Channel'), \
                    self.DictionaryParameter(CHP, 'Cluster')
                self.Print((
                    'Warning!!! \'TimeWindowPlotMode\'-parameter with '
                    F'\'{RS}\'-value will be ignored for {CH}.{CL} ',
                    F'channel because \'{RS}\' response data is None.\n'))

            OOPP, TW = self.DictionaryParameter(
                CHP, 'OnOffPlotParameters')[Count], None
            if TWPM == RS:
                if OOPP is not None:
                    TW = self.DictionaryParameter(OOPP, 'TimeWindow')
                else:
                    PrintError(CHP, RS)
            return OOPP, TW

        for CHP in self.GetParameters_Channels():
            TWPM = self.DictionaryParameter(CHP, 'TimeWindowPlotMode').lower()
            if TWPM == 'on' or TWPM == 'off':
                OOPP1, TW1 = GetTW(CHP, TWPM, 0, 'on')
                OOPP2, TW2 = GetTW(CHP, TWPM, 1, 'off')
                List = [[TW1, OOPP2], [TW2, OOPP1]]
                for TW, OOPP in List:
                    if TW is not None and OOPP is not None:
                        OOPP['TimeWindow'] = TW

    def CheckAndModifyParameters(self):
        self.CheckAndModifyOnOffParameters()
        self.CheckAndModifyOnOffTimeWindows()

    def Processing(self, PLXFileName):

        def PrintError():
            self.Print((
                'Error!!! DynamicChannelData->Processing.\n'
                'A called function is not realized in derived class.\n'))

        self.LoadData_CheckSpikesAndTTLEvents(PLXFileName)
        self.CheckStimuli()
        if self.StatusFlag:
            self.StimulusTimeStampCorrection()
            self.PrintSpikeNumbers()
            for CD in self.ChannelDataList:
                try:
                    CD.ReverseCorrelationOnOff(self.OnOffStimuliList)
                    CD.FilterThresholdData()
                    CD.UpdateMaximums()
                    CD.MaximumPositions_HotSpotMaximum_Strength()
                except NotImplementedError:
                    self.StatusFlag = False
                    PrintError()
                    break

#
# Print - functions
#

    def PrintHotSpotMaximums(self):
        self.ShowError_DCFR()
        raise NotImplementedError()

    def PrintStrengths(self):
        self.ShowError_DCFR()
        raise NotImplementedError()

    def PrintSignIndexes(self):
        self.ShowError_DCFR()
        raise NotImplementedError()

    def PrintCDParameters(
           self, CDFunc, ScaleIndex=None,
            HeaderList2=['Maximum', 'Strength']):
        for CD in self.ChannelDataList:
            CDFunc(CD, ScaleIndex)

    def PrintSignIndexes_Slice(
            self, FuncList, ScaleIndex=None,
            HeaderList=['Maximum', 'Strength']):

        def Func(v1, v2):
            RV = (v1 - v2) / (v1 + v2)
            return self.FMT(RV)

        def PrintError(Count, CCS):
            self.Print((
                F'Error!!! Can not provide sign index for {CCS} channel ',
                F'N{Count + 1} in channel list. ',
                'Please check configuration file.'))

        self.Print(('On/Off sign indexes'))
        HeaderList, VList, ErrorList = ['CH.CL'] + HeaderList, [], []

        for Count, CD in enumerate(self.ChannelDataList):
            CCS = '{0}.{1}'.format(
                CD.GetCHP_Channel(), CD.GetCHP_Cluster())
            List = [CCS]

            for F in FuncList:
                v1, v2 = F(CD, ScaleIndex)
                if v1 and v2:
                    V = Func(v1, v2)
                else:
                    V = 'None'
                    PrintError(Count, CCS)
                List.append(V)

            VList.append(List)

        T = tabulate.tabulate(VList, headers=HeaderList)
        self.Print((T, '\n'))

#
# Show/Present-funstions
#

    def ShowPlots(self, PlotGraphicFileName):
        self.ShowError_DCFR()
        raise NotImplementedError()

    def ShowFittingCurves_HotSpotMaximum(
            self, FittingCurves_HotSpotMaximum_GraphicFileName):
        self.ShowError_DCFR()
        raise NotImplementedError()

    def ShowPlots_Slice(
            self, PlotGraphicFileName, ScaleIndex=None, TitleExtention=None):

        fig, aspect, pad_fraction, Title = plt.figure(), 20, 1.0, \
            'Plots and dynamics plots'
        if TitleExtention is not None:
            Title += ' - ' + TitleExtention
        fig.suptitle(Title, fontsize=20)
        fig.subplots_adjust(
            top=self.GetParameters_PlotAdjustment_Top(),
            left=self.GetParameters_PlotAdjustment_Left(),
            right=self.GetParameters_PlotAdjustment_Right(),
            bottom=self.GetParameters_PlotAdjustment_Bottom(),
            hspace=self.GetParameters_PlotAdjustment_HSpace(),
            wspace=self.GetParameters_PlotAdjustment_WSpace())
        RowNumber, JL = 0, range(len(self.TimeWindowList))
        if self.GetParameters_PlotStep() and\
                self.GetParameters_PlotStep() > 1:
            JL = [
                PlotCount for PlotCount in JL
                if PlotCount % self.GetParameters_PlotStep() == 0]

        for CD in self.ChannelDataList:
            MPDS, DPDS = CD.GetMPD_DPDL_Slices(ScaleIndex)
            for i in range(2):
                if MPDS[i] is not None and DPDS[i] is not None:
                    RowNumber += 1

        GridSize, RowCount, RL = (RowNumber, len(JL) + 1), 0, {'On', 'Off'}

        def PPFunc(CD, PD, RS, SNA, ColumnCount=0):
            self.PresentPlot(
                CD, PD, GridSize, aspect, pad_fraction, RS, SMA,
                RowCount, ColumnCount)

        for CD in self.ChannelDataList:
            MPDS, DPDS = CD.GetMPD_DPDL_Slices(ScaleIndex)
            for i, RS in enumerate(RL):
                MPD, DPDL = MPDS[i], DPDS[i]
                if MPD is not None and DPDL is not None:
                    SMA = MPD.GetStrengthMaskArray()

                    # main plot data
                    PPFunc(CD, MPD, RS, SMA)
                    # dynamics plot datas
                    for ColumnCount, DPDCount in enumerate(JL):
                        PPFunc(CD, DPDL[DPDCount], RS, SMA, ColumnCount + 1)

                    RowCount += 1

        if self.Parameters['TightLayoutFlag']:
            plt.tight_layout()
        plt.get_current_fig_manager().window.showMaximized()
        plt.show()
        self.SaveFigure(
            fig, PlotGraphicFileName,
            self.GetParameters_PlotGraphicFileDPI(), 'plot')

    def PresentPlot(
            self, CD, PD, GridSize, aspect, pad_fraction, RS, SMA,
            RowCount, ColumnCount):

        def PrintError(CD, PD):
            self.Print((
                F'{CD.GetCHP_Channel()}.{CD.GetCHP_Cluster()}-{RS} '
                F'{PD.GetTimeWindowPoint(0)} - {PD.GetTimeWindowPoint(1)} ms: '
                'strength mask is empty.'))

        def Func(PSPSF, Shape, Count, s1, s2):

            def LF(s, V):
                return PSPSF[s] if PSPSF[s] > 0 and\
                    PSPSF[s] < Shape[Count] else V

            v1, v2 = LF(s1, 0), LF(s2, Shape[Count])
            if v1 > v2:
                v1, v2 = v2, v1
            return v1, v2

        ax, PSPSF, NPA, Extent = plt.subplot2grid(
            GridSize, (RowCount, ColumnCount)), \
            self.GetParameters_ShowPlotSubfield(), None, None
        I1, I2, J1, J2 = PSPSF['I1'], PSPSF['I2'], PSPSF['J1'], PSPSF['J2']

        if I1 is None and I2 is None and J1 is None and J2 is None:
            NPA = PD.GetFilteredThresholdData2DNPArray()
        else:
            Shape = PD.GetFilteredThresholdData2DNPArray_Shape()
            (I1, I2), (J1, J2) = Func(PSPSF, Shape, 0, 'I1', 'I2'), \
                Func(PSPSF, Shape, 1, 'J1', 'J2')
            NPA = PD.GetFilteredThresholdData2DNPArray()[I1: I2, J1: J2]
            Extent = [J1-0.5, J2-0.5, I2-0.5, I1-0.5]

        Title = (
            F'{CD.GetCHP_Channel()}.{CD.GetCHP_Cluster()}-{RS}: '
            F'{PD.GetTimeWindowPoint(0)} - {PD.GetTimeWindowPoint(1)} ms\n'
            F'Spikes: {CD.GetSpikeListLength()}')
        ax.set_title(Title)

        im = ax.imshow(
            NPA, interpolation=self.GetParameters_Interpolation(),
            cmap='jet', vmax=PD.Maximum, extent=Extent)
        divider, width = make_axes_locatable(ax), \
            axes_size.AxesY(ax, aspect=1/aspect)
        pad = axes_size.Fraction(pad_fraction, width)
        cax = divider.append_axes('right', size=width, pad=pad)
        plt.colorbar(im, cax=cax)
        if self.GetParameters_MaximumFieldShow():
            if SMA:
                for SMA_I, SMA_J in SMA:
                    ax.add_patch(
                        mpl.patches.Rectangle(
                            (SMA_J-0.5, SMA_I-0.5), 1, 1, hatch='//////',
                            fill=False, linewidth=0.3, edgecolor='black'))
            else:
                PrintError(CD, PD)

    def ShowCurves(
            self, TitleString, YLabel, GraphicFileName, Colors,
            GraphicFileDPI, Func, ScaleIndex=None):

        fig, MTSDL, TWL = self.ShowCurves_Init(TitleString), [], \
            [TW for j, TW in enumerate(self.TimeWindowList)]
        for RowCount, CD in enumerate(self.ChannelDataList):
            ax, OnOffData, FittingCurveEnableFlags = \
                self.ShowCurves_Middle(CD, YLabel, RowCount), [[], []], \
                [False, False]
            for Count, OOPP in enumerate(CD.GetCHP_OOPP()):
                if OOPP is not None:
                    FittingCurveEnableFlags[Count] = \
                        self.DictionaryParameter(OOPP, 'FittingCurveEnable')
                    if FittingCurveEnableFlags[Count]:
                        OnOffData[Count] = [
                            Func(CD, ScaleIndex, Count, j)
                            for j, TW in enumerate(self.TimeWindowList)]

            MTSD = self.PlotFittingCurves(
                Data=OnOffData, SubPlot=ax, TimeWindowList=TWL, Colors=Colors,
                FittingCurve=CD.GetCHP_FittingCurve(),
                PolyfitDegree=CD.GetCHP_PolyfitDegree(),
                MinMaxXFit=CD.GetCHP_MinMaxXFit(),
                OnOffBaseLineIntervals=CD.GetBaseLineIntervals(ScaleIndex),
                StdLatencyTimeDelay=CD.GetCHP_StdLatencyTimeDelay())
            MTSDL.append(MTSD)

        self.ShowCurves_End(fig, GraphicFileName, GraphicFileDPI)
        return MTSDL

    def ShowCurves_Init(self, TitleString):
        fig = plt.figure()
        fig.suptitle(TitleString, fontsize=20)
        fig.subplots_adjust(
            top=self.GetParameters_CurveAdjustment_Top(),
            left=self.GetParameters_CurveAdjustment_Left(),
            right=self.GetParameters_CurveAdjustment_Right(),
            bottom=self.GetParameters_CurveAdjustment_Bottom(),
            hspace=self.GetParameters_CurveAdjustment_HSpace(),
            wspace=self.GetParameters_CurveAdjustment_WSpace())
        return fig

    def ShowCurves_Middle(self, CD, YLabel, RowCount):
        ax = plt.subplot2grid((len(self.ChannelDataList), 1), (RowCount, 0))
        Title = F'{CD.GetCHP_Channel()}.{CD.GetCHP_Cluster()} On-Off'
        ax.set_title(Title, fontsize=14)
        ax.set_xlabel('Time Windows (ms)', fontsize=14)
        ax.set_ylabel(YLabel, fontsize=14)
        return ax

    def ShowCurves_End(self, fig, GraphicFileName, GraphicFileDPI):
        if self.GetParameters_CW_Maximized():
            plt.get_current_fig_manager().window.showMaximized()
        else:
            fig.set_size_inches(
                self.GetParameters_CW_Width(),
                self.GetParameters_CW_Height(), forward=True)
        plt.show()
        self.SaveFigure(fig, GraphicFileName, GraphicFileDPI, 'fitting curve')

    def PlotFittingCurves(self, **KWArgs):

        def GetBaseLine(X, Y, BLI):
            return [y for (x, y) in zip(X, Y) if x >= BLI[0] and x <= BLI[1]]

        ax, OnOffData, TWL, Colors, FittingCurve, PolyfitDegree, MinMaxXFit, \
            StdLatencyTimeDelay, BaseLineIntervals = KWArgs['SubPlot'], \
            KWArgs['Data'], KWArgs['TimeWindowList'], KWArgs['Colors'], \
            KWArgs['FittingCurve'], KWArgs['PolyfitDegree'], \
            KWArgs['MinMaxXFit'], KWArgs['StdLatencyTimeDelay'], \
            KWArgs['OnOffBaseLineIntervals']

        # show x-ticks
        XD, BaseLines, OutputTuple = self.GetXTickPositions(TWL), \
            [None, None], [None, None]
        XPositions, MinX, MaxX = XD['XPositions'], XD['MinX'], XD['MaxX']
        ax.set_xticks(XD['XTickPositions'])
        for Count, (BLI, OOD) in enumerate(zip(BaseLineIntervals, OnOffData)):
            if BLI is not None and OOD is not None:
                BaseLines[Count] = GetBaseLine(XD['XPositions'], OOD, BLI)

        # show curves
        for Count, (OOD, Color, BL) in enumerate(zip(
                OnOffData, Colors, BaseLines)):
            if BL is not None:
                OutputTuple[Count] = self.PlotFittingCurveFunc(
                    ax, Color, XPositions, OOD, MinX, MaxX, FittingCurve,
                    PolyfitDegree, self.GetParameters_LSPN(), MinMaxXFit, BL,
                    StdLatencyTimeDelay)
        return OutputTuple

    def PlotFittingCurveFunc(
            self, ax, Color, X, Y, MinX, MaxX, FittingCurve, PolyfitDegree,
            LinspacePointNumber, MinMaxXFit, BaseLine,
            StdLatencyTimeDelay):
        ax.plot(X, Y, 'X', color=Color)
        Flag, XX, YY = self.CutXYData(X, Y, MinX, MaxX, MinMaxXFit)
        STD_Threshold, BaseLineY, ReturValue = None, None, None
        if BaseLine is not None:
            STD_Threshold = np.std(BaseLine) * self.GetParameters_STC()
            BaseLineY = np.mean(BaseLine)
            STD_Threshold += BaseLineY
        if Flag and XX.size and YY.size:
            if FittingCurve == 'polyfit':
                ReturValue = self.PlotFittingCurve_Polyfit(
                    ax, Color, X, Y, XX, YY, MinX, MaxX, LinspacePointNumber,
                    PolyfitDegree, BaseLineY, STD_Threshold,
                    StdLatencyTimeDelay)
            if FittingCurve == 'gaussian':
                ReturValue = self.PlotFittingCurve_Gaussian(
                    ax, Color, X, Y, XX, YY, MinX, MaxX, LinspacePointNumber,
                    BaseLineY, STD_Threshold, StdLatencyTimeDelay)
            if FittingCurve == 'two_exponential':
                ReturValue = self.PlotFittingCurve_TwoExponential(
                    ax, Color, X, Y, XX, YY, MinX, MaxX, LinspacePointNumber,
                    BaseLineY, STD_Threshold, StdLatencyTimeDelay)
        return ReturValue

    def PlotFittingCurve_Polyfit(
            self, ax, Color, X, Y, XX, YY, MinX, MaxX, LinspacePointNumber,
            PolyfitDegree, BaseLineY, STD_Threshold, StdLatencyTimeDelay):
        P = np.polyfit(XX, YY, PolyfitDegree)
        DerminationCoefficient = r2_score(Y, np.polyval(P, X))
        XX = ar(np.linspace(MinX, MaxX, LinspacePointNumber))
        YY = np.polyval(P, XX)
        PeakParameterDict = self.PlotMaximumAndMedianLines(
            ax, Color, XX, YY, BaseLineY, STD_Threshold, StdLatencyTimeDelay)
        PeakParameterDict['R2'] = '{:0.2f}'.format(
            100 * DerminationCoefficient)
        return PeakParameterDict

    def PlotFittingCurve_Gaussian(
            self, ax, Color, X, Y, XX, YY, MinX, MaxX, LinspacePointNumber,
            BaseLineY, STD_Threshold, StdLatencyTimeDelay):
        Mean = sum(XX * YY) / sum(YY)
        Sigma = sum(YY * (XX - Mean) ** 2) / sum(YY)
        if Sigma != 0:
            P0, Fit, Success = [max(YY), Mean, Sigma], None, None
            try:
                Fit, Success = curve_fit(self.Gaussian, XX, YY, p0=P0)
            except RuntimeError:
                self.PrintError(('Error!!! Optimal parameters are not found.'))
            if Success is not None:
                DerminationCoefficient = r2_score(Y, self.Gaussian(X, *Fit))
                XX = ar(np.linspace(MinX, MaxX, LinspacePointNumber))
                YY = self.Gaussian(XX, *Fit)
                PeakParameterDict = self.PlotMaximumAndMedianLines(
                    ax, Color, XX, YY, BaseLineY, STD_Threshold,
                    StdLatencyTimeDelay)
                PeakParameterDict['R2'] = '{:0.2f}'.format(
                    100 * DerminationCoefficient)
                return PeakParameterDict
            else:
                self.PrintError(('Error!!! Division by zero. Sigma == 0.'))
        return None

    def PlotFittingCurve_TwoExponential(
            self, ax, Color, X, Y, XX, YY, MinX, MaxX, LinspacePointNumber,
            BaseLineY, STD_Threshold, StdLatencyTimeDelay):
        P0, Fit, Success = [max(YY), 55, 0.1, 0.1], None, None
        try:
            Fit, Success = curve_fit(self.TwoExponential, XX, YY, p=P0)
        except RuntimeError:
            self.PrintError(('Error: Optimal parameters are not found.'))
        if Success is not None:
            DerminationCoefficient = r2_score(Y, self.TwoExponential(X, *Fit))
            XX = ar(np.l/inspace(MinX, MaxX, LinspacePointNumber))
            YY = self.TwoExponential(XX, *Fit)
            PeakParameterDict = self.PlotMaximumAndMedianLines(
                ax, Color, XX, YY, BaseLineY, STD_Threshold,
                StdLatencyTimeDelay)
            PeakParameterDict['R2'] = F'{100 * DerminationCoefficient:0.2f}'
            return PeakParameterDict
        return None

    def PlotMaximumAndMedianLines(
            self, ax, Color, XX, YY, BaseLineY, STD_Threshold,
            StdLatencyTimeDelay):

        def DrawLine(verts):
            path = Path(verts, [Path.MOVETO, Path.LINETO])
            ax.add_patch(patches.PathPatch(path, color=Color, lw=2))

        ax.plot(XX, YY, color=Color)
        ax.axhline(color='black', lw=1)
        MaxY = max(YY)
        XM = [
            X for i, (X, Y) in enumerate(zip(XX, YY)) if Y == MaxY][0]
        DrawLine([(XM, 0), (XM, MaxY)])
        SL = MaxY / 2 if BaseLineY is None else\
            BaseLineY + (MaxY - BaseLineY) / 2
        XL, X1, X2 = [
            X for i, (X, Y) in enumerate(zip(XX, YY))
            if Y >= SL], None, None
        if XL:
            X1, X2 = XL[0], XL[-1]
            DrawLine([(XL[0], SL), (XL[-1], SL)])
        Dict = {
            'XPeak': XM, 'T0': X1, 'Sigma': X2 - X1 if X1 and X2 else None,
            'StdThresholdTime': None}
        if STD_Threshold:
            XT = [
                X for i, (X, Y) in enumerate(zip(XX, YY))
                if Y >= STD_Threshold and (
                    True if StdLatencyTimeDelay is None
                    else X >= StdLatencyTimeDelay)]
            if XT:
                Dict['StdThresholdTime'] = XT[0]
                ax.axvline(x=XT[0], color=Color, ls='dashed', lw=2)
        return Dict

    def ShowMaximumAndPowerMTSD(self):
        self.ShowError_DCFR()
        raise NotImplementedError()

    def ShowPeakParameters(self, MTSDL, Title):
        self.Print((Title))
        ResponseList, VList, HeaderList = ['On', 'Off'], [], [
            'CH.CL', 'On/Off', 'Peak T', 'Sigma', 'T0', 'R2 %', 'Latency']
        for CD, MTSD in zip(self.ChannelDataList, MTSDL):
            for Count, (R, R_MTSD) in enumerate(zip(ResponseList, MTSD)):
                if R_MTSD is not None:
                    VList.append([
                        '{0}.{1}'.format(
                            CD.GetCHP_Channel(), CD.GetCHP_Cluster()),
                        R, self.FMT(R_MTSD['XPeak']),
                        self.FMT(R_MTSD['Sigma']), self.FMT(R_MTSD['T0']),
                        self.FMT(R_MTSD['R2']),
                        self.FMT(R_MTSD['StdThresholdTime'])])
        T = tabulate.tabulate(VList, headers=HeaderList)
        self.Print((T, '\n'))

#
# Other
#

    def GetXTickPositions(self, TWL):
        XPositions = [TW[0] + (TW[1]-TW[0])/2 for TW in TWL]
        MinX, MaxX = min(XPositions), max(XPositions)
        Bin = self.GetParameters_FCHB()
        X1 = (MinX // Bin) * Bin
        if MinX - X1 < Bin / 2:
            X1 -= Bin
        X2 = (MaxX // Bin + 1) * Bin
        XTickPositions = [
            x for x in range(int(X1), int(X2+1))
            if x % Bin == 0]
        return {'XTickPositions': XTickPositions, 'XPositions': XPositions,
                'MinX': MinX, 'MaxX': MaxX}

    def CutXYData(self, X, Y, MinX, MaxX, MinMaxXFit):
        MinXFit, MaxXFit = MinMaxXFit
        ReturnValue = (False, None, None)
        if MinXFit is None and MaxXFit is None:
            ReturnValue = (True, ar(X), ar(Y))
        else:
            if MinXFit is None:
                MinXFit = MinX
            if MaxXFit is None:
                MaxXFit = MaxX
            if MinXFit < MaxXFit:
                XX = ar([xv for xv in X if xv >= MinXFit and xv <= MaxXFit])
                YY = ar([
                    yv for (xv, yv) in zip(X, Y)
                    if xv >= MinXFit and xv <= MaxXFit])
                if XX.size and YY.size:
                    ReturnValue = (True, XX, YY)
                else:
                    self.PrintError((
                        'Error!!! There is no data in MinX: MaxX interval.'))
        return ReturnValue

    def Gaussian(self, X, *P):
        Amp, Mean, Sigma = P
        return Amp * np.exp(-(X - Mean) ** 2 / (2 * Sigma ** 2))

    def TwoExponential(self, X, *P):
        # A - peak area
        # B - elution time
        # C - width of gaussian
        # D - exponential damping term
        A, B, C, D = P
        return (
            A / 2 / D * np.exp(C ** 2 / 2.0 / D ** 2 + (B - X) / D) *
            (erf((X - B) / (np.sqrt(2.0) * C) - C / np.sqrt(2.0) / D) + 1.0))
