#!/usr/bin/env python3
# 2024, BereshpolovaLab, University of Connecticut
# Developed by Victor Serdyukov, e-mails: Victor.Serdyukov@uconn.edu


import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import tabulate

from matplotlib.path import Path
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size

from Common import ErrorObject as C_EO
from Common import Filters as C_F
from Common import PlexonHeaders as C_PH
from Common import StimulusEvent as C_SE


###########################################################
# CBasicPlexonData
###########################################################

class CBasicPlexonData(C_EO.CErrorObject):
    PL_SingleWFType, PL_ExtEventType, PL_StrobedExtChannel, \
        MAX_CHANNEL_NUMBER, MAX_CLUSTER_NUMBER, LogFile, ResponseStringTuple, \
        GridSize, ScreenHeight = 1, 4, 257, 32, 4, None, ('On', 'Off'), \
        (30, 22), 768

    def __init__(self, Parameters, LogFile):
        super().__init__(LogFile)

        self.StatusFlag, self.Parameters, self.TTLList, self.StimuliList, \
            self.ChannelDataList, self.OnOffStimuliList = False, Parameters, \
            [], [], [], [[], []]

    def __del__(self):
        ST, FT = ('Well done, my friend!'), \
            ('Something is wrong here. This is not good work.')
        self.Print(ST if self.StatusFlag else FT)

#
# GetParameters - functions
#

    def GetParameters_ProtocolNumber(self):
        return self.DictionaryParameter(self.Parameters, 'ProtocolNumber')

    def GetParameters_VRR(self):
        return self.DictionaryParameter(self.Parameters, 'VRR')

    def GetParameters_VRR_ADF(self):
        return self.DictionaryParameter(self.Parameters, 'VRR_AutoDetectFlag')

    def GetParameters_FileIntervalPoint(self, Count):
        return self.DictionaryParameter(self.Parameters, 'FileInterval')[Count]

    def GetParameters_Filter(self):
        return self.DictionaryParameter(self.Parameters, 'Filter')

    def GetParameters_GaussianSigma(self):
        return self.DictionaryParameter(self.Parameters, 'GaussianSigma')

    def GetParameters_Interpolation(self):
        return self.DictionaryParameter(self.Parameters, 'Interpolation')

    def GetParameters_Channels(self):
        return self.DictionaryParameter(self.Parameters, 'Channels')

    def GetParameters_PlotGridLineVisible(self):
        return self.DictionaryParameter(self.Parameters, 'PlotGridLineVisible')

    def GetParameters_PlotGridLineColor(self):
        return self.DictionaryParameter(self.Parameters, 'PlotGridLineColor')

    def GetParameters_PlotGridLineStyle(self):
        return self.DictionaryParameter(self.Parameters, 'PlotGridLineStyle')

#
# Get - functions
#

    def GetStatusFlag(self): return self.StatusFlag

#
# Set - functions
#

    def SetStatusFlag(self, StatusFlag): self.StatusFlag = StatusFlag

#
# Print - functions
#

    def PrintError(self, List):
        self.Print(List)
        self.StatusFlag = False

    def PrintSpikeNumbers(self):
        HeaderList, Data2DList = ['CH.CL', 'Numeber of spikes'], [[
            F'{CD.GetCHP_Channel()}.{CD.GetCHP_Cluster()}',
            CD.GetSpikeListLength()] for CD in self.ChannelDataList]
        Tabulate = tabulate.tabulate(Data2DList, headers=HeaderList)
        self.Print((F'{Tabulate}\n\n'))

#
# Analysis - functions
#

    def LoadDataFromPLXFile(self, PLXFileName, TTLFlag=True):
        f = None
        try:
            f = open(PLXFileName, 'rb')
        except IOError:
            self.PrintError((
                'Error!!! Cannot open PLX-file! Something is wrong.'))
        finally:
            if f:
                Data = f.read()
                f.close()
                FileHeader = C_PH.CPL_FileHeader(Data)
                N1 = FileHeader.GetSize()
                if FileHeader.GetNumDSPCH() != 0:
                    ChanHeader = C_PH.CPL_ChanHeader()
                    N1 += FileHeader.GetNumDSPCH() * ChanHeader.GetSize()
                if FileHeader.GetNumEventCH() != 0:
                    EventHeader = C_PH.CPL_EventHeader()
                    N1 += FileHeader.GetNumEventCH() * EventHeader.GetSize()
                if FileHeader.GetNumSlowCH() != 0:
                    SlowChannelHeader = C_PH.CPL_SlowChannelHeader()
                    N1 += FileHeader.GetNumSlowCH() * \
                        SlowChannelHeader.GetSize()
                DBH, T1, T2 = C_PH.CPL_DataBlockHeader(), \
                    self.GetParameters_FileIntervalPoint(0), \
                    self.GetParameters_FileIntervalPoint(1)
                while N1 > 0:
                    N1 = DBH.ExtractData(Data, N1)
                    if N1 > 0:
                        TimeStamp = DBH.GetTimeStamp() * \
                            FileHeader.GetTimeStampPeriod()
                        if (T1 is None or
                                (T1 is not None and TimeStamp >= T1 * 1000)) \
                                and (T2 is None or (T2 is not None and
                                     TimeStamp <= T2 * 1000)):
                            # spike-events
                            if DBH.GetType() == self.PL_SingleWFType and \
                                    DBH.GetChannel() >= 1 and \
                                    DBH.GetChannel() <= \
                                    self.MAX_CHANNEL_NUMBER and \
                                    DBH.GetUnit() >= 1 and DBH.GetUnit() <= \
                                    self.MAX_CLUSTER_NUMBER and \
                                    self.ChannelDataList:
                                for CD in self.ChannelDataList:
                                    CD.AddSpike(
                                       DBH.GetChannel(), DBH.GetUnit(),
                                       TimeStamp)
                            # TTL-events
                            elif TTLFlag and \
                                    DBH.GetType() == self.PL_ExtEventType and \
                                    DBH.GetChannel() == \
                                    self.PL_StrobedExtChannel and \
                                    (DBH.GetUnit() & 0x8000) >> 15 == \
                                    self.Parameters['PlexonDigitalCard']:
                                TTLWord = np.int16(DBH.Unit & 0x7FFF)
                                self.TTLList.append([TimeStamp, TTLWord])
                for CD in self.ChannelDataList:
                    CD.SortSpikes()
                if (TTLFlag and self.TTLList) or not TTLFlag:
                    self.StatusFlag = True

    def ConvertTTLEventsToStimuli(self):
        if self.GetParameters_ProtocolNumber() == 11:
            FBPNPattern, BurstNumber = int(0x6C00), 7
            for Count, GTE in enumerate(self.TTLList):
                if (int(GTE[1]) & int(0xFC00)) == FBPNPattern and \
                        Count <= len(self.TTLList) - BurstNumber:
                    self.AddToStimuliListAndOnOffStimuliList(
                        Count, BurstNumber, GTE[0])
        elif self.GetParameters_ProtocolNumber() == 0:
            BurstNumber = 6
            for Count in range(0, len(self.TTLList), 6):
                if Count <= len(self.TTLList) - BurstNumber:
                    self.AddToStimuliListAndOnOffStimuliList(
                        Count, BurstNumber, self.TTLList[Count][0])

    def AddToStimuliListAndOnOffStimuliList(self, Count, BurstNumber, T):
        BurstWordList = [
            LTE[1] for LTE in self.TTLList[Count:Count + BurstNumber]]
        SE = C_SE.CStimulusEvent(T, BurstWordList)
        if SE.IsStimulusPhase():
            self.StimuliList.append(SE)
            OnOffCount = SE.GetOnOffCount()
            if OnOffCount is not None and SE.GetGridJ() >= 0 and \
                    SE.GetGridJ() < self.GridSize[0] and \
                    SE.GetGridI() >= 0 and SE.GetGridI() < self.GridSize[1]:
                self.OnOffStimuliList[OnOffCount].append(SE)

    def StimulusTimeStampCorrection(self, OOSL=None, PrintFlag=True):
        if not OOSL:
            OOSL = self.OnOffStimuliList
        VRR = self.GetParameters_VRR()
        str = F'Vertical Refresh Rate = {VRR}'
        if self.GetParameters_VRR_ADF():
            str = None
            if len(self.StimuliList) < 2:
                str = 'VRR autodetect flag is on ' \
                    'but number of stimuli is 1. \n' \
                    'Preset VRR value will be used for timestamp correction'
            else:
                VRR = self.IdentifyVRR()
                str = F'AutoDetect Vertical Refresh Rate = {VRR}'
        if PrintFlag:
            self.Print((F'{str}\n'))

        for StimuliList in OOSL:
            for SE in StimuliList:
                SE.CorrectTimeStamp(self.ScreenHeight, VRR)

    def IdentifyVRR(self):
        PeriodList = [
            (SE2.GetTimeStamp() - SE1.GetTimeStamp()) /
            (SE1.GetStimulusDuration() + SE1.GetRestDuration())
            for SE1, SE2 in zip(self.StimuliList[:-1], self.StimuliList[1:])]
        Frequency = round(1000.0 / C_F.MedianFilter(PeriodList))
        # 100 or 160 Hz
        Frequency = 160.0 if Frequency >= 130.0 else 100.0
        return Frequency

    def LoadData_CheckSpikesAndTTLEvents(self, PLXFileName, TTLFlag=True):
        self.LoadDataFromPLXFile(PLXFileName, TTLFlag)
        if self.StatusFlag:
            for CD in self.ChannelDataList:
                if not CD.GetStatus():
                    self.PrintError((
                        'Error!!! There are no spikes in file for CH.CL: '
                        F'{CD.GetCHP_Channel()}.{CD.GetCHP_Cluster()}'))
                    self.StatusFlag = False
            if TTLFlag:
                if self.StatusFlag and not self.TTLList:
                    self.PrintError((
                        'Error!!! There are no TTL-events in file.'))
                    self.StatusFlag = False
                if self.StatusFlag:
                    self.ConvertTTLEventsToStimuli()

    def CheckStimuli(self):
        if self.StatusFlag:
            for StimuliList, ResponseString in zip(
                    self.OnOffStimuliList, self.ResponseStringTuple):
                if not StimuliList:
                    self.PrintError((
                        F'There are no {ResponseString}-stimuli'
                        ' in this file. '
                        'It is possible that something is wrong.'))
                    self.StatusFlag = False

#
# Draw - functions
#

    def PresentPlot(
            self, CD, PD, aspect, pad_fraction, SubplotGridSizeTuple, IJTuple,
            Maximum=None):
        ax = plt.subplot2grid(SubplotGridSizeTuple, IJTuple)
        if Maximum is None:
            Maximum = PD.GetMaximum()
        im = ax.imshow(
            PD.GetFilteredThresholdData2DNPArray(),
            interpolation=self.GetParameters_Interpolation(), cmap='jet',
            vmin=0, vmax=Maximum)
        self.DrawTicksTitleGridLines(
            ax, im, aspect, pad_fraction, CD, PD)
        return ax

    def DrawTicksTitleGridLines(self, ax, im, aspect, pad_fraction, CD, PD):

        def GetXYAndLabels(Index):
            V = [V - 0.5 for V in range(0, self.GridSize[Index], 5)]
            VLabels = [V for V in range(0, self.GridSize[Index], 5)]
            return V, VLabels

        X, XLabels = GetXYAndLabels(0)
        Y, YLabels = GetXYAndLabels(1)
        plt.xticks(X, XLabels)
        plt.yticks(Y, YLabels)

        divider, width = make_axes_locatable(ax), \
            axes_size.AxesY(ax, aspect=1/aspect)
        pad = axes_size.Fraction(pad_fraction, width)
        cax = divider.append_axes('right', size=width, pad=pad)
        plt.colorbar(im, cax=cax)

        TitleString = (
            F'{CD.GetCHP_Channel()}.{CD.GetCHP_Cluster()}'
            F'-{PD.GetPlotParameters_Response()}')
        TimeWindow = PD.GetTimeWindow()
        if TimeWindow is not None:
            TitleString += F': {TimeWindow[0]} - {TimeWindow[1]} ms'
        TitleString += F'\nSpikes: {CD.GetSpikeListLength()}'
        ax.set_title(TitleString)
        if self.GetParameters_PlotGridLineVisible():
            self.DrawGridLines(
                ax, self.GetParameters_PlotGridLineColor(),
                self.GetParameters_PlotGridLineStyle(), Delta=0.5)

    def DrawGridLines(self, ax, LineColor, LineStyle, Delta=0):
        codes, VertList = [Path.MOVETO, Path.LINETO], []
        for X in range(self.GridSize[0]):
            VertList.append(
                [(X+Delta, -Delta), (X+Delta, self.GridSize[1]-Delta)])
        for Y in range(self.GridSize[1]):
            VertList.append(
                [(-Delta, Y+Delta), (self.GridSize[0]-Delta, Y+Delta)])
        for Vert in VertList:
            path = Path(Vert, codes)
            ax.add_patch(
                patches.PathPatch(path, color=LineColor, linestyle=LineStyle))

#
# Other
#

    def ShowAndSave(self, fig, FileName, Parameter, PE):
        plt.get_current_fig_manager().window.showMaximized()
        plt.show()
        self.SaveFigure(fig, FileName, self.Parameters[Parameter], PE)

    def SaveFigure(self, fig, GraphicFileName, GraphicFileDPI, TS):
        if GraphicFileName:
            try:
                fig.savefig(
                    GraphicFileName,
                    dpi=GraphicFileDPI if GraphicFileDPI else fig.dpi)
            except IOError:
                self.Print((
                    F'Cannot write the {TS} graphics-file! '
                    'Something is wrong.'))
