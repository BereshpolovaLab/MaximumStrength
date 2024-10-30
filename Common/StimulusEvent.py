#!/usr/bin/env python3
# 24, BereshpolovaLab, University of Connecticut


import numpy as np


###########################################################
# CStimulusEvent
###########################################################

class CStimulusEvent():

    def __init__(self, TimeStamp, BurstWordList):
        self.TimeStamp, BWLLLength = TimeStamp, len(BurstWordList)

        # 1 - word: background index, color
        self.BackgroundIndex, self.Color = \
            BurstWordList[0] >> 8 & 0x03, BurstWordList[0] & 0xFF
        # 2 - word: X-position of the grid
        self.XPosition = self.ExtractPosition(BurstWordList[1])
        # 3 - word: Y-position of the grid
        self.YPosition = self.ExtractPosition(BurstWordList[2])
        # 4 - word: horizontal scale and stimulus duration
        self.HorizontalScale, self.StimulusDuration = \
            self.ExtractScale(BurstWordList[3]), BurstWordList[3] & 0xFF
        # 5 - word: vertical scale and rest duration
        self.VerticalScale, self.RestDuration = \
            self.ExtractScale(BurstWordList[4]), BurstWordList[4] & 0xFF
        # 6 - word: andgle, I, J
        self.Angle = (BurstWordList[5] >> 10 & 0x0F) * 15
        self.GridI, self.GridJ = BurstWordList[5] >> 5 & 0x1F, \
            BurstWordList[5] & 0x1F
        # Protocol N11
        if BWLLLength == 7:
            self.CustomizedBackgroundColor, self.SquareSize = \
                BurstWordList[6] >> 6 & 0xFF, BurstWordList[6] & 0x3F
        # Protocol N0
        if BWLLLength == 6:
            self.CustomizedBackgroundColor = None
            SquareSizeLowPart = BurstWordList[1] >> 12 & 0x07
            SquareSizeHiPart = BurstWordList[2] >> 9 & 0x38
            self.SquareSize = SquareSizeHiPart | SquareSizeLowPart
        self.BackgroundColor = {0: 0, 1: 127, 2: 255}.get(
            self.BackgroundIndex, self.CustomizedBackgroundColor)

#
# Get - functions
#

    def GetTimeStamp(self): return self.TimeStamp

    def GetYPosition(self): return self.YPosition

    def GetHorizontalScale(self): return self.HorizontalScale

    def GetVerticalScale(self): return self.VerticalScale

    def GetStimulusDuration(self): return self.StimulusDuration

    def GetRestDuration(self): return self.RestDuration

    def GetGridI(self): return self.GridI

    def GetGridJ(self): return self.GridJ

    def GetSquareSize(self): return self.SquareSize

#
# Set - functions
#

    def SetGridI(self, GridI): self.GridI, self.GridJ = GridI, 0

#
# Service - functions
#

    def CorrectTimeStamp(self, ScreenHeight, VRR):
        TimeDelta = 1000.0 * (
            self.GetYPosition() + (0.5 + self.GetGridI()) *
            self.GetSquareSize()) / (ScreenHeight * VRR)
        self.TimeStamp += TimeDelta

    def ExtractPosition(self, Word):
        Position = np.int16(Word & 0xFFF)
        if Position & 0x800:
            Position = np.int16(Position | 0xF000)
        return int(Position)

    def ExtractScale(self, Word):
        W = Word >> 8 & 0x3F
        return {0x21: 0.1, 0x22: 0.25, 0x23: 0.5}.get(W, float(W))

    def GetOnOffCount(self):
        OnOffCount = None
        if self.BackgroundColor:
            if self.Color > self.BackgroundColor:
                OnOffCount = 0
            if self.Color < self.BackgroundColor:
                OnOffCount = 1
        return OnOffCount

    def IsStimulusPhase(self):
        return True if self.Color != self.BackgroundColor else False

    def IsRestPhase(self): return ~self.IsStimulusPhase()
