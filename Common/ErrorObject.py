#!/usr/bin/env python3
# 2024, BereshpolovaLab, University of Connecticut
# Developed by Victor Serdyukov, e-mails: Victor.Serdyukov@uconn.edu


import inspect
import sys

from Common import PrintObject as C_PO


###########################################################
# CErrorObject
###########################################################

class CErrorObject(C_PO.CPrintObject):

    def __init__(self, LogFile):
        super().__init__(LogFile)

#
# Get - functions
#

    def GetModuleAndFunctionNames(self, FL, SL):
        return sys._getframe(FL).f_globals['__name__'], \
            inspect.stack()[SL][0].f_code.co_name

#
# Error - functions
#

    def ShowError_DCFR(self, FL=2, SL=2):
        ModuleName, FuncName = self.GetModuleAndFunctionNames(FL, SL)
        self.Print((
            F'Errror!!! Module: {ModuleName}. '
            F'Function {FuncName} is not realized in derived class.\n'))

    def ShowError_NoneIndex(self, FL=2, SL=2):
        ModuleName, FuncName = self.GetModuleAndFunctionNames(FL, SL)
        self.Print((
            F'Errror!!! Module: {ModuleName}. '
            F'Function: {FuncName}. Index is None.\n'))

    def CheckDictionaryKeyAndPrintError(self, Dictionary, Key, FL=3, SL=3):
        if Key not in Dictionary:
            ModuleName, FuncName = self.GetModuleAndFunctionNames(FL, SL)
            self.Print((
                F'Errror!!! Module: {ModuleName}. Function: {FuncName}. '
                F'The key \'{Key}\' does not exist in dictionary.\n'))

    def CheckDictionaryKeyInDictionaryAndPrintError(
            self, Dictionary, Key1, Key2):
        if Key1 in Dictionary:
            self.CheckDictionaryKeyAndPrintError(Dictionary[Key1], Key2)

    def DictionaryParameter(self, Dictionary, Key):
        self.CheckDictionaryKeyAndPrintError(Dictionary, Key)
        return Dictionary[Key]

    def DictInDictParameter(self, Dictionary, Key1, Key2):
        self.CheckDictionaryKeyAndPrintError(Dictionary, Key1)
        self.CheckDictionaryKeyInDictionaryAndPrintError(
            Dictionary, Key1, Key2)
        return Dictionary[Key1][Key2]
