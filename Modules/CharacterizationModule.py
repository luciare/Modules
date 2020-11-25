# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 08:17:10 2020

@author: Lucia
"""

import numpy as np
import pickle
import datetime
from scipy.signal import welch

import pyqtgraph.parametertree.parameterTypes as pTypes

from PyQt5 import Qt
from PyQt5.QtCore import QObject
from PyQt5.QtWidgets import QFileDialog

import PyGFETdb.DataStructures as PyData
import PyGFETdb.PlotDataClass as PyFETpl
import PyqtTools.PlotModule as PltBuffer2D

################PARAMETER TREE################################################

ConfigSweepsParams = {'name': 'SweepsConfig',
                      'type': 'group',
                    'children': ({'name': 'Start/Stop Sweep',
                                  # 'title': 'Start Sweep',
                                  'type': 'action', },
                                 {'name': 'Pause/ReStart Sweep',
                                  # 'title': 'Start Sweep',
                                  'type': 'action', },
                                 {'name': 'ACenable',
                                  'title': 'AC Characterization',
                                  'type': 'bool',
                                  'value': True, },
                                 {'name': 'StabCriteria',
                                  'type': 'list',
                                  'values': ['All channels', 'One Channel', 'Mean'],
                                  'value': 'Mean',
                                  'visible': True},
                                 {'name': 'VgSweep',
                                          'type': 'group',
                                          'children': ({'name': 'Vinit',
                                                        'type': 'float',
                                                        'value': 0,
                                                        'siPrefix': True,
                                                        'suffix': 'V'},
                                                       {'name': 'Vfinal',
                                                        'type': 'float',
                                                        'value': -0.4,
                                                        'siPrefix': True,
                                                        'suffix': 'V'},
                                                       {'name': 'NSweeps',
                                                        'type': 'int',
                                                        'value': 1,
                                                        'siPrefix': True,
                                                        'suffix': 'Sweeps'},
                                                       )},
                                 {'name': 'VdSweep',
                                  'type': 'group',
                                  'children': ({'name': 'Vinit',
                                                'type': 'float',
                                                'value': 0.02,
                                                'siPrefix': True,
                                                'suffix': 'VRMS'},
                                               {'name': 'Vfinal',
                                                'type': 'float',
                                                'value': 0.2,
                                                'siPrefix': True,
                                                'suffix': 'VRMS'},
                                               {'name': 'NSweeps',
                                                'type': 'int',
                                                'value': 1,
                                                'siPrefix': True,
                                                'suffix': 'Sweeps'},
                                               )},
                                 {'name': 'MaxSlope',
                                  'title': 'Maximum Slope',
                                  'type': 'float',
                                  'value': 1e-10,
                                  'siPrefix': True,
                                  'suffix': 'A/s'},
                                 {'name': 'TimeOut',
                                  'title': 'Max Time for Stabilization',
                                  'type': 'int',
                                  'value': 10,
                                  'siPrefix': True,
                                  'suffix': 's'},
                                 {'name': 'TimeBuffer',
                                  'title': 'Buffer Time for Stabilization',
                                  'type': 'int',
                                  'value': 1,
                                  'siPrefix': True,
                                  'suffix': 's'},
                                 {'name': 'DelayTime',
                                  'title': 'Time to wait for acquisition',
                                  'type': 'int',
                                  'value': 1,
                                  'siPrefix': True,
                                  'suffix': 's'},
                                 )}

SaveSweepsParams = ({'name': 'SaveSweepConfig',
                    'type': 'group',
                    'children': ({'name': 'Save File',
                                  'type': 'action'},
                                 {'name': 'Folder',
                                  'type': 'str',
                                  'value': ''},
                                 {'name': 'Oblea',
                                  'type': 'str',
                                  'value': ''},
                                 {'name': 'Disp',
                                  'type': 'str',
                                  'value': ''},
                                 {'name': 'Name',
                                  'type': 'str',
                                  'value': ''},
                                 {'name': 'Cycle',
                                  'type': 'int',
                                  'value': 0},
                                 )
                    })

class SweepsConfig(pTypes.GroupParameter):
    def __init__(self, QTparent, **kwargs):
        pTypes.GroupParameter.__init__(self, **kwargs)
        self.QTparent = QTparent
        
        self.addChild(ConfigSweepsParams)
        self.SwConfig = self.param('SweepsConfig')

        self.VgParams = self.SwConfig.param('VgSweep')
        self.VdParams = self.SwConfig.param('VdSweep')

        self.VgParams.sigTreeStateChanged.connect(self.on_Sweeps_Changed)
        self.VdParams.sigTreeStateChanged.connect(self.on_Sweeps_Changed)
        self.on_Sweeps_Changed()
        
        self.addChild(SaveSweepsParams)
        self.SvSwParams = self.param('SaveSweepConfig')
        self.SvSwParams.param('Save File').sigActivated.connect(self.FileDialog)

    def on_Sweeps_Changed(self):
        '''
        Creates the numpy arrays used for Vgs and Vds sweeps during the
        characterization process
        '''
        self.VgSweepVals = np.linspace(self.VgParams.param('Vinit').value(),
                                       self.VgParams.param('Vfinal').value(),
                                       self.VgParams.param('NSweeps').value())

        self.VdSweepVals = np.linspace(self.VdParams.param('Vinit').value(),
                                       self.VdParams.param('Vfinal').value(),
                                       self.VdParams.param('NSweeps').value())

    def FileDialog(self):
        RecordFile = QFileDialog.getExistingDirectory(self.QTparent,
                                                      "Select Directory",
                                                      )
        if RecordFile:
            self.SvSwParams.param('Folder').setValue(RecordFile)

    def FilePath(self):
        return self.param('Folder').value()
    
    def GetConfigSweepsParams(self):
        '''Returns parameters to do the sweeps
           SwConfig={'AEnable': True,
                     'StabCriteria': Mean,
                     'VgSweep': array([ 0. , -0.1, -0.2, -0.3]),
                     'VdSweep': array([0.1]),
                     'MaxSlope': 1e-10,
                     'TimeOut': 10,
                     'TimeBuffer': 1,
                     'DelayTime': 1,
                     }
        '''
        SwConfig = {}
        for Config in self.SwConfig.children():
            if Config.name() == 'Start/Stop Sweep':
                continue
            if Config.name() == 'Pause/ReStart Sweep':
                continue
            if Config.name() == 'VgSweep':
                SwConfig[Config.name()] = self.VgSweepVals
                continue
            if Config.name() == 'VdSweep':
                SwConfig[Config.name()] = self.VdSweepVals
                continue
            SwConfig[Config.name()] = Config.value()
        return SwConfig
    
    def GetSaveSweepsParams(self):
        '''Returns de parameters to save the caracterization
           Config={'Folder': 'C:/Users/Lucia/Dropbox (ICN2 AEMD - GAB GBIO)/
                              TeamFolderLMU/FreqMux/Lucia/DAQTests/SweepTests
                              /18_12_19',
                   'Oblea': 'testPyCont',
                   'Disp': 'Test',
                   'Name': 'Test',
                   'Cycle': 0
                   }
        '''
        Config = {}
        for Conf in self.SvSwParams.children():
            if Conf.name() == 'Save File':
                continue
            Config[Conf.name()] = Conf.value()

        return Config  

################CHARACTERIZATION THREAD#######################################
        
class StbDetThread(Qt.QThread):
    DataStab = Qt.pyqtSignal()
    NextVg = Qt.pyqtSignal()
    NextVd = Qt.pyqtSignal()
    CharactEnd = Qt.pyqtSignal()

    def __init__(self, ACenable, StabCriteria, VdSweep, VgSweep, MaxSlope, TimeOut, 
                 TimeBuffer, DelayTime, nChannels, ChnName, PlotterDemodKwargs, 
                 **kwargs):
        '''
        Thread used to do the characterization process of the transistors.

        Parameters
        ----------
        ACenable : bool
            Defines if AC characterization is done.
        StabCriteria : str
            Defines the criteria used to determine the stabilization of the 
            acquired data. It can be:
                - All channels: All channels must have a slope lower than 
                                MaxSlope
                - One Channel: it is enough if only one channel has a slope
                                lower than MaxSlope
                - Mean: the arithmetic mean of the slopes of all the channels
                        must be lower than MaxSlope
        VdSweep : numpy array
            Contains the values used for Vds sweep
        VgSweep : numpy array
            Contains the values used for Vgs sweep
        MaxSlope : float
            The maximum slope permited to consider the data stable.
        TimeOut : int
            Maximum time that the system waits for the stabilization of the 
            acquired data. After TimeOut seconds, the data is processed even
            though the data has not reached the MaxSlope value
        TimeBuffer : int
            Specifies the time laps to acquire the data (buffersize) 
        DelayTime : int
            Indicates an amount of time, at the begining of the recording,
            during which data is not acquired
        nChannels : int
            Number of channels enabled
        ChnName : Dictionary
            Contains the names from each  channel and its index
            {'Ch04': 0, 'Ch05': 1, 'Ch06': 2}
        PlotterDemodKwargs : dictionary
            Ploter kwargs
        **kwargs : kwargs
            Parameters received but not used

        Returns
        -------
        None.

        '''
        
        super(StbDetThread, self).__init__()
        # Init threads and flags
        self.threadCalcPSD = None
        self.ToStabData = None
        self.Wait = True
        self.Stable = False
        self.StabTimeOut = False
        # Define global variables
        self.ACenable = ACenable
        self.StabCriteria = StabCriteria
        self.MaxSlope = MaxSlope
        self.TimeOut = TimeOut
        self.DelayTime = DelayTime
        self.ElapsedTime = 0
        self.FsDemod = PlotterDemodKwargs['Fs']
        # Define global variables for Vg and Vd sweep
        self.VgIndex = 0
        self.VdIndex = 0
        self.VgSweepVals = VgSweep
        self.VdSweepVals = VdSweep
        self.NextVgs = self.VgSweepVals[self.VgIndex]
        self.NextVds =self.VdSweepVals[self.VdIndex]
        # init the timer
        self.Timer = Qt.QTimer()
        # Define the buffer size
        self.Buffer = PltBuffer2D.Buffer2D(self.FsDemod,
                                           nChannels,
                                           TimeBuffer)
        #Define DC and AC dictionaries
        self.SaveDCAC = SaveDicts(ACenable=self.ACenable,
                                  SwVdsVals=VdSweep,
                                  SwVgsVals=VgSweep,
                                  Channels=ChnName,
                                  nFFT=int(PlotterDemodKwargs['nFFT']),
                                  FsDemod=self.FsDemod
                                  )
        # init the PSD characterization plot
        if self.ACenable:
            self.PSDPlotVars = ('PSD',)
            self.threadCalcPSD = CalcPSD(**PlotterDemodKwargs)
            self.threadCalcPSD.PSDDone.connect(self.on_PSDDone)
            self.SaveDCAC.PSDSaved.connect(self.on_NextVgs)
            self.PlotSwAC = PyFETpl.PyFETPlot()
            self.PlotSwAC.AddAxes(self.PSDPlotVars)   
            
        else:
            self.SaveDCAC.DCSaved.connect(self.on_NextVgs)
        # Define the DC characterization plots   
        self.DCPlotVars = ('Ids', 'Rds', 'Gm', 'Ig')
        self.PlotSwDC = PyFETpl.PyFETPlot()
        self.PlotSwDC.AddAxes(self.DCPlotVars)

    def run(self):
        '''
        This function is executed after thread.run() is done and is the main
        function of the thread. It is executed constantly until the thread
        is stopped
        '''
        while True:
            # Waits until the buffer is filled with the acquired data
            if self.Buffer.IsFilled():
                # The function CalcSlope is called to determine the stability
                self.CalcSlope()
                if self.Stable:
                    # If the data is stable, a signal is emited
                    self.DataStab.emit()
                    # And the Timer is stop as TimeOut is not necessary
                    self.Timer.stop()
                    self.Timer.deleteLater()
                    print('IsStable')
                    # if AC characterization is wanted
                    if self.ACenable:
                        # The PSD of the data is executed
                        self.threadCalcPSD.start()
                    # And the DC data is saved before the sweep continuous
                    self.SaveDCAC.SaveDCDict(Ids=self.DCIds,
                                             Dev=self.Dev,
                                             SwVgsInd=self.VgIndex,
                                             SwVdsInd=self.VdIndex)    
                self.Buffer.Reset()

            else:
                Qt.QThread.msleep(10)

    def AddData(self, NewData):
        '''
        Function used to add continuously data to the buffer.

        Parameters
        ----------
        NewData : numpy array
            array with the DC or AC data acquired by the system

        Returns
        -------
        None.

        '''
        # if the data is not stable
        if self.Stable is False:
            while self.Buffer.IsFilled():
                # if the buffer is filled the data is ignored until the 
                # characterization process finishes and the buffer is empty
                continue
            # if wait is true
            if self.Wait:
                # the time elapsed since the begining is calculated
                self.ElapsedTime = self.ElapsedTime+len(NewData[:,0])*(1/self.FsDemod)
                Diff = self.DelayTime-self.ElapsedTime
                # when the delay time is finished
                if Diff <= 0:
                    print('Delay Time finished')
                    self.Wait = False
                    self.ElapsedTime = 0
                    # The TimeOut timer is set and started
                    self.Timer = Qt.QTimer()
                    self.Timer.timeout.connect(self.printTime)
                    self.Timer.setSingleShot(True)
                    self.Timer.start(self.TimeOut*1000)
            # if wait is false
            else:            
                # the acquired data is added to the buffer
                self.Buffer.AddData(NewData)
        # if AC characterization is required
        if self.ACenable:
            # and the acquired data is stable
            if self.Stable is True:
                # the data is added directly processed to calculated the PSD
                self.threadCalcPSD.AddData(NewData)

    def printTime(self):
        '''
        This function is executed when TimeOut time is over
        '''
        print('TimeOut')
        # The timer is stopped
        self.Timer.stop()
        self.Timer.deleteLater()
        # And the final slope is calculated
        self.CalcSlope()
        # Stable flag is set to True
        self.Stable = True
        # and the DataStab signal is emited
        self.DataStab.emit()
        # If AC characterization is required
        if self.ACenable:
            # PSD thread is started
            self.threadCalcPSD.start()
        # And the DC data is saved
        self.SaveDCAC.SaveDCDict(Ids=self.DCIds,
                                 Dev=self.Dev,
                                 SwVgsInd=self.VgIndex,
                                 SwVdsInd=self.VdIndex)  
        self.Buffer.Reset()

    def CalcSlope(self):
        '''
        this Function is used to calculate the slope of the acquired data.
        the slope is used to determine if the data in the buffer is stable.
        '''
        # a Deviation and an Ids numpy arrays are created
        self.Dev = np.ndarray((self.Buffer.shape[1],))
        self.DCIds = np.ndarray((self.Buffer.shape[1], 1))
        # the buffer data is separate into the different channels
        for ChnInd, dat in enumerate(self.Buffer.transpose()):
            r = len(dat)
            t = np.arange(0, (1/self.FsDemod)*r, (1/self.FsDemod))
            # with the time and the data, the value of Ids
            # and the deviation of the curve are obtained
            mm, oo = np.polyfit(t, dat, 1)
            self.Dev[ChnInd] = np.abs(np.mean(mm)) #slope (uA/s)
            self.DCIds[ChnInd] = oo
        print('Dev',self.Dev)
        # Stab value is set to 0
        Stab = 0
        # Depending on the stabilization criteria Stab value changes
        if self.StabCriteria == 'All channels':
            # the deviations of all channels must be considered
            for slope in self.Dev:
                # If one channel has an slope higher than the MaxSlope
                if slope > self.MaxSlope:
                    #Stab is set to -1
                    Stab = -1
                    # And the loop finishes
                    break
            if Stab == -1:
                # And Stable variable is set to False
                self.Stable = False
            else:
                # If none of the channels has a slope higher than MaxSlope
                # Stable variable is set to True
                self.Stable = True
        # If the stabilizaton criteria is that only one channels has to be stable
        elif self.StabCriteria == 'One Channel':
            # The deviations are checked one after the other
            for slope in self.Dev:
                # if one slope has a lower value than MaxSlope
                if slope < self.MaxSlope:
                    # Stable variable is set to True
                    self.Stable = True
                    # and the loop finishes
                    break
        # if stabilization criteria used is the mean of all slopes
        elif self.StabCriteria == 'Mean':
            # the mean of Deviation is calculated
            slope = np.mean(self.Dev)
            # If the value is lower than MaxSlope
            if slope < self.MaxSlope:
                # Stable value is set to True
                self.Stable = True
                
    def on_PSDDone(self):
        '''
        Function used for saving AC characterization data

        Returns
        -------
        None.

        '''
        # frequencies and data from the PSD thread are obtained
        self.freqs = self.threadCalcPSD.ff
        self.PSDdata = self.threadCalcPSD.psd
        # The PSD thread is stopped
        self.threadCalcPSD.stop()
        # And the AC characterization is saved
        self.SaveDCAC.SaveACDict(psd=self.PSDdata,
                                 ff=self.freqs,
                                 SwVgsInd=self.VgIndex,
                                 SwVdsInd=self.VdIndex
                                 )
        # Finally the PSD characterization plot is updated
        self.UpdateAcPlots(self.SaveDCAC.DevACVals)

    def on_NextVgs(self):
        '''
        Function executed to change the Vgs value
        '''
        # The buffer is reset
        self.Buffer.Reset()
        # Stable is set to False
        self.Stable = False
        # The index of Vg is incremented by 1
        self.VgIndex += 1
        # If the index value is lower than the Vgs sweep array  length
        if self.VgIndex < len(self.VgSweepVals):
            # The value of Vgs is changed to the next one of the array
            self.NextVgs = self.VgSweepVals[self.VgIndex]
            # Wait is set to True
            self.Wait = True
            print(self.VgIndex)
            # and the DC characterization plots are updated
            self.UpdateSweepDcPlots(self.SaveDCAC.DevDCVals)
            # Finally a signal is emit to notify the change of Vgs
            self.NextVg.emit()
        # if the sweep of Vgs has finished
        else:
            # the Index is initialized to 0
            self.VgIndex = 0
            # and the Vgs value is set to the first value of the array
            self.NextVgs = self.VgSweepVals[self.VgIndex]
            # DC characterization plots are updated
            self.UpdateSweepDcPlots(self.SaveDCAC.DevDCVals)
            # finally on_NextVds function is called
            self.on_NextVds()

    def on_NextVds(self):
        '''
        Function executed to change the Vds value
        '''
        # The index of Vd is incremented by 1
        self.VdIndex += 1
        # IF the index value is lower than the Vds sweep array  length
        if self.VdIndex < len(self.VdSweepVals):
            # Next value of the array is set as Vds
            self.NextVds = self.VdSweepVals[self.VdIndex]
            # Wait is set to True
            self.Wait = True
            print(self.VdIndex)
            # And a signal is emitted to notify the change of Vds
            self.NextVd.emit()
        # If the sweep of Vds has finished
        else:
            # The index is initialized to 0
            self.VdIndex = 0
            # and the value of Vds is set to the first value of the sweep array
            self.NextVds = self.VdSweepVals[self.VdIndex]
            # DC characterization is saved
            self.DCDict = self.SaveDCAC.DevDCVals
            # If AC characterization is required
            if self.ACenable:
                # AC characterization is saved
                self.ACDict = self.SaveDCAC.DevACVals
            else:
                # otherwise, Ac characterization is saved as None
                self.ACDict = None
            # A signal notifying the end of the characterization is emited
            self.CharactEnd.emit()
           
    def UpdateSweepDcPlots(self, Dcdict):
        '''
        this function is used to uptade the DC characterization plots
        '''
        
        if self.PlotSwDC:
            self.PlotSwDC.ClearAxes()
            self.PlotSwDC.PlotDataCh(Data=Dcdict)
            self.PlotSwDC.AddLegend()
            self.PlotSwDC.Fig.canvas.draw()  
            
    def UpdateAcPlots(self, Acdict):
        '''
        this function is used to uptade the AC characterization plots
        '''
        
        if self.PlotSwAC:
            self.PlotSwAC.ClearAxes()
            self.PlotSwAC.PlotDataCh(Data=Acdict)
            self.PlotSwAC.Fig.canvas.draw()
            
    def stop(self):
        '''
        Function used to disconnect the characterization thread emitting 
        signals, and stop the PSD and the characterization threads
        '''
        
        if self.threadCalcPSD is not None:
            self.SaveDCAC.PSDSaved.disconnect()
            self.threadCalcPSD.PSDDone.disconnect()
            self.threadCalcPSD.stop()
        self.terminate()
        
################CALC PSD THREAD###############################################
        
class CalcPSD(Qt.QThread):
    PSDDone = Qt.pyqtSignal()
    def __init__(self, Fs, nFFT, nAvg, nChannels, scaling, **kwargs):
        '''
        Initialization of the thread that is used to calculate the PSD

        Parameters
        ----------
        Fs : float
            Sampling frequency of the buffered data
        nFFT : float
            nperseg of signal.welch. 
        nAvg : integer
            DESCRIPTION.
        nChannels : integer
            Number of acquisition channels active.
        scaling : str
            Two options, Density or Spectrum
        **kwargs : kwargs
            Parameters received but not used

        Returns
        -------
        None.

        '''

        super(CalcPSD, self).__init__()
        # class variables are assigned
        self.scaling = scaling
        self.nFFT = 2**nFFT
        self.nChannels = nChannels
        self.Fs = Fs
        # Buffer size is calculated
        self.BufferSize = self.nFFT * nAvg
        # A buffer is created
        self.Buffer = PltBuffer2D.Buffer2D(self.Fs, self.nChannels,
                                           self.BufferSize/self.Fs)

    def run(self, *args, **kwargs):
        while True:
            # When the PSD buffer is filled
            if self.Buffer.IsFilled():
                # the welch is donde to obtain the frequencies and values
                # of the PSD obtained from the data
                self.ff, self.psd = welch(self.Buffer,
                                          fs=self.Fs,
                                          nperseg=self.nFFT,
                                          scaling=self.scaling,
                                          axis=0)
                # The buffer is reset
                self.Buffer.Reset()
                # And a signal notifying that PSD has been done is emitted
                self.PSDDone.emit()

            else:
                Qt.QCoreApplication.processEvents()
                Qt.QThread.msleep(200)

    def AddData(self, NewData):
        '''
        Function used to add the new data to the PSD buffer

        Parameters
        ----------
        NewData : numpy array
            Data acquired that has been considered stable and it is wanted
            to obtain its PSD

        Returns
        -------
        None.

        '''
        
        self.Buffer.AddData(NewData)

    def stop(self):
        self.Buffer.Reset()
        self.terminate()

################SAVE CHARACTERIZATION DICTs###################################
        
class SaveDicts(QObject):
    PSDSaved = Qt.pyqtSignal()
    DCSaved = Qt.pyqtSignal()

    def __init__(self, SwVdsVals, SwVgsVals, Channels,
                 nFFT, FsDemod, Gate=False, ACenable=True):
        '''Initialize the Dictionaries to Save the Characterization
           SwVdsVals: array. Contains the values for the Vd sweep
                             [0.1, 0.2]
           SwVgsVals: array. Contains the values for the Vg sweep
                             [ 0.,  -0.1, -0.2, -0.3]
           Channels: dictionary. Contains the names from each demodulated
                                 channel and column and its index
                                 {'Ch04Col1': 0, 'Ch05Col1': 1, 'Ch06Col1': 2}
           nFFT: int.
                   8
           FsDemod: float. Sampling Frequency used in the Demodulation Process
                           5000.0
        '''
        super(SaveDicts, self).__init__()

        self.ChNamesList = sorted(Channels)
        self.ChannelIndex = Channels
        self.DevDCVals = self.InitDCRecord(nVds=SwVdsVals,
                                           nVgs=SwVgsVals,
                                           ChNames=self.ChNamesList,
                                           Gate=Gate)
        # AC dictionaries
        if ACenable:
            self.PSDnFFT = 2**nFFT
            self.PSDFs = FsDemod
    
            Fpsd = np.fft.rfftfreq(self.PSDnFFT, 1/self.PSDFs)
            nFgm = np.array([])
    
            self.DevACVals = self.InitACRecord(nVds=SwVdsVals,
                                               nVgs=SwVgsVals,
                                               nFgm=nFgm,
                                               nFpsd=Fpsd,
                                               ChNames=self.ChNamesList)
        
    def InitDCRecord(self, nVds, nVgs, ChNames, Gate):

        Time = datetime.datetime.now()
        DevDCVals={}
        for Ch in ChNames:
            DCVals={'Ids':np.ones((len(nVgs),len(nVds)))*np.NaN,
                    'Dev':np.ones((len(nVgs),len(nVds)))*np.NaN,
                    'Vds':nVds,
                    'Vgs':nVgs,
                    'ChName':Ch,
                    'Name':Ch,
                    'DateTime':Time}
            DevDCVals[Ch]=DCVals
    
        if Gate:
            GateDCVals = {'Ig':np.ones((len(nVgs),len(nVds)))*np.NaN,
                        'Vds':nVds,
                        'Vgs':nVgs,
                        'ChName':'Gate',
                        'Name':'Gate',
                        'DateTime':Time}
            DevDCVals['Gate']=GateDCVals
    
        return DevDCVals

    def InitACRecord(self, nVds, nVgs, nFgm, nFpsd, ChNames):
    
        Time = datetime.datetime.now()
        DevACVals={}
        for Ch in ChNames:
            noise = {}
            gm = {}
            for i in range(nVds.size):
                noise['Vd{}'.format(i)] = np.ones((len(nVgs),nFpsd.size))*np.NaN
                gm['Vd{}'.format(i)] = np.ones((len(nVgs),nFgm.size))*np.NaN*np.complex(1)
    
            ACVals={'PSD':noise,
                    'gm':gm,
                    'Vgs':nVgs,
                    'Vds':nVds,
                    'Fpsd':nFpsd,
                    'Fgm':nFgm,
                    'ChName':Ch,
                    'Name':Ch,
                    'DateTime':Time}
            DevACVals[Ch]=ACVals
    
        return DevACVals
    
    def SaveDCDict(self, Ids, Dev, SwVgsInd, SwVdsInd):
        '''Function that Saves Ids Data in the Dc Dict in the appropiate form
           for database
           Ids: array. Contains all the data to be saved in the DC dictionary
           SwVgsInd: int. Is the index of the actual Vg Sweep Iteration
           SwVdsInd: int. Is the Index of the actual Vd Sweep iteration
        '''
        for chn, inds in self.ChannelIndex.items():
            self.DevDCVals[chn]['Ids'][SwVgsInd,
                                       SwVdsInd] = Ids[inds]
            self.DevDCVals[chn]['Dev'][SwVgsInd,
                                       SwVdsInd] = Dev[inds]
        self.DCSaved.emit()

        # print('DCSaved')

    def SaveACDict(self, psd, ff, SwVgsInd, SwVdsInd):
        '''Function that Saves PSD Data in the AC Dict in the appropiate form
           for database
           psd: array(matrix). Contains all the PSD data to be saved in the AC
                               dictionary
           ff: array. Contains the Frequencies of the PSD to be saved in the AC
                      dictionary
           SwVgsInd: int. Is the index of the actual Vg Sweep Iteration
           SwVdsInd: int. Is the Index of the actual Vd Sweep iteration
        '''
        for chn, inds in self.ChannelIndex.items():
            self.DevACVals[chn]['PSD']['Vd{}'.format(SwVdsInd)][
                    SwVgsInd] = psd[:, inds].flatten()
            self.DevACVals[chn]['Fpsd'] = ff
       
        self.PSDSaved.emit()

    def SaveDicts(self, Dcdict, Folder, Oblea, Disp, Name, Cycle, Acdict=None):
        '''Creates the appropiate Folder NAme to be upload to the database
           Dcdict: dictionary. Dictionary with DC characterization that has
                               the structure to be read and save correctly
                               in the database
                               {'Ch04Col1': {'Ids': array([[1.94019351e-02],
                                                           [5.66072141e-08],
                                                           [5.66067698e-08],
                                                           [5.65991858e-08]
                                                           ]),
                                             'Dev': array([])
                                             'Vds': array([0.1]),
                                             'Vgs': array([ 0. , -0.1,
                                                           -0.2, -0.3]),
                                             'ChName': 'Ch04Col1',
                                             'Name': 'Ch04Col1',
                                             'DateTime': datetime.datetime
                                                         (2019, 12, 19, 16, 20,
                                                         59, 52661)
                                             },
  
           Acdict: dictionary. Dictionary with AC characterization that has the
                               structure to be read and save correctly in the
                               database
                               {'Ch04Col1': {'PSD': {'Vd0': array([
                                                            [4.67586928e-26,
                                                            1.61193712e-25],
                                                            ...
                                                            [5.64154950e-26,
                                                            2.10064857e-25]
                                                                   ])
                                                    },
                                             'gm': {'Vd0': array([],
                                                                 shape=(4, 0),
                                                                 dtype=complex128
                                                                 )
                                                    },
                                             'Vgs': array([ 0. , -0.1,
                                                           -0.2, -0.3]),
                                             'Vds': array([0.1]),
                                             'Fpsd': array([   0., 19.53125,
                                                            ...  2500.]),
                                             'Fgm': array([], dtype=float64),
                                             'ChName': 'Ch04Col1',
                                             'Name': 'Ch04Col1',
                                             'DateTime': datetime.datetime
                                                         (2019, 12, 19, 16, 20,
                                                         59, 52661)
                                             },
                               
                            }
           Folder, Oblea, Disp, Name, Cycle: str.
        '''
        self.FileName = '{}/{}-{}-{}-Cy{}.h5'.format(Folder,
                                                     Oblea,
                                                     Disp,
                                                     Name,
                                                     Cycle)
#        print(self.FileName)
        with open(self.FileName, "wb") as f:
            pickle.dump(Dcdict, f)
            if Acdict is not None:
                pickle.dump(Acdict, f)
        print('Saved')