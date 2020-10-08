# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 12:45:22 2019

@author: Lucia
"""
from PyQt5 import Qt
import pyqtgraph.parametertree.parameterTypes as pTypes

from scipy import signal
import numpy as np

DemodulParams = ({'name': 'DemodConfig',
                  'type': 'group',
                  'children': ({'name': 'DemEnable',
                                'title': 'On/Off',
                                'type': 'bool',
                                'value': True}, 
                               {'name': 'FsDemod',
                                'type': 'float',
                                'value': 2e6,
                                'readonly': True,
                                'siPrefix': True,
                                'suffix': 'Hz'},
                               {'name': 'DSFs',
                                'title': 'DownSampling Fs',
                                'type': 'float',
                                'readonly': True,
                                'value': 10e3,
                                'siPrefix': True,
                                'suffix': 'Hz'},
                               {'name': 'DSFact',
                                'title': 'DownSampling Factor',
                                'type': 'int',
                                'value': 100},
                               {'name': 'FiltOrder',
                                'title': 'Filter Order',
                                'type': 'int',
                                'value': 2},
                               {'name': 'OutType',
                                'title': 'Output Var Type',
                                'type': 'list',
                                'values': ['Real', 'Imag', 'Angle', 'Abs'],
                                'value': 'Abs'},
                               )
                  })


class DemodParameters(pTypes.GroupParameter):
    def __init__(self, **kwargs):
        pTypes.GroupParameter.__init__(self, **kwargs)

        self.addChild(DemodulParams)
        self.DemConfig = self.param('DemodConfig')
        self.DemEnable = self.DemConfig.param('DemEnable')
        self.FsDem = self.DemConfig.param('FsDemod')
        self.DSFs = self.DemConfig.param('DSFs')
        self.DSFact = self.DemConfig.param('DSFact')
        self.on_DSFact_changed()
        self.FsDem.sigValueChanged.connect(self.on_FsDem_changed)
        self.FiltOrder = self.DemConfig.param('FiltOrder')
        self.OutType = self.DemConfig.param('OutType')

    def ReCalc_DSFact(self, BufferSize):
        while BufferSize % self.DSFact.value() != 0:
            self.DSFact.setValue(self.DSFact.value()+1)
        self.on_DSFact_changed()
        print('DSFactChangedTo'+str(self.DSFact.value()))

    def on_FsDem_changed(self):
        self.on_DSFact_changed()

    def on_DSFact_changed(self):
        DSFs = self.FsDem.value()/self.DSFact.value()
        self.DSFs.setValue(DSFs)

    def GetParams(self):
        '''Functions that return the Parameters of the Demodulation Process
           Returns a dictionary:
           {'FsDemod': 500000.0,
            'DSFact': 100,
            'FiltOrder': 2,
            'OutType': 'Abs'}
        '''
        Demod = {}
        for Config in self.DemConfig.children():
            if Config.name() == 'DemEnable':
                continue
            if Config.name() == 'DSFs':
                continue
            Demod[Config.name()] = Config.value()
        return Demod

    def GetChannels(self, Rows, Fcs):
        '''Function that returns a dictionary with the names of demodulation
           channels and indexes
            {'Ch01Col1':0,
             'Ch02Col1':1,
             'Ch03Col1':2,
             'Ch04Col1':3,
             'Ch05Col1':4,
             'Ch06Col1':5,
             'Ch07Col1':6,
             'Ch08Col1':7,
            }
        '''
        DemChnNames = {}
        i = 0
        for r in Rows:
            for col, f in Fcs.items():
                DemChnNames[r+col] = i
                i = i + 1
        return DemChnNames
    
    def GetChannelsNames(self, Rows, Fcs):
        '''Function that returns an array with the names of demodulation
           channels in dtype S10
            ['Ch01Col1',
             'Ch02Col1',
             'Ch03Col1',
             'Ch04Col1',
             'Ch05Col1',
             'Ch06Col1',
             'Ch07Col1',
             'Ch08Col1',
            ]
        '''
        DemChnNames = []
        for r in Rows:
            for col, f in Fcs.items():
                DemChnNames.append(r+col)
        DemChns = np.array(DemChnNames, dtype='S10')
        return DemChns


class Filter():
    def __init__(self, Fs, Freqs, btype, Order):
        freqs = np.array(Freqs)/(0.5*Fs)
        self.b, self.a = signal.butter(Order,
                                       freqs,
                                       btype,
                                       )
        self.zi = signal.lfilter_zi(self.b,
                                    self.a,
                                    )

    def Apply(self, Sig):
        sigout, self.zi = signal.lfilter(b=self.b,
                                        a=self.a,
                                        x=Sig,
                                        axis=0,
                                        zi=self.zi
                                        )
        #probar con filtfilt
        # sigout = signal.filtfilt(b=self.b,
        #                          a=self.a,
        #                          x=Sig,
        #                          axis=0,
        #                          padtype='constant',
        #                          # method='gust'
        #                         )
        return sigout


class Demod():
    def __init__(self, Fc, FetchSize, Fs, DownFact, Order, Signal):
        ''' Demodulation Class, applies the filters and the resampling process.
            Fc: float. Frequency of the Carrier used for Modulation
            FetchSize: int. Defines the number of samples of the buffer of
                            acquisition
            Fs: float. Sampling Frequency used for acquisition process
            DownFact: int. Down Sampling Factor to calculate Sampling
                           Frequency of the demodulation process
            Order: int. Order of the internal filter of the process
            Signal: array. Contains the values that forms the carrier signal
                           used in Modulation
        '''
        self.Fs = Fs
        self.Fc = Fc
        self.DownFact = DownFact
        self.FsOut = Fs/DownFact

        self.FiltR = Filter(Fs, self.FsOut/2, 'lp', Order)
        self.FiltI = Filter(Fs, self.FsOut/2, 'lp', Order)
        self.FiltH = Filter(Fs, 2*5e3, 'highpass', 4)
        # self.FiltL = Filter(Fs, 200e3, 'lp', 4)
        
        # self.FiltR2 = Filter(Fs, self.FsOut/2, 'lp', Order)
        # self.FiltI2 = Filter(Fs, self.FsOut/2, 'lp', Order)

        self.vcoi = Signal

    def Apply(self, SigInput):
        SigIn = self.FiltH.Apply(SigInput)
        # SigIn = self.FiltH.Apply(SigInHigh)
        rdem = np.real(self.vcoi*SigIn)
        idem = np.imag(self.vcoi*SigIn)

        FilterRPart = self.FiltR.Apply(rdem)
        FilterIPart = self.FiltI.Apply(idem)
        
        # FilterRPart = self.FiltR2.Apply(FilterRPart)
        # FilterIPart = self.FiltI2.Apply(FilterIPart)

        sObject = slice(None, None, self.DownFact)

        RSrdem = FilterRPart[sObject]
        RSidem = FilterIPart[sObject]

        complexDem = RSrdem + (RSidem*1j)
        # complexDem = FilterRPart + (FilterIPart*1j)

        return complexDem


class DemodThread(Qt.QThread):
    NewData = Qt.pyqtSignal()

    def __init__(self, Carriers, ScopeChannels, BufferSize, FsDemod, DSFact,
                 FiltOrder, Signal, Gain, **Keywards):
        '''Initialization of Demodulation Process Thread
           Fcs: dictionary. returns the name of the columns with its carrier
                            frequency
                            {'Col1': 30000.0}
           RowList: array. Contains the names of the rows that are being
                           acquired
                           ['Ch04', 'Ch05', 'Ch06']
           FetchSize: int. Defines the number of samples to fill the buffer.
           FsDemod: float. Specifies the Sampling Frequency of the Acquisition
                           process
           DSFacti: int. Specifies de DownSampling Factor to reduce sampling
                         frequency
           FiltOrder: int. Defines the order of the internal filter of the
                           demodulation process
           Signal: array. Contains the values that forms the carrier signal
                          used in Modulation process
           Keywords: dictionary. Contains the output Type of demodulation,
                                 absolut, real, imaginary or angle
                                 {'OutType': 'Abs'}
        '''
        super(DemodThread, self).__init__()
        self.ToDemData = None

        self.Gain = Gain
        self.DemOutputs = []
        self.NamesForDict = []
        for Row in ScopeChannels:
            DemOut = []
            for Cols, Freq in Carriers.items():
                Dem = Demod(Freq, BufferSize, FsDemod, DSFact,
                            FiltOrder, Signal)
                DemOut.append(Dem)
                self.NamesForDict.append(str(Row+Cols))
            self.DemOutputs.append(DemOut)
        self.OutDemodData = np.ndarray((round(BufferSize/DSFact),
                                        round(len(ScopeChannels)*len(Carriers.keys()))),
                                        dtype=complex)

    def run(self):
        while True:
            if self.ToDemData is not None:
                ind = 0
                for ir, rows in enumerate(self.DemOutputs):
                    for instance in rows:
                        data = instance.Apply(self.ToDemData[:, ir])
                        #factor 2 a causa de la demodulación ya que el resultado
                        #es (1/2)*Vin*Vcoi y dividido por la ganancia para tener
                        #corriente
                        self.OutDemodData[:, ind] = 2*data/self.Gain
                        ind = ind + 1

                self.NewData.emit()
                self.ToDemData = None
            else:
                Qt.QThread.msleep(10)
#        #multiprocessing

    def AddData(self, NewData):
        if self.ToDemData is not None:
            print('Error Demod !!!!')
        self.ToDemData = NewData

    def stop(self):
        self.terminate()
