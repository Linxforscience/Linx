#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Linx.py
Lock-in numerique

Created by X. FABREGES on 2012-02-13.
Modified: 2014-01-23

"""

import sys
import time
import platform

os_version = platform.system()
if os_version.find('Windows')>= 0:
    from ui_linx_win import Ui_MainWindow
else:
    from ui_linx import Ui_MainWindow
del os_version

import array as array_test
from numpy import frombuffer, savetxt
from scipy import zeros, cos, pi, sin, fft, real, mean, optimize,\
fromfile, floor, exp, argmax, sqrt, diff, float32
from scipy.integrate import trapz, cumtrapz
from scipy.fftpack import fftfreq

from scipy.signal import buttord, ellipord, cheb1ord, ellip, cheby1, filtfilt
from modified_filter import butter # error free scipy butter

import PyQt4.Qwt5 as Qwt 
from PyQt4.QtCore import Qt,  SIGNAL
from PyQt4.QtGui import QApplication,  QMainWindow,  QPen,  QFileDialog

class switch(object):
    ''' Switch/case like class '''
    def __init__(self, value):
        self.value  =  value
        self.fall  =  False

    def __iter__(self):
        yield self.match
        raise StopIteration
    
    def match(self, *args):
        if self.fall or not args:
            return True
        elif self.value in args:
            self.fall  =  True
            return True
        else:
            return False
        
class Numlockin_Window(QMainWindow, Ui_MainWindow):
    ''' Global definition of window and functions '''
    def __init__(self):
        ''' Initialization of GUI '''
        QMainWindow.__init__(self)
        Ui_MainWindow.__init__(self)
        self.setupUi(self)
        
        # Global variables
        self.data = 0
        self.data_loaded = 0
        self.label_40.setText('Version: 20140224')
        self.outputdir=''

        # Initialization of fields to default values
        # Acquisition fields
        self.lineEdit.setText("1000")
        self.cutoff = 1000
        self.label_15.setText('us')
        self.lineEdit_2.setText("0")
        self.lineEdit_3.setText("0")
        self.dephase = float(self.lineEdit_3.text())
        
        # Input file fields
        self.lE_ncol.setText("3")
        self.nbcol = int(self.lE_ncol.text())
        self.phasemin = zeros(int(self.nbcol)-2)
        self.lE_dBdtcol.setText("1")
        self.coldbdt = int(self.lE_dBdtcol.text())
        self.lE_refcol.setText("2")
        tempstring = self.lE_refcol.text()
        tempstring = tempstring.split(', ')
        self.colrefi = zeros(len(tempstring))
        for i in range(0, len(tempstring)):
            self.colrefi[i] = int(tempstring[i])            
        self.lE_datcol.setText("3")
        tempstring = self.lE_datcol.text()
        tempstring = tempstring.split(', ')
        self.coldatai = zeros(len(tempstring))
        for i in range(0, len(tempstring)):
            self.coldatai[i] = int(tempstring[i])
        
        # Experimental parameters
        self.lE_srate.setText("500000")
        self.sample_rate = float(self.lE_srate.text())
        self.lE_current.setText("10")
        self.current = float(self.lE_current.text())
        self.lE_PUarea.setText("4.666e-4")
        self.pu_area = float(self.lE_PUarea.text())
        self.lE_gain.setText("200")
        self.gain = float(self.lE_gain.text())
        self.current = float(self.lE_current.text())
        self.lineEdit_13.setText("10")
        self.checkBox_4.setChecked(True)
        
        # Signal processing parameters
        self.lineEdit_16.setText("auto")
        self.comboBox_2.addItem('ButterWorth')
        self.comboBox_2.addItem('Elliptic')
        self.comboBox_2.addItem('Chebyshev')
        self.lineEdit_17.setText("3")
        self.lineEdit_18.setText("24")
        self.lineEdit_19.setText("2")
        self.lineEdit_20.setText("auto")
        self.lineEdit_21.setText("2")
        self.zp = 2
        self.UI_cutoffHz()
        
        self.lineEdit_10.setText("0")
        self.checkBox.setChecked(False)
        
        self.checkBox_2.setChecked(True)
        self.checkBox_3.setChecked(True)
        self.radioButton.setChecked(True)
        
        self.checkBox_6.setChecked(False)
        self.dbdted = 0
        
        self.checkBox_7.setChecked(True)
        self.checkBox_8.setChecked(True)
        
        self.comboBox_3.addItem('Sinusoidal')
        self.comboBox_3.addItem('Triangular')
        
        self.checkBox_5.setChecked(True)
        
        # Initialization of Qwt;Plot areas
        self.qwtPlot.setTitle("Raw Signal")
        self.qwtPlot.setAxisTitle(Qwt.QwtPlot.yLeft, "Raw Voltage (V)")
        self.qwtPlot.setAxisTitle(Qwt.QwtPlot.xBottom, "time (s)")
        self.qwtPlot_2.setTitle("Signal")
        self.qwtPlot_2.setAxisTitle(Qwt.QwtPlot.yLeft, "Resistance (Ohms)")
        self.qwtPlot_2.setAxisTitle(Qwt.QwtPlot.xBottom, "Field (T)")
        self.qwtPlot_3.setTitle("Ref. Signal")
        self.qwtPlot_3.setAxisTitle(Qwt.QwtPlot.yLeft, "Ref. Voltage (V)")
        self.qwtPlot_3.setAxisTitle(Qwt.QwtPlot.xBottom, "time (s)")
        self.qwtPlot_4.setTitle("Magnetic field")
        self.qwtPlot_4.setAxisTitle(Qwt.QwtPlot.yLeft, "B (T), PU (V)")
        self.qwtPlot_4.setAxisTitle(Qwt.QwtPlot.xBottom, "time (s)")
        
        # Connect each field to corresponding funtion
        self.connect(self.lineEdit, SIGNAL("returnPressed ()"), self.UI_cutoffHz)
        self.connect(self.lineEdit_4, SIGNAL("returnPressed ()"), self.UI_cutoffms)
        
        self.connect(self.lE_srate, SIGNAL("returnPressed ()"), self.UI_config)
        self.connect(self.lineEdit_3, SIGNAL("returnPressed ()"), self.UI_config)
        self.connect(self.lE_PUarea, SIGNAL("returnPressed ()"), self.UI_config)
        self.connect(self.lE_current, SIGNAL("returnPressed ()"), self.UI_config)
        self.connect(self.lE_gain, SIGNAL("returnPressed ()"), self.UI_config)
        self.connect(self.checkBox_4, SIGNAL("clicked()"), self.plot_lpass)
        
        self.connect(self.lE_ncol, SIGNAL("editingFinished ()"), self.UI_config_raw)
        self.connect(self.lE_dBdtcol, SIGNAL("editingFinished ()"), self.UI_config_raw)
        self.connect(self.lE_refcol, SIGNAL("editingFinished ()"), self.UI_config_raw)
        self.connect(self.lE_datcol, SIGNAL("editingFinished ()"), self.UI_config_raw)
        
        self.connect(self.pushButton_5, SIGNAL("clicked()"), self.f_loadconfig) #Chargement d'un fichier de config
        self.connect(self.pushButton_4, SIGNAL("clicked()"), self.f_saveconfig) #Sauvegarde d'un fichier de config
        self.connect(self.pushButton_6, SIGNAL("clicked()"), self.f_batch)      #Traitement en série de données
        self.connect(self.pushButton_7, SIGNAL("clicked()"), self.f_loadfile)   #Chargement d'un fichier de données
        self.connect(self.pushButton_8, SIGNAL("clicked()"), self.f_savefile)   #Sauvegarde d'un fichier chargé
        
        self.connect(self.pushButton_9, SIGNAL("clicked()"), sys.exit) 
        
        self.connect(self.comboBox_2, SIGNAL("activated(int)"), self.plot_lpass) 
        self.connect(self.lineEdit_16, SIGNAL("returnPressed ()"), self.plot_lpass)
        self.connect(self.lineEdit_17, SIGNAL("returnPressed ()"), self.plot_lpass)
        self.connect(self.lineEdit_18, SIGNAL("returnPressed ()"), self.plot_lpass)
        self.connect(self.lineEdit_19, SIGNAL("returnPressed ()"), self.plot_lpass)
        self.connect(self.lineEdit_20, SIGNAL("returnPressed ()"), self.data_artefact)
        self.connect(self.lineEdit_21, SIGNAL("returnPressed ()"), self.plot_lpass)
        
        self.connect(self.radioButton, SIGNAL("clicked()"), self.plot_combobox)
        self.connect(self.radioButton_2, SIGNAL("clicked()"), self.plot_combobox)
        
        self.connect(self.checkBox, SIGNAL("clicked()"), self.data_artefact)
        self.connect(self.lineEdit_10, SIGNAL("returnPressed ()"), self.data_artefact)
        
        self.connect(self.comboBox, SIGNAL("activated(int)"), self.plot_combobox)
        self.connect(self.checkBox_2, SIGNAL("clicked()"), self.plot_combobox)
        self.connect(self.checkBox_3, SIGNAL("clicked()"), self.plot_combobox)
        
        self.connect(self.tabWidget_2, SIGNAL("currentChanged(int)"), self.plot_combobox)
        
        self.connect(self.checkBox_6, SIGNAL("clicked()"), self.data_artefact)
        
        self.connect(self.checkBox_7, SIGNAL("clicked()"), self.plot_combobox)
        self.connect(self.checkBox_8, SIGNAL("clicked()"), self.plot_combobox)

        self.tabWidget.setCurrentIndex(0)
        self.tabWidget_2.setCurrentIndex(0)
        self.data_loaded = 0
        
        # Initialization of Qwt.Plot zoomer
        self.zoomer  =  Qwt.QwtPlotZoomer(self.qwtPlot.xBottom, 
                        self.qwtPlot.yLeft, 
                        Qwt.QwtPicker.DragSelection, 
                        Qwt.QwtPicker.AlwaysOff, 
                        self.qwtPlot.canvas())
        self.zoomer.setRubberBandPen(QPen(Qt.black))
        self.zoomer_2  =  Qwt.QwtPlotZoomer(self.qwtPlot_2.xBottom, 
                        self.qwtPlot_2.yLeft, 
                        Qwt.QwtPicker.DragSelection, 
                        Qwt.QwtPicker.AlwaysOff, 
                        self.qwtPlot_2.canvas())
        self.zoomer_2.setRubberBandPen(QPen(Qt.black))
        self.zoomer_3  =  Qwt.QwtPlotZoomer(self.qwtPlot_3.xBottom, 
                        self.qwtPlot_3.yLeft, 
                        Qwt.QwtPicker.DragSelection, 
                        Qwt.QwtPicker.AlwaysOff, 
                        self.qwtPlot_3.canvas())
        self.zoomer_3.setRubberBandPen(QPen(Qt.black))
        self.zoomer_4  =  Qwt.QwtPlotZoomer(self.qwtPlot_4.xBottom, 
                        self.qwtPlot_4.yLeft, 
                        Qwt.QwtPicker.DragSelection, 
                        Qwt.QwtPicker.AlwaysOff, 
                        self.qwtPlot_4.canvas())
        self.zoomer_4.setRubberBandPen(QPen(Qt.black))

    def f_loadconfig(self):
        ''' Open configuration file, modify every field accordingly '''
        inputfile = QFileDialog.getOpenFileName(self, "Open Configuration File", self.outputdir)
        if inputfile!= '':
            finput = open(str(inputfile), 'r')
        
            for tempstring in finput.read().split('\n'):
                tempstring = tempstring.split()
                for i in range(0, len(tempstring)):
                    tempstring[i] = tempstring[i].strip(' = ')

                if len(tempstring)>0:
                    for case in switch(tempstring[0]):
                        if case('#'):
                            break
                        if case('inputfile'):
                            self.inputdir=tempstring[2]
                            break
                        if case('outputfile'):
                            self.outputdir=tempstring[2]
                            break
                        if case('cutoff'):
                            self.lineEdit.setText(tempstring[2])
                            break
                        if case('phase'):
                            self.lineEdit_3.setText(tempstring[2])
                            break
                        if case('pickup'):
                            self.lE_PUarea.setText(tempstring[2])
                            break
                        if case('sample_rate'):
                            self.lE_srate.setText(tempstring[2])
                            break
                        if case('nb_col'):
                            self.lE_ncol.setText(tempstring[2])
                            break
                        if case('col_dbdt'):
                            self.lE_dBdtcol.setText(tempstring[2])
                            break
                        if case('col_ref'):
                            self.lE_refcol.setText(tempstring[2])
                            break
                        if case('current'):
                            self.lE_current.setText(tempstring[2])
                            break
                        if case('gain'):
                            self.lE_gain.setText(tempstring[2])
                            break
                        if case('thinning'):
                            self.lineEdit_13.setText(tempstring[2])
                            break
                        if case('col_data'):
                            self.lE_datcol.setText(tempstring[2])
                            break
                        if case('spikes'):
                            if float(tempstring[2])>0:
                                self.lineEdit_10.setText(tempstring[2])
                                self.checkBox.setChecked(True)
                            break
                        if case('plot_phase'):
                            if float(tempstring[2]) == 1:
                                self.checkBox_2.setChecked(True)
                                self.checkBox_3.setChecked(False)
                            elif float(tempstring[2]) == 2:
                                self.checkBox_2.setChecked(True)
                                self.checkBox_3.setChecked(False)
                            elif float(tempstring[2]) == 3:
                                self.checkBox_2.setChecked(True)
                                self.checkBox_3.setChecked(True)
                            else:
                                self.checkBox_2.setChecked(False)
                                self.checkBox_3.setChecked(False)
                            break
                        if case('binary'):
                            if float(tempstring[2]) == 1:
                                self.checkBox_5.setChecked(True)
                            if float(tempstring[2]) == 0:
                                self.checkBox_5.setChecked(False)    
                            break
                        if case('sub_dbdt'):
                            if float(tempstring[2]) == 1:
                                self.checkBox_6.setChecked(True)
                            break
                        if case('RMS'):
                            if float(tempstring[2]) == 1:
                                self.checkBox_4.setChecked(True)
                            break
            del inputfile
            self.UI_cutoffHz()
             
    def f_saveconfig(self,*args):
        ''' Save configuration file based on current GUI status '''
        if len(args)==0:
            inputfile = QFileDialog.getSaveFileName(self, "Save Configuration File",self.outputdir)
        else:
            inputfile=args[0]
            
        if inputfile!= '':
            finput = open(str(inputfile), 'w')
            # Ajout des commentaires d'aide pour l'utilisateur
            finput.write('# Linx configuration file\n')
            if self.data_loaded==1:
                finput.write('# Files\n')
                finput.write('inputfile =  '+str(self.inputdir)+'\n')
                finput.write('outputfile =  '+str(self.outputdir)+'\n')
            finput.write('# Cutoff frequency (Hz)\n')
            finput.write('cutoff =  '+str(self.lineEdit.text())+'\n')
            finput.write('# Manual phase (deg)\n')
            finput.write('phase =  '+str(self.lineEdit_3.text())+'\n')
            finput.write('# Pickup area (m2)\n')
            finput.write('pickup =  '+str(self.lE_PUarea.text())+'\n')
            finput.write('# Sample rate (Hz)\n')
            finput.write('sample_rate =  '+str(self.lE_srate.text())+'\n')
            finput.write('# Number of columns in raw file\n')
            finput.write('nb_col =  '+str(self.lE_ncol.text())+'\n')
            finput.write('# dBdt col. number\n')
            finput.write('col_dbdt =  '+str(self.lE_dBdtcol.text())+'\n')
            finput.write('# Signal reference col. number\n')
            finput.write('col_ref =  '+str(self.lE_refcol.text())+'\n')
            finput.write('# Data from sample col. number\n')
            finput.write('col_data =  '+str(self.lE_datcol.text())+'\n')
            finput.write('# Injected current (mA)\n')
            finput.write('current =  '+str(self.lE_current.text())+'\n')
            finput.write('# Gain/Amp.\n')
            finput.write('gain =  '+str(self.lE_gain.text())+'\n')
            finput.write('# Thinning of save files\n')
            finput.write('thinning =  '+str(self.lineEdit_13.text())+'\n')
            finput.write('# Automatic spikes suppression (experimental)\n')
            finput.write('spikes =  '+str(self.lineEdit_10.text())+'\n')
            finput.write('# Signal displayed/saved (in and/or out of phase)\n')
            finput.write('# 3 = in+out, 2 = out, 1 = in\n')
            if self.checkBox_2.isChecked() and self.checkBox_3.isChecked():
                finput.write('plot_phase =  3\n')
            elif self.checkBox_2.isChecked():
                finput.write('plot_phase =  1\n')
            elif self.checkBox_3.isChecked():
                finput.write('plot_phase =  2\n')
            else:
                finput.write('plot_phase =  1\n')
            finput.write('# up/down filed displayed/saved\n')
            finput.write('# 3 = up+down, 2 = dwn, 1 = up\n')
            if self.checkBox_7.isChecked() and self.checkBox_8.isChecked():
                finput.write('plot_updwn =  3\n')
            elif self.checkBox_8.isChecked():
                finput.write('plot_updwn =  2\n')
            elif self.checkBox_7.isChecked():
                finput.write('plot_updwn =  1\n')
            else:
                finput.write('plot_updwn =  2\n')
            finput.write('# RMS current (1=yes, 0=no)\n')
            if self.checkBox_4.isChecked():
                finput.write('RMS =  1\n')
            else:
                finput.write('RMS =  0\n')
            finput.write('# Bianry file (1=yes, 0=no)\n')
            if self.checkBox_5.isChecked():
                finput.write('binary =  1\n')
            else:
                finput.write('binary =  0\n')
            finput.write('# Remove dBdt from raw data (1=yes, 0=no)\n')
            if self.checkBox_6.isChecked():
                finput.write('sub_dbdt =  1\n')
            else:
                finput.write('sub_dbdt =  0\n')
            finput.close()
            del inputfile

    def f_load_data(self, inputfile):
        ''' Load data from user file '''
        if inputfile!= '':
            self.inputdir = inputfile
            if self.outputdir=='':
                self.outputdir=self.inputdir
            
            self.UI_config_raw()
            self.UI_config_light()
            
            self.label_23.setText('<font color = Red>Loading data</font>')
            app.processEvents()
            time.sleep(0.01)
            
            if self.checkBox_5.isChecked():
                inputtest = open(str(inputfile), 'rb')
                self.temp = fromfile(inputtest, float32, -1, "")
            else:
                inputtest = open(str(inputfile), 'r')
                self.temp = fromfile(inputtest, float, -1, " ")
                
            inputtest.close()
            
            os_version = platform.system()
            if os_version.find('Windows')>= 0:
                if len(self.temp)>524288:
                    if self.checkBox_5.isChecked():
                        taille = len(self.temp)
                        inputest = open(str(inputfile), 'rb')
                        
                        self.temp = array_test.array('f')
                        self.temp.fromfile(inputest, taille)
                        self.temp = frombuffer(self.temp, dtype = float32)

                        inputest.close()
            
            self.temp.resize(len(self.temp)/self.nbcol, self.nbcol)
            
            # On crée le tableau de données reformaté
            del self.data
            self.data = zeros((len(self.temp[:, 1]), len(self.coldatai)+len(self.colrefi)+3))
            # La première colonne contient les infos de temps
            for i in range(0, len(self.temp[:, 1])):
                self.data[i, 0] = i/self.sample_rate
        
            # La seconde colonne le dbdt
            self.data[:, 1] = self.temp[:, self.coldbdt-1]
            self.data[:, 1] = self.data[:, 1]-mean(self.data[:, 1])
            
            # Les len(self.colref) colonnes suivantes les données de références
            # prise en charge de plusieurs lock-in
            self.colref = self.colrefi
            self.freq = zeros(len(self.colref))
            for i in range(0, len(self.colref)):
                self.data[:, i+2] = self.temp[:, self.colref[i]-1]-mean(self.temp[:, self.colref[i]-1])
                self.colref[i] = i+2
            
                # On calcule les fréquences de modulation de chaque lock-in
                # On approche la valeur avec une FFT
                data_fourier = abs(real(fft(self.data[:, i+2])))
                freq = fftfreq(len(self.data[:, i+2]), 1/self.sample_rate)
                
                posj = argmax(data_fourier[0:round(len(data_fourier)/2)+1])
                
                self.freq[i] = freq[posj]
                
                # On effectue une DFT pour affiner le calcul
                deltaf = (freq[2]-freq[1])/1000
                fourier_lenght = 16084
                
                F = [0]*3
                expon = -2j*pi*self.data[0:fourier_lenght, 0]
                
                F[1] = abs(trapz(self.data[0:fourier_lenght, self.colref[i]]*exp(expon*self.freq[i]), self.data[0:fourier_lenght, 0]))
                F[2] = abs(trapz(self.data[0:fourier_lenght, self.colref[i]]*exp(expon*(self.freq[i]+deltaf)), self.data[0:fourier_lenght, 0]))
                F[0] = abs(trapz(self.data[0:fourier_lenght, self.colref[i]]*exp(expon*(self.freq[i]-deltaf)), self.data[0:fourier_lenght, 0]))
                
                if F[2]>F[1]:
                    esssaimax = F[1]
                    
                    while abs(deltaf)>0.0002:
                        F[2] = abs(trapz(self.data[0:fourier_lenght, self.colref[i]]*exp(expon*(self.freq[i]+deltaf)), self.data[0:fourier_lenght, 0]))
                        if F[2]>essaimax:
                            essaimax = F[2]
                            self.freq[i] = self.freq[i]+deltaf
                        else:
                            deltaf = -deltaf/10
                elif F[0]>F[1]:
                    deltaf = -deltaf
                    essaimax = F[1]
                    
                    while abs(deltaf)>0.0002:
                        F[0] = abs(trapz(self.data[0:fourier_lenght, self.colref[i]]*exp(expon*(self.freq[i]+deltaf)), self.data[0:fourier_lenght, 0]))
                        if F[0]>essaimax:
                            essaimax = F[0]
                            self.freq[i] = self.freq[i]+deltaf
                        else:
                            deltaf = -deltaf/10

            # Les len(self.coldata)    colonnes suivantes contiennent les données        
            self.coldata = zeros(len(self.coldatai))
            for i in range(0, len(self.coldata)):
                self.data[:, i+2+len(self.colref)] = self.temp[:, self.coldatai[i]-1]
                self.coldata[i] = i+2+len(self.colref)
            self.phasemin = zeros(len(self.coldata))
            self.phasemini = zeros(len(self.coldata))
        
            # Enfin la dernière colonne contient le champ intégré
            self.data[1:, 2+len(self.coldata)+len(self.colref)] = cumtrapz(self.data[:, 1], self.data[:, 0])
            self.data[1:, 2+len(self.coldata)+len(self.colref)] = abs(self.data[1:, 2+len(self.coldata)+len(self.colref)])

            self.f_max = abs(self.data[1:, 2+len(self.coldata)+len(self.colref)]).argmax()
            B_max = abs(self.data[self.f_max, 2+len(self.coldata)+len(self.colref)]/self.pu_area)
            if B_max<200:
                self.f_start = 1
                self.f_stop = len(self.data[:, 0])
            else:
                for i in range(0, int(floor(len(self.data[:, 1])/100))):
                    if abs(mean(self.data[i*100:i*100+99, 2+len(self.coldata)+len(self.colref)]))>1e-3*self.data[self.f_max, 2+len(self.coldata)+len(self.colref)]:
                        self.f_start = (i-1)*100
                        break
                for i in range(int(self.f_max/100+10), int(floor(len(self.data[:, 1])/100))):
                    if abs(mean(self.data[i*100:i*100+99, 2+len(self.coldata)+len(self.colref)]))<1e-3*self.data[self.f_max, 2+len(self.coldata)+len(self.colref)]:
                        self.f_stop = (i+3)*100
                        break
            
            # On affiche par défaut les infos relatives au premier échantillon
            self.lineEdit_2.setText(str(float(int(self.freq[0]*100))/100))
            app.processEvents()
            time.sleep(0.01)
        
            self.comboBox.clear()
            for i in range(0, len(self.coldata)):
                self.comboBox.addItem(str(i+1))
        
            self.checkBox.setChecked(False)
            self.spiked = 0
            self.checkBox_6.setChecked(False)
            self.dbdted = 0
            self.data_loaded = 1
        
            
            #cProfile.runctx('self.data_analyze()', globals(), locals(), 'essai.lprof')
            self.data_analyze()
            self.plot_analyzed()
            self.plot_raw()
            if B_max<0.5:
                self.label_23.setText('Low field detected, check dBdt column or PU area.')
            else:
                self.label_23.setText(str(inputfile))
        else:
            self.label_23.setText('Please select a data file')

    def f_loadfile(self):
        if self.data_loaded == 0:
            inputfile = QFileDialog.getOpenFileName(self, "Open Raw Data File")
        else:
            inputfile = QFileDialog.getOpenFileName(self, "Open Raw Data File", self.inputdir)
            
        self.f_load_data(inputfile)

    def f_save_data(self, outputfile):
        ''' Save treated data to ASCII '''
        if outputfile!= '':
            self.outputdir = outputfile
            if self.inputdir=='':
                self.inputdir=self.outputdir
            
            thinout = int(self.lineEdit_13.text())
            
            # Sauvegarde de la montee et de la descente du champ
            if self.checkBox_7.isChecked() and self.checkBox_8.isChecked():
                out_data = zeros((floor(len(self.data[0::thinout, 0])), 3+2*+len(self.coldata)))
                out_data[:, 0] = self.data[0::thinout, 0]
                out_data[:, 1] = self.data[0::thinout, 2+len(self.colref)+len(self.coldata)]/self.pu_area
                out_data[:, 2] = self.data[0::thinout, 1]
                for j in range(0, len(self.coldata)):
                    out_data[:, 2*j+3] = self.sig_out[0::thinout, 2*j]/self.intgain*1e3
                    out_data[:, 2*j+4] = self.sig_out[0::thinout, 2*j+1]/self.intgain*1e3
            
            # Sauvegarde de la montee du champ uniquement
            elif self.checkBox_7.isChecked():
                out_data = zeros((floor(len(self.data[0:self.f_max:thinout, 0])), 3+2*+len(self.coldata)))
                out_data[:, 0] = self.data[0:self.f_max:thinout, 0]
                out_data[:, 1] = self.data[0:self.f_max:thinout, 2+len(self.colref)+len(self.coldata)]/self.pu_area
                out_data[:, 2] = self.data[0:self.f_max:thinout, 1]
                for j in range(0, len(self.coldata)):
                    out_data[:, 2*j+3] = self.sig_out[0:self.f_max:thinout, 2*j]/self.intgain*1e3
                    out_data[:, 2*j+4] = self.sig_out[0:self.f_max:thinout, 2*j+1]/self.intgain*1e3
                
            # Sauvegarde de la descente du champ uniquement
            elif self.checkBox_8.isChecked():
                out_data = zeros((floor(len(self.data[self.f_max::thinout, 0])), 3+2*+len(self.coldata)))
                out_data[:, 0] = self.data[self.f_max::thinout, 0]
                out_data[:, 1] = self.data[self.f_max::thinout, 2+len(self.colref)+len(self.coldata)]/self.pu_area
                out_data[:, 2] = self.data[self.f_max::thinout, 1]
                for j in range(0, len(self.coldata)):
                    out_data[:, 2*j+3] = self.sig_out[self.f_max::thinout, 2*j]/self.intgain*1e3
                    out_data[:, 2*j+4] = self.sig_out[self.f_max::thinout, 2*j+1]/self.intgain*1e3
                
            # Sinon on sauvegarde tout
            else:
                out_data = zeros((floor(len(self.data[0::thinout, 0])), 3+2*+len(self.coldata)))
                out_data[:, 0] = self.data[0::thinout, 0]
                out_data[:, 1] = self.data[0::thinout, 2+len(self.colref)+len(self.coldata)]/self.pu_area
                out_data[:, 2] = self.data[0::thinout, 1]
                for j in range(0, len(self.coldata)):
                    out_data[:, 2*j+3] = self.sig_out[0::thinout, 2*j]/self.intgain*1e3
                    out_data[:, 2*j+4] = self.sig_out[0::thinout, 2*j+1]/self.intgain*1e3
                
            f_handle  =  file(str(outputfile), 'w')
            f_handle.write('#time\tB\tdBdt\tin_phase\tout_phase\n')
            savetxt(f_handle, out_data[0:len(out_data[:, 0])-2, :], fmt = '%10g', delimiter = '\t')
            f_handle.close()
            self.label_23.setText('Data saved')
    
    def f_savefile(self):
        if self.outputdir=='':
            outputfile = QFileDialog.getSaveFileName(self, "Save Data File")
        else:
            outputfile = QFileDialog.getSaveFileName(self, "Save Data File", self.outputdir)
        self.f_save_data(outputfile)
        outputfile_mod=outputfile.split('.')
        self.f_saveconfig(outputfile_mod[0]+'.conf')
            
    def f_batch(self):
        ''' Batch processing ''' 
        inputfile = QFileDialog.getOpenFileNames(self, "Open Raw Data Files")
        #Si des fichiers ont été choisi on lace le traitement
        if inputfile!= '':
            self.UI_config_raw()
            self.data_loaded = 1
        
            self.label_24.setText('<font color = Red>Batch Treatment In Progress</font>')
            app.processEvents()
            time.sleep(0.01)
            
            # Chargement de chacun des fichiers et sauvegarde
            for input_i in range(0, len(inputfile)):
                self.f_load_data(str(inputfile[input_i]))
                self.f_save_data(str(inputfile[input_i])+'_out.dat')
                self.f_saveconfig(str(inputfile[input_i])+'_out.conf')

                self.label_22.setText(str(input_i+1)+'/'+str(len(inputfile)))
                app.processEvents()
                time.sleep(.01)
        else:
            self.label_23.setText('Please select a data file')
                
        self.label_24.setText('<font color = Green>Batch Treatment Done</font>')
        
    def UI_cutoffHz(self):
        ''' Adjust time constant based on cuttof frequency '''
        self.cutoff = float(self.lineEdit.text())
        if self.cutoff == 0:
            self.lineEdit_4.setText('0')
        else:
            self.lineEdit_4.setText(str(int(1e9/(2*pi*self.cutoff))/1000))
        
        if self.data_loaded == 1:
            self.data_analyzefast()
            self.plot_analyzed()    
        
    def UI_cutoffms(self):
        ''' Adjust cutoff frequency based on time constant '''
        self.time_cste = float(self.lineEdit_4.text())
        if self.time_cste == 0:
            self.lineEdit.setText('0')
            self.cutoff = 0
        else:
            self.lineEdit.setText(str(int(1e9/(2*pi*self.time_cste))/1000))
            self.cutoff = float(self.lineEdit.text())
            
        if self.data_loaded == 1:
            self.data_analyzefast()
            self.plot_analyzed()
    
    def UI_config_light(self):
        ''' Get constant from field values '''
        self.sample_rate = float(self.lE_srate.text())
        self.dephase = float(self.lineEdit_3.text())
        self.pu_area = float(self.lE_PUarea.text())
        self.gain = float(self.lE_gain.text())
        self.current = float(self.lE_current.text())
        
    def UI_config(self):
        ''' Get constant from field values and analyze data accordingly '''
        self.sample_rate = float(self.lE_srate.text())
        self.dephase = float(self.lineEdit_3.text())
        self.pu_area = float(self.lE_PUarea.text())
        self.gain = float(self.lE_gain.text())
        self.current = float(self.lE_current.text())
        if self.data_loaded == 1:
            self.data_analyzefast()
            self.plot_analyzed()
            self.plot_raw()
            
    def UI_config_raw(self):
        ''' Get raw file constant from GUI '''
        self.nbcol = int(self.lE_ncol.text())
        self.coldbdt = int(self.lE_dBdtcol.text())
        
        tempstring = self.lE_datcol.text()
        tempstring = tempstring.split(', ')
        self.coldatai = zeros(len(tempstring))
        for i in range(0, len(tempstring)):
            self.coldatai[i] = int(tempstring[i])
        
        del tempstring
        tempstring = self.lE_refcol.text()
        tempstring = tempstring.split(', ')
        self.colrefi = zeros(len(tempstring))
        for i in range(0, len(tempstring)):
            self.colrefi[i] = tempstring[i]
    
    def plot_lpass(self):
        ''' If data already loaded, execute analysis '''
        if self.data_loaded == 1:
            self.data_analyze()
            self.plot_analyzed()    
    
    def plot_combobox(self):
        ''' If data loaded, refresh plots '''
        if self.data_loaded == 1:
            self.plot_analyzed()
            self.plot_raw()
        
    def plot_analyzed(self):
        ''' Plot analyzed data '''
        # On affiche la bonne fréquence de modulation
        act_sample = int(self.comboBox.currentIndex())
        if len(self.colref) == 1:
            thefreq = self.freq[0]
        else:
            thefreq = self.freq[act_sample]
            
        self.lineEdit_2.setText(str(float(int(thefreq*100))/100))
        app.processEvents()
        time.sleep(0.01)
        
        # On choisi les abscisses en fonction du choix utilisateur
        if self.radioButton.isChecked():
            if self.checkBox_7.isChecked() and self.checkBox_8.isChecked():
                x_space = abs(self.data[self.f_start:self.f_stop:10, 2+len(self.colref)+len(self.coldata)]/self.pu_area)
            elif self.checkBox_7.isChecked():
                x_space = abs(self.data[self.f_start:self.f_max:10, 2+len(self.colref)+len(self.coldata)]/self.pu_area)
            else:
                x_space = abs(self.data[self.f_max:self.f_stop:10, 2+len(self.colref)+len(self.coldata)]/self.pu_area)
            
            self.qwtPlot_2.setAxisTitle(Qwt.QwtPlot.xBottom, "Field (T)")
        else:
            x_space = self.data[self.f_start:self.f_stop:10, 0]
            self.qwtPlot_2.setAxisTitle(Qwt.QwtPlot.xBottom, "time (s)")
            
            
        self.qwtPlot_2.setTitle("Signal")
        self.qwtPlot_2.setAxisTitle(Qwt.QwtPlot.yLeft, "Resistance (Ohms)")
        # On sélectionne les données en fonction de l'échantillon
        # On prend en compte le gain et le courant d'excitation
        if self.checkBox_4.isChecked():
            self.intgain = 0.5*self.amplitude[act_sample]*(self.gain*self.current*sqrt(2))
        else:
            self.intgain = 0.5*self.amplitude[act_sample]*(self.gain*self.current)            
        
        if self.checkBox_7.isChecked() and self.checkBox_8.isChecked():
            y_space = self.sig_out[self.f_start:self.f_stop:10, 2*(act_sample)]/self.intgain*1e3
            y_space_2 = self.sig_out[self.f_start:self.f_stop:10, 2*(act_sample)+1]/self.intgain*1e3
        elif self.checkBox_7.isChecked():
            y_space = self.sig_out[self.f_start:self.f_max:10, 2*(act_sample)]/self.intgain*1e3
            y_space_2 = self.sig_out[self.f_start:self.f_max:10, 2*(act_sample)+1]/self.intgain*1e3
        else:
            y_space = self.sig_out[self.f_max:self.f_stop:10, 2*(act_sample)]/self.intgain*1e3
            y_space_2 = self.sig_out[self.f_max:self.f_stop:10, 2*(act_sample)+1]/self.intgain*1e3
    
        
        # On efface l'ancien plot
        self.qwtPlot_2.autoDelete()
        self.qwtPlot_2.detachItems()
        
        # On dimensionne les axes à la main
        self.qwtPlot_2.setAxisScale(Qwt.QwtPlot.xBottom, min(x_space), max(x_space))
        y_space_bound = [min(y_space), min(y_space_2), max(y_space), max(y_space_2)]
        if self.checkBox_2.isChecked() and self.checkBox_3.isChecked():
            self.qwtPlot_2.setAxisScale(Qwt.QwtPlot.yLeft, min(y_space_bound), max(y_space_bound))
        elif self.checkBox_3.isChecked():
            self.qwtPlot_2.setAxisScale(Qwt.QwtPlot.yLeft, y_space_bound[1], y_space_bound[3])
        elif self.checkBox_2.isChecked():
            self.qwtPlot_2.setAxisScale(Qwt.QwtPlot.yLeft, y_space_bound[0], y_space_bound[2])        
        
        #Initialisation de la base de zoom (zoom out max)
        self.zoomer_2.setZoomBase()
        
        # On réaffiche le nouveau
        if self.checkBox_2.isChecked():
            self.qwtPlot_2.curve  =  Qwt.QwtPlotCurve("In phase")
            self.qwtPlot_2.curve.detach()
            self.qwtPlot_2.curve.attach(self.qwtPlot_2)
            self.qwtPlot_2.curve.setRenderHint(Qwt.QwtPlotItem.RenderAntialiased)
            self.qwtPlot_2.curve.setPen(QPen(Qt.red))
            self.qwtPlot_2.curve.setData(x_space, y_space)
        if self.checkBox_3.isChecked():
            self.qwtPlot_2.curve2  =  Qwt.QwtPlotCurve("Out of phase")
            self.qwtPlot_2.curve2.detach()
            self.qwtPlot_2.curve2.attach(self.qwtPlot_2)
            self.qwtPlot_2.curve2.setRenderHint(Qwt.QwtPlotItem.RenderAntialiased)
            self.qwtPlot_2.curve2.setPen(QPen(Qt.black))
            self.qwtPlot_2.curve2.setData(x_space, y_space_2)

        self.qwtPlot_2.replot()
        
    def plot_raw(self):
        ''' Plot raw data '''
        self.qwtPlot.autoDelete()
        self.qwtPlot.detachItems()
        self.qwtPlot_3.autoDelete()
        self.qwtPlot_3.detachItems()
        self.qwtPlot_4.autoDelete()
        self.qwtPlot_4.detachItems()
        
        longueur = len(self.data[:, 0])
        
        act_sample = int(self.comboBox.currentIndex())
        y_space = self.data[1:longueur:4, 2+len(self.colref)+act_sample]
        if len(self.colref) == 1:
            y_space_2 = self.data[1:longueur:10, self.colref[0]]
        else:
            y_space_2 = self.data[1:longueur:10, self.colref[act_sample]]

        # Champ integre et pickup
        y_space_3 = self.data[1:longueur:10, 2+len(self.coldata)+len(self.colref)]/self.pu_area
        y_space_4 = max(y_space_3)*self.data[1:longueur:10, 1]/max(self.data[1:longueur:10, 1])
        
        self.qwtPlot.setTitle("Raw Signal")
        self.qwtPlot.setAxisTitle(Qwt.QwtPlot.yLeft, "Raw Voltage (V)")
        self.qwtPlot.setAxisTitle(Qwt.QwtPlot.xBottom, "time (s)")
        self.qwtPlot.curve  =  Qwt.QwtPlotCurve("Raw Voltage")
        self.qwtPlot.curve.attach(self.qwtPlot)
        self.qwtPlot.curve.setPen(QPen(Qt.red))
        self.qwtPlot.curve.setData(self.data[1:longueur:4, 0], y_space)
        self.zoomer.setZoomBase()
        self.qwtPlot.replot()
        
        self.qwtPlot_3.setTitle("Ref. Signal")
        self.qwtPlot_3.setAxisTitle(Qwt.QwtPlot.yLeft, "Ref. Voltage (V)")
        self.qwtPlot_3.setAxisTitle(Qwt.QwtPlot.xBottom, "time (s)")
        self.qwtPlot_3.curve  =  Qwt.QwtPlotCurve("Ref. Voltage")
        self.qwtPlot_3.curve.attach(self.qwtPlot_3)
        self.qwtPlot_3.curve.setPen(QPen(Qt.red))
        self.qwtPlot_3.curve.setData(self.data[1:longueur:10, 0], y_space_2)
        self.zoomer_3.setZoomBase()
        self.qwtPlot_3.replot()
        
        #Initialisation de la base de zoom (zoom out max)
        self.qwtPlot_4.setTitle("B(t)")
        self.qwtPlot_4.setAxisTitle(Qwt.QwtPlot.yLeft, "B (T), PU (V)")
        self.qwtPlot_4.setAxisTitle(Qwt.QwtPlot.xBottom, "time (s)")
        self.qwtPlot_4.curve  =  Qwt.QwtPlotCurve("Field")
        self.qwtPlot_4.curve.attach(self.qwtPlot_4)
        self.qwtPlot_4.curve.setRenderHint(Qwt.QwtPlotItem.RenderAntialiased)    
        self.qwtPlot_4.curve.setPen(QPen(Qt.red))
        self.qwtPlot_4.curve.setData(self.data[0::10, 0], y_space_3)
        self.qwtPlot_4.curve_2  =  Qwt.QwtPlotCurve("Pick-up")
        self.qwtPlot_4.curve_2.attach(self.qwtPlot_4)
        self.qwtPlot_4.curve_2.setRenderHint(Qwt.QwtPlotItem.RenderAntialiased)    
        self.qwtPlot_4.curve_2.setPen(QPen(Qt.black))
        self.qwtPlot_4.curve_2.setData(self.data[0::10, 0], y_space_4)
        self.zoomer_4.setZoomBase()
        self.qwtPlot_4.replot()
    
    def data_analyze(self):
        ''' Treat raw data '''
        self.sig_out = zeros((len(self.data[:, 2]), len(self.coldata)*2))
        phase_ref = zeros(len(self.data[:, 2]))
        antiphase_ref = zeros(len(self.data[:, 2]))
        self.amplitude = zeros(len(self.coldata))
        self.zp = int(self.lineEdit_21.text())
        
        for j in range(0, len(self.coldata)):
            if len(self.colref) == 1:
                posref = self.colref
                thefreq = self.freq[0]
            else:
                posref = self.colref[j]
                thefreq = self.freq[j]
            
            # Phase detector
            self.amplitude[j] = (max(self.data[1:10000, int(posref)])-min(self.data[1:10000, int(posref)]))/2
            phase_ref[:] = self.amplitude[j]*cos(2*pi*thefreq*self.data[:, 0])
            
            # Phase detector using sinusoidal fit
            fitfunc       = lambda p, x: self.amplitude[j]*cos(2*pi*thefreq*x+p[0])
            errfunc       = lambda p, x, y: fitfunc(p, x) - y
            p0            = [0]
            fit_le        = int(self.zp*round(self.sample_rate/thefreq))
            p1, success   = optimize.leastsq(errfunc, p0[:], args=(self.data[0:fit_le,0], self.data[0:fit_le, int(posref)]))
            self.phasemin = p1*180/pi
            
            phase_ref[:] = self.amplitude[j]*cos(2*pi*thefreq*self.data[:, 0]+self.phasemin[j]*pi/180+self.dephase*pi/180)
            self.sig_out[:, 2*j] = phase_ref[:]*self.data[:, self.coldata[j]]
            self.sig_out[:, 2*j] = self.data_lpass(self.sig_out[:, 2*j], self.cutoff, self.sample_rate)
            
            antiphase_ref[:] = -self.amplitude[j]*sin(2*pi*thefreq*self.data[:, 0]+self.phasemin[j]*pi/180+self.dephase*pi/180)
            self.sig_out[:, 2*j+1] = antiphase_ref[:]*self.data[:, self.coldata[j]]
            self.sig_out[:, 2*j+1] = self.data_lpass(self.sig_out[:, 2*j+1], self.cutoff, self.sample_rate)

    def data_analyzefast(self):
        ''' Same as data_analyze without phase detection '''
        phase_ref = zeros(len(self.data[:, 2]))
        antiphase_ref = zeros(len(self.data[:, 2]))
        self.sig_out = zeros((len(self.data[:, 2]), len(self.coldata)*2))
        
        for j in xrange(0, len(self.coldata)):
            if len(self.colref) == 1:
                thefreq = self.freq[0]
                posref = self.colref
            else:
                thefreq = self.freq[j]
                posref = self.colref[j]            
            
            self.amplitude[j] = (max(self.data[1:1000, int(posref)])-min(self.data[1:1000, int(posref)]))/2
            antiphase_ref[:] = -self.amplitude[j]*sin(2*pi*thefreq*self.data[:, 0]+self.phasemin[j]*pi/180+self.dephase*pi/180)    
            phase_ref[:] = self.amplitude[j]*cos(2*pi*thefreq*self.data[:, 0]+self.phasemin[j]*pi/180+self.dephase*pi/180)
            self.sig_out[:, 2*j+1] = antiphase_ref[:]*self.data[:, self.coldata[j]]
            self.sig_out[:, 2*j+1] = self.data_lpass(self.sig_out[:, 2*j+1], self.cutoff, self.sample_rate)
            self.sig_out[:, 2*j] = phase_ref[:]*self.data[:, self.coldata[j]]
            self.sig_out[:, 2*j] = self.data_lpass(self.sig_out[:, 2*j], self.cutoff, self.sample_rate)
    
    def data_lpass(self, x, Wp, srate):
        ''' Low-pass filter using various filter type '''
        tempstring = self.lineEdit_16.text()
        
        if tempstring == 'auto':
            Wp = float(Wp*2/srate)
            Ws = Wp*float(self.lineEdit_19.text())
            Rp = float(self.lineEdit_17.text())
            Rs = float(self.lineEdit_18.text())
            
            if self.comboBox_2.currentIndex() == 0:
                (norder, Wn) = buttord(Wp, Ws, Rp, Rs)
            elif self.comboBox_2.currentIndex() == 1:
                (norder, Wn) = ellipord(Wp, Ws, Rp, Rs)
            else:
                (norder, Wn) = cheb1ord(Wp, Ws, Rp, Rs)

        else:
            norder = float(tempstring)
            Wp = float(Wp*2/srate)
            Ws = Wp*2
            self.lineEdit_19.setText(str(Ws/Wp))
            Rp = 3
            self.lineEdit_17.setText(str(Rp))
            Rs = 0.3*norder*20
            self.lineEdit_18.setText(str(Rs))
            
            if self.comboBox_2.currentIndex() == 0:
                (norder, Wn) = buttord(Wp, Ws, Rp, Rs)
            elif self.comboBox_2.currentIndex() == 1:
                (norder, Wn) = ellipord(Wp, Ws, Rp, Rs)
            else:
                (norder, Wn) = cheb1ord(Wp, Ws, Rp, Rs)
            
        if self.comboBox_2.currentIndex() == 0:
            (b, a)  =  butter(norder, Wn)
        elif self.comboBox_2.currentIndex() == 1:
            (b, a)  =  ellip(norder, Rp, Rs, Wn)
        else:
            (b, a)  =  cheby1(norder, Rp, Wn)
            
        
        y  =  filtfilt(b, a, x)
        
        return(y)
        
    def data_hpass(self, x, Wp, srate):
        ''' High-pass filter '''
        Wp = float(Wp*2/srate)
        Ws = Wp*float(self.lineEdit_19.text())
        Rp = float(self.lineEdit_17.text())
        Rs = float(self.lineEdit_18.text())

        tempstring = self.lineEdit_16.text()
        if tempstring == 'auto':
            if self.comboBox_2.currentIndex() == 0:
                (norder, Wn) = buttord(Wp, Ws, Rp, Rs)
            elif self.comboBox_2.currentIndex() == 1:
                (norder, Wn) = ellipord(Wp, Ws, Rp, Rs)
            else:
                (norder, Wn) = cheb1ord(Wp, Ws, Rp, Rs)
        else:
            norder = float(tempstring)
            Wn = Wp

        if self.comboBox_2.currentIndex() == 0:
            (b, a)  =  butter(norder, Wn, btype = 'high')
        elif self.comboBox_2.currentIndex() == 1:
            (b, a)  =  ellip(norder, Rp, Rs, Wn)
        else:
            (b, a)  =  cheby1(norder, Rp, Wn)


        y  =  filtfilt(b, a, x)

        return(y)
                    
    def data_artefact(self):
        ''' Supress dB/dt from raw signal '''
        ''' Supress spikes from raw signal ... experimental '''
        if self.data_loaded == 1:
            if self.dbdted == 0 and self.spiked == 0:
                self.data_raw = zeros((len(self.data[:, 0]), len(self.data[0, :])))
                self.data_raw[:, :] = self.data[:, :]

            if self.checkBox_6.isChecked():
                # dbdt removal : high-pass filter on raw data
                self.dbdted = 1
                tempstring = self.lineEdit_20.text()
                if tempstring == 'auto':
                    dbdt_cutoff = min(self.freq)/20
                else:
                    dbdt_cutoff = float(self.lineEdit_20.text())
                
                self.label_24.setText('<font color = Red>Calculations In Progress</font>')
                app.processEvents()
                time.sleep(0.01)
                for i in range(0, len(self.coldata)):
                    self.data[:, self.coldata[i]] = self.data_hpass(self.data_raw[:, self.coldata[i]], dbdt_cutoff, self.sample_rate)

            if (self.checkBox.isChecked()) and float(self.lineEdit_10.text())>0:
                self.spiked = 1
                datai = zeros((len(self.data[:, 0]), len(self.data[0, :])))
                self.label_24.setText('<font color = Red>Calculations In Progress</font>')
                app.processEvents()
                time.sleep(0.01)
                for i in range(0, len(self.coldata)):
                    if self.dbdted == 1:
                        datai[:, self.coldata[i]] = self.data[:, self.coldata[i]]
                    else:
                        datai[:, self.coldata[i]] = self.data_raw[:, self.coldata[i]]

                    datai[1:, self.coldata[i]] = diff(datai[:, self.coldata[i]])
                    
                    spikes_power = 10./float(self.lineEdit_10.text())
                    
                    for j in range(20, len(datai[:, 0])-20):
                        if abs(datai[j, self.coldata[i]])>= 1.0001*spikes_power:
                            datai[j, self.coldata[i]] = (datai[j-20, self.coldata[i]]+datai[j+20, self.coldata[i]])/2
                        self.data[j, self.coldata[i]] = self.data[j-1, self.coldata[i]]+datai[j, self.coldata[i]]
                        
            if (self.checkBox.isChecked() == False) and (self.checkBox_6.isChecked() == False):
                self.data[:, :] = self.data_raw[:, :]
            
            self.data_analyze()
            self.plot_analyzed()
            self.plot_raw()

        self.label_24.setText('')

    def check_update(self):
        pass
        
if __name__  ==  '__main__':
    app  =  QApplication(sys.argv)
    window  =  Numlockin_Window()
    window.show()
    sys.exit(app.exec_())
