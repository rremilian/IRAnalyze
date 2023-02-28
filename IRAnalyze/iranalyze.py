import sys
import os
import numpy as np
import scipy as sp
import csv

"""
IRAnalyze - A Python module to analyze experimental and theoretical spectra
@github: https://github.com/rremilian
"""

class IRAnalyze:
    def __init__(self):
        self.__version__ = "1.0.0"

class File:
    
    def __init__(self, location):
        self.location = location
        if os.name == "nt":
            self.name = self.location.split("\\")[-1]
        elif os.name == "posix":
            self.name = self.location.split("/")[-1]
    
class LogFile(File):
    
    def __init__(self, location, raman=False):
        File.__init__(self,location)
        self.freqlist = []
        self.irintens = []
        self.ramanact = []
         
        with open(self.location, "r") as f:
            self.lines = f.read()
            self.lines = self.lines.split("\n")
            for line in self.lines:
                if "Frequencies ---" in line:
                    line = line.split()
                    line = " ".join(line)
                    freq = line.split(" ")[2:7]
                    
                    for f in freq:
                        f = float(f)
                        self.freqlist.append(f)
                
                elif "IR Intensities ---" in line:
                    line = line.split()
                    line = " ".join(line)
                    intens = line.split(" ")[3:8]
                    
                    for i in intens:
                        i = float(i)
                        self.irintens.append(i)
                
                elif raman == True and "Raman Activities ---" in line:
                    line = line.split()
                    line = " ".join(line)
                    intens = line.split(" ")[3:8]
                    
                    for i in intens:
                        i = float(i)
                        self.ramanact.append(i)
                    
            self.nfreq = len(self.freqlist)
            self.freqlist = np.array(self.freqlist)
            self.irintes = np.array(self.irintens)
            self.ramanact = np.array(self.ramanact)

    def calcspectrum(self, emin, emax, step, sigma, scale=1.00, raman=False, mode="lorentz"):
        if raman == False:
            return TheoreticalSpectrum(self.freqlist, self.irintens, emin, emax, sigma, step, scale, mode="lorentz")
        else:
            return TheoreticalSpectrum(self.freqlist, self.ramanact, emin, emax, sigma, step, scale, mode="lorentz")

class TheoreticalSpectrum:
    
    def __init__(self, freqlist, intlist, emin, emax, sigma, step, scale, mode="lorentz"):
        self.freqlist = freqlist
        self.intensities = intlist
        self.sigma = sigma
        self.emax = emax
        self.emin = emin
        self.step = step
        self.mode = mode
        self.scale = scale
        self.nstep = int(round((self.emax - self.emin)/self.step)+1)
        self.temp = np.empty((self.nstep,len(self.freqlist)))
        self.frequencies = []
        
        if self.scale != 1.00:
            self.freqlist *= self.scale
            
        if self.emax < self.emin:
            raise Exception("The max frequency is lower than the min frequency")
            sys.exit()
        
        for i in range(self.nstep):
            self.frequencies.append(self.emin + i * self.step)
        self.frequencies = np.array(self.frequencies)
        
        if self.mode == "gauss":
            gauss = lambda params, x: (1.0/(params[1]*np.sqrt(2*np.pi))*np.exp(-0.5*((x-params[0])/params[1])**2))
        elif self.mode == "lorentz":
            lorentz  = lambda params, x: (1.0/(np.pi*params[1])*(params[1]**2)/((x-params[0])**2+params[1]**2))
            
        for i in range(len(self.freqlist)):
            params = [self.freqlist[i],sigma]
            if self.mode == "lorentz":
                self.temp[:,i] = lorentz(params,self.frequencies) * self.intensities[i]
            if self.mode == "gauss":
                self.temp[:,i] = gauss(params, self.frequencies) * self.intensities[i]
        
        self.intensities = np.sum(self.temp, axis=1)
        
    def resetscale(self):
        if self.scale != 1.00:
            self.freqlist /= self.scale
    
    def normalize(self):
        max_intensity = self.intensities.max()
        self.intensities /= max_intensity
    
class ExperimentalSpectrum(File):
    
    def __init__(self, location):
        File.__init__(self, location)
        self.frequencies = []
        self.intensities = []
        
        with open(location, "r") as expspectrum:
            self.spectrum = csv.reader(expspectrum)
            self.spectrum = list(self.spectrum)
        
        for i in range(len(self.spectrum)):
            self.frequencies.append(float(self.spectrum[i][0]))
            self.intensities.append(float(self.spectrum[i][1]))
        
        self.frequencies = np.array(self.frequencies)
        self.intensities = np.array(self.intensities)
    
    def normalize(self):
        max_intensity = self.intensities.max()
        self.intensities /= max_intensity
    
    def interpolate(self, thspectrum):
        frequencies_reversed = self.frequencies[::-1]
        intensities_reversed = self.intensities[::-1]
        self.interpolated_intensities = np.interp(thspectrum.frequencies, frequencies_reversed, intensities_reversed)

    def compare(self, thspectrum, mode="pearson"):
        self.interpolate(thspectrum)
        if mode == "pearson":
            p = sp.stats.pearsonr(self.interpolated_intensities, thspectrum.intensities)[0]
            thspectrum.resetscale()
            return p
        elif mode == "spearman":
            s = sp.stats.spearmanr(self.interpolated_intensities, thspectrum.intensities)[0]
            thspectrum.resetscale()
            return s
        
iranalyze = IRAnalyze()