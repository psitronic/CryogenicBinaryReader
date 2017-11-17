"""!/usr/bin/python
Created by A.A. Sidorenko aka spintronic

Class and functions to extract wave-forms V(x) from the binary data files .dat
created by a SQUID magnetometer (Cryogenic Ltd.)

It is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

It is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with CryogenicBinaryReader.py. If not, see <http://www.gnu.org/licenses/>.

"""
from struct import unpack
import numpy as np
from fitting import *

class DataS700x():
    """
    Class to extract wave-froms V(x) from binary .dat files
    """
    def __init__(self,file):
        self.endian = '>' # Big endian for LabView
        self.index = 4
        self.Nr = 0
        self.Nc = 0
        infile = open(file,"rb")
        self.data = infile.read()
        infile.close()
        self.comments = ''
        self.sequenceInfo = ''
        self.acqisitionInfo = ''
        self.sampleInfo = float()
        self.calibrationInfo = float()
        self.recordInfoBefore = float()
        self.recordInfoAfter = float()
        self.DCAnalysis = float()
        self.scansDC = float()
        self.nmeasurements = 0
    
    def getArraySize(self, offset):
        datatype_str = self.endian + str(self.index/2) + 'l'
        Nrow, Ncol = unpack(datatype_str,self.data[0+offset:2*self.index+offset])
        return Nrow, Ncol
        
    def getArray(self,offset,dt_str):
        self.Nr, self.Nc = self.getArraySize(offset)
        if self.Nr == 0:
            IndexSize = self.Nc
        else:
            IndexSize = self.Nc * self.Nr
        StartPos = 2*self.index + offset
        StopPos = self.index * (IndexSize + 2) + offset
        datatype_str = self.endian + str(IndexSize) + dt_str
        Index = unpack(datatype_str, self.data[StartPos:StopPos])        
        return Index
        
    def nonZeroArray(self,offset,dt_str):
        Nrecord = 0
        Index = self.getArray(offset,dt_str)
        for ind in Index:
            if ind==0:
                break
            else:
                Nrecord+=1
        nonZeroIndex = Index[0:Nrecord]
        return nonZeroIndex, Nrecord
        
    def getPrimaryIndex(self):
        dt_str = 'l'
        data = self.nonZeroArray(0,dt_str)[:][:]
        pindex = data[0]
        nm = data[1]
        return pindex, nm
        
    def getSecondaryIndex(self,Pindex):
        dt_str = 'l'
        return self.nonZeroArray(Pindex, dt_str)[0]
        
    def getStringArray(self,position,count):
        dtype = self.endian + str(count) + 'c'
        return unpack(dtype,self.data[position:position+count])
   
    def getDataArray(self,position):
        dt_str = 'f'
        return self.getArray(position,dt_str)   
        
    def convert1Dto2D(self,scans):
        Nrows, Ncols = self.Nr, self.Nc
        slices = np.zeros(Nrows*Ncols).reshape(Nrows,Ncols)
        for i in range(0,Nrows):
            slices[i,:] = scans[i * Ncols: (i + 1) * Ncols]
        return slices
                
    def getData(self, npoint):
        pindexArray, self.nmeasurements = self.getPrimaryIndex()
        pindex = pindexArray[npoint]
        sindex = self.getSecondaryIndex(pindex)
        nextP = pindexArray[npoint+1] #
        max_index = len(self.getSecondaryIndex(pindex))
        DCscans1DArray = []
                
        for i in range(0,max_index/2):
            shift = 2 * i
            position = sindex[shift]
            nextElement = sindex[shift+1]
            representation = 255 & nextElement
            datatype = nextElement >> 8            
           
            if (shift+2) < max_index:
                count = sindex[shift+2] - position
            else:
                count = nextP - position
              
            if npoint == 0:
                if   datatype == 1 and representation == 0:
                    self.comments = self.getStringArray(position,count)
                elif datatype == 1 and representation == 1:
                    self.sampleInfo = self.getDataArray(position)
                elif datatype == 2 and representation == 1:
                    self.calibrationInfo = self.getDataArray(position)
                elif datatype == 3 and representation == 0:
                    self.sequenceInfo = self.getStringArray(position,count)
            else:
                if   datatype == 4 and representation == 1:
                    self.recordInfoBefore = self.getDataArray(position)
                elif datatype == 5 and representation == 1:
                    self.recordInfoAfter = self.getDataArray(position)
                elif datatype == 6 and representation == 2:
                    DCscans1DArray = self.getDataArray(position)
                    self.scansDC = self.convert1Dto2D(DCscans1DArray)                    
                elif datatype == 8 and representation == 1:
                    self.DCAnalysis = self.getDataArray(position)
                elif datatype == 6 and representation == 0:
                    self.acqisitionInfo = self.getStringArray(position,count)
                    
def averageScans(pts):
    avrScansUp = 0
    avrScansDown = 0
    positionUp = []
    positionDown = []
    position = []
    avrScansUp = []
    avrScansDown = []
    avrScans = []
    nscans = len(pts[0]) -1        
    npoints = len(pts)
    for i in range(0,npoints/2):
        positionUp = np.append(positionUp,pts[i,0])
        positionDown = np.append(positionDown,pts[npoints-1-i,0])
        sumScansUp = 0
        sumScansDown = 0
        for j in range(1,nscans+1): # to start from scan #2 change the beginning 1 to 2 and in avrScansUp and avrScansDown nscans to (nscans-1)
            sumScansUp += pts[i,j]
            sumScansDown +=  pts[npoints-1-i,j]
        avrScansUp = np.append(avrScansUp,sumScansUp/(nscans-0))
        avrScansDown = np.append(avrScansDown,sumScansDown/(nscans-0))
        
    avrScans = (avrScansUp + avrScansDown)/2
    position = (positionUp + positionDown)/2
    return position, avrScans

                     
def fit(x_data,y_data):
    """
    Function to fit the wave form V(x) measured by SQUID
    Takes two arrays: position and voltage
    Returns: position array, function value for optimal fitting parameters, voltage at the sample position
    """
    A1 = 0.
    A2 = 0.
    A3 = 0.
    A4 = 0.

    p0 = [A1, A2, A3, A4]

    f = lambda x_data, A1, A2, A3, A4: A1+A2*x_data+\
        A3*(2.*pow(1.25**2 + (x_data + A4)**2,-3./2.) - pow(1.25**2 +\
        (0.70 + (x_data + A4))**2,-3./2.) - pow(1.25**2 + (-0.70 + (x_data + A4))**2,-3./2.))

    popt, punc, rc, d = general_fit(f, x_data, y_data, p0)

    A1, A2, A3, A4 = popt


    print('optimal parameters: ', popt)
    print('uncertainties of parameters: ')
    A1err, A2err, A3err, A4err = punc

    print('A1 = %12.7e +/- %7.4e' % (A1, A1err))
    print('A2 = %10.4e +/- %7.4e' % (A2, A2err))
    print('A3 = %10.4e +/- %7.4e' % (A3, A3err))
    print('A4 = %10.4e +/- %7.4e' % (A4, A4err))

    return x_data, f(x_data,A1, A2, A3, A4), A3                    

       
dataSample = DataS700x("") # Sample and background data file
dataBckg = DataS700x("") # Background data file
dataSample.getData(0)
dataBckg.getData(0)

nmax = dataSample.nmeasurements-2

ofile = open("") # Save final data to file

for i in range(2,nmax+1):
    
    dataSample.getData(i)
    pos_sample, value_sample = averageScans(dataSample.scansDC)
    
    dataBckg.getData(i)
    pos_bckg, value_bckg = averageScans(dataBckg.scansDC)
    
    x, y = pos_sample, value_sample - value_bckg # to subtract background replace 0 with value_bckg and vice versa
    
    moment = 0.3326379003 * fit(x,y)[2] * -dataSample.calibrationInfo[0] # The prefactor is from calibration
    
    ofile.write('%12.5f %12.5f %12.5e\n' % ((dataSample.recordInfoBefore[3]+dataSample.recordInfoAfter[3])/2.,dataSample.recordInfoBefore[5],moment))
    
ofile.close()

