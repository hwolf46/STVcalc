# -*- coding: utf-8 -*-
"""
Created on Fri Sep 15 16:30:21 2017

@author: HWolf
"""
import sys
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui
import numpy as np
import pymsgbox

#import pyqtgraph.parametertree.parameterTypes as pTypes
#from pyqtgraph.parametertree import Parameter, ParameterTree, ParameterItem, registerParameterType
mApp = QtGui.QApplication([])

def ReadBinary():
#read data from binary MOSOS file with .ctg (new version)
    global mFileName
    global mXdata
    global mFHR, mQual
    global mUtP
    global mIntV
    try:
        File = open(mFileName, 'br')
    except FileNotFoundError:
        pymsgbox.alert('File not opened. Program will shut down.')
        raise SystemExit
    mTel = 0
    mB = 0
    mSingle = -1
    mData = []
    mFHR = []
    mXdata = []
    mIntV = []
    mQual = []
    mUtP = []
    with open(mFileName, 'br') as File:
        byte_content = File.read()
#        mQ = byte_content
        mData = [byte_content[i + 1] << 8 | byte_content[i] for i in range(0, len(byte_content), 2)]
    File.close()
    mTel = 0
    while mTel < len(mData) and mSingle < 0:
        mB = bin(mData[mTel])
        mB1 = '0b' + ('0000000000000000' + mB[2:])[-16:]
        mB = bin(mData[mTel + 4])
        mB2 = '0b' + ('0000000000000000' + mB[2:])[-16:]
        if (int(mB1[-10:-2], 2)) > 0:
            mSingle = 0
        elif (int(mB2[-10:-2], 2)) > 0:
            mSingle = 4
        mTel += 12
    mTel -= 12
    while mTel < len(mData):
        mB = bin(mData[mTel + mSingle])
        mB = '0b' + ('0000000000000000' + mB[2:])[-16:]
        mFHR.append(int(mB[-10:-2], 2))
        mQual.append(int(mB[-13:-10], 2))
        mB = bin(mData[mTel + 1  + mSingle])
        mB = '0b' + ('0000000000000000' + mB[2:])[-16:]
        mFHR.append(int(mB[-10:-2], 2))
        mQual.append(int(mB[-13:-10], 2))
        mB = bin(mData[mTel + 2  + mSingle])
        mB = '0b' + ('0000000000000000' + mB[2:])[-16:]
        mFHR.append(int(mB[-10:-2], 2))
        mQual.append(int(mB[-13:-10], 2))
        mB = bin(mData[mTel + 3  + mSingle])
        mB = '0b' + ('0000000000000000' + mB[2:])[-16:]
        mFHR.append(int(mB[-10:-2], 2))
        mQual.append(int(mB[-13:-10], 2))
# uterus tonus
        mB = bin(mData[mTel+10])
        mB = '0b' + ('0000000000000000' + mB[2:])[-16:]
        mUtP.append(int(mB[-8:], 2))
        mUtP.append(int(mB[-8:], 2))
        mUtP.append(int(mB[-16:-8], 2))
        mUtP.append(int(mB[-16:-8], 2))
        mTel += 12
    mTel = 0
    while mFHR[mTel] == 0:
        del mFHR[mTel]
        del mUtP[mTel]
        mTel += 1
    mTel = 0
    while mTel < len(mFHR):
        mXdata.append(mTel/240)
        mTel += 1
    mTel = 0
    while mTel < len(mFHR):
        if mFHR[mTel] > 0:
            if mFHR[mTel] < 30 or mFHR[mTel] > 200:
                mFHR[mTel] = 0
                mIntV.append(0)
            else:
                mIntV.append(60000/mFHR[mTel])
        else:
            mIntV.append(0)
        mTel += 1
    RejectData()
    EpochCalc()
    MakePlot1()

def ReadDataFile():
#read daat from prn file (heart rate values in single column at 4 Hz)
    global mFHR, mQual
#    global mUtP
    global mIntV
    global mXdata
    global mFileName
    mTel = 0.00
    mFr = 0
    mI = 0
    mQ = 0
    mU = 0
    mXdata = []
    mFHR = []
    mQual = []
    mUtP = []
    mIntV = []
    mLLine = []
    try:
        f = open(mFileName, 'r')
    except FileNotFoundError:
        pymsgbox.alert('File not opened. Program will shut down.')
        raise SystemExit
    for line in f:
        i = 0
        b = 0
        mLLine = []
        mLine = line
        while i < len(mLine):
            if mLine[i] == ';':
                mLLine.append(int(mLine[b:i]))
                b = i + 1
            i += 1
            if i == len(mLine) and i > 1:
                mLLine.append(int(mLine[b:]))
        mFr = int(mLLine[0])
        if len(mLLine) > 1:
            mQ = int(mLLine[1])
        if len(mLLine) > 2:
            mU = int(mLLine[2])
        if (mFr < 30) or (mFr > 200): mFr = 0
        if mFr == 0:
            mI = 0
            mQ = 0
        else:
            mI = 60000/mFr
        mFHR.append(mFr)
        mQual.append(mQ)
        mUtP.append(mU)
        mIntV.append(mI)
        mXdata.append(mTel)
        mTel += 0.25/60
    f.close()
    RejectData()
    EpochCalc()
    MakePlot1()
        
def RejectData():
    # rejection algorithm for data very far from meean
    global mIntV, mFHR
    mI1 = 0
    mI2 = 0
    mI3 = 0
    mI4 = 0
    mTel = 0
    mIMean = 0
    mIList = []
    mIMeanTot = np.median(list(filter((0).__ne__, mIntV)))
    mTel = 0
    while mTel < len(mIntV):
        mI1 = mI2
        mI2 = mI3
        mI3 = mI4
        mI4 = mIntV[mTel]
        if mI4 > 0:
            if mI1 > 0 or mI2 > 0 or mI3 > 0:
                mIList = [mI1, mI2, mI3]
                mIMean = np.mean(list(filter((0).__ne__, mIList)))
            else:
                mIMean = mIMeanTot
            if (mI4 <= 0.6 * mIMean) or (mI4 >= 1.55 * mIMean):
                mI4 = 0
                mIntV[mTel] = 0
                mFHR[mTel] = 0
        mTel += 1
    mTel = 0

def EpochCalc():
#epoch calculation
    global mIntV, mIntEp, mFHRcode, mFHRXep, mFHRep
    global mVersion, Lbl1
    mIntEp = []
    mFHRep = []
    mFHRcode = []
    mFHRXep = []
    mTel = 0
    mCode = 0
    mIntVmean = 0
    while mTel < len(mIntV):
        if np.bincount(mIntV[mTel:mTel+15])[0] < 15 and len(mIntV[mTel:mTel+15]) == 15:
            mIntVmean = np.mean(list(filter((0).__ne__, mIntV[mTel:mTel+15])))
            mCode = 1
        else:
            mIntVmean = 0
            mCode = 9
        mIntEp.append(mIntVmean)
        if mIntVmean > 0:
            mFHRep.append(60000/mIntVmean)
        else:
            mFHRep.append(0)
        mFHRcode.append(mCode)
        mFHRXep.append(mTel/(15*16))
        mTel += 15

def CalcBaseline():
    global mFHRep, mIntEp
    global mFHRbase
    global mFHRXep
    global mSTVmin
    mSTVmin = []
    mDist = []
    mRefHRe = 0
    mRefTe = 0
    mStartHRe = 0
    mStartTe = 0
# Calculation reference point according to Dawes
    mTel = 90
    mLimit = 0
    mCount = 0
    mDist = np.bincount(mFHRep)
    mMaxDist = np.argmax(mDist[90:]) + 90
    while mTel < len(mDist) - 5:
        mLimit += mDist[mTel]
        if mLimit > np.sum(mDist) / 6 :
            if (mDist[mTel] > mDist[mTel+1] and mDist[mTel] > mDist[mTel+2] and
                    mDist[mTel] > mDist[mTel+3] and mDist[mTel] > mDist[mTel+4] and mDist[mTel] > mDist[mTel+5]):
                mRefHRe = mTel
                mCount = mDist[mTel]
                mTel = 200
        mTel += 1
#    print(mCount, np.sum(mDist))
    if mRefHRe == 0:
        mRefHRe = mMaxDist
    elif np.absolute(60000/mRefHRe - 60000/mMaxDist) > 30 and (mCount / np.sum(mDist)) < 0.02:
            mRefHRe = mMaxDist
    if mRefHRe > 0:
        mRefTe = 60000 / mRefHRe        
# Calculation starting point
    mTel = 0
    mLimit = 10
    mTe = 0
    i = 0
    while mLimit <= 40:
        for i in mIntEp[0:64]:
            if (i > mRefTe - mLimit) and (i < mRefTe + mLimit):
                mTe += i
                mTel += 1
        if mTel > 10 and mTe > 0:
            mStartTe = mTe / mTel
            mStartHRe = 60000 / mStartTe
            mLimit = 50
        mLimit += 10
    if mStartHRe == 0:
        mStartHRe = mRefHRe
        mStartTe = mRefTe
# write and smooth baseline first pass
    mLimit = 40
    BaselineWrite(mStartTe, mRefTe, mStartHRe, mLimit)
# check for baseline difference of more than 10 minutes of more than 150 msec is missing yet
    mTel = 0
    i1 = 0
    i2 = 0
    d = 0
    while d < 6:
        while mTel < len(mFHRbase):
            if mFHRep[mTel] > mFHRbase[mTel]:
                i1 += 1
            elif mFHRep[mTel] < mFHRbase[mTel]:
                i2 += 1
            else:
                i1 = 0
                i2 = 0
            if i1 > 160 or i2 > 160:
                mTel = len(mFHRbase)
                mLimit += 20
                BaselineWrite(mStartTe, mRefTe, mStartHRe, mLimit)
            else:
                d = 6
            mTel += 1  
        d += 1

def BaselineWrite(mStartTe, mRefTe, mStartHRe, mLimit):
    global mFHRbase, mTbase
    mTel = 0
    mMeanIntRec = mStartTe
    mFHRbase = []
    mFHRbaseR = []
    mTbase = []
# Set invalid epochs to mStartHre
    while mTel < len(mIntEp):
        if mTel >= 8:
            mMeanIntRec = np.mean(mTbase[mTel-8:mTel-1])
        if mTel < 8 and (mIntEp[mTel] == 0 or mIntEp[mTel] < (mStartTe - mLimit) or mIntEp[mTel] > (mStartTe + mLimit)):
            mTbase.append(mStartTe)
            mFHRbase.append(mStartHRe)
        elif mIntEp[mTel] == 0 or mIntEp[mTel] < (mMeanIntRec - mLimit) or mIntEp[mTel] > (mMeanIntRec + mLimit):
            mTbase.append(mMeanIntRec)
            mFHRbase.append(60000/mMeanIntRec)
        elif mIntEp[mTel] == 0 or mIntEp[mTel] < (mRefTe - mLimit) or mIntEp[mTel] > (mRefTe + mLimit):
            mTbase.append(mMeanIntRec)
            mFHRbase.append(60000/mMeanIntRec)            
        else:
            mTbase.append(mIntEp[mTel])
            mFHRbase.append(mFHRep[mTel])
        mTel += 1
# Smoothing algorithm FHRN baseline
    mTel = 1
    mFHRbaseR.append(mFHRbase[0])
    while mTel < len(mFHRbase):
        mFHRbaseR.append(mFHRbase[mTel]*0.05 + mFHRbaseR[mTel-1]*0.95)
        mTel += 1
    mTel = len(mFHRbaseR)-1
    mFHRbase[mTel] = mFHRbaseR[mTel]
    mTel -= 1
    while mTel >= 0:
        mFHRbase[mTel] = mFHRbaseR[mTel] * 0.05 + mFHRbase[mTel + 1] * 0.95
        mTel -= 1

def CalcSTV():
    global mFHRX, mFHRbase, mFHRep, mIntEp, mFHRcode, mAccelLoc, mDecelLoc
    global mReport, mExport, mFileName, mReportFileName, mVersion
    global mSTVmin, mLTVmin
    CalcBaseline()
    mSTVmin = []
    mSTVX = []
    mAccelLoc = []
    mDecelLoc = []
    mSTVtot = 0
    mSTV = 0
    mTel = 0
    mFHRmean = 0
# Coderen hartfrequenties 10 of 20 slagen boven / onder basislijn of 75ms van basis interval
    mTel = 0
    while mTel < len(mFHRep):
        if mIntEp[mTel] < (60000/mFHRbase[mTel]) - 75 and mFHRep[mTel] > 0:
            mFHRcode[mTel] = 180
        elif mFHRep[mTel] > mFHRbase[mTel] + 10:
            mFHRcode[mTel] = 175
        elif mIntEp[mTel] > (60000/mFHRbase[mTel]) + 75:
            mFHRcode[mTel] = 165
        elif mFHRep[mTel] <= mFHRbase[mTel] - 20 and mFHRep[mTel] > 0:
            mFHRcode[mTel] = 167
        elif mFHRep[mTel] <= mFHRbase[mTel] - 10 and mFHRep[mTel] > 0:
            mFHRcode[mTel] = 170
        mTel += 1            
# Terugzetten code van FHR meer of minder dan <20 van baseline die niet in deceleratie zitten en
#        alle accelleraties.
# Exclusie blijft alle deceleraties en van mFHRep meer dan +-75 msec van baseline, die niet in accel. of decel. zit,
# dit is code 172 (decelleratie), en 165 en 180 (uitbijters)
    mTel = 0
    while mTel < len(mFHRcode):
#        if mFHRcode[mTel] == 170 or mFHRcode[mTel] == 175:
#        if mFHRcode[mTel] == 165 or mFHRcode[mTel] == 170 or mFHRcode[mTel] == 175 or mFHRcode[mTel] == 180:
        if mFHRcode[mTel] == 167 or mFHRcode[mTel] == 170 or mFHRcode[mTel] == 175 or mFHRcode[mTel] == 177:
            mFHRcode[mTel] = 1
        mTel += 1
# Naastliggende waarden van uitbijters ook excluderen
    mTel = 1
    while mTel < len(mFHRcode) - 1:
        if mFHRcode[mTel] == 165 or mFHRcode[mTel] == 180:
            mFHRcode[mTel-1] = mFHRcode[mTel]
            if not (mFHRcode[mTel+1] == 165 or mFHRcode[mTel+1] == 180):
                mFHRcode[mTel+1] = mFHRcode[mTel]
                mTel += 1
        mTel += 1        
# Tellen accelleraties
    mTel = 0
    while mTel < len(mFHRep):
        mAccelLoc.append(0)
        mDecelLoc.append(0)
        mTel += 1
    mAccell = 0
    mAccell60 = 0    
    mMax = 0
    mTel = 0
    i = 0
    while mTel < len(mFHRep):
        if mFHRep[mTel] >= mFHRbase[mTel] and mFHRcode[mTel] != 165:
            while mTel < len(mFHRep) and mFHRep[mTel] > mFHRbase[mTel] and mFHRcode[mTel] != 165:
                if mMax < mFHRep[mTel] - mFHRbase[mTel]:
                    mMax = mFHRep[mTel] - mFHRbase[mTel]
                i += 1
                mTel += 1
        else:
            mMax = 0
            i = 0
            mTel += 1
        if i > 4  and mMax > 10:
            mAccell += 1
            if mTel <= 960:
                mAccell60 += 1
            while i > 0:
                mAccelLoc[mTel - i] = 177
                if mFHRcode[mTel - i] == 180:
                    mFHRcode[mTel - i] = 1
                i -= 1
            i = 0
# Tellen decelleraties
    mDecell = 0
    mDecell60 = 0
    mTel = 0
    i = 0
    mMin = 0
    while mTel < len(mFHRep):
        if (mFHRep[mTel] < mFHRbase[mTel]) and (mFHRep[mTel] > 10):
            while mTel < len(mFHRep) and mFHRep[mTel] < mFHRbase[mTel]:
                if mFHRep[mTel] > 10:
                    if mMin < mFHRbase[mTel] - mFHRep[mTel]:
                        mMin = mFHRbase[mTel] - mFHRep[mTel]
                    i += 1
                mTel += 1
        else:
            mMin = 0
            i = 0
            mTel += 1
        if (i > 15  and mMin > 10) or (i > 7 and mMin > 20):
            mDecell += 1
            if mTel <= 960:
                mDecell60 += 1
            while i > 0:
                mFHRcode[mTel - i] = 172
                mDecelLoc[mTel - i] = 172
                i -= 1
            i = 0
# Berekenen STV over minuten met alle epochs valide (mFHRcode = 1), dus zonder uitbijters of decelleratie
# Alleen complete minuten
#    mTel = 0
#    i = 0
#    while mTel < len(mFHRep)-16:
#        if np.bincount(mFHRcode[mTel:mTel+16])[1] < 16:
#            mSTVmin.append(99)
#        else:
#            i = 1
#            mTe = 0
#            while i < 16:
#                mTe += np.absolute(mIntEp[mTel + i] - mIntEp[mTel + i - 1])
#                i += 1
#            mSTVmin.append(mTe/15)
#        mSTVX.append(mTel/16)
#        mTel += 16
# Berekenen STV over minuten met < 50% dataverlies en geen deceleratie (begin of einde),
# zonder uitbijters, verschil naast elkaar gelegen paren.
    mTel = 0
    while mTel < len(mFHRep)-16:
        i = 0
        mTe = 0
        n = 1
#        d = 1
        while n < 16 and mTe < 9000:
            if ((mFHRcode[mTel + n] == 172 or mFHRcode[mTel + n - 1] == 172) or 
                (mFHRcode[mTel + n] == 180 or mFHRcode[mTel + n - 1] == 180) or (mFHRcode[mTel + n] == 165 or mFHRcode[mTel + n - 1] == 165)):
                mTe = 9999
            elif mFHRcode[mTel + n] == 1 and mFHRcode[mTel + n - 1] == 1:
                mTe += np.absolute(mIntEp[mTel + n] - mIntEp[mTel + n - 1])
                i += 1
#            elif mFHRcode[mTel + n] == 1:
#                d = 1
#                while n - d > 1:
#                    if mFHRcode[mTel + n - d] == 1:
#                        mTe += np.absolute(mIntEp[mTel + n] - mIntEp[mTel + n - d])
#                        i += 1
#                    d += 1
            n += 1
        if i >= 8 and mTe > 0 and mTe < 9000:
            mSTVmin.append(mTe / i)
        else:
            mSTVmin.append(99)
        mSTVX.append(mTel/16)
        mTel += 16
# Berekenen van Longterm variation
    mTel = 0
    mRangeMin = []
    while mTel < len(mIntEp) - 16:
        n = 0
        i = 0
        mMax = 0
        mMin = 1000
        while n < 16:
            if mFHRcode[mTel + n] == 1:
                if mMax < mIntEp[mTel+n]:
                    mMax = mIntEp[mTel+n]
                if mMin > mIntEp[mTel+n]:
                    mMin = mIntEp[mTel+n]
                    if mMin > mTbase[mTel + n]:
                        mMin = mTbase[mTel + n]
                i += 1
            n += 1
        if i > 8:
            mRangeMin.append(mMax - mMin)
        else:
            mRangeMin.append(999)
        mTel += 16
    mTel = 0
    mLTVmin = []
    while mTel < len(mRangeMin):
        mLTVmin.append(0)
        mTel += 1
    mTel = 0
    n = 0
    i = 0
    while mTel < len(mRangeMin) - 6:
        n = 0
        i = 0
        while n < 6:
            if mRangeMin[mTel + n] > 31 and mRangeMin[mTel + n] < 900:
                i += 1
            n += 1
        if i >= 5:
            n = 0
            while n < 6:
                mLTVmin[mTel + n] = 175
                n += 1
        mTel += 1
    mTel = 0
    i = 0
    while mTel < len(mLTVmin):
        if mLTVmin[mTel] == 175:
            i += 1
        mTel += 1
    mHighLTVtime = i
    mLTV = np.mean(list(filter((999).__ne__, mRangeMin)))
# Berekenen FHR baseline

    i = 0
    mTel = 0
    while mTel < len(mFHRbase):
        if mFHRbase[mTel] > 0:
            i += 1
        mTel += 1
    mFHRmean = sum(mFHRbase) / i
# Berekenen gemiddelde minuut-STV over hele registratie
    mTel = 0
    i = 0
    i2 = 0
    mSTV = 0
    mSTVtot = 0
    mSTV60 = 0
    mSTV60tot = 0
    mSTVmin[0] = 99
    while mTel < len(mSTVmin):
        if mSTVmin[mTel] > 0 and mSTVmin[mTel] < 90:
            if mTel <= 60:
                mSTV60tot = mSTV60tot + mSTVmin[mTel]
                i2 += 1
            mSTVtot = mSTVtot + mSTVmin[mTel]
            i += 1
        mTel += 1
    mSTV = mSTVtot / i
    if i2 > 0:
        mSTV60 = mSTV60tot / i2
# Rapportage
    mDuration = round(len(mFHRep)/16)
    mReport = ('Base FHR : ' + str(round(mFHRmean)) + ', STV : ' + str(round(mSTV, 2))
               + ', STV 60min: ' + str(round(mSTV60, 2)) + ', LTV: '  + str(round(mLTV, 2)) +
               ', High LTV minutes:' + str(mHighLTVtime) + ', Duration CTG : ' +
               str(mDuration) + ' minutes, Valid minutes : ' + str(i) + ' Lost minutes : ' + str(round(((mDuration - i)/mDuration) * 100)) +
               '%, \n' + 'Accelleraties all : ' + str(mAccell) + ', <60 min.' + str(mAccell60) + 
               ', Decelleraties all : ' + str(mDecell) + ', <60 min.' + str(mDecell60) + ', Type : ' +
               str(CBox._chosenText))
    mExport = (mFileName + '; ' + str(round(mFHRmean)) + '; ' + str(round(mSTV, 2)) + '; ' + 
               str(round(mSTV60, 2)) + '; ' + str(round(mLTV, 2)) + '; ' + str(mHighLTVtime) + ";"  +
               str(round(mDuration)) + '; ' + str(i) + '; ' + str(mAccell) + '; ' + str(mAccell60) + '; ' + 
               str(mDecell) + '; ' + str(mDecell60) + '; ' + str(CBox._chosenText))
    mF = open(mReportFileName, 'a')
    mF.write(mExport + '\n')
    mF.close    
#    mF = open(mFileName[0:-3] + 'epo', 'w')
#    mTel = 0
#    while mTel < len(mFHRep):
#        mF.write(str(round(mFHRep[mTel], 2)) + '\n')
#        mTel += 1
#    mF.close
    if mVersion == 1:
        MakePlot1()

def GetFile():
    global mFileName, mFolderName
    global mReport, mRange
    global mSTVmin, mFHRbase, mFHRcode, mXdata, mFHR, mFHRX, mAccelLoc, mDecelLoc
    mRange = []
    mSTVmin = []
    mFHRbase = []
    mFHRcode = []
    mXdata = []
    mFHR = []
    mFHRX = []
    mAccelLoc = []
    mDecelLoc = []
    mReport = 'Report'
    try:
        mFileName = (QtGui.QFileDialog.getOpenFileName
                     (None, 'Open file for STV calculation', 
                      mFolderName, '(*.prn *.ctg *.txt)'))
    except FileNotFoundError:
        pymsgbox.alert(text='You did not read a CTG data file. \nThe program will cancel.', 
                       title='Pas op', button='Gezien')
        raise SystemExit
    mFileName = mFileName[0]
    if (mFileName[-3:] == 'prn') or (mFileName[-3:] == 'txt'):
        ReadDataFile()
    elif mFileName[-3:] == 'ctg':
        ReadBinary()
    else:
        pymsgbox.alert(text='You did not read a CTG data file. \nOther formats are not acceptable. Try again.', 
                       title='Pas op', button='Gezien')


# functions for individual CTG and plot
def SelectOut():
    global mXdata, mFHR, mIntV, mSTVmin, mMouseAction
    global mRange
#    CalcRange()
    mMouseAction = 0
    Lbl1.setText("Not calculated")
    vb.setMouseMode(vb.PanMode)
    mTel = 0
    mSt = mRange[0][0]
    mEnd = mRange[0][1]
    while mTel < len(mXdata):
        if mXdata[mTel] >= mSt and mXdata[mTel] < mEnd:
            mFHR[mTel] = 0
            mIntV[mTel] = 0
        mTel += 1
    mSTVmin = []
    EpochCalc()
    CalcSTV()

def SelectAll():
    global mFHRbase, mSTVmin
    global mFileName
    mReport = "Not calculated"
    Lbl1.setText(mReport)
    mSTVmin = []
    mFHRbase = []
    if mFileName[-3:] == 'prn' or mFileName[-3:] == 'txt':
        ReadDataFile()
    elif mFileName[-3:] == 'ctg':
        ReadBinary()
    CalcSTV()

def MakePlot1():
## Create plotwindow
    global mSTVmin
    global mFileName, mReport
    global mFHRep, mFHRXep, mFHRbase, mFHRcode, mAccelLoc, mDecelLoc, mLTVmin
#    global mUtP, mXdata
    global pl1, pl2
    mTel = 0
    while mTel < len(mFHRcode):
        if mFHRcode[mTel] == 172 and mDecelLoc[mTel] == 172:
            mFHRcode[mTel] = 0
        mTel += 1    
    mFHRepN = []
    mTel = 0
    while mTel < len(mFHRep):
        if mFHRep[mTel] > 0:
            mFHRepN.append(mFHRep[mTel])
        else:
            mFHRepN.append(np.nan)
        mTel += 1
#    mFHRN = []
#    mTel = 0
#    while mTel < len(mFHR):
#        if mFHR[mTel] > 0:
#            mFHRN.append(mFHR[mTel])
#        else:
#            mFHRN.append(np.nan)
#        mTel += 1
    font = QtGui.QFont()
    font.setPixelSize(20)
    plt.clear()
    pl2.clear()
    pl1.getAxis('right').setLabel('STV', color='#0000ff', font='20pt')
    pl1.getAxis('left').setLabel('FHR', color='#0000ff', font='20pt')
    pl1.plot(mFHRXep, mFHRepN, connect="finite", pen='k', color='b')               
#    pl1.plot(mFHRXep, mFHRepN, connect="finite", pen='k', color='b', symbol='o', symbolSize=3, symbolBrush='b')               
#    pl1.plot(mXdata, mFHRN, connect="finite", pen='k', color='b')
    pl1.setLabel('bottom', text='Minutes', color='#0000ff', font='20pt')
    pl1.setYRange(80, 200)
    pl1.getAxis('bottom').tickFont = font
    pl1.getAxis('left').tickFont = font
    pl1.getAxis('right').linkToView(pl2)
    pl1.getAxis('right').tickFont = font
    if np.sum(mSTVmin) > 0:
        pl2.addItem((pg.ScatterPlotItem(range(0, len(mSTVmin)), mSTVmin,
                    symbol='o', brush='b', color='b', size=12)))
        pl2.setYRange(0, 40)
    if np.sum(mFHRcode) > 0:
        pl1.plot(mFHRXep, mFHRcode, symbol='s', pen=None, symbolBrush='r', symbolPen='r', symbolSize=10)
    if np.sum(mFHRbase) > 0:
        pl1.plot(mFHRXep, mFHRbase, pen='k')
#    if len(mUtP) > 0:
#        pl1.plot(mXdata, mUtP, pen='k')
    if np.sum(mDecelLoc) > 0:
        pl1.plot(mFHRXep, mDecelLoc, pen=None, symbol='s', symbolBrush='m', symbolPen=None, symbolSize=12)        
    if np.sum(mAccelLoc) > 0:
        pl1.plot(mFHRXep, mAccelLoc, pen=None, symbol='s', symbolPen='g', symbolBrush='g', symbolSize=12)
    if np.sum(mLTVmin) > 0:
        pl1.plot(range(0, len(mLTVmin)), mLTVmin, pen=None, symbol ='s', symbolPen=(0,128,0), symbolSize=12)
    plt.showGrid(x=True, y=True, alpha=1)
    plt.setXRange(0, 60)
    plt.setBackground('#e6e6e6')
    plt.setTitle(' CTG file  ' + mFileName)
    Lbl1.setText(mReport)
    updateViews()

def updateViews():
    ## view has resized; update auxiliary views to match
    global pl1, pl2
    pl2.setGeometry(pl1.vb.sceneBoundingRect())
    pl2.linkedViewChanged(pl1.vb, pl2.XAxis)

class CustomViewBox(pg.ViewBox):
    def __init__(self, *args, **kwds):
        pg.ViewBox.__init__(self, *args, **kwds)
    def mouseClickEvent(self, ev):
        global mMouseAction
 #       if ev.button() == QtCore.Qt.RightButton and ev.button() == QtCore.Qt.Double:
        if ev.button() == QtCore.Qt.RightButton:
            self.setMouseMode(self.PanMode)
            mMouseAction = 0
        if ev.double():
            mMouseAction = 1
            self.setMouseMode(self.RectMode)
    def mouseDragEvent(self, ev):
        global mMouseAction, mRange
        if ev.button() == QtCore.Qt.RightButton:
            ev.ignore()
        elif mMouseAction == 1:
            pg.ViewBox.mouseDragEvent(self, ev)
            ev.accept()
            mRange = pl1.viewRange()
            if ev.isFinish():
                SelectOut()

#Start of program.......................

mFileName = ''
mFolderName = ''
mOutFileName = ''
mReportFileName = ''
mReport = 'Report'
mRange = []
mSTVmin = []
mLTVmin = []
mFHRbase = []
mTbase = []
mFHRcode = []
mFHRX = []
mXdata = []
mFHR = []
mQual = []
mUtP = []
mFHRep = []
mIntEp = []
mFHRXep = []
mIntV = []
mIntEp = []
mAccelLoc = []
mDecelLoc = []
mMouseAction = 0
mVersion = 1

vb = CustomViewBox()
#Btn1 = QtGui.QPushButton(text='  Select to remove  ')
#Btn1.setToolTip('With this button you can delete the marked CTG part')
#Btn1.clicked.connect(SelectOut)
Btn2 = QtGui.QPushButton(text='  Restore complete registration  ')
Btn2.clicked.connect(SelectAll)
Btn3 = QtGui.QPushButton(text='  Get new file  ')
Btn3.clicked.connect(GetFile)
#Btn4 = QtGui.QPushButton(text='  Calculate baseline  ')
#Btn4.clicked.connect(CalcBaseline)
Btn5 = QtGui.QPushButton(text='  Calculate STV  ')
Btn5.clicked.connect(CalcSTV)
CBox = pg.ComboBox(items={'No decellerations' :1, 'Variable decelerations' : 2, 'Recurrent decelerations' : 3})
CBox.setValue(1)
CBox.currentIndexChanged.connect(CalcSTV)
Lbl1 = QtGui.QLabel("""Results""")
Lbl1.setText("Report")
Lbl2 = QtGui.QLabel("""Informatie""")
Lbl2.setText("Double-click with left mouse button to select and remove a part of the CTG for analysis.")
plt = pg.PlotWidget(viewBox=vb)
pl1 = plt.plotItem
pl1.showAxis('right')
pl1.setMouseEnabled(x=True, y=False)
pl1.clear()
pl2 = pg.ViewBox()
pl2.setXLink(pl1)
pl1.scene().addItem(pl2)
pl2.setMouseEnabled(x=True, y=False)

win = QtGui.QWidget()
layout = QtGui.QGridLayout()
layout.setRowStretch(0, 1)
layout.setRowStretch(1, 1)
layout.setRowStretch(2, 20)
#layout.setRowStretch(3, 4)
win.setLayout(layout)
layout.addWidget(Lbl1, 0, 1, 1, 3)
layout.addWidget(Lbl2, 0, 0, 1, 1)
layout.addWidget(Btn2, 1, 0, 1, 1)
layout.addWidget(Btn3, 1, 1, 1, 1)
layout.addWidget(Btn5, 1, 2, 1, 1)
layout.addWidget(CBox, 1, 3, 1, 1)
layout.addWidget(plt, 2, 0, 1, 4)
#layout.addWidget(pl2, 3, 0, 3, 3)

mFolderName = sys.argv[0]
i = len(mFolderName) - 1
while i >= 0:
    if mFolderName[i] == '/':
        mFolderName = mFolderName[0:i+1]
        i = 0
    i -= 1
mFolderName =  "C://USER//CTGcomp//CTGfiles//"
mReportFileName = "C://USER//CTGcomp//CTGfiles//STVexport.txt"
if len(sys.argv) > 1:
    mFileName = sys.argv[1]
    if len(sys.argv) > 2:
        mOutFileName = sys.argv[2]
    elif len(mFileName) > 0:
        i = len(mFileName) - 1
        while i >= 0:
            if mFileName[i] == '/':
                mReportFileName = mFileName[0:i+1] + 'STVexport.txt'
                i = 0
            i -= 1
    if (mFileName[-3:] == 'prn') or (mFileName[-3:] == 'txt'):
        ReadDataFile()
    elif mFileName[-3:] == 'ctg':
        ReadBinary()
    else:
        GetFile()
else:
    GetFile()


#pl1.vb.sigResized.connect(updateViews)
#win.show()
win.showMaximized()

## Start Qt event loop unless running in interactive mode or using pyside.
if __name__ == '__main__':
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        mApp.instance().exec_()
