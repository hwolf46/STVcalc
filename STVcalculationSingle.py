# -*- coding: utf-8 -*-
"""
Created on 03-28-2017
Python 3.6
Author: HWolf
Basic program for cCTG analysis
CTG files should be in same folder as this application
"""
import sys
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui
import numpy as np
import pymsgbox

mApp = QtGui.QApplication([])
sys.setrecursionlimit(5000)

def GetFile():
    global mFileName, mFolderName, mTableName, mPlotTitle
    global mReport, mExport, mRange
    global mSTVmin, mFHRbase, mFHRcode, mXdata, mFHR, mFHRX, mAccelLoc, mDecelLoc
    global mFileCounter, mCTGtooShort
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
    mExport = ''
    mFileName = ''
    mFileNameR = ''
    mFileCounter = 0
    try:
        mFileNameR = (QtGui.QFileDialog.getOpenFileName
                      (None, 'Open file for STV calculation',
                      mFolderName, '(*.prn *.ctg *.txt)'))
    except FileNotFoundError:
        pymsgbox.alert(text='You did not read a CTG data file. \nThe program will cancel.',
                       title='Pas op', button='Gezien')
        raise SystemExit
    mFileName = mFileNameR[0]
    mPlotTitle = mFileName
    if (mFileName[-3:] == 'prn') or (mFileName[-3:] == 'txt'):
        ReadDataFile()
    elif mFileName[-3:] == 'ctg':
        ReadBinary()
    else:
        pymsgbox.alert(text='You did not read a CTG data file. \nOther formats are not acceptable. Try again.',
                       title='Pas op', button='Gezien')
    if np.sum(mFHR) > 0:
        RejectData()
        if mCTGtooShort == 1:
            mExport = mFileName + '; 999; 999; 999; 999; 999; 999; 999; 999; 999'
            mF = open(mReportFileName, 'a')
            mF.write(mExport + '\n')
            mF.close
        else:
            EpochCalc()
            CalcSTV()
            CBox.setValue(1)
            MakePlot1()

def ReadBinary():
#read data from binary MOSOS file with .ctg (new version)
    global mFileName
    global mXdata
    global mFHR, mQual
    global mUtP, mFetMov
    global mIntV
    global mKindNr
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
    mFetMov = []
    mXdata = []
    mIntV = []
    mQual = []
    mUtP = []
    with open(mFileName, 'br') as File:
        byte_content = File.read()
        mData = [byte_content[i + 1] << 8 | byte_content[i] for i in range(0, len(byte_content), 2)]
    File.close()
    mTel = 0
# To detect if channel 1 or 2 was connected to the US probe (should always be 1 in singletons, but soemtimes a mistake is made)
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
    if mSingle == -1:
        pymsgbox.alert('File contains insufficient data. Program will shut down.')
        raise SystemExit
    mTel -= 12
    if mKindNr == 1:
        mSingle = 0
    elif mKindNr == 2:
        mSingle = 4
    while mTel < len(mData):
        mB = bin(mData[mTel + mSingle])
        mB = '0b' + ('0000000000000000' + mB[2:])[-16:]
        mFHR.append(int(mB[-10:-2], 2))
        mQual.append(int(mB[-13:-10], 2))
        mFetMov.append(int(mB[-15:-13], 2))
        mB = bin(mData[mTel + 1  + mSingle])
        mB = '0b' + ('0000000000000000' + mB[2:])[-16:]
        mFHR.append(int(mB[-10:-2], 2))
        mQual.append(int(mB[-13:-10], 2))
        mFetMov.append(int(mB[-15:-13], 2))
        mB = bin(mData[mTel + 2  + mSingle])
        mB = '0b' + ('0000000000000000' + mB[2:])[-16:]
        mFHR.append(int(mB[-10:-2], 2))
        mQual.append(int(mB[-13:-10], 2))
        mFetMov.append(int(mB[-15:-13], 2))
        mB = bin(mData[mTel + 3  + mSingle])
        mB = '0b' + ('0000000000000000' + mB[2:])[-16:]
        mFHR.append(int(mB[-10:-2], 2))
        mQual.append(int(mB[-13:-10], 2))
        mFetMov.append(int(mB[-15:-13], 2))
# uterus tonus
        mB = bin(mData[mTel+10])
        mB = '0b' + ('0000000000000000' + mB[2:])[-16:]
        mUtP.append(int(mB[-8:], 2))
        mUtP.append(int(mB[-8:], 2))
        mUtP.append(int(mB[-16:-8], 2))
        mUtP.append(int(mB[-16:-8], 2))
        mTel += 12
#Delete zero values at the start of the registrtation, usually caused by connecting the machine before fetal signal is stable
    mTel = 0
    while mFHR[mTel] == 0:
        del mFHR[mTel]
        del mUtP[mTel]
        mTel += 1
# Setting the X-axis values
    mTel = 0
    while mTel < len(mFHR):
        mXdata.append(mTel/240)
        mTel += 1
# Calculation of interval from FHR and exclusion of extremes
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
    if np.count_nonzero(mFHR)< 2400:
        pymsgbox.alert('File contains less than 10 minutes CTG data. Program will shut down.')
        raise SystemExit
#  Make a prn file
#    mF = open(mFileName[0:-3] + 'prn', 'w')
#    mTel = 0
#    while mTel < len(mFHR):
#        mF.write(str(mFHR[mTel]) + ', ' + str(mUtP[mTel] + ', ' + str(mQual[mTel] + '\n')
#        mTel += 1
#    mF.close

def ReadDataFile():
#Read data from prn file (heart rate, uterine pressure, and  quality bit values in single column at 4 Hz)
#If the file contains only FHR then these data will be read and uterine pressure remains 0
    global mFHR, mQual
    global mUtP
    global mIntV
    global mXdata
    global mFileName
    mTel = 0.00
    mFr = 0
    mI = 0
    mQ = 0
    mU = 0
    i = 0
    b = 0
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
            mU = int(mLLine[1])
        if len(mLLine) > 2:
            mQ = int(mLLine[2])
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
    if np.count_nonzero(mFHR)< 2400:
        pymsgbox.alert('File contains less than 10 minutes CTG data. Program will shut down.')
        raise SystemExit

def RejectData():
    # rejection algorithm for data very far from meean
    global mIntV, mFHR, mCTGtooShort
    mI1 = 0
    mI2 = 0
    mI3 = 0
    mI4 = 0
    mTel = 0
    mIMean = 0
    mIList = []
    mCTGtooShort = 0
    mIMeanTot = np.median(list(filter((0).__ne__, mIntV)))
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
    if len(mFHR) < 2400:
        mCTGtooShort = 1

def EpochCalc():
#epoch calculation
    global mIntV, mIntEp, mFHRcode, mFHRXep, mFHRep, mUtEp, mUtP
    global Lbl1
    mIntEp = []
    mFHRep = []
    mFHRcode = []
    mFHRXep = []
    mUtEp = []
    mTel = 0
    mCode = 0
    mIntVmean = 0
    mUTmean = 0
    mListInt = []
    while mTel < len(mIntV):
        if np.bincount(mIntV[mTel:mTel+15])[0] < 13 and len(mIntV[mTel:mTel+15]) == 15:
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
        mListInt = [num for num in mUtP[mTel:mTel+15] if mUtP[num] > 0]
        if len(mListInt) > 0:
            mUTmean = (np.mean(mListInt))*100/256
        else:
            mUTmean = 0
        mUtEp.append(mUTmean)
        mTel += 15
    mTel = 1
# wissen geisoleerde epochs
    while mTel < len(mIntEp) - 1:
        if mIntEp[mTel] > 0 and mIntEp[mTel-1] == 0 and mIntEp[mTel+1] == 0:
            mIntEp[mTel] = 0
            mFHRep[mTel] = 0
            mFHRcode[mTel] = 210
        if ((mTel > 2 and mTel < len(mIntEp) - 5) and (mIntEp[mTel] > 0 and mIntEp[mTel+1] > 0) and
            (np.sum(mIntEp[mTel-3:mTel]) == 0 and np.sum(mIntEp[mTel+2:mTel+5]) == 0)):
            mIntEp[mTel] = 0
            mFHRep[mTel] = 0
            mFHRcode[mTel] = 210
            mIntEp[mTel+1] = 0
            mFHRep[mTel+1] = 0
            mFHRcode[mTel+1] = 210
        mTel += 1
#histogram of IBI epochs
#    n=0
#    bins=0
#    patches=''
#    L = len(np.bincount(mIntEp))
#    n, bins, patches = plt.hist(mIntEp, L, normed=1, facecolor='blue', alpha=0.85)
#    plt.grid(True)
#    plt.axis([300, 550, 0, 0.07])
#    plt.xlabel('IBI')
#    plt.ylabel('Frequency')
#    plt.savefig(mFileName[:-3] + 'png')
#    plt.close()


def CalcBaseline():
    global mFHRep, mIntEp, mFHRcode
    global mFHRbase
    global mFHRXep
    global mSTVmin
    global mRefHRe, mStartHRe, mRefTe
    mSTVmin = []
    mDist = []
    mRefHRe = 0
    mRefTe = 0
    mStartHRe = 0
    mStartTe = 0
# Calculation reference point according to Dawes with interval
    mLimit = 0
    mCount = 0
    mDist = np.bincount(mIntEp)
    mMaxDist = np.argmax(mDist[300:]) + 300
    mLimitTot = np.sum(mDist) - mDist[0]
    mTel = len(mDist) - 1
    while mTel > 300:
        mLimit += mDist[mTel]
        if mLimit > mLimitTot / 5:
            if (mDist[mTel] > mDist[mTel-1] and mDist[mTel] > mDist[mTel-2] and
                    mDist[mTel] > mDist[mTel-3] and mDist[mTel] > mDist[mTel-4] and mDist[mTel] > mDist[mTel-5]):
                mRefTe = mTel
                mCount = mDist[mTel]
#                if (mCount / mLimitTot) > 0.005:
                if (mCount / mLimitTot) > 0.01:
                    mTel = 0
        mTel -= 1
    if mRefTe == 0:
        mRefTe = mMaxDist
    elif np.absolute(mRefTe - mMaxDist) > 30 and (mCount / mLimitTot) < 0.005:
        mRefTe = mMaxDist
    if mRefTe > 0:
        mRefHRe = 60000 / mRefTe
# Calculation starting point
    mTel = 0
    mLimit = 20
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
        mTel = 0
        mLimit = 20
        mTe = 0
        i = 0
        while mLimit <= 40:
            for i in mIntEp[64:128]:
                if (i > mRefTe - mLimit) and (i < mRefTe + mLimit):
                    mTe += i
                    mTel += 1
            if mTel > 10 and mTe > 0:
                mStartTe = mTe / mTel
                mStartHRe = 60000 / mStartTe
                mLimit = 50
            mLimit += 10
    if mStartHRe == 0 or np.absolute(mStartHRe - mRefHRe) > 20:
        mStartHRe = mRefHRe
        mStartTe = mRefTe
# exclude incidental jumps of 35/min or more
    mTel = 1
    while mTel < len(mFHRep)-2:
        if mFHRep[mTel] > 0 and mFHRep[mTel-1] > 0:
            if ((np.absolute(mFHRep[mTel] - mFHRep[mTel-1]) > 35) and
               (np.absolute(mFHRep[mTel+1] - mFHRep[mTel-1]) < 35 or mFHRep[mTel+1] == 0)):
                if mFHRep[mTel] - mFHRep[mTel-1] > 0:
                    mFHRcode[mTel] = 200
                else:
                    mFHRcode[mTel] = 185
                mFHRep[mTel] = 0
                mIntEp[mTel] = 0
        elif mFHRep[mTel] > 0 and mFHRep[mTel-1] == 0:
            if ((np.absolute(mFHRep[mTel] - mRefHRe) > 35) and
                    (np.absolute(mFHRep[mTel+1] - mRefHRe) < 35 or mFHRep[mTel+1] == 0)):
                if mFHRep[mTel] - mRefHRe > 0:
                    mFHRcode[mTel] = 200
                else:
                    mFHRcode[mTel] = 185
                mFHRep[mTel] = 0
                mIntEp[mTel] = 0
        mTel += 1
# write and smooth baseline first pass
    mLimit = 30
    BaselineWrite(mStartTe, mRefTe, mStartHRe, mLimit)
# check for baseline difference of more than 10 minutes of more than 150 msec is missing yet
    d = 0
    while d < 7:
        mTel = 0
        i1 = 0
        i2 = 0
        mMark1 = 0
        mMark2 = 0
        mMark3 = 0
        while mTel < len(mTbase):
            if mTbase[mTel] > mIntEp[mTel]:
                i1 += 1
                if i2 > 20 and mMark2 < 3:
                    mMark2 += 1
                else:
                    i2 = 0
                    mMark2 = 0
            elif mTbase[mTel] < mIntEp[mTel]:
                i2 += 1
                if i1 > 20 and mMark3 < 3:
                    mMark3 += 1
                else:
                    i1 = 0
                    mMark3 = 0
            else:
                i1 = 0
                i2 = 0
            if i1 > 160 or i2 > 160:
                mTel = len(mTbase)
                mLimit += 10
                BaselineWrite(mStartTe, mRefTe, mStartHRe, mLimit)
                mMark1 = 1
            mTel += 1
        if mMark1 == 0:
            d = 7
        else:
            d += 1

def BaselineWrite(mStartTe, mRefTe, mStartHRe, mLimit):
    global mFHRbase, mTbase, mIntEp
    mTel = 0
    mLapseT = 8
    mLapseB = 0
    mLapseE = 0
    mListInt = []
    mMeanIntRec = mStartTe
    mFHRbase = []
    mFHRbaseR = []
    mTbase = []
# Set mean values for a minute 2 x mLapseT) that are within the constraints set by mLimit to baseline
    while mTel < len(mIntEp):
        if mTel < mLapseT:
            mLapseB = mTel
        else:
            mLapseB = mLapseT
        if mTel > len(mIntEp) - mLapseT:
            mLapseE = len(mIntEp) - mTel
        else:
            mLapseE = mLapseT
        mListInt = ([num for num in mIntEp[mTel-mLapseB:mTel+mLapseE]
                     if num > mRefTe-mLimit and num < mRefTe + mLimit])
        if len(mListInt) > 2:
            mMeanIntRec = np.mean(mListInt)
        mTbase.append(mMeanIntRec)
        mFHRbase.append(60000/mMeanIntRec)
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
    global mFHRX, mFHRbase, mTbase, mFHRep, mIntEp, mFHRcode, mAccelLoc, mDecelLoc
    global mReport, mExport, mFileName, mReportFileName, mCheckPlot
    global mStrUpdate, mStrInsert, mZis, mDat, mKindNr
    global mSTVmin, mLTVmin
    CalcBaseline()
    mSTVmin = []
    mSTVX = []
    mAccelLoc = []
    mDecelLoc = []
    mSTV = 0
    mTel = 0
    mFHRmean = 0
    mValid = 0
# Coderen hartfrequenties 10 of 20 slagen boven / onder basislijn of 75ms van basis interval
    mTel = 0
    while mTel < len(mFHRep):
        if mIntEp[mTel] < (60000/mFHRbase[mTel]) - 75 and mFHRep[mTel] > 0:
            mFHRcode[mTel] = 200
        elif mIntEp[mTel] > (60000/mFHRbase[mTel]) + 75:
            mFHRcode[mTel] = 185
        mTel += 1
# Terugzetten code van FHR meer of minder dan <20 van baseline die niet in deceleratie zitten en
#        alle accelleraties.
# Exclusie blijft alle deceleraties en van mFHRep meer dan +-75 msec van baseline, die niet in accel. of decel. zit,
# dit is code 192 (decelleratie), en 185 en 200 (uitbijters)
    mTel = 0
    while mTel < len(mFHRcode):
        if (mFHRcode[mTel] == 187 or mFHRcode[mTel] == 190
                or mFHRcode[mTel] == 195 or mFHRcode[mTel] == 197):
            mFHRcode[mTel] = 1
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
    mTot = len(mFHRep)
    while mTel < mTot:
        if mFHRep[mTel] >= mFHRbase[mTel] and mFHRcode[mTel] != 185:
            while(mTel < mTot and mFHRep[mTel] > mFHRbase[mTel] and mFHRcode[mTel] != 185):
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
                mAccelLoc[mTel - i] = 197
                if mFHRcode[mTel - i] == 200:
                    mFHRcode[mTel - i] = 1
                i -= 1
            i = 0
# Tellen decelleraties
    mDecell = 0
    mDecell60 = 0
    mTel = 0
    i = 0
    mMin = 0
    while mTel < mTot:
        if (mFHRep[mTel] < mFHRbase[mTel]) and (mFHRep[mTel] > 10):
            while ((mTel < mTot and mFHRep[mTel] < mFHRbase[mTel] and  mFHRep[mTel] > 10) or
                   (mTel < mTot-1 and mFHRep[mTel] < 10 and mFHRep[mTel+1] < mFHRbase[mTel+1] and mFHRep[mTel+1] > 10) or
                   (mTel < mTot-2 and mFHRep[mTel] < 10 and mFHRep[mTel+1] < 10 and mFHRep[mTel+2] < mFHRbase[mTel+2] and mFHRep[mTel+2] > 10)):
                if (mMin < mFHRbase[mTel] - mFHRep[mTel]) and mFHRep[mTel] > 10:
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
                mFHRcode[mTel - i] = 192
                mDecelLoc[mTel - i] = 192
                i -= 1
            i = 0
# Naastliggende waarden van uitbijters ook excluderen
    mTel = 1
    while mTel < len(mFHRcode) - 1:
        if mFHRcode[mTel] == 185 or mFHRcode[mTel] == 200:
            mFHRcode[mTel-1] = mFHRcode[mTel]
            if not (mFHRcode[mTel+1] == 185 or mFHRcode[mTel+1] == 200):
                mFHRcode[mTel+1] = mFHRcode[mTel]
                mTel += 1
        mTel += 1
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
        while n < 16 and mTe < 9000:
            if (mFHRcode[mTel + n] == 192 or mFHRcode[mTel + n - 1] == 192):
                mTe = 9999
            elif mFHRcode[mTel + n] == 1 and mFHRcode[mTel + n - 1] == 1:
                mTe += np.absolute(mIntEp[mTel + n] - mIntEp[mTel + n - 1])
                i += 1
            n += 1
        if i >= 8 and mTe > 0 and mTe < 9000:
            mSTVmin.append(mTe / i)
        else:
            mSTVmin.append(99)
        mSTVX.append(mTel/16)
        mTel += 16
# Berekenen van Longterm variation
    mHighLTVtime = 0
    mLowLTVtime = 0
    mHighLTVtime60 = 0
    mLowLTVtime60 = 0
    mRangeMin = []
    mLTVmin = []
    mTel = 0
    while mTel < len(mIntEp) - 16:
        n = 0
        i = 0
        mMax = 0
        mMin = 0
        while n < 16:
            if mFHRcode[mTel + n] == 1:
                if mMax < mIntEp[mTel+n] - mTbase[mTel + n]:
                    mMax = mIntEp[mTel+n] - mTbase[mTel + n]
                if mMin < mTbase[mTel + n] - mIntEp[mTel+n]:
                    mMin = mTbase[mTel+n] - mIntEp[mTel+n]
                i += 1
            n += 1
        if i > 8:
            mRangeMin.append(mMax + mMin)
        else:
            mRangeMin.append(999)
        mTel += 16
    mTel = 0
    while mTel < len(mRangeMin):
        mLTVmin.append(0)
        mTel += 1
# Calculate high LTV minutes
    mTel = 0
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
# Calculate low LTV minutes
    while mTel < len(mRangeMin) - 6:
        n = 0
        i = 0
        while n < 6:
            if mRangeMin[mTel + n] < 30 and mRangeMin[mTel + n] > 0:
                i += 1
            n += 1
        if i >= 5:
            n = 0
            while n < 6:
                mLTVmin[mTel + n] = 170
                n += 1
        mTel += 1
    mTel = 0
    i = 0
    j = 0
    while mTel < len(mLTVmin):
        if mLTVmin[mTel] == 175:
            i += 1
        elif mLTVmin[mTel] == 170:
            j += 1
        if mTel < 60:
            mHighLTVtime60 = i
            mLowLTVtime60 = j
        mTel += 1
    mHighLTVtime = i
    mLowLTVtime = j
    mLTV = np.mean(list(filter((999).__ne__, mRangeMin)))
# Berekenen FHR baseline
    i = 0
    mTel = 0
    while mTel < len(mFHRbase):
        if mFHRbase[mTel] > 0:
            i += 1
        mTel += 1
    mFHRmean = sum(mFHRbase) / i
# Berekenen gemiddelde minuut-STV over hele registratie en exclusie minute STV that is > 3 x MoM
    mTel = 0
    mSTV = 0
    mSTV60 = 0
    mSTVmin[0] = 99
    mSTV = np.mean([num for num in mSTVmin[1:] if num > 0 and num < 90])
    while mTel < len(mSTVmin):
        if mSTVmin[mTel] > 4 * mSTV:
            mSTVmin[mTel] = 99
        mTel += 1
    mSTV = np.mean([num for num in mSTVmin[1:] if num > 0 and num < 90])
    mSTV60 = np.mean([num for num in mSTVmin[1:61] if num > 0 and num < 90])
# Rapportage
    mDuration = round(len(mFHRep)/16)
    mValid = len([num for num in mSTVmin[1:] if num > 0 and num < 90])
    mReport = ('Base FHR : ' + str(round(mFHRmean)) + ', STV : ' + str(round(mSTV, 2))
               + ', STV 60min: ' + str(round(mSTV60, 2)) + ', LTV: '  + str(round(mLTV, 2)) +
               ', High LTV minutes:' + str(mHighLTVtime) + ', High LTV minutes 60 min:' + str(mHighLTVtime60) +
               ', Low  LTV minutes:' + str(mLowLTVtime) + ', Low  LTV minutes 60 MIN:' + str(mLowLTVtime60) +
               '\n Duration CTG : ' + str(mDuration) + ' minutes, Valid minutes : ' + str(mValid) +
               ', Lost minutes : ' + str(round((1 -(mValid/mDuration)) * 100)) +
               ', Accellerations all : ' + str(mAccell) + ', <60 min. :' + str(mAccell60) +
               ', Decellerations all : ' + str(mDecell) + ', <60 min. :' + str(mDecell60) + ', Type : ' +
               str(CBox._chosenText))
    mExport = (mFileName + '; ' + str(round(mFHRmean)) + '; ' + str(round(mSTV, 2)) + '; ' +
               str(round(mSTV60, 2)) + '; ' + str(round(mLTV, 2)) + '; ' + str(mHighLTVtime)
               + '; ' + str(mHighLTVtime60) + '; ' + str(mLowLTVtime) + '; ' + str(mLowLTVtime60) + '; '  +
               str(round(mDuration)) + '; ' + str(mValid) + '; ' + str(mAccell) + '; ' + str(mAccell60) + '; ' +
               str(mDecell) + '; ' + str(mDecell60) + '; ' + str(CBox.value()))
    if mSTV > 0 and mFHRmean > 0:
        mStrUpdate = ('FHR = ' + str(round(mFHRmean)) + ', STV = ' + str(round(mSTV, 2)) + ', STV60 = ' +
                      str(round(mSTV60, 2)) + ', LTV = ' + str(round(mLTV, 2)) + ', HighMin = ' +
                      str(mHighLTVtime) + ', HighMin60 = ' + str(mHighLTVtime60) + ', LowMin = ' +
                      str(mLowLTVtime) + ', LowMin60 = ' + str(mLowLTVtime60) + ', Duration = '  +
                      str(round(mDuration)) + ', Valid = ' + str(mValid) + ', Accel = ' + str(mAccell) +
                      ', Decel = ' + str(mDecell) + ', DecelType = '  + str(CBox.value()))
        mStrInsert = ("(Regis_Key, Zis_Nr, KindNr, FileName, CTGdate, FHR, STV, STV60, LTV, HighMin, " +
                      "HighMin60, LowMin, LowMin60, Duration, Valid, Accel, Decel, DecelType) VALUES(" +
                      str(mRegKey) + ", '" + mZis + "', " + str(mKindNr) + ", '" + mFileName +
                      "', #" + mDat + "#, " + str(round(mFHRmean)) + ", " + str(round(mSTV, 2)) + ", " +
                      str(round(mSTV60, 2)) + ", " + str(round(mLTV, 2)) + ", " + str(mHighLTVtime) +
                      ", " + str(mHighLTVtime60) + ", " + str(mLowLTVtime) + ", " + str(mLowLTVtime60) +
                      ", "  + str(round(mDuration)) + ", " + str(mValid) +
                      ", " + str(mAccell) + ", " + str(mDecell) + ", "  + str(CBox.value()) + ")")
    else:
        mStrUpdate = ''
        mStrInsert = ''
    if len(mExport) > 0:
        mF = open(mReportFileName, 'a')
        mF.write(mExport + '\n')
        mF.close
#  Make an epoch file
#    mF = open(mFileName[0:-3] + 'epo', 'w')
#    mTel = 0
#    while mTel < len(mFHRep):
#        mF.write(str(round(mFHRep[mTel], 2)) + '\n')
#        mTel += 1
#    mF.close


def CBoxCheckChanged():
    if np.sum(mFHR) > 0:
        CalcSTV()
        MakePlot1()

def Recalculate():
    CalcSTV()
    MakePlot1()

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
    MakePlot1()

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
    RejectData()
    EpochCalc()
    CalcSTV()
    MakePlot1()

def MakePlot1():
## Create plotwindow
    global mSTVmin
    global mFileName, mReport, mPlotTitle
    global mFHRep, mFHRXep, mFHRbase, mFHRcode, mAccelLoc, mDecelLoc, mLTVmin, mUtEp
    global mUtP, mXdata
    global pl1, pl2, pl3, win
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
    Plot1.clear()
    pl1.clear()
    pl2.clear()
    Plot2.clear()
    pl3.clear()
    font = QtGui.QFont()
    font.setPixelSize(20)
    Plot1.showGrid(x=True, y=False, alpha=100)
    Plot1.setBackground('#e6e6e6')
    Plot2.showGrid(x=True, y=True, alpha=100)
    Plot2.setBackground('#e6e6e6')
    Plot1.setTitle(mPlotTitle)
    pl1.setYRange(80, 200)
    Plot1.setXRange(0, 60)
    pl1.getAxis('bottom').tickFont = font
    pl1.getAxis('left').tickFont = font
    pl1.getAxis('left').setGrid(100)
    pl1.getAxis('left').setLabel('FHR', color='#0000ff', font='20pt')
    pl1.getAxis('right').setLabel('STV', color='#0000ff', font='20pt')
    pl1.getAxis('right').tickFont = font
    pl1.showAxis('right')
    pl2.setXLink(pl1)
    pl2.setYRange(0, 30)
    pl2.setXRange(0, 60)
    pl3.setLabel('bottom', text='Minutes', color='#0000ff', font='20pt')
    pl3.setYRange(0, 100)
    pl3.setXRange(0, 60)
    pl3.setXLink(pl1)
    pl3.showAxis('right')
    pl3.getAxis('left').setLabel('Uterus', color='#0000ff', font='20pt')
    pl3.getAxis('right').setLabel('Pressure', color='#0000ff', font='20pt')
    pl3.getAxis('bottom').tickFont = font
    pl3.getAxis('left').tickFont = font
    pl3.getAxis('right').tickFont = font
    pl3.setXRange(0, 60)
    pl1.vb.sigResized.connect(updateViews)
    pl3.vb.sigResized.connect(updateViews)
    pl1.plot(mFHRXep, mFHRepN, connect="finite", pen=pg.mkPen('k', width=2))
#    pl1.plot(mFHRXep, mFHRepN, connect="finite", pen='k', color='b', symbol='o', symbolSize=3, symbolBrush='b')
#    pl1.plot(mXdata, mFHRN, connect="finite", pen='k', color='b')
    pl2.addItem(pg.ScatterPlotItem(range(0, len(mSTVmin)), mSTVmin, symbol='o', brush='b', color='b', size=12))
    if np.sum(mFHRcode) > 0:
        pl1.plot(mFHRXep, mFHRcode, symbol='s', pen=None, symbolBrush='r', symbolPen='r', symbolSize=10)
    if np.sum(mFHRbase) > 0:
        pl1.plot(mFHRXep, mFHRbase, pen='k')
    if np.sum(mDecelLoc) > 0:
        pl1.plot(mFHRXep, mDecelLoc, pen=None, symbol='s', symbolBrush='m', symbolPen=None, symbolSize=12)
    if np.sum(mAccelLoc) > 0:
        pl1.plot(mFHRXep, mAccelLoc, pen=None, symbol='s', symbolPen='g', symbolBrush='g', symbolSize=12)
    if np.sum(mLTVmin) > 0:
        pl1.plot(range(0, len(mLTVmin)), mLTVmin, pen=None, symbol='s', symbolPen=(0,128,0), symbolSize=12)
    Lbl1.setText(mReport)
    pl3.plot(mFHRXep, mUtEp, pen=pg.mkPen('k', width=2))
    updateViews()
    win.activateWindow()
    win.showMaximized()

def updateViews():
    ## view has resized; update auxiliary views to match
    global pl1, pl2
    pl2.setGeometry(pl1.vb.sceneBoundingRect())
    pl2.linkedViewChanged(pl1.vb, pl2.XAxis)

class CustomViewBox(pg.ViewBox):
    def __init__(self, *args, **kwds):
        pg.ViewBox.__init__(self, *args, **kwds)
        self.setMouseMode(self.PanMode)
    def mouseClickEvent(self, ev):
        global mMouseAction
        if ev.button() == QtCore.Qt.RightButton:
            self.setMouseMode(self.PanMode)
            mMouseAction = 0
        elif ev.button() == QtCore.Qt.LeftButton:
            self.setMouseMode(self.PanMode)
            mMouseAction = 0
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
def DeSelect():
    global mMouseAction, vb
    mMouseAction = 1
    vb.setMouseMode(vb.RectMode)

def SaveExit():
    sys.exit("done")

#Start of program.......................
mFileName = ''
mFolderName = ''
mTableName = ''
mOutFileName = ''
mReportFileName = ''
mPlotTitle = ''
mStrInsert = ''
mStrUpdate = ''
mReport = 'Report'
mExport = ''
mZis = ''
mDat = ''
mKindNr = 0
mRegKey = 0
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
mFetMov = []
mFHRep = []
mIntEp = []
mFHRXep = []
mIntV = []
mIntEp = []
mAccelLoc = []
mDecelLoc = []
mMouseAction = 0
mVersion = 1
mRefHRe = 0
mStartHRe = 0
mCTGtooShort = 0

vb = CustomViewBox()

Btn1 = QtGui.QPushButton(text=' Exit ')
Btn1.setToolTip('Exit program')
Btn1.clicked.connect(SaveExit)
Btn2 = QtGui.QPushButton(text=' New File ')
Btn2.clicked.connect(GetFile)
Btn2.setToolTip("Click to open the file selection dialog")
Btn3 = QtGui.QPushButton(text='  Restore complete registration  ')
Btn3.clicked.connect(SelectAll)
Btn3.setToolTip("Click to restore the complete CTG after deselection of a part of the tracing")
Btn4 = QtGui.QPushButton(text='  De-select part of CTG ')
Btn4.clicked.connect(DeSelect)
Btn4.setToolTip("Click this button to activate mouse drag for selection of a part of the CTG that should be removed")
Btn5 = QtGui.QPushButton(text='  Calculate STV  ')
Btn5.clicked.connect(Recalculate)
Btn5.setToolTip("Click to re-calculate STV")
CBox = pg.ComboBox(items={'No decellerations' :1,
                          'Variable decelerations 1/hr' : 2,
                          'Variable decelerations 2/hr' : 3,
                          'Variable decelerations >=3/hr' : 4,
                          'Recurrent decelerations' : 5,
                          'Signal insufficient' : 6})
CBox.setValue(1)
CBox.setToolTip("Select the most appropriate setting describing decelerations in the CTG")
CBox.currentIndexChanged.connect(CBoxCheckChanged)
Lbl1 = QtGui.QLabel("""Results""")
Lbl1.setText("Report")
Lbl1.setFont(QtGui.QFont('SansSerif', 14))
Plot1 = pg.PlotWidget(viewBox=vb)
Plot2 = pg.PlotWidget(viewbox=vb)
pl1 = Plot1.plotItem
pl2 = pg.ViewBox()
pl1.getAxis('right').linkToView(pl2)
pl1.scene().addItem(pl2)
pl1.setMouseEnabled(x=True, y=False)
pl2.setMouseEnabled(x=True, y=False)
pl3 = Plot2.plotItem
pl3.setMouseEnabled(x=True, y=False)
win = QtGui.QWidget()
layout = QtGui.QGridLayout()
layout.setRowStretch(0, 1)
layout.setRowStretch(1, 1)
layout.setRowStretch(2, 20)
layout.setRowStretch(3, 5)
win.setLayout(layout)
layout.addWidget(Lbl1, 0, 0, 1, 6)
layout.addWidget(Btn2, 1, 0, 1, 1)
layout.addWidget(Btn3, 1, 1, 1, 1)
layout.addWidget(Btn4, 1, 2, 1, 1)
layout.addWidget(Btn5, 1, 3, 1, 1)
layout.addWidget(CBox, 1, 4, 1, 1)
layout.addWidget(Btn1, 1, 5, 1, 1)
layout.addWidget(Plot1, 2, 0, 1, 6)
layout.addWidget(Plot2, 3, 0, 1, 6)


mArgs = sys.argv
mFolderName = mArgs[0]
i = len(mFolderName) - 1
while i >= 0:
    if mFolderName[i] == "/":
        mFolderName = mFolderName[0:i+1]
        i = 0
    i -= 1
if len(mArgs) > 1:
    mFileName = mArgs[1]
else:
    mFileName = ''
if len(mArgs) > 2:
    mPlotTitle = mArgs[2]
    mZis = mArgs[2][:7]
    mDat = mArgs[2][9:]
if len(mArgs) > 3:
    mRegKey = mArgs[3]

mReportFileName = mFolderName + "STVexport.txt"

if len(mFileName) > 0:
    ReadBinary()
    RejectData()
    EpochCalc()
    CalcSTV()
    CBox.setValue(1)
    MakePlot1()
else:
    GetFile()



win.showMaximized()

## Start Qt event loop unless running in interactive mode or using pyside.
if __name__ == '__main__':
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        mApp.instance().exec_()
