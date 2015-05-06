# -*- coding: cp936 -*-
"""
Tool Name:  EOF
Source Name: EOF.py
Version: ArcGIS 10.1
Author: WangHuachen of Hohai University


"""

################### Imports ########################
import os as OS
import sys as SYS
import collections as COLL
import numpy as NUM
import math as MATH
import numpy.random as RAND
import arcpy as ARCPY
import arcpy.management as DM
import arcpy.analysis as ANA
import arcpy.da as DA
import ErrorUtils as ERROR
import SSUtilities as UTILS
import locale as LOCALE
import copy
import pylab as PL

LOCALE.setlocale(LOCALE.LC_ALL, '')

################### GUI Interface ###################
def setupEOF():
    """SDFSDFSDFSDFSFDSDFSDF."""

    inputFFCC = ARCPY.GetParameterAsText(0)
    selectFld = ARCPY.GetParameterAsText(1).upper()
    selectFld = selectFld.split(";")
    way = ARCPY.GetParameterAsText(2)
    outputNum = ARCPY.GetParameterAsText(3)
    displayTime = ARCPY.GetParameter(4)
    code = ARCPY.GetParameter(5)
    outputControl = ARCPY.GetParameter(6)
    outputPath = ARCPY.GetParameterAsText(7)
    

    ARCPY.AddMessage(way)
    
    
    if way == "OriginData": 
        way = 1
    elif way == "Standardization":
        way = 3
    elif way == "RowBalance":
        way = 2
    else:
        way = 4

    ARCPY.AddMessage(way)
    
    

    #######正式开始调用EOF函数#########

    e = EOF(inputFFCC, selectFld, way, outputNum, displayTime, outputControl, outputPath)

    

class EOF(object):
    """aaaa
    """

    def __init__(self, inputFFCC, selectFld, way, outputNum, displayTime, outputControl, outputPath):

        #### 设置初始属性 ####
        UTILS.assignClassAttr(self, locals())
       

        #### 初始化 ####
        self.initialize() 

        #### 标准化或者矩平 ####
        self.transf()

        #### 求协方差矩阵 ####
        self.forma()

        #### 雅克比分解 ####
        self.jcb()
        #### arrange ####
        self.arrang()

        #### tcoeff ####
        self.tcoeff()

        self.createOutput()


    def initialize(self):
        

        #######
        ##result = []
        ##originData = []
        originI = 0
        originJ = len(self.selectFld)

        with DA.SearchCursor(self.inputFFCC,["OBJECTID"]) as b:
            for row in b:
                originI += 1


        originX = NUM.zeros((originJ,originI))
        ##ARCPY.AddMessage(originI)
        ##ARCPY.AddMessage(originJ)
        

        for i in xrange(originJ):
            tmpfld = self.selectFld
            tmpresult = []
            ##ARCPY.AddMessage(tmpfld[i])
            with DA.SearchCursor(self.inputFFCC,tmpfld[i]) as b:
                for row in b:
                    tmpresult.append(row[0])
            ##ARCPY.AddMessage(tmpresult)
            originX[i] = tmpresult
            
            
            del tmpfld
            del tmpresult
        

        
                

        originX = originX.T

        


                
        
        
        minIJ = min(originI,originJ)
        

        self.originI = originI
        self.originJ = originJ
        self.originX = originX
        self.minIJ = minIJ


    def transf(self):
        """Sets the study area for the k-function."""

        originI = self.originI
        originJ = self.originJ
        originX = self.originX
        way = self.way

        ###### 矩平 #########
        if (way == 2):
            for i in xrange(0,originI):
                AVF = 0.0
                for j in xrange(0,originJ):
                    AVF = AVF + originX[i,j]

                AVF = AVF/float(originJ)
                for j in xrange(0,originJ):
                    originX[i,j] = originX[i,j] - AVF

            del AVF
            ARCPY.AddMessage("矩平了")

        ##### 标准化 #####
        elif (way == 3):
            
            for i in xrange(0,originI):
                DF = 0.0

                for j in xrange(0,originJ):
                    DF=DF+float(originX[i,j]*originX[i,j])

                DF = MATH.sqrt(DF/float(originJ))

                if (DF == 0):
                    DF = 200.0 #### 测试用 ###

                for j in xrange(0,originJ):
                    originX[i,j]=originX[i,j]/DF

            del DF
            ARCPY.AddMessage("标准化了")

        ##### 矩平标准化 #####
        elif (way == 4):
            for i in xrange(0,originI):
                AVF = 0.0
                for j in xrange(0,originJ):
                    AVF = AVF + originX[i,j]

                AVF = AVF/float(originJ)
                for j in xrange(0,originJ):
                    originX[i,j] = originX[i,j] - AVF

            del AVF

            for i in xrange(0,originI):
                DF = 0.0

                for j in xrange(0,originJ):
                    DF=DF+float(originX[i,j]*originX[i,j])

                DF = MATH.sqrt(DF/float(originJ))

                if (DF == 0):
                    DF = 200.0 #### 测试用 ###

                for j in xrange(0,originJ):
                    originX[i,j]=originX[i,j]/DF

            del DF

            ARCPY.AddMessage("矩平标准化了")

            
            

        self.originX = originX


    def forma(self):
        """Finalizes original GA Table."""
        ###XXt或者XtX###
        originI = self.originI
        originJ = self.originJ
        P = copy.copy(self.originX)
        P = NUM.matrix(P)###矩阵运算 ###
        Pt = P.T
    
        if (originI > originJ):
            P = Pt*P 
        else:
            P = P*Pt
    
        self.P = P

        

    def jcb(self):

        originI = self.originI
        originJ = self.originJ
        minIJ = self.minIJ
        originX = copy.copy(self.originX)
        originX = NUM.matrix(originX)
        P = copy.copy(self.P)

        Q = NUM.linalg.eig(P)
        TZZ = Q[0]
        TZXL = NUM.matrix(Q[1].T)##因为numpy输出的特征向量格式问题要加T##

        

    
        if (originI > originJ):		##空间点数大于时间数##
            temptzxl = []
            TZXLL = NUM.zeros((minIJ,originI))##因为向量横过来存##
            TZXLL = NUM.matrix(TZXLL)
            for i in xrange(minIJ):
                temptzxlt = copy.copy(TZXL[i].T)	##横的V调成竖的V##
                temptzxl = originX*temptzxlt##V=X.Vr##
                TZXLt = temptzxl.T
                TZXLL[i] = copy.copy(TZXLt/MATH.sqrt(TZXLt*temptzxl))
            TZXL = TZXLL
            	
        self.TZZ = TZZ
        self.TZXL = TZXL


    def arrang(self):
        """Returns Ripley's Edge Correction value."""

        TZZ = copy.copy(self.TZZ)
        ##TZZ = map(float,TZZ)
        TZXL = copy.copy(self.TZXL)
        ##TZXL = map(float,TZXL)
        minIJ = self.minIJ
    
        #### ####
        TZlist = []
        for i in xrange(0,minIJ):
            TZlist.append([i,TZZ[i],TZXL[i]])
        
        #### 排序 ####
        TZlist = sorted(TZlist,key = lambda TZlist:TZlist[1])		###从小到大排列####
        TZlist.reverse()	####从大到小排列####
        
        self.TZlist = TZlist

    def tcoeff(self):
        """ """
        TZlist =copy.copy(self.TZlist)
        originI = self.originI
        originJ = self.originJ
        minIJ = self.minIJ
        originX = copy.copy(self.originX)
        originX = NUM.matrix(originX)

        
        ### 时间系数矩阵####
        Time=NUM.matrix(NUM.zeros((minIJ,originJ)))
        TZXL=NUM.matrix(NUM.zeros((minIJ,originI)))
        TZZ=[]
        ID=[]
        for i in xrange(minIJ):
            temp = TZlist[i]
            TZXL[i]=temp[2]
            TZZ.append(temp[1])
            ID.append(temp[0])
        TZXL=NUM.matrix(TZXL)
        for i in xrange(0,minIJ):
            temp = TZXL[i]*originX
            Time[i]=temp

    
        
        ### 方差贡献 ###
        sumTZZ = 0.0
        FCGX = []
        LJGX = []
        for i in xrange(0,minIJ):
            sumTZZ = sumTZZ + TZZ[i]
    
        ###方差贡献####
        for i in xrange(0,minIJ):
            GX = TZZ[i]/sumTZZ
            FCGX.append([ID[i],GX])
    
        ###累计方差贡献#####
        tempSUM = 0.0
        for i in xrange(0,minIJ):
            tempSUM = tempSUM + TZZ[i]
            GX = tempSUM/sumTZZ
            LJGX.append([ID[i],GX])
        
        self.Time = Time
        self.FCGX = FCGX
        self.LJGX = LJGX
        self.TZXLL = TZXL
        


    def createOutput(self):
        ""
        TZlist=self.TZlist

        TZXL = self.TZXLL
        Time = self.Time
        FCGX = self.FCGX
        LJGX = self.LJGX
        originI= self.originI
        originJ = self.originJ
        minIJ = self.minIJ
        nnn = "EOFys"


        if (self.way == 1):
            nnn = "EOFys"
        elif (self.way == 2):
            nnn = "EOFjp"
        elif (self.way == 3):
            nnn = "EOFbzh"
        elif (self.way == 4):
            nnn = "EOFjpbzh"

        for i in xrange (int(self.outputNum)):
            tmpname = nnn+str(i+1)
            try:
                ARCPY.AddField_management(self.inputFFCC, tmpname, "DOUBLE")
            except:
                ARCPY.AddMessage("【ERROR！！】:添加eof字段出错！！")
                raise SystemExit()
            

        ##try:
            ##ARCPY.AddField_management(self.inputFC, "EOF", "DOUBLE")
        ##except:
            ##ARCPY.AddMessage("【ERROR！！】:添加eof字段出错！！")
            ##raise SystemExit()

        for i in xrange(int(self.outputNum)):
            
            ttt = nnn+str(i+1)
            ii = 0
            ##ARCPY.AddMessage("【】【】【】【】【】【】【】【】")
            ##ARCPY.AddMessage(TZlist)
            tempp = TZlist[0+i][2]
            tempp1 = NUM.array(tempp, dtype = float)
            ##ARCPY.AddMessage(tempp1)
            ##ARCPY.AddMessage(tempp1.T)
            temp = map(float,tempp1.T)
        
            with DA.UpdateCursor(self.inputFFCC,(ttt)) as cursor:
                for row in cursor:
                    row = temp[ii]
                    try:
                        cursor.updateRow([row])
                    except:
                        ARCPY.AddMessage("【ERROR！！】:更新eof出错！！")
                        raise SystemExit()
                    ii+=1
            del ii


        if (self.displayTime == True):
            for i in xrange(int(self.outputNum)):
                y = map(float,Time[i].T)
                x = [1980, 1985, 1990, 1995, 2000, 2005, 2010]
                PL.figure(i+1)
                PL.plot(x,y)
                ##ARCPY.AddMessage(FCGX[i][1])
                if (i == 0):
                    aaa = "The 1st Mod With Variance Devoting Rates:"
                elif (i == 1):
                    aaa = "The 2nd Mod With Variance Devoting Rates:"
                elif (i == 2):
                    aaa = "The 3rd Mod With Variance Devoting Rates:"
                else :
                    aaa = "The " +str(i+1)+"th Mod With Variance Devoting Rates:"

                aaa = aaa + str((FCGX[i][1]*100000//10)/100) + "%"
                if (self.way == 1):
                    aaa += " (OriginData)"
                elif (self.way == 2):
                    aaa += " (RowBalance)"
                elif (self.way == 3):
                    aaa += " (Standardization)"
                elif (self.way == 4):
                    aaa += " (RowBalance & Standardization)"
                
                PL.title(aaa)
            PL.show()


        if (self.outputControl == True):
            listfile = OS.listdir(self.outputPath)
            tail = len(listfile)+1

            outFile = open(self.outputPath+"\\"+str(tail)+".txt",'w')

            for testIter in xrange(int(self.outputNum)):
                ldVal = FCGX[i][1]
                outLine = str(ldVal) + '\n'
                outFile.write(outLine)
        
               
            outFile.close()

        
            


            

    
        ##ARCPY.AddMessage(Time)
        self.code = True



        

        


        


if __name__ == "__main__":
    setupEOF()
