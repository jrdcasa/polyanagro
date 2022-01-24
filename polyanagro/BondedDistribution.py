import numpy as np
import polyanagro as pag


class BondedDistributions(object):

########################################################################
    def __init__(self):

        #Histogram distribution characteristics
        self._deltaBond = 0.005
        self._maxbinBond = 2000

        self._deltaAngle  = 0.5
        self._maxbinAngle = 360

        self._deltaDih  = 1
        self._maxbinDih = 360
#        self._deltaDih  = 1.0
#        self._maxbinDih = 360

        self._deltaImp  = 1.0
        self._maxbinImp = 360

        #Distribution files
        # self._fileBondDist  = open("Bond_Dist.dat",'w')
        # self._fileAngleDist = open("Angle_Dist.dat",'w')
        # self._fileDihDist   = open("Dih_Dist.dat",'w')
        # self._fileImpDist   = open("Imp_Dist.dat",'w')
        # self._fileDihDistFlory   = open("Dih_Dist_Flory.dat",'w')

        #Bin Histograms
        self._bdist = np.zeros([self._maxbinBond],dtype=np.float64)
        self._adist = np.zeros([self._maxbinAngle],dtype=np.float64)
        self._ddist = np.zeros([self._maxbinDih],dtype=np.float64)
        self._ddistFlory = np.zeros([self._maxbinDih],dtype=np.float64)
        self._idist = np.zeros([self._maxbinImp],dtype=np.float64)
        self._tacticitydist = np.zeros([self._maxbinDih],dtype=np.float64)

        self._bondHist = np.zeros([self._maxbinBond],dtype=np.int32)
        self._angleHist = np.zeros([self._maxbinAngle],dtype=np.int32)
        self._dihHist = np.zeros([self._maxbinDih],dtype=np.int32)
        self._dihHistFlory = np.zeros([self._maxbinDih],dtype=np.int32)
        self._impHist = np.zeros([self._maxbinImp],dtype=np.int32)
        self._tactHist = np.zeros([self._maxbinDih],dtype=np.int32)

        #Setup histograms
        self._setupHistograms()


########################################################################
    def _setupHistograms(self):

        pag.setup_hist_bondC(self._deltaBond, self._maxbinBond, self._bdist)
        pag.setup_hist_angleC(self._deltaAngle, self._maxbinAngle, self._adist)
        pag.setup_hist_dihC(self._deltaDih, self._maxbinDih, self._ddist)
        pag.setup_hist_dihC(self._deltaDih, self._maxbinDih, self._ddistFlory)
        pag.setup_hist_dihC(self._deltaImp, self._maxbinImp, self._idist)
        pag.setup_hist_dihC(self._deltaDih, self._maxbinDih, self._tacticitydist)

########################################################################
    def bondDist(self,bond2DArray,X1,Y1,Z1):

        iserror = 1
        if len(bond2DArray) != 0:
            iserror = pag.bondDistC(bond2DArray, X1, Y1, Z1, self._bondHist)
        return iserror

# ########################################################################
#     def angleDist(self,angle2DArray,X1,Y1,Z1):
#
#         iserror = 1
#         if len(angle2DArray) != 0:
#             iserror = distC.angleDistC(angle2DArray, X1, Y1, Z1, self._angleHist)
#         return iserror
#
# ########################################################################
#     def dihDist(self,dih2DArray,X1,Y1,Z1):
#
#         iserror = 1
#         if len(dih2DArray) != 0:
#             iserror = distC.dihDistC(dih2DArray, X1, Y1, Z1, self._dihHist)
#         return iserror
#
#
# ########################################################################
#     def dihDistFlory(self,dih2DArray,X1,Y1,Z1):
#
#         iserror = 1
#         if len(dih2DArray) != 0:
#             iserror = distC.dihDistFloryC(dih2DArray, X1, Y1, Z1, self._dihHistFlory)
#         return iserror
#
# ########################################################################
#     def impDist(self,imp2DArray,X1,Y1,Z1):
#
#         iserror = 1
#         if len(imp2DArray) != 0:
#             iserror = distC.dihDistC(imp2DArray, X1, Y1, Z1, self._impHist)
#         return iserror
#
# ########################################################################
#     def tacticity(self,tact2DArray,X1,Y1,Z1):
#
#         iserror = 1
#         if len(tact2DArray) != 0:
#             iserror = distC.tactDistC(tact2DArray, X1, Y1, Z1, self._tactHist)
#         return iserror
#
# ########################################################################
#     @staticmethod
#     def write_tacticity_hist(self):
#
#         """
#                                    (at4)
#                                  D1  |    D3
#            (at7)---(at5)---(at2)---(at1) --- (at3)---(at6)
#
#            D1 (R1) --> 3-1-2-5
#            D2 (R2) --> 1-2-5-7
#            D3 (R3) --> 6-3-1-2
#         """
#
#         # Calculate the histogram for tacticity
#         f = open('tacticity_D.dat', 'r')
#         valdih1 = []
#         valdih2 = []
#         valdih3 = []
#         valdih1_3 = []
#         while 1:
#             line = f.readline()
#             if not line: break
#             d1 = line.split()[0]
#             d2 = line.split()[1]
#             d3 = line.split()[2]
#             valdih1.append(float(d1))
#             valdih2.append(float(d2))
#             valdih3.append(float(d3))
#             valdih1_3.append(float(d1))
#             valdih1_3.append(float(d3))
#         f.close()
#
#         # hist = np.histogram(valdih,bins=1,range=[-180,180],normed=True)
#         #freq, bins = np.histogram(valdih, bins=360, range=[-180, 180], normed=True)
#         freq1, bins1 = np.histogram(valdih1, bins=360, range=[0, 360], density=True)
#         data1 = zip(bins1, freq1)
#         freq2, bins2 = np.histogram(valdih2, bins=360, range=[0, 360], density=True)
#         data2 = zip(bins2, freq2)
#         freq3, bins3 = np.histogram(valdih3, bins=360, range=[0, 360], density=True)
#         data3 = zip(bins3, freq3)
#         freq1_3, bins1_3 = np.histogram(valdih1_3, bins=360, range=[0, 360], density=True)
#         data1_3 = zip(bins1_3, freq1_3)
#
#         f = open('Dih_Dist_tactD1.dat', 'w')
#         for item in data1:
#             line= str(item[0])+" "+str(item[1])+"\n"
#             f.writelines(line)
#         f.close()
#         f = open('Dih_Dist_tactD2.dat', 'w')
#         for item in data2:
#             line= str(item[0])+" "+str(item[1])+"\n"
#             f.writelines(line)
#         f.close()
#
#         f = open('Dih_Dist_tactD3.dat', 'w')
#         for item in data3:
#             line= str(item[0])+" "+str(item[1])+"\n"
#             f.writelines(line)
#         f.close()
#
#         f = open('Dih_Dist_tactD1_3.dat', 'w')
#         for item in data1_3:
#             line= str(item[0])+" "+str(item[1])+"\n"
#             f.writelines(line)
#         f.close()
#
#         # Calculate the histogram for tacticity
#         f = open('tacticity_R.dat', 'r')
#         valdih1 = []
#         valdih2 = []
#         valdih3 = []
#         valdih1_3 = []
#         while 1:
#              line = f.readline()
#              if not line: break
#              d1 = line.split()[0]
#              d2 = line.split()[1]
#              d3 = line.split()[2]
#              valdih1.append(float(d1))
#              valdih2.append(float(d2))
#              valdih3.append(float(d3))
#              valdih1_3.append(float(d1))
#              valdih1_3.append(float(d3))
#
#         f.close()
#
#         # hist = np.histogram(valdih,bins=1,range=[-180,180],normed=True)
#         #freq, bins = np.histogram(valdih, bins=360, range=[-180, 180], normed=True)
#         freq1, bins1 = np.histogram(valdih1, bins=360, range=[0, 360], density=True)
#         data1 = zip(bins1, freq1)
#         freq2, bins2 = np.histogram(valdih2, bins=360, range=[0, 360], density=True)
#         data2 = zip(bins2, freq2)
#         freq3, bins3 = np.histogram(valdih3, bins=360, range=[0, 360], density=True)
#         data3 = zip(bins3, freq3)
#         freq1_3, bins1_3 = np.histogram(valdih1_3, bins=360, range=[0, 360], density=True)
#         data1_3 = zip(bins1_3, freq1_3)
#
#         f = open('Dih_Dist_tactR1.dat', 'w')
#         for item in data1:
#             line= str(item[0])+" "+str(item[1])+"\n"
#             f.writelines(line)
#         f.close()
#
#         f = open('Dih_Dist_tactR2.dat', 'w')
#         for item in data2:
#             line= str(item[0])+" "+str(item[1])+"\n"
#             f.writelines(line)
#         f.close()
#
#         f = open('Dih_Dist_tactR3.dat', 'w')
#         for item in data3:
#             line= str(item[0])+" "+str(item[1])+"\n"
#             f.writelines(line)
#         f.close()
#
#         f = open('Dih_Dist_tactR1_3.dat', 'w')
#         for item in data1_3:
#             line= str(item[0])+" "+str(item[1])+"\n"
#             f.writelines(line)
#         f.close()
#
#
# ########################################################################
#     def writeDist(self):
#
#         for i in range(0,self._maxbinBond):
#             if np.sum(self._bondHist) == 0: break
#             self._fileBondDist.writelines(str(self._bdist[i])+" "+str(self._bondHist[i])+" "+\
#                                           str(self._bondHist[i]*(1.0/np.sum(self._bondHist))*(1.0/self._deltaBond))+"\n")
#         for i in range(0,self._maxbinAngle):
#             if np.sum(self._angleHist) == 0: break
#             self._fileAngleDist.writelines(str(self._adist[i])+" "+str(self._angleHist[i])+" "+\
#                                            str(self._angleHist[i]*(1.0/np.sum(self._angleHist))*(1.0/self._deltaAngle))+"\n")
#         for i in range(0,self._maxbinDih):
#             if np.sum(self._dihHist) != 0:
#                 self._fileDihDist.writelines(str(self._ddist[i])+" "+str(self._dihHist[i])+" "+\
#                                              str(self._dihHist[i]*(1.0/np.sum(self._dihHist))*(1.0/self._deltaDih))+"\n")
#         for i in range(0,self._maxbinImp):
#             if np.sum(self._impHist) != 0:
#                 self._fileImpDist.writelines(str(self._idist[i])+" "+str(self._impHist[i])+" "+\
#                                              str(self._impHist[i]*(1.0/np.sum(self._impHist))*(1.0/self._deltaImp))+"\n")
#
#         for i in range(0,self._maxbinDih):
#             if np.sum(self._dihHistFlory) != 0:
#                 self._fileDihDistFlory.writelines(str(self._ddistFlory[i])+" "+str(self._dihHistFlory[i])+" "+\
#                                                   str(self._dihHistFlory[i]*(1.0/np.sum(self._dihHistFlory))*(1.0/self._deltaDih))+"\n")
#
#
#
#         self._fileBondDist.close()
#         self._fileAngleDist.close()
#         self._fileDihDist.close()
#         self._fileDihDistFlory.close()
#         self._fileImpDist.close()
