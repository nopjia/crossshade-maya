import maya.cmds as cmds
import maya.mel as mel
import numpy as np
from scipy.optimize import minimize

csN = 0       # cross section count number
csC = None    # cross section curves names
chP = None    # cross hair positions
chT = None    # cross hair tangents

def readCrossSections():
  global csN
  global csC 
  global chP  
  global chT  

  mel.eval("layerEditorSelectObjects layerCross;")
  csC = mel.eval("ls -sl -type transform")

  # no curve drawn
  if not csC:
    return

  csN = len(csC)
  chP = [[None]*csN for x in range(csN, 0, -1)]  
  chT = [[None]*csN for x in range(csN, 0, -1)]

  for ci in range(csN):
    for cj in range(ci+1, csN):

      rawIntersects = cmds.curveIntersect(csC[ci], csC[cj], useDirection=True, direction=(0,0,1))
      
      if rawIntersects:
        intersects = [float(i) for i in rawIntersects.split()]
        intT_i = intersects[0]
        intT_j = intersects[1]
        
        pi = cmds.pointOnCurve(csC[ci], pr=intT_i, p=True)        
        chP[ci][cj] = pi;  # save
        cmds.spaceLocator(p=pi)   # place locator        
        
        t_ij = cmds.pointOnCurve(csC[ci], pr=intT_i, nt=True)
        chT[ci][cj] = t_ij;
        cmds.spaceLocator(p=[ p+t for p,t in zip(pi,t_ij) ])
        t_ji = cmds.pointOnCurve(csC[cj], pr=intT_j, nt=True)
        chT[cj][ci] = t_ji;
        cmds.spaceLocator(p=[ p+t for p,t in zip(pi,t_ji) ])
        
        print "(%s,%s) processed" % (ci, cj)
        
      else:
        print "(%s,%s) no intersect" % (ci, cj)
  
  # clear all selection
  mel.eval("select -cl")
        
def printCrossSectionData1():
  global csN
  global csC
  global chP
  global chT
  
  for i in range(csN):
    print "%s : %s" % (i, csC[i])
  
  for i in range(csN):
    for j in range(csN):
      print "x_(%s,%s) : %s" % (i, j, chP[i][j])
      
  for i in range(csN):
    for j in range(csN):
      print "t_(%s,%s) : %s" % (i, j, chT[i][j])

consList = None
def minOptimize(): 
  global csN
  global csC
  global chP
  global chT
  
  EPSILON = 0.1
  
  # construct constraints
  global consList
  consList = []
  interIdx = 0
  for i in range(csN):
    for j in range(i+1, csN):
      if (chT[i][j]):
        
        # -e < n_i . n_j < e
        consList.append({
          "type": "ineq",
          "fun" : lambda x: (x[i*2]*x[j*2] + x[i*2+1]*x[j*2+1] + 1.0) + e,
          "test" : "(x[%s]*x[%s] + x[%s]*x[%s] + 1.0) + e" % (i*2,j*2,i*2+1,j*2+1)
        })
        consList.append({
          "type": "ineq",
          "fun" : lambda x: -(x[i*2]*x[j*2] + x[i*2+1]*x[j*2+1] + 1.0) + e,
          "test" : "-(x[%s]*x[%s] + x[%s]*x[%s] + 1.0) + e" % (i*2,j*2,i*2+1,j*2+1)
        })
        
        # -e < t_ij . t_ji < e
        consList.append({
          "type": "ineq",
          "fun" : lambda x: 
            (chT[i][j][0]*chT[j][i][0] + chT[i][j][0]*chT[j][i][0] + x[csN*2+interIdx]*x[csN*2+1+interIdx]) + e,
          "test": "(chT[i][j][0]*chT[j][i][0] + chT[i][j][0]*chT[j][i][0] + x[%s]*x[%s]) + e" % (csN*2+interIdx, csN*2+1+interIdx)
            
        })
        consList.append({
          "type": "ineq",
          "fun" : lambda x: 
            -(chT[i][j][0]*chT[j][i][0] + chT[i][j][0]*chT[j][i][0] + x[csN*2+interIdx]*x[csN*2+1+interIdx]) + e,
          "test": "-(chT[i][j][0]*chT[j][i][0] + chT[i][j][0]*chT[j][i][0] + x[%s]*x[%s]) + e" % (csN*2+interIdx, csN*2+1+interIdx)
        })
        
        print "cons for (%s,%s)" % (i, j)
        
        interIdx+=1
      
def runCrossShade():
  readCrossSections()
  printCrossSectionData1()
  
runCrossShade()