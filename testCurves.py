import maya.cmds as cmds
import maya.mel as mel
import numpy as np
from scipy.optimize import minimize

csN = 0       # cross section count number
csC = None    # cross section curves names
chN = 0       # cross hair count number
chP = None    # cross hair positions
chT = None    # cross hair tangents

def readCrossSections():
  global csN
  global csC
  global chN
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
        
        chN += 1
        
      else:
        print "(%s,%s) no intersect" % (ci, cj)
  
  # clear all selection
  mel.eval("select -cl")
        
def printCrossSectionData1():
  global csN
  global csC
  global chN
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

def minOptimize(): 
  global csN
  global csC
  global chN
  global chP
  global chT
  
  EPSILON = 0.1
  
  # construct constraints
  consList = []
  tanPairIdx = 0
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
            (chT[i][j][0]*chT[j][i][0] + chT[i][j][0]*chT[j][i][0] + x[csN*2+tanPairIdx]*x[csN*2+1+tanPairIdx]) + e,
          "test": "(chT[i][j][0]*chT[j][i][0] + chT[i][j][0]*chT[j][i][0] + x[%s]*x[%s]) + e" % (csN*2+tanPairIdx, csN*2+1+tanPairIdx)            
        })
        consList.append({
          "type": "ineq",
          "fun" : lambda x: 
            -(chT[i][j][0]*chT[j][i][0] + chT[i][j][0]*chT[j][i][0] + x[csN*2+tanPairIdx]*x[csN*2+1+tanPairIdx]) + e,
          "test": "-(chT[i][j][0]*chT[j][i][0] + chT[i][j][0]*chT[j][i][0] + x[%s]*x[%s]) + e" % (csN*2+tanPairIdx, csN*2+1+tanPairIdx)
        })
        
        # t_ij . n_i = 0
        consList.append({
          "type": "eq",
          "fun" : lambda x: (chT[i][j][0]*x[i*2] + chT[i][j][1]*x[i*2+1] + x[csN*2+tanPairIdx]),
          "test": "(chT[i][j][0]*x[i*2] + chT[i][j][1]*x[i*2+1] + x[%s])" % (csN*2+tanPairIdx)
        })
        # t_ji . n_j = 0
        consList.append({
          "type": "eq",
          "fun" : lambda x: (chT[j][i][0]*x[j*2] + chT[j][i][1]*x[j*2+1] + x[csN*2+tanPairIdx+1]),
          "test": "(chT[i][j][0]*x[i*2] + chT[i][j][1]*x[i*2+1] + x[%s])" % (csN*2+tanPairIdx+1)
        })
        
        print "cons for (%s,%s)" % (i, j)
        
        tanPairIdx += 2
    
  # print constraints
  print "CONSTRAINTS:"
  for i in consList:
    print i["test"]
  
  # construct cost function
  funcString = "0"
  tanPairIdx = 0
  for i in range(csN):
    for j in range(i+1, csN):
      if (chT[i][j]):
        ni = "np.array([x["+str(i*2)+"], x["+str(i*2+1)+"], 1])"
        nj = "np.array([x["+str(j*2)+"], x["+str(j*2+1)+"], 1])"
        tij = "np.array(["+str(chT[i][j][0])+", "+str(chT[i][j][1])+", x["+str(csN*2+tanPairIdx)+"]])"
        tji = "np.array(["+str(chT[j][i][0])+", "+str(chT[j][i][1])+", x["+str(csN*2+tanPairIdx+1)+"]])"
        
        sum = "((np.power(np.linalg.norm(np.cross("+tji+","+ni+")),2)) + (np.power(np.linalg.norm(np.cross("+tij+","+nj+")),2)) + np.power(x["+str(csN*2+tanPairIdx)+"],2) + np.power(x["+str(csN*2+tanPairIdx+1)+"],2))"
        
        funcString += "+"+sum
        
        tanPairIdx+=2
  
  print "FUNCTION:"
  print funcString
  
  # define function
  funcDefString = ("""def cmFunction(x):
  return %s""" % (funcString))
  exec(funcDefString)
  
  print cmFunction
  
  # initial guesses
  x0 = [-1 for x in range(csN*2+chN*2)]
  #res = minimize(cmFunction, x0, method='SLSQP', constraints=cons, options={'xtol': 1e-8, 'disp': True})
  
      
def runCrossShade():
  readCrossSections()
  printCrossSectionData1()
  minOptimize()
  
runCrossShade()