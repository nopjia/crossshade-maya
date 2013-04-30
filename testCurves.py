import maya.cmds as cmds
import maya.mel as mel
import numpy as np
import scipy as sp

def drange(start, stop, step):
  r = start
  while r < stop:
  	yield r
  	r += step

def normalize(v):
  return v/np.linalg.norm(v)

# angle t
# axis u
def rotation(t, u):
  return [
    [np.cos(t)+u[0]*u[0], u[0]*u[1]*(1-np.cos(t))-u[2]*np.sin(t), u[0]*u[2]*(1-np.cos(t))+u[1]*np.sin(t)],
    [u[1]*u[0]*(1-np.cos(t))+u[2]*np.sin(t), np.cos(t)+u[1]*u[1], u[1]*u[2]*(1-np.cos(t))-u[0]*np.sin(t)],
    [u[2]*u[0]*(1-np.cos(t))-u[1]*np.sin(t), u[2]*u[1]*(1-np.cos(t))+u[0]*np.sin(t), np.cos(t)+u[2]*u[2]]
  ]

class CrossSection:
  def __init__(self):
    self.name = None
    self.nor = None
    
class CrossHair:
  def __init__(self):
    self.t = -1
    self.pos = None
    self.tan = None
    self.nor = None
  
csNum = 0       # cross section count number
cs = None    # cross section curves names
chNum = 0       # cross hair count number
ch = None    # cross hair t-parameter

def readCrossSections():
  global csNum
  global cs
  global chNum
  global ch

  mel.eval("layerEditorSelectObjects layerCross;")
  curveNames = mel.eval("ls -sl -type transform")

  # no curve drawn
  if not curveNames:
    return
  
  # init globals
  
  csNum = len(curveNames)
  cs = [CrossSection() for x in range(csNum)]
  for i in range(csNum):
    cs[i].name = curveNames[i]
  
  chNum = 0
  ch = [[] for x in range(csNum)]
  for i in range(csNum):
    ch[i] = [CrossHair() for x in range(csNum)]

  
  # process intersections
  
  for i in range(csNum):
    for j in range(i+1, csNum):

      rawIntersects = cmds.curveIntersect(cs[i].name, cs[j].name, useDirection=True, direction=(0,0,1))
      
      if rawIntersects:
        intersects = [float(k) for k in rawIntersects.split()]
        ch[i][j].t = intersects[0]
        ch[j][i].t = intersects[1]
        
        pi = cmds.pointOnCurve(csNam[i], pr=ch[i][j].t, p=True)        
        ch[i][j].pos = np.array(pi)  # save
        ch[j][i].pos = np.array(pi)
        #cmds.spaceLocator(p=pi)   # place locator        
        
        t_ij = cmds.pointOnCurve(csNam[i], pr=ch[i][j].t, nt=True)
        ch[i][j].tan = np.array(t_ij)
        
        t_ji = cmds.pointOnCurve(csNam[j], pr=ch[j][i].t, nt=True)
        ch[j][i].tan = np.array(t_ji)
        
        
        print "(%s,%s) processed" % (i, j)
        
        chNum += 1
        
      else:
        print "(%s,%s) no intersect" % (i, j)
  
  # clear all selection
  mel.eval("select -cl")
        
def printCrossSectionData1():
  global csNum
  global cs
  global chNum
  global ch
  
  for i in range(csNum):
    print "%s : %s" % (i, cs[i].name)
  
  for i in range(csNum):
    for j in range(csNum):
      print "x_(%s,%s) : at %s : %s" % (i, j, ch[i][j].t, ch[i][j].pos)
      
  for i in range(csNum):
    for j in range(csNum):
      print "t_(%s,%s) : %s" % (i, j, ch[i][j].tan)
      
readCrossSections()
printCrossSectionData1()

#consList
def minOptimize(): 
  global csNum
  global csNam
  global csNor
  global chNum
  global chTee
  global chPos
  global chTan
  global chNor
  
  EPSILON = 0.1
  
  # construct constraints
  global consList
  consList = []
  tanPairIdx = csNum*2
  for i in range(csNum):
    for j in range(i+1, csNum):
      if (chTan[i][j] is not None):
        
        # -e < n_i . n_j < e
        consList.append({
          "type": "ineq",
          "fun" : (lambda i,j: lambda x: (x[i*2]*x[j*2] + x[i*2+1]*x[j*2+1] + 1.0) + EPSILON)(i,j)
        })
        consList.append({
          "type": "ineq",
          "fun" : (lambda i,j: lambda x: -(x[i*2]*x[j*2] + x[i*2+1]*x[j*2+1] + 1.0) + EPSILON)(i,j)
        })
        
        # -e < t_ij . t_ji < e
        consList.append({
          "type": "ineq",
          "fun" : (lambda tijx, tijy, tjix, tjiy, tanPairIdx:
            lambda x: (tijx*tjix + tijy*tjiy + x[tanPairIdx]*x[tanPairIdx+1]) + EPSILON)
            (chTan[i][j][0], chTan[i][j][1], chTan[j][i][0], chTan[j][i][1], tanPairIdx)
        })
        consList.append({
          "type": "ineq",
          "fun" : (lambda tijx, tijy, tjix, tjiy, tanPairIdx:
            lambda x: -(tijx*tjix + tijy*tjiy + x[tanPairIdx]*x[tanPairIdx+1]) + EPSILON)
            (chTan[i][j][0], chTan[i][j][1], chTan[j][i][0], chTan[j][i][1], tanPairIdx)
        })
        
        # t_ij . n_i = 0
        consList.append({
          "type": "eq",
          "fun" : (lambda tijx, tijy, i, tanPairIdx:
            lambda x: (tijx*x[i*2] + tijy*x[i*2+1] + x[tanPairIdx]) )
            (chTan[i][j][0], chTan[i][j][1], i, tanPairIdx)
        })
        # t_ji . n_j = 0
        consList.append({
          "type": "eq",
          "fun" : (lambda tjix, tjiy, j, tanPairIdx:
            lambda x: (tjix*x[j*2] + tjiy*x[j*2+1] + x[tanPairIdx+1]) )
            (chTan[j][i][0], chTan[j][i][1], j, tanPairIdx)
        })
        
        # print "constraints for (%s,%s)" % (i, j)
        
        tanPairIdx += 2
  
  # construct cost function
  funcString = "0"
  tanPairIdx = csNum*2
  for i in range(csNum):
    for j in range(i+1, csNum):
      if (chTan[i][j] is not None):
        ni = "np.array([x["+str(i*2)+"], x["+str(i*2+1)+"], 1])"
        nj = "np.array([x["+str(j*2)+"], x["+str(j*2+1)+"], 1])"
        tij = "np.array(["+str(chTan[i][j][0])+", "+str(chTan[i][j][1])+", x["+str(tanPairIdx)+"]])"
        tji = "np.array(["+str(chTan[j][i][0])+", "+str(chTan[j][i][1])+", x["+str(tanPairIdx+1)+"]])"
        
        sum = "((np.power(np.linalg.norm(np.cross("+tji+","+ni+")),2)) + (np.power(np.linalg.norm(np.cross("+tij+","+nj+")),2)) + np.power(x["+str(tanPairIdx)+"],2) + np.power(x["+str(tanPairIdx+1)+"],2))"
        
        funcString += "+"+sum
        
        tanPairIdx += 2
  
  print "FUNCTION:"
  print funcString
  
  # define function
  exec ("""def cmFunction(x):
  return %s""" % (funcString)) in globals(), locals()
  
  # run optimization
  x0 = [-1 for x in range(csNum*2+chNum*2)]  # initial guess
  res = sp.optimize.minimize(cmFunction, x0, method='SLSQP', constraints=tuple(consList), options={'disp': True})
  
  # store cross section plane normals
  for i in range(csNum):
    csNor[i] = normalize( np.array([res.x[i*2], res.x[i*2+1], 1]) )    
    print "n_%s : %s" % (i, csNor[i])
  
  # store t_ij_z's
  tanPairIdx = csNum*2
  for i in range(csNum):
    for j in range(i+1, csNum):
      if (chTan[i][j] is not None):
        chTan[i][j][2] = res.x[tanPairIdx]
        chTan[j][i][2] = res.x[tanPairIdx+1]
        
        chTan[i][j] = normalize(chTan[i][j])
        chTan[j][i] = normalize(chTan[j][i])
        
        print "t_(%s,%s) : %s" % (i, j, chTan[i][j])
        print "t_(%s,%s) : %s" % (j, i, chTan[j][i])
        
        tanPairIdx += 2
        
  # solve for n_ij
  for i in range(csNum):
    for j in range(i+1, csNum):
      if (chTan[i][j] is not None):
        nor = np.cross(chTan[i][j], chTan[j][i])
        
        # normal flip hack
        if nor[2] < 0: nor = -nor
        
        chNor[i][j] = normalize(nor)
        print "n_(%s,%s) : %s" % (i, j, chNor[i][j])
        
        cmds.spaceLocator( p=(chPos[i][j]+chNor[i][j]).tolist() )

# get interpolated normal of ch_ij along curve i at t
def getCHNormAtT(chI, chJ, tparam):
  global csNum
  global csNam
  global csNor
  global chNum
  global chTee
  global chPos
  global chTan
  global chNor

  origN = chNor[chI][chJ]
  origT1 = chTan[chI][chJ]
  origT2 = chTan[chJ][chI]
  targetT = cmds.pointOnCurve(csNam[chI], pr=tparam, nt=True)   # z comp incorrect
  targetT[2] = -(origT2[0]*targetT[0]+origT2[1]*targetT[1])/origT2[2]   # solve t_z, t2 dot t = 0
    
  crossProd = np.cross(origT1, targetT)
  axis = normalize(crossProd)
  angle = np.arcsin( np.linalg.norm(crossProd) / (np.linalg.norm(origT1)*np.linalg.norm(targetT)) )
  rot = rotation(angle,axis)
  
  return normalize(np.dot(rot, origN))
        
def propagateCurve():
  global csNum
  global csNam
  global csNor
  global chNum
  global chTee
  global chPos
  global chTan
  global chNor
  
  T_STEP = .1

def propagatePatch():
  print "propagate coons patch"
        
#----------------------------------------------------------
# RUN COMMAND
#----------------------------------------------------------
        
def runCrossShade():
  readCrossSections()
  printCrossSectionData1()
  minOptimize()
  propagateCurve()
  propagatePatch()
  
runCrossShade()