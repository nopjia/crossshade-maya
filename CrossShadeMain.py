import sys

import maya.OpenMaya as OpenMaya
import maya.OpenMayaMPx as OpenMayaMPx
import maya.cmds as cmds
import maya.mel as mel

import numpy as np
from scipy.optimize import minimize


#--------------------------------------------------------------------
# ALGORITHM
#--------------------------------------------------------------------

def normalize(v):
  return v/np.linalg.norm(v)

csNum = 0       # cross section count number
csNam = None    # cross section curves names
csNor = None    # cross section normals
chNum = 0       # cross hair count number
chPos = None    # cross hair positions
chTan = None    # cross hair tangents
chNor = None    # cross hair normals

def readCrossSections():
  global csNum
  global csNam
  global csNor
  global chNum
  global chPos
  global chTan
  global chNor

  mel.eval("layerEditorSelectObjects layerCross;")
  csNam = mel.eval("ls -sl -type transform")

  # no curve drawn
  if not csNam:
    return
  
  # init globals
  
  csNum = len(csNam)
  csNor = [None for x in range(csNum)]
  chNum = 0
  chPos = [[None]*csNum for x in range(csNum, 0, -1)]  
  chTan = [[None]*csNum for x in range(csNum, 0, -1)]
  chNor = [[None]*csNum for x in range(csNum, 0, -1)]

  
  # process intersections
  
  for ci in range(csNum):
    for cj in range(ci+1, csNum):

      rawIntersects = cmds.curveIntersect(csNam[ci], csNam[cj], useDirection=True, direction=(0,0,1))
      
      if rawIntersects:
        intersects = [float(i) for i in rawIntersects.split()]
        intT_i = intersects[0]
        intT_j = intersects[1]
        
        pi = cmds.pointOnCurve(csNam[ci], pr=intT_i, p=True)        
        chPos[ci][cj] = np.array(pi)  # save
        #cmds.spaceLocator(p=pi)   # place locator        
        
        t_ij = cmds.pointOnCurve(csNam[ci], pr=intT_i, nt=True)
        chTan[ci][cj] = np.array(t_ij)
        #cmds.spaceLocator( p=(chPos[ci][cj]+chTan[ci][cj]).tolist() )
        t_ji = cmds.pointOnCurve(csNam[cj], pr=intT_j, nt=True)
        chTan[cj][ci] = np.array(t_ji)
        #cmds.spaceLocator( p=(chPos[ci][cj]+chTan[cj][ci]).tolist() )
        
        print "(%s,%s) processed" % (ci, cj)
        
        chNum += 1
        
      else:
        print "(%s,%s) no intersect" % (ci, cj)
  
  # clear all selection
  mel.eval("select -cl")
        
def printCrossSectionData1():
  global csNum
  global csNam
  global csNor
  global chNum
  global chPos
  global chTan
  global chNor
  
  for i in range(csNum):
    print "%s : %s" % (i, csNam[i])
  
  for i in range(csNum):
    for j in range(csNum):
      print "x_(%s,%s) : %s" % (i, j, chPos[i][j])
      
  for i in range(csNum):
    for j in range(csNum):
      print "t_(%s,%s) : %s" % (i, j, chTan[i][j])

def minOptimize(): 
  global csNum
  global csNam
  global csNor
  global chNum
  global chPos
  global chTan
  global chNor
  
  EPSILON = 0.1
  
  # construct constraints
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
        
        print "cons for (%s,%s)" % (i, j)
        
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
  res = minimize(cmFunction, x0, method='SLSQP', constraints=tuple(consList), options={'disp': True})
  
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


#--------------------------------------------------------------------
# COMMAND
#--------------------------------------------------------------------

kPluginCmdName = "CrossShadeCmd"

# command
class scriptedCommand(OpenMayaMPx.MPxCommand):
  def __init__(self):
    OpenMayaMPx.MPxCommand.__init__(self)
    
  def doIt(self, argList):
    readCrossSections()
    printCrossSectionData1()
    minOptimize()

# Creator
def cmdCreator():
  return OpenMayaMPx.asMPxPtr( scriptedCommand() )
  
# Initialize the script plug-in
def initializePlugin(mobject):
  mplugin = OpenMayaMPx.MFnPlugin(mobject)
  
  # call UI MEL file
  mel.eval( "source \"" + mplugin.loadPath() + "/ui/mainUI.mel\"" )
  mel.eval( "createCrossShadeUI" )
  
  try:
    mplugin.registerCommand( kPluginCmdName, cmdCreator )
  except:
    sys.stderr.write( "Failed to register command: %s\n" % kPluginCmdName )
    raise

# Uninitialize the script plug-in
def uninitializePlugin(mobject):
  mplugin = OpenMayaMPx.MFnPlugin(mobject)
  
  mel.eval( "deleteCrossShadeUI" )
  
  try:
    mplugin.deregisterCommand( kPluginCmdName )
  except:
    sys.stderr.write( "Failed to unregister command: %s\n" % kPluginCmdName )
    raise
