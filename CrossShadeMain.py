import sys

import maya.OpenMaya as OpenMaya
import maya.OpenMayaMPx as OpenMayaMPx
import maya.cmds as cmds
import maya.mel as mel

import numpy as np
import scipy as sp
import scipy.optimize


#--------------------------------------------------------------------
# ALGORITHM
#--------------------------------------------------------------------

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
  def __init__(self, i):
    self.id = i
    self.name = None
    self.nor = None
    self.ch = []   # sorted in t order along curve
    
  def __str__(self):
    return "cs(%s) %s" % (self.id, self.name)
    
class CrossHair:
  def __init__(self,i,j):
    self.i = i
    self.j = j
    self.t = -1
    self.pos = None
    self.tan = None
    self.nor = None
    
  def __str__(self):
    return "ch(%s,%s) @ t=%s" % (self.i, self.j, self.t)
  
csNum = 0       # cross section count number
cs = None    # cross section curves names
chNum = 0       # cross hair count number
ch = None    # cross hair t-parameter

def readCrossSections():
  global csNum
  global cs
  global chNum
  global ch

  #mel.eval("layerEditorSelectObjects layerCross;")
  #mel.eval("SelectAllNURBSCurves;")
  curveNames = mel.eval("ls -sl -type transform")

  # no curve drawn
  if not curveNames:
    return
  
  # init globals
  
  csNum = len(curveNames)
  cs = [None for x in range(csNum)]
  for i in range(csNum):
    cs[i] = CrossSection(i)
    cs[i].name = curveNames[i]
  
  chNum = 0
  ch = [[None]*csNum for x in range(csNum)]
  
  
  # process intersections
  
  for i in range(csNum):
    for j in range(i+1, csNum):

      rawIntersects = cmds.curveIntersect(cs[i].name, cs[j].name, useDirection=True, direction=(0,0,1))
      
      if rawIntersects:
        ch[i][j] = CrossHair(i,j)
        ch[j][i] = CrossHair(j,i)
      
        intersects = [float(k) for k in rawIntersects.split()]
        ch[i][j].t = intersects[0]
        ch[j][i].t = intersects[1]
        
        pi = cmds.pointOnCurve(cs[i].name, pr=ch[i][j].t, p=True)        
        ch[i][j].pos = np.array(pi)  # save
        ch[j][i].pos = np.array(pi)
        #cmds.spaceLocator(p=pi)   # place locator        
        
        t_ij = cmds.pointOnCurve(cs[i].name, pr=ch[i][j].t, nt=True)
        ch[i][j].tan = np.array(t_ij)
        
        t_ji = cmds.pointOnCurve(cs[j].name, pr=ch[j][i].t, nt=True)
        ch[j][i].tan = np.array(t_ji)
        
        
        print "(%s,%s) processed" % (i, j)
        
        chNum += 1
        
      else:
        print "(%s,%s) no intersect" % (i, j)
  
  # clear all selection
  mel.eval("select -cl")

  # store sorted list of ch for each cs
  for i in range(csNum):
    for j in range(csNum):
      if ch[i][j] is not None:
        cs[i].ch.append(ch[i][j])
        
    cs[i].ch = sorted(cs[i].ch, key=lambda ch: ch.t)
  
def printCSData():
  global csNum
  global cs
  global chNum
  global ch
  
  for i in range(csNum):
    print "%s : %s" % (i, cs[i].name)
  
  for i in range(csNum):
    for j in range(csNum):
      if ch[i][j] is not None:
        print "x_(%s,%s) : at %s : %s" % (i, j, ch[i][j].t, ch[i][j].pos)
      
  for i in range(csNum):
    for j in range(csNum):
      if ch[i][j] is not None:
        print "t_(%s,%s) : %s" % (i, j, ch[i][j].tan)

#consList
def minOptimize(): 
  global csNum
  global cs
  global chNum
  global ch
  
  EPSILON = 0.1
  
  # construct constraints
  global consList
  consList = []
  tanPairIdx = csNum*2
  for i in range(csNum):
    for j in range(i+1, csNum):
      if (ch[i][j] is not None):
        
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
            (ch[i][j].tan[0], ch[i][j].tan[1], ch[j][i].tan[0], ch[j][i].tan[1], tanPairIdx)
        })
        consList.append({
          "type": "ineq",
          "fun" : (lambda tijx, tijy, tjix, tjiy, tanPairIdx:
            lambda x: -(tijx*tjix + tijy*tjiy + x[tanPairIdx]*x[tanPairIdx+1]) + EPSILON)
            (ch[i][j].tan[0], ch[i][j].tan[1], ch[j][i].tan[0], ch[j][i].tan[1], tanPairIdx)
        })
        
        # t_ij . n_i = 0
        consList.append({
          "type": "eq",
          "fun" : (lambda tijx, tijy, i, tanPairIdx:
            lambda x: (tijx*x[i*2] + tijy*x[i*2+1] + x[tanPairIdx]) )
            (ch[i][j].tan[0], ch[i][j].tan[1], i, tanPairIdx)
        })
        # t_ji . n_j = 0
        consList.append({
          "type": "eq",
          "fun" : (lambda tjix, tjiy, j, tanPairIdx:
            lambda x: (tjix*x[j*2] + tjiy*x[j*2+1] + x[tanPairIdx+1]) )
            (ch[j][i].tan[0], ch[j][i].tan[1], j, tanPairIdx)
        })
        
        # print "constraints for (%s,%s)" % (i, j)
        
        tanPairIdx += 2
  
  # construct cost function
  funcString = "0"
  tanPairIdx = csNum*2
  for i in range(csNum):
    for j in range(i+1, csNum):
      if (ch[i][j] is not None):
        ni = "np.array([x["+str(i*2)+"], x["+str(i*2+1)+"], 1])"
        nj = "np.array([x["+str(j*2)+"], x["+str(j*2+1)+"], 1])"
        tij = "np.array(["+str(ch[i][j].tan[0])+", "+str(ch[i][j].tan[1])+", x["+str(tanPairIdx)+"]])"
        tji = "np.array(["+str(ch[j][i].tan[0])+", "+str(ch[j][i].tan[1])+", x["+str(tanPairIdx+1)+"]])"
        
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
    cs[i].nor = normalize( np.array([res.x[i*2], res.x[i*2+1], 1]) )    
    print "n_%s : %s" % (i, cs[i].nor)
  
  # store t_ij_z's
  tanPairIdx = csNum*2
  for i in range(csNum):
    for j in range(i+1, csNum):
      if (ch[i][j] is not None):
        ch[i][j].tan[2] = res.x[tanPairIdx]
        ch[j][i].tan[2] = res.x[tanPairIdx+1]
        
        ch[i][j].tan = normalize(ch[i][j].tan)
        ch[j][i].tan = normalize(ch[j][i].tan)
        
        #print "t_(%s,%s) : %s" % (i, j, ch[i][j].tan)
        #print "t_(%s,%s) : %s" % (j, i, ch[j][i].tan)
        
        tanPairIdx += 2
        
  # solve for n_ij
  for i in range(csNum):
    for j in range(i+1, csNum):
      if (ch[i][j] is not None):
        nor = np.cross(ch[i][j].tan, ch[j][i].tan)
        
        # normal flip hack
        if nor[2] < 0: nor = -nor
        
        ch[i][j].nor = normalize(nor)
        ch[j][i].nor = nor
        #print "n_(%s,%s) : %s" % (i, j, ch[i][j].nor)        
        #cmds.spaceLocator( p=(ch[i][j].pos+ch[i][j].nor).tolist() )

# get interpolated normal of ch_ij along curve i at t
def getCHNormAtT(chI, chJ, tparam):
  global csNum
  global cs
  global chNum
  global ch

  origN =  ch[chI][chJ].nor
  origT1 = ch[chI][chJ].tan
  origT2 = ch[chJ][chI].tan
  targetT = cmds.pointOnCurve(cs[chI].name, pr=tparam, nt=True)   # z comp incorrect
  targetT[2] = -(origT2[0]*targetT[0]+origT2[1]*targetT[1])/origT2[2]   # solve t_z, t2 dot t = 0
    
  crossProd = np.cross(origT1, targetT)
  axis = normalize(crossProd)
  angle = np.arcsin( np.linalg.norm(crossProd) / (np.linalg.norm(origT1)*np.linalg.norm(targetT)) )
  rot = rotation(angle,axis)
  
  return normalize(np.dot(rot, origN))
  
def createPatchMesh(vertices, normals):
  width = len(vertices)
  height = len(vertices[0])
  
  cmds.polyPlane(n="tp", sx=width-1, sy=height-1, ax=[0, 0, 1])
  
  for j in range(height):
    for i in range(width):
      cmds.select("tp.vtx[%s]" % (j*width+i), r=True)
      cmds.move(vertices[i][j][0], vertices[i][j][1], vertices[i][j][2], a=True)
      cmds.polyNormalPerVertex(xyz=normals[i][j].tolist())

  normalZ = cmds.polyInfo("tp.f[0]", fn=True)[0].split()[-1]
  if normalZ[0]=="-":
    cmds.polyNormal("tp", normalMode=0, userNormalMode=0, ch=1);
  
  cmds.rename("tp", "patch1")
  cmds.select(cl=True)

vertices = None
normals = None
# input 4 pairs (i,j) of crosshair intersections, connecting indices
def createCoonsPatch(cpairs):
  global csNum
  global cs
  global chNum
  global ch
  
  global vertices
  global normals

  # square patch dimension T_STEP by T_STEP
  T_STEPS = 10  
  
  vertices = [[None]*(T_STEPS) for x in range(T_STEPS)]
  normals = [[None]*(T_STEPS) for x in range(T_STEPS)]
  
  # ALONG CURVE
  
  # go from each ch to ch
  for p in range(4):
    cStart = ch[ cpairs[p][1] ][ cpairs[p][0] ]
    cEnd = None
    if p<3:
      cEnd = ch[ cpairs[p+1][0] ][ cpairs[p+1][1] ]
    else:      
      cEnd = ch[ cpairs[0][0] ][ cpairs[0][1] ]
    tStep = (cEnd.t-cStart.t)/(T_STEPS-1)
    
    # print "%s to %s : %s" % (cStart.t, cEnd.t, tStep)
    
    # go down curve
    
    t = cStart.t
    for step in range(T_STEPS-1):
      
      if step == 0:
        pos = cStart.pos
        nor = cStart.nor
      else:
        # get position
        pos = np.array(cmds.pointOnCurve(cs[cStart.i].name, pr=t, p=True))
        
        # get normal
        nor1 = getCHNormAtT(cStart.i, cStart.j, t)
        nor2 = getCHNormAtT(cEnd.i, cEnd.j, t)      
        blendT = (t-cStart.t) / (cEnd.t-cStart.t)
        nor = blendT*nor2+(1-blendT)*nor1
      
      if p == 0:
        coord = (step,0)
      elif p == 1:
        coord = (T_STEPS-1,step)
      elif p == 2:
        coord = (T_STEPS-1-step,T_STEPS-1)
      elif p == 3:
        coord = (0, T_STEPS-1-step)
      
      #print "(%s,%s)" % (coord[0], coord[1])
      vertices[coord[0]][coord[1]] = pos
      normals[coord[0]][coord[1]] = nor
      
      #cmds.spaceLocator( p=(pos+nor).tolist() )
      
      t = t + tStep
      
  
  # ALONG PATCH
  
  n = T_STEPS-1
  
  for i in range(1, T_STEPS-1):
    for j in range(1, T_STEPS-1):      
      
      vertices[i][j] = np.array([0.0,0.0,0.0])
      normals[i][j] = np.array([0.0,0.0,0.0])
      
      fi = float(i)
      fj = float(j)
      
      for k in range(3):
          
        vertices[i][j][k] = (
          (1.0-fi/n)*vertices[0][j][k] + fi/n*vertices[n][j][k] +
          (1.0-fj/n)*vertices[i][0][k] + fj/n*vertices[i][n][k] -
          (
            vertices[0][0][k]+(vertices[n][0][k]-vertices[0][0][k])*(fi/n) + 
            (
              (vertices[0][n][k]+(vertices[n][n][k]-vertices[0][n][k])*(fi/n)) - 
              (vertices[0][0][k]+(vertices[n][0][k]-vertices[0][0][k])*(fi/n))
            ) * (fj/n)
          )
        )
        
        normals[i][j][k] = (
          (1.0-fi/n)*normals[0][j][k] + fi/n*normals[n][j][k] +
          (1.0-fj/n)*normals[i][0][k] + fj/n*normals[i][n][k] -
          (
            normals[0][0][k]+(normals[n][0][k]-normals[0][0][k])*(fi/n) + 
            (
              (normals[0][n][k]+(normals[n][n][k]-normals[0][n][k])*(fi/n)) - 
              (normals[0][0][k]+(normals[n][0][k]-normals[0][0][k])*(fi/n))
            ) * (fj/n)
          )
        )

      # cmds.spaceLocator( p=(vertices[i][j]+normals[i][j]).tolist() )
    
  createPatchMesh(vertices, normals)

pairsList = None
def propagateCurve():
  global csNum
  global cs
  global chNum
  global ch
  
  global pairsList

  # CONSTRUCT GRID TOPOLOGY
  
  listX = [c.j for c in cs[0].ch]
  listXr = list(listX)
  listXr.reverse()
  listY = None

  # find two orthogonal axis
  for i in range(len(cs)):
    l = [c.j for c in cs[i].ch]
    if not (l==listX or l==listXr):
      listY = l
      break;

  width = len(listX)
  height = len(listY)

  """
  for j in range(height):
    line = ""
    for i in range(width):
      line = line + "%s " % ( [listA[i], listB[j]] )
    print line
  """
  
  # construct list of quad patches
  pairsList = []
  for j in range(height-1):
    for i in range(width-1):
      pairs = 4*[None]
      pairs[0] = [listX[i],   listY[j]]
      pairs[1] = [listX[i],   listY[j+1]]
      pairs[2] = [listX[i+1], listY[j+1]]
      pairs[3] = [listX[i+1], listY[j]]
      
      # order to connect
      # rearrange first pair
      if (pairs[0][0]==pairs[1][0] or pairs[0][0]==pairs[1][1]):
        pairs[0].reverse()
      # rearrange rest of chain
      for k in range(1,4):
        if (pairs[k][0] != pairs[k-1][1]):
          pairs[k].reverse()
      
      pairsList.append(pairs)
  
  for pairs in pairsList:
    print "patch %s" % (pairs)
    createCoonsPatch(pairs)
    

#--------------------------------------------------------------------
# COMMAND
#--------------------------------------------------------------------

kPluginCmdName = "CrossShadeCmd"

# command
class scriptedCommand(OpenMayaMPx.MPxCommand):
  def __init__(self):
    OpenMayaMPx.MPxCommand.__init__(self)
    
  def doIt(self, argList):
    global chNum

    readCrossSections()
    
    if chNum <= 0: 
      cmds.warning("No crosshairs found.")
      return      
    
    printCSData()
    minOptimize()
    propagateCurve()

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
