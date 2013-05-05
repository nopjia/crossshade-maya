import sys

import maya.OpenMaya as OpenMaya
import maya.OpenMayaMPx as OpenMayaMPx
import maya.cmds as cmds
import maya.mel as mel

import numpy as np
import scipy as sp
import scipy.optimize


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

vertices = None # vertex locations

class CrossHairSingle:
  def __init__(self):
    self.i = -1
    self.iName = None
    self.j = -1
    self.jName = None
    self.iT = -1
    self.jT = -1
    self.pos = None
    self.nor = None
    self.dist = 999

  def __str__(self):
    return "cs(%s,%s)" % (self.i, self.j)

class Region:
  def __init__(self,i):
    self.ID = i
    self.corners = None

  def __str__(self):
    return "cs(%s)" % (self.ID)

chs = None
regions = None
numReg = 0
grid = None

def readCrossSections():
  global csNum
  global cs
  global chNum
  global ch

  global chs
  global regions

  mel.eval("layerEditorSelectObjects layerCross;")
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

  chs = [None for x in range (csNum*csNum)]
  
  
  # process intersections
  
  for i in range(csNum):
    for j in range(i+1, csNum):

      rawIntersects = cmds.curveIntersect(cs[i].name, cs[j].name, useDirection=True, direction=(0,0,1))
      
      if rawIntersects:
        ch[i][j] = CrossHair(i,j)
        ch[j][i] = CrossHair(j,i)
 
        # Build 1D cross hair array
        chs[chNum] = CrossHairSingle()
        chs[chNum].i = i
        chs[chNum].j = j
        chs[chNum].iName = cs[i].name
        chs[chNum].jName = cs[j].name
      
        intersects = [float(k) for k in rawIntersects.split()]
        ch[i][j].t = intersects[0]
        ch[j][i].t = intersects[1]
        chs[chNum].iT = intersects[0]
        chs[chNum].jT = intersects[1]
        
        pi = cmds.pointOnCurve(cs[i].name, pr=ch[i][j].t, p=True)   
        ch[i][j].pos = np.array(pi)  # save
        ch[j][i].pos = np.array(pi)
        chs[chNum].pos = np.array(pi)
        cmds.spaceLocator(p=pi)   # place locator        
        
        t_ij = cmds.pointOnCurve(cs[i].name, pr=ch[i][j].t, nt=True)
        ch[i][j].tan = np.array(t_ij)
        
        t_ji = cmds.pointOnCurve(cs[j].name, pr=ch[j][i].t, nt=True)
        ch[j][i].tan = np.array(t_ji)
        
        
        print "[%s,%s] processed at (%s,%s)" % (i, j, pi[0], pi[1])  
        
        chNum += 1
        
      else:
        print "[%s,%s] no intersect" % (i, j)
  
  # clear all selection
  mel.eval("select -cl")

  # store sorted list of ch for each cs
  for i in range(csNum):
    for j in range(csNum):
      if ch[i][j] is not None:
        cs[i].ch.append(ch[i][j])
        
    cs[i].ch = sorted(cs[i].ch, key=lambda ch: ch.t)

  for i in range(csNum):
    for j in range(csNum):
      if ch[i][j]:
        print "[%s,%s] at (%s, %s)" % (ch[i][j].i, ch[i][j].j, ch[i][j].pos[0], ch[i][j].pos[1])

  # Print all cross hairs in chs
  for x in range(chNum):
    print "Cross Hair Single [%s,%s] at (%s, %s)" % (chs[x].i, chs[x].j, chs[x].pos[0], chs[x].pos[1])
  
def printCrossSectionData1():
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
        
        print "t_(%s,%s) : %s" % (i, j, ch[i][j].tan)
        print "t_(%s,%s) : %s" % (j, i, ch[j][i].tan)
        
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
        print "n_(%s,%s) : %s" % (i, j, ch[i][j].nor)
        
        # cmds.spaceLocator( p=(ch[i][j].pos+ch[i][j].nor).tolist() )

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

# input np.arrays
def createPatchMesh(vertices, normals):
  width = len(vertices)
  height = len(vertices[0])

  numFaces = (width-1)*(height-1)
  numVertices = width*height
  numFaceConnects = 4*numFaces

  # enumerate x first for each y
  mVertices = OpenMaya.MFloatPointArray()  
  for j in range(height):
    for i in range(width):
      mVertices.append( OpenMaya.MFloatPoint(vertices[i][j][0], vertices[i][j][1], vertices[i][j][2]) )

  toIdx = lambda x,y: y*width + x

  # construct face connects
  faceConnects = []
  for i in range(width-1):
    for j in range(height-1):
      faceConnects.append(toIdx( i  , j  ))
      faceConnects.append(toIdx( i+1, j  ))
      faceConnects.append(toIdx( i+1, j+1))
      faceConnects.append(toIdx( i  , j+1))
  mFaceConnects = OpenMaya.MIntArray()
  scriptUtil = OpenMaya.MScriptUtil()
  scriptUtil.createIntArrayFromList( faceConnects,  mFaceConnects )

  mFaceCounts = OpenMaya.MIntArray(numFaces, 4)

  # create mesh
  meshFS = OpenMaya.MFnMesh()
  meshFS.create(numVertices, numFaces, mVertices, mFaceCounts, mFaceConnects)
  
  # set normals
  
  mNormals = OpenMaya.MFloatVectorArray()
  for j in range(height):
    for i in range(width):
      mNormals.append( OpenMaya.MFloatVector(normals[i][j][0], normals[i][j][1], normals[i][j][2]) )
  meshFS.setNormals(mNormals)

  # set face normals

  # assign lambert shader
  mel.eval("string $myLambert = `shadingNode -asShader lambert`;")
  mel.eval("setAttr lambert2.color 0.279 0.275 0.5;")
  poly = cmds.ls(type='mesh')
  cmds.select(poly)
  mel.eval("hyperShade -assign $myLambert;")

  # reverse normals
  cmds.select(poly)
  mel.eval("ReversePolygonNormals;")
  #mel.eval("polySetToFaceNormal;")
  #mel.eval("polySoftEdge -a 180 -ch 1 polySurface1;")

def createLight():
  # create light
  mel.eval("pointLight;")
  mel.eval("setAttr pointLight1.translateZ 5.0")
  mel.eval("setAttr pointLight1.translateY 5.0")


def testCreatePatch():
  width = 6
  height = 4
  
  vertices = [[None]*height for x in range(width)]
  for i in range(width):
    for j in range(height):
      vertices[i][j] = np.array([i*0.75, j*0.75, 0.0])
      
  normals = [[None]*height for x in range(width)]
  for i in range(width):
    for j in range(height):
      normals[i][j] = normalize( np.array([ 0.0, float(height-j)/height, 1.0-float(height-j)/height]) )

  createPatchMesh(vertices, normals)

def getMaxI(x):
  mi1 = cmds.pointOnCurve(chs[x].iName, pr=cmds.getAttr(chs[x].iName + '.maxValue'), p=True)
  mi2 = cmds.pointOnCurve(chs[x].iName, pr=cmds.getAttr(chs[x].iName + '.minValue'), p=True)
  if mi1[0] > mi2[0]:
    maxI = cmds.pointOnCurve(chs[x].iName, pr=cmds.getAttr(chs[x].iName + '.maxValue'), p=True)   
  else:
    maxI = cmds.pointOnCurve(chs[x].iName, pr=cmds.getAttr(chs[x].iName + '.minValue'), p=True) 

  return maxI 


def getMaxJ(x):
  mj1 = cmds.pointOnCurve(chs[x].jName, pr=cmds.getAttr(chs[x].jName + '.maxValue'), p=True)  
  mj2 = cmds.pointOnCurve(chs[x].jName, pr=cmds.getAttr(chs[x].jName + '.minValue'), p=True)  
  if mj1[0] > mj2[0]:
    maxJ = cmds.pointOnCurve(chs[x].jName, pr=cmds.getAttr(chs[x].jName + '.maxValue'), p=True)  
  else:
    maxJ = cmds.pointOnCurve(chs[x].jName, pr=cmds.getAttr(chs[x].jName + '.minValue'), p=True) 

  return maxJ

def getMinI(x):
  mi1 = cmds.pointOnCurve(chs[x].iName, pr=cmds.getAttr(chs[x].iName + '.maxValue'), p=True)
  mi2 = cmds.pointOnCurve(chs[x].iName, pr=cmds.getAttr(chs[x].iName + '.minValue'), p=True)
  if mi1[0] < mi2[0]:
    maxI = cmds.pointOnCurve(chs[x].iName, pr=cmds.getAttr(chs[x].iName + '.maxValue'), p=True)   
  else:
    maxI = cmds.pointOnCurve(chs[x].iName, pr=cmds.getAttr(chs[x].iName + '.minValue'), p=True) 

  return maxI 


def getMinJ(x):
  mj1 = cmds.pointOnCurve(chs[x].jName, pr=cmds.getAttr(chs[x].jName + '.maxValue'), p=True)  
  mj2 = cmds.pointOnCurve(chs[x].jName, pr=cmds.getAttr(chs[x].jName + '.minValue'), p=True)  
  if mj1[0] < mj2[0]:
    maxJ = cmds.pointOnCurve(chs[x].jName, pr=cmds.getAttr(chs[x].jName + '.maxValue'), p=True)  
  else:
    maxJ = cmds.pointOnCurve(chs[x].jName, pr=cmds.getAttr(chs[x].jName + '.minValue'), p=True) 

  return maxJ


def createRegions():
  global regions
  global ch
  global chs
  global chNum
  global numReg

  # Pre-determine num of regions
  numReg = 0
  currentReg = 0

  if chNum < 4:
    return
  elif chNum == 4:
    numReg = 1
  else:
    #numReg = np.ceil((chNum - 4)/2.0) + 1
    numReg = np.floor((chNum - 4)/2.0) + 2

  regions = [None for x in range (numReg)]
  print "Num regions: %s" % (numReg)


  # Build regions
  for x in range(chNum):
    cpairs = []
    print "[%s,%s] at (%s, %s)" % (chs[x].i, chs[x].j, chs[x].pos[0], chs[x].pos[1]) # print current cross hair
    closest = None
    closestD = 999
    ################## second neighbor search ################
    # search all other cross hairs
    firstN = x
    secN = -1
    for y in range(chNum):
      # if another cross hair has the same first term as the current's second
      # if another cross hair is on the same curve as this one
      # then it is a potential closest neighbor
      if x != y and (chs[y].i == chs[x].j or chs[y].j == chs[x].j or chs[y].j == chs[x].i \
        or chs[y].i == chs[x].i):
        print "Checking [%s,%s] with [%s,%s]" % (chs[x].i,chs[x].j,chs[y].i,chs[y].j)
        xd = chs[y].pos[0] - chs[x].pos[0]
        print xd
        yd = chs[y].pos[1] - chs[x].pos[1]
        #dist = np.sqrt(xd*xd + yd*yd)
        dist = np.sqrt(np.fabs(xd) + yd*yd)
        # store closest neighbor based on distance
        if closest == None and xd > 0:
          closest = y
          closestD = dist
          secN = y
        elif dist < closestD and xd > 0:
          closest = y
          closestD = dist
          secN = y
        print closestD  
        print secN
    if secN > -1:
      print "secN = [%s,%s]" % (chs[secN].i,chs[secN].j)
      #cpairs[0] = [chs[x].i, chs[x].j]   
      maxI = getMaxI(x)
      maxJ = getMaxJ(x)
      if maxI[0] > maxJ[0]:
        cpairs.append([chs[x].j, chs[x].i])
      else:
        cpairs.append([chs[x].i, chs[x].j])
      #cpairs[1] = [chs[secN].i, chs[secN].j]
      maxI = getMinI(secN)
      maxJ = getMinJ(secN)
      if maxI[0] > maxJ[0]:
        cpairs.append([chs[secN].j, chs[secN].i])
      else:
        cpairs.append([chs[secN].i, chs[secN].j])

      ################## third neighbor search ################
      thirdN = -1
      closest2 = None
      closestD2 = 999
      for z in range(chNum):
        # if another cross hair has the same first term as the current's second
        # then it is a potential closest neighbor
        if secN != z and x != z and (chs[z].i == chs[secN].j or chs[z].j == chs[secN].j or \
          chs[z].j == chs[secN].i or chs[z].i == chs[secN].i):
          print "Checking [%s,%s] with [%s,%s]" % (chs[secN].i,chs[secN].j,chs[z].i,chs[z].j)
          xd = chs[z].pos[0] - chs[secN].pos[0]
          yd = chs[z].pos[1] - chs[secN].pos[1]
          #dist = np.sqrt(xd*xd + yd*yd)
          dist = np.sqrt(xd*xd + np.fabs(yd))
          # store closest neighbor based on distance
          if closest2 == None and yd < 0:
            closest2 = z
            closestD2 = dist
            thirdN = z
          elif dist < closestD2 and yd < 0:
            closest2 = z
            closestD2 = dist
            thirdN = z
          print closestD2
          print thirdN
      if thirdN > -1:
        print "thirdN = [%s,%s]" % (chs[thirdN].i,chs[thirdN].j)
        maxI = getMinI(thirdN)
        maxJ = getMinJ(thirdN)
        if maxI[0] < maxJ[0]:  
          cpairs.append([chs[thirdN].j, chs[thirdN].i])
        else:
          cpairs.append([chs[thirdN].i, chs[thirdN].j])
        ################## fourth neighbor search ################
        fourthN = -1
        closest3 = None
        closestD3 = 999
        for w in range(chNum):
          # if another cross hair has the same first term as the current's second
          # then it is a potential closest neighbor
          if thirdN != w and w != x and w != z and (chs[w].i == chs[thirdN].j or \
            chs[w].j == chs[thirdN].j or chs[w].j == chs[thirdN].i or chs[w].i == chs[thirdN].i):
            xd = chs[w].pos[0] - chs[thirdN].pos[0]
            yd = chs[w].pos[1] - chs[thirdN].pos[1]
            #dist = np.sqrt(xd*xd + yd*yd)
            dist = np.sqrt(np.fabs(xd) + yd*yd)
            # store closest neighbor based on distance
            if closest3 == None and xd < 0:
              closest3 = w
              closestD3 = dist
              fourthN = w
            elif dist < closestD3 and xd < 0:
              closest3 = w
              closestD3 = dist
              fourthN = w
            print closestD3
        if fourthN > -1:
          #cpairs[3] = [chs[fourthN].i, chs[fourthN].j]
          print "fourthN = [%s,%s]" % (chs[fourthN].i,chs[fourthN].j)
          maxI = getMaxI(fourthN)
          maxJ = getMaxJ(fourthN)
          if maxI[0] < maxJ[0]:
            cpairs.append([chs[fourthN].j, chs[fourthN].i])
          else:
            cpairs.append([chs[fourthN].i, chs[fourthN].j])

    # add region if cpairs is full
    if len(cpairs) == 4 and currentReg < numReg:
      regions[currentReg] = Region(currentReg)
      regions[currentReg].corners = cpairs
      currentReg = currentReg + 1  

  # Print all regions
  for x in range(numReg):
    if regions[x]:
      print "Region %s" % (regions[x].ID)
      print "[%s, %s]" % (regions[x].corners[0][0], regions[x].corners[0][1])
      print "[%s, %s]" % (regions[x].corners[1][0], regions[x].corners[1][1])
      print "[%s, %s]" % (regions[x].corners[2][0], regions[x].corners[2][1])
      print "[%s, %s]" % (regions[x].corners[3][0], regions[x].corners[3][1])


def propagateCurve():
  global csNum
  global cs
  global chNum
  global ch

  global vertices
  global regions
  global numReg

  # TO DO: generalize
  """regions = [None for x in range (0,4)]
  regions[0] = Region(0)
  regions[0].corners = [
    [0,3],
    [3,1],
    [1,2],
    [2,0]
  ]
  regions[1] = Region(1)
  regions[1].corners = [
    [0,4],
    [4,1],
    [1,3],
    [3,0]
  ]
  regions[2] = Region(2)
  regions[2].corners = [
    [1,3],
    [3,5],
    [5,2],
    [2,1]
  ]
  regions[3] = Region(3)
  regions[3].corners = [
    [1,4],
    [4,5],
    [5,3],
    [3,1]
  ]"""

  for y in range(numReg):
    #region null check
    if regions[y]:


      # square patch dimension T_STEP by T_STEP
      T_STEPS = 10

     #given 4 corner points, CW
      """cpairs = [
        [0,3],
        [3,1],
        [1,2],
        [2,0]
      ]"""
      cpairs = regions[y].corners


      vertices = [[None]*(T_STEPS) for x in range(T_STEPS)]
      normals = [[None]*(T_STEPS) for x in range(T_STEPS)]
        
      # ALONG CURVE
        
      # go from each ch to ch
      for p in range(4):
        cStart = ch[ cpairs[p][1] ][ cpairs[p][0] ]
        if p<3:
          cEnd = ch[ cpairs[p+1][0] ][ cpairs[p+1][1] ]
        else:      
          cEnd = ch[ cpairs[0][0] ][ cpairs[0][1] ]
        tStep = (cEnd.t-cStart.t)/(T_STEPS-1)
          
        print "%s to %s : %s" % (cStart.t, cEnd.t, tStep)
      
        # go down curve
          
        t = cStart.t
        #cmds.spaceLocator( p=cmds.pointOnCurve(cs[cStart.i].name, pr=t, p=True))

        for step in range(T_STEPS-1):    
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
            
          print "(%s,%s)" % (coord[0], coord[1])
          vertices[coord[0]][coord[1]] = pos
          normals[coord[0]][coord[1]] = nor
          
          #cmds.spaceLocator( p=(pos+nor).tolist() )
          
          t = t + tStep
      
      # test 2d array
      for i in range(T_STEPS):
        line = ""
        for j in range(T_STEPS):
          if vertices[i][j] is not None:
            line = line + ". "
          else:
            line = line + "  "
        print line


      """halfway = (T_STEPS)/2.0 #5.0
      inc = 0.5
      incF = 0.5
      print halfway
      # Set z values
      for i in range(1, T_STEPS - 1):
        vertices[i][0][2] = vertices[i][0][2] + inc
        vertices[0][i][2] = vertices[0][i][2] + inc
        vertices[i][T_STEPS-1][2] = vertices[i][T_STEPS-1][2] + inc
        vertices[T_STEPS-1][i][2] = vertices[T_STEPS-1][i][2] + inc
        if i < halfway:
          incF = incF/2.0
          inc = inc + incF
        elif i > halfway:
          incF = incF*2.0
          inc = inc - incF
        print inc"""
            
        #print vertices[i][0][2]
        #print vertices[0][i][2]

        
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
      #for x in range(1, n):
        #for y in range(1, n):
          #cmds.spaceLocator(p=vertices[x][y])
      createPatchMesh(vertices, normals)
      #print "Y VALUE"
      #print y


  #createPatchMesh(vertices, normals)
  createLight()
        
#----------------------------------------------------------
# RUN COMMAND
#----------------------------------------------------------
        
def runCrossShade():
  readCrossSections()
  printCrossSectionData1()
  minOptimize()
  createRegions()
  propagateCurve()
  
runCrossShade()