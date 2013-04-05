import maya.cmds as cmds
import maya.mel as mel

numCrossSections = 0
crossSections = None
crossHairPos  = None
crossHairTan  = None

def processCrossSections():
  global numCrossSections
  global crossSections 
  global crossHairPos  
  global crossHairTan  

  mel.eval('layerEditorSelectObjects layerCross;')
  crossSections = mel.eval('ls -sl -type transform')

  # no curve drawn
  if not crossSections:
    return

  numCrossSections = len(crossSections)
  crossHairPos = [[0]*numCrossSections for x in range(numCrossSections, 0, -1)]  
  crossHairTan = [[0]*numCrossSections for x in range(numCrossSections, 0, -1)]

  for ci in range(numCrossSections):
    for cj in range(ci+1, numCrossSections):

      rawIntersects = cmds.curveIntersect(crossSections[ci], crossSections[cj], useDirection=True, direction=(0,0,1))
      
      if rawIntersects:
        intersects = [float(i) for i in rawIntersects.split()]
        intT_i = intersects[0]
        intT_j = intersects[1]
        
        pi = cmds.pointOnCurve(crossSections[ci], pr=intT_i, p=True)        
        crossHairPos[ci][cj] = pi;  # save
        cmds.spaceLocator(p=pi)   # place locator        
        
        t_ij = cmds.pointOnCurve(crossSections[ci], pr=intT_i, nt=True)
        crossHairTan[ci][cj] = t_ij;
        cmds.spaceLocator(p=[ p+t for p,t in zip(pi,t_ij) ])
        t_ji = cmds.pointOnCurve(crossSections[cj], pr=intT_j, nt=True)
        crossHairTan[cj][ci] = t_ji;
        cmds.spaceLocator(p=[ p+t for p,t in zip(pi,t_ji) ])
        
        print 'process (%s,%s)' % (ci, cj)
        
      else:
        print 'no curve intersect'
        
        
def printCrossSectionData():
  global numCrossSections
  global crossSections
  global crossHairPos
  global crossHairTan
  
  for i in range(numCrossSections):
    print "%s : %s" % (i, crossSections[i])
  
  for i in range(numCrossSections):
    for j in range(numCrossSections):
      print "x_(%s,%s) : %s" % (i, j, crossHairPos[i][j])
      
  for i in range(numCrossSections):
    for j in range(numCrossSections):
      print "t_(%s,%s) : %s" % (i, j, crossHairTan[i][j])