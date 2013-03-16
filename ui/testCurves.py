import maya.cmds as cmds
import maya.mel as mel

crossHairPos
crossSections

def processCrossSections():

  mel.eval("layerEditorSelectObjects layerCross;")
  crossSections = mel.eval('ls -sl -type "transform"')

  # no curve drawn
  if not crossSections:
    return

  numCrossSections = len(crossSections)
  crossHairPos = [[0]*numCrossSections for x in range(numCrossSections, 0, -1)]

  for ci in range(numCrossSections):
    for cj in range(ci+1, numCrossSections):

      rawIntersects = cmds.curveIntersect(crossSections[ci], crossSections[cj], useDirection=True, direction=(0,0,1))
      
      if rawIntersects:
        intersects = [float(i) for i in rawIntersects.split()]
        
        pi = cmds.pointOnCurve(crossSections[ci], pr=intersects[0], p=True)

        # place locator
        cmds.spaceLocator(p=pi)
        
        # save
        crossHairPos[ci][cj] = pi;
        
        print "(%s,%s)" % (ci, cj)
        
      else:
        print "no curve intersect"