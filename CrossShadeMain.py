import sys
import maya.OpenMaya as OpenMaya
import maya.OpenMayaMPx as OpenMayaMPx
import maya.cmds as cmds
import maya.mel as mel

kPluginCmdName = "CrossShadeCmd"

numCrossSections = 0
crossSections = None
crossHairPos  = None
crossHairTan  = None

def readCrossSections():
  global numCrossSections
  global crossSections 
  global crossHairPos  
  global crossHairTan  

  mel.eval("layerEditorSelectObjects layerCross;")
  crossSections = mel.eval("ls -sl -type transform")

  # no curve drawn
  if not crossSections:
    return

  numCrossSections = len(crossSections)
  crossHairPos = [[None]*numCrossSections for x in range(numCrossSections, 0, -1)]  
  crossHairTan = [[None]*numCrossSections for x in range(numCrossSections, 0, -1)]

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
        
        print "(%s,%s) processed" % (ci, cj)
        
      else:
        print "(%s,%s) no intersect" % (ci, cj)

  # clear all selection
  mel.eval("select -cl")
        
def printCrossSectionData1():
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

      
# command
class scriptedCommand(OpenMayaMPx.MPxCommand):
  def __init__(self):
    OpenMayaMPx.MPxCommand.__init__(self)
    
  def doIt(self, argList):
    readCrossSections()
    printCrossSectionData1()

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
