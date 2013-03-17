import sys
import maya.OpenMaya as OpenMaya
import maya.OpenMayaMPx as OpenMayaMPx
import maya.cmds as cmds
import maya.mel as mel

kPluginCmdName = "CrossShadeCmd"


def processCrossSections():

  mel.eval('layerEditorSelectObjects layerCross;')
  crossSections = mel.eval('ls -sl -type transform')

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
        
        # print '(%s,%s)' % (ci, cj)
        
      else:
        print 'no curve intersect'

# command
class scriptedCommand(OpenMayaMPx.MPxCommand):
  def __init__(self):
    OpenMayaMPx.MPxCommand.__init__(self)
    
  def doIt(self, argList):
    processCrossSections()

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
