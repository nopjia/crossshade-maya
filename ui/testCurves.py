import maya.cmds as cmds
import maya.mel as mel

curves = mel.eval("ls -sl")
intersects = [float(i) for i in cmds.curveIntersect(curves[0], curves[1], useDirection=True, direction=(0,0,1)).split()]
p1 = cmds.pointOnCurve(curves[0], pr=intersects[0], p=True)
p2 = cmds.pointOnCurve(curves[1], pr=intersects[1], p=True)
n1 = cmds.pointOnCurve(curves[0], pr=intersects[0], nn=True)
n2 = cmds.pointOnCurve(curves[1], pr=intersects[1], nn=True)
