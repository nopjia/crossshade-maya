from sympy import *
import maya.cmds as cmds
import maya.mel as mel

curves = mel.eval("ls -sl")
intersects = [float(i) for i in cmds.curveIntersect(curves[0], curves[1], useDirection=True, direction=(0,0,1)).split()]
pi = cmds.pointOnCurve(curves[0], pr=intersects[0], p=True)
pj = cmds.pointOnCurve(curves[1], pr=intersects[1], p=True)
ni = cmds.pointOnCurve(curves[0], pr=intersects[0], nn=True)
nj = cmds.pointOnCurve(curves[1], pr=intersects[1], nn=True)
ti = cmds.pointOnCurve(curves[0], pr=intersects[0], t=True)
tj = cmds.pointOnCurve(curves[1], pr=intersects[1], t=True)

cmds.spaceLocator(p=pi)

niz = Symbol("niz")
njz = Symbol("njz")
tiz = Symbol("tiz")
tjz = Symbol("tjz")
solution = solve(
  [
    ni[0]*nj[0] + ni[1]*nj[1] + niz*njz,
    ti[0]*tj[0] + ti[1]*tj[1] + tiz*tjz,
    ti[0]*ni[0] + ti[1]*ni[1] + tiz*niz,
    tj[0]*nj[0] + tj[1]*nj[1] + tjz*njz
  ], 
  niz, njz, tiz, tjz)