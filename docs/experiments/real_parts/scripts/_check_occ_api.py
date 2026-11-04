from OCC.Core.Geom import Geom_BezierSurface, Geom_BSplineSurface
print("Bezier setters:", [m for m in dir(Geom_BezierSurface) if m.startswith("Set")])
print("Bezier pole/weight:", [m for m in dir(Geom_BezierSurface) if "Pole" in m or "Weight" in m or "Rational" in m])
print("BSpline ctor available:", Geom_BSplineSurface is not None)
from OCC.Core.GeomAPI import GeomAPI_IntSS
print("IntSS:", [m for m in dir(GeomAPI_IntSS) if not m.startswith("_")])
