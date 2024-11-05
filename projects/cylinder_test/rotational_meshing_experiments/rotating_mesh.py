import gmsh

gmsh.initialize()

gmsh.model.occ.synchronize()

radius_cylinder = 1
airgap_width = 0.1
radius_airgap = radius_cylinder + airgap_width
lc_airgap = 0.1
print("radius airgap", radius_airgap)

gmsh.model.occ.addPoint(0, 0, 0, 1, 51) 

gmsh.model.occ.addPoint(0, radius_airgap, 0, lc_airgap , 43)  
gmsh.model.occ.addPoint(radius_airgap, 0, 0, lc_airgap , 44)  
gmsh.model.occ.addPoint(0, -radius_airgap, 0, lc_airgap, 45) 
gmsh.model.occ.addPoint(-radius_airgap, 0, 0, lc_airgap, 46) 

gmsh.model.occ.addCircleArc(43, 51, 44, 53)

gmsh.model.occ.addCircleArc(44, 51, 45, 54)
gmsh.model.occ.addCircleArc(45, 51, 46, 55)
gmsh.model.occ.addCircleArc(46, 51, 43, 56)

gmsh.model.occ.addCurveLoop([53, 54, 55, 56], 17)

radius_airgap_rotor = 1.3

gmsh.model.occ.addPoint(0, radius_airgap_rotor, 0, lc_airgap , 47)  
gmsh.model.occ.addPoint(radius_airgap_rotor, 0, 0, lc_airgap , 48)  
gmsh.model.occ.addPoint(0, -radius_airgap_rotor, 0, lc_airgap, 49) 
gmsh.model.occ.addPoint(-radius_airgap_rotor, 0, 0, lc_airgap, 50) 

gmsh.model.occ.addCircleArc(47, 51, 48, 57)
gmsh.model.occ.addCircleArc(48, 51, 49, 58)
gmsh.model.occ.addCircleArc(49, 51, 50, 59)
gmsh.model.occ.addCircleArc(50, 51, 47, 60)

gmsh.model.occ.addCurveLoop([57, 58, 59, 60], 19)

gmsh.model.occ.addPlaneSurface([17, 19], 10)

gmsh.model.occ.synchronize()

gmsh.model.mesh.generate(2)

gmsh.merge("stator.msh")
gmsh.merge("rotor.msh")


gmsh.fltk.run()










