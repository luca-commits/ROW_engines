import gmsh

gmsh.initialize()
gmsh.merge("rotor.msh")

# Extract node tags and coordinates from the specified physical group
nodeTags, nodeCoords = gmsh.model.mesh.getNodesForPhysicalGroup(1, 1)

