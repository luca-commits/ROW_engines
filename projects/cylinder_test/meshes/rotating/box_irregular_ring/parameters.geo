cm = 1e-1;

radius_airgap = 1.3*cm;  // Radius of the conducting cylinder
radius_airbox = 1*cm;
radius_cylinder = 1*cm;

airgap_stator_width = 0.1*cm;
airgap_rotor_width = 0.13*cm; 

airgap_width = 0.05*cm;
radius_ring_internal = radius_airgap + airgap_width;
ring_width = 0.3*cm;
radius_ring_external = radius_ring_internal + ring_width;
radius_airgap_stator = radius_cylinder + airgap_stator_width;
radius_airgap_rotor = radius_cylinder + airgap_stator_width + airgap_rotor_width;


lc_airbox = 0.5*cm;
lc_ring_external = 0.2*cm;
lc_ring_internal = 0.05*cm;
lc_airgap = 0.07*cm;
lc_cylinder = 0.5*cm;

box_side = 1.52*cm;
lc_box = 0.03*cm;
lc_cube = 0.11*cm; 
