radius_airgap = 1.3;  // Radius of the conducting cylinder
radius_airbox = 3;
radius_cylinder = 1;

airgap_stator_width = 0.1;
airgap_rotor_width = 0.1; 

airgap_width = 0.2;
radius_ring_internal = radius_airgap + airgap_width;
ring_width = 0.3;
radius_ring_external = radius_ring_internal + ring_width;
radius_airgap_stator = radius_cylinder + airgap_stator_width;
radius_airgap_rotor = radius_cylinder + airgap_stator_width + airgap_rotor_width;


lc_airbox = 0.5;
lc_ring_external = 0.01;
lc_ring_internal = 0.01;
lc_airgap = 0.005;
lc_cylinder = 0.1;

