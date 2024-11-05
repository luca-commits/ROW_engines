#include "rotate_mesh.hpp"
#include "remesh_airgap.hpp"
#include <cmath> 
#include <string> 
#include <iostream>

int main(int argc, char *argv[]){
  std::string current_exec_name = argv[0]; // Name of the current exec program
  std::vector<std::string> all_args;
  if (argc > 1) {
      // Convert each argument to std::string and add it to all_args
      for (int i = 1; i < argc; ++i) {
          all_args.push_back(std::string(argv[i]));
      }
  }
  double number_of_meshes;
  sscanf(argv[1],"%lf",&number_of_meshes);
  std::string rotated_airgap_rotor_boundary_file = "airgap_rotor_boundary_rotated.msh";
  std::string rotor_file = "rotor.msh";
  std::string rotated_rotor_file = "rotated_rotor.msh";
  double angle_step = M_PI / 180 ; 

//   remeshAirgap("stator.msh", "rotor.msh", "airgap.msh",  0);
  std::cout << "C ------" << std::flush;
//   add_offset_to_tags("airgap.msh", 10000);
//   rewrite_tags("airgap.msh", "rotor.msh");
// //   std::cout << "D -------" << std::endl;
//   rewrite_tags("airgap.msh", "stator.msh");
// //   std::cout << "E -------" << std::endl;

  for (unsigned i = 0; i < number_of_meshes; ++i){
    double rel_angle = angle_step * i;
    remeshAirgap("stator.msh", "rotor.msh", "airgap.msh",  M_PI * rel_angle);
    rotateAllNodes_alt(rotor_file, rotated_rotor_file, M_PI * rel_angle);
    mergeEverything("stator.msh", rotated_rotor_file, "airgap.msh", "motor_" + std::to_string(0) + ".msh");
  }

//   std::map<int, int> target_to_source {{10000, 100}};
//   change_tags_in_file("rotor.msh", target_to_source); 

  return 0;
}