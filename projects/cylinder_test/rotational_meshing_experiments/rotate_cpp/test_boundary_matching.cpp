#include <gmsh.h>
#include <Eigen/Dense>
#include <iostream>
#include <set>
#include <vector>
#include <map>
#include <cmath>

struct Node {
    double x, y;
    std::size_t tag;
    
    bool operator==(const Node& other) const {
        return std::abs(x - other.x) < 1e-7 && std::abs(y - other.y) < 1e-7;
    }
    
    bool operator<(const Node& other) const {
        if (std::abs(x - other.x) < 1e-7) {
            return y < other.y;
        }
        return x < other.x;
    }
};

std::vector<Node> getNodesFromPhysicalGroup(const std::string& filename, int physicalTag) {
    gmsh::clear();
    gmsh::merge(filename);
    
    // Changed dimension from 2 to 1 for curves
    std::vector<int> entityTags;
    gmsh::model::getEntitiesForPhysicalGroup(1, physicalTag, entityTags);
    
    std::set<std::size_t> nodeTags;
    for (const auto& tag : entityTags) {
        std::vector<std::size_t> elementTags, nodeTagsForElement;
        // Changed type from 2 to 1 (line elements)
        gmsh::model::mesh::getElementsByType(1, elementTags, nodeTagsForElement, tag);
        nodeTags.insert(nodeTagsForElement.begin(), nodeTagsForElement.end());
    }
    
    std::vector<std::size_t> allNodeTags;
    std::vector<double> coords;
    std::vector<double> paramCoords;
    gmsh::model::mesh::getNodes(allNodeTags, coords, paramCoords);
    
    std::vector<Node> nodes;
    for (const auto& tag : nodeTags) {
        auto it = std::find(allNodeTags.begin(), allNodeTags.end(), tag);
        if (it != allNodeTags.end()) {
            size_t idx = std::distance(allNodeTags.begin(), it);
            nodes.push_back({coords[idx * 3], coords[idx * 3 + 1], tag});
        }
    }
    
    // Print debug info
    std::cout << "Found " << nodes.size() << " nodes for physical group " << physicalTag 
              << " in file " << filename << std::endl;
    
    return nodes;
}

bool testBoundaryMatching(const std::string& rotatedRotorFile, 
                         const std::string& airgapFile,
                         const std::string& statorFile,
                         double angle,
                         double tolerance = 1e-7) {
    gmsh::initialize();
    
    // Get boundary nodes from each mesh
    std::vector<Node> rotorNodes = getNodesFromPhysicalGroup(rotatedRotorFile, 34); // rotor-airgap boundary
    std::vector<Node> airgapNodes = getNodesFromPhysicalGroup(airgapFile, 34);      // airgap-rotor boundary
    std::vector<Node> statorNodes = getNodesFromPhysicalGroup(statorFile, 9);        // stator-airgap boundary
    std::vector<Node> airgapStatorNodes = getNodesFromPhysicalGroup(airgapFile, 9);  // airgap-stator boundary
    
    bool allMatched = true;
    int unmatched_count = 0;
    
    // Test rotor-airgap interface
    std::cout << "\nTesting rotor-airgap interface..." << std::endl;
    for (const auto& rotorNode : rotorNodes) {
        bool found = false;
        for (const auto& airgapNode : airgapNodes) {
            if (std::abs(rotorNode.x - airgapNode.x) < tolerance && 
                std::abs(rotorNode.y - airgapNode.y) < tolerance) {
                found = true;
                break;
            }
        }
        if (!found) {
            if (unmatched_count < 5) { // Limit the number of printed mismatches
                std::cout << "Unmatched rotor node at (" << rotorNode.x << ", " << rotorNode.y 
                         << ") tag: " << rotorNode.tag << std::endl;
            }
            unmatched_count++;
            allMatched = false;
        }
    }
    if (unmatched_count > 5) {
        std::cout << "... and " << (unmatched_count - 5) << " more unmatched nodes" << std::endl;
    }
    
    // Test stator-airgap interface
    std::cout << "\nTesting stator-airgap interface..." << std::endl;
    unmatched_count = 0;
    for (const auto& statorNode : statorNodes) {
        bool found = false;
        for (const auto& airgapNode : airgapStatorNodes) {
            if (std::abs(statorNode.x - airgapNode.x) < tolerance && 
                std::abs(statorNode.y - airgapNode.y) < tolerance) {
                found = true;
                break;
            }
        }
        if (!found) {
            if (unmatched_count < 5) {
                std::cout << "Unmatched stator node at (" << statorNode.x << ", " << statorNode.y 
                         << ") tag: " << statorNode.tag << std::endl;
            }
            unmatched_count++;
            allMatched = false;
        }
    }
    if (unmatched_count > 5) {
        std::cout << "... and " << (unmatched_count - 5) << " more unmatched nodes" << std::endl;
    }
    
    // Print statistics
    std::cout << "\nBoundary Statistics:" << std::endl;
    std::cout << "Rotor boundary nodes: " << rotorNodes.size() << std::endl;
    std::cout << "Airgap-rotor boundary nodes: " << airgapNodes.size() << std::endl;
    std::cout << "Stator boundary nodes: " << statorNodes.size() << std::endl;
    std::cout << "Airgap-stator boundary nodes: " << airgapStatorNodes.size() << std::endl;
    std::cout << "Current rotation angle: " << angle << " radians" << std::endl;
    
    gmsh::finalize();
    return allMatched;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " <number_of_meshes>" << std::endl;
        return 1;
    }
    
    double number_of_meshes;
    sscanf(argv[1], "%lf", &number_of_meshes);
    double angle_step = M_PI / 180;
    
    for (unsigned i = 0; i < number_of_meshes; ++i) {
        double rel_angle = angle_step * i;
        std::cout << "\nTesting mesh " << i << " at angle " << rel_angle << " radians" << std::endl;
        
        if (!testBoundaryMatching("rotated_rotor.msh", 
                                 "airgap.msh",
                                 "stator.msh",
                                 rel_angle)) {
            std::cout << "Boundary matching test failed for mesh " << i << std::endl;
            return 1;
        }
    }
    
    std::cout << "All boundary matching tests passed!" << std::endl;
    return 0;
}