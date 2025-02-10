// #include <igl/edges.h>
// #include <ipc/ipc.hpp>
// #include <ipc/barrier/adaptive_stiffness.hpp>

// #include "../include/FEM3D.h"

// void FEM3D::buildIPCRestData()
// {
//     ipc_vertices = surface_vertices;
//     ipc_faces = surface_indices;
//     igl::edges(ipc_faces, ipc_edges);
    
//     // if (ipc::has_intersections(ipc_vertices, ipc_edges, ipc_faces))
//     //     std::cout << "ipc has ixn in rest state! have to debug the input" << std::endl;
//     collision_mesh = ipc::CollisionMesh(ipc_vertices,ipc_edges,ipc_faces);

//     TV min_corner, max_corner;
//     computeBoundingBox(min_corner, max_corner);
//     T bb_diag = (max_corner - min_corner).norm();
//     VectorXT dedx(num_nodes * 3), dbdx(num_nodes * 3);
//     dedx.setZero(); dbdx.setZero();
//     ipc_barrier_weight = 1.0;
//     addIPCForceEntries(dbdx); dbdx *= -1.0;
//     computeResidual(dedx); dedx *= -1.0; dedx -= dbdx;
//     ipc_barrier_weight = ipc::initial_barrier_stiffness(bb_diag, ipc_barrier_distance, 1.0, dedx, dbdx, max_barrier_weight);
//     std::cout << "barrier weight " <<  ipc_barrier_weight << " max_barrier_weight " << max_barrier_weight << std::endl;
// }

// T FEM3D::addIPCEnergy()
// {
//     Eigen::MatrixXd ipc_vertices_deformed = ipc_vertices;
    
//     for (int i = 0; i < num_nodes; i++)
//     {
//         ipc_vertices_deformed.row(i).segment<3>(0) = deformed.segment<3>(i * 3);
//     }

//     ipc::CollisionConstraints ipc_constraints;
//     ipc_constraints.set_use_convergent_formulation(true);
//     ipc_constraints.build(collision_mesh, ipc_vertices_deformed, ipc_barrier_distance);
    
//     T energy_ipc = ipc_barrier_weight * ipc_constraints.compute_potential(collision_mesh, ipc_vertices_deformed, ipc_barrier_distance);

//     return energy_ipc;
// }


// void FEM3D::addIPCForceEntries(VectorXT& residual)
// {
//     Eigen::MatrixXd ipc_vertices_deformed = ipc_vertices;
    
//     for (int i = 0; i < num_nodes; i++)
//     {
//         ipc_vertices_deformed.row(i).segment<3>(0) = deformed.segment<3>(i * 3);
//     }
    
//     ipc::CollisionConstraints ipc_constraints;
//     ipc_constraints.set_use_convergent_formulation(true);
//     ipc_constraints.build(collision_mesh, ipc_vertices_deformed, ipc_barrier_distance);

//     VectorXT contact_gradient = ipc_barrier_weight * ipc_constraints.compute_potential_gradient(collision_mesh, ipc_vertices_deformed, ipc_barrier_distance);    

//     residual -= contact_gradient;
// }

// void FEM3D::addIPCHessianEntries(std::vector<Entry>& entries)
// {
//     Eigen::MatrixXd ipc_vertices_deformed = ipc_vertices;
    
//      for (int i = 0; i < num_nodes; i++)
//     {
//         ipc_vertices_deformed.row(i).segment<3>(0) = deformed.segment<3>(i * 3);
//     }
    
//     ipc::CollisionConstraints ipc_constraints;
//     ipc_constraints.set_use_convergent_formulation(true);
//     ipc_constraints.build(collision_mesh, ipc_vertices_deformed, ipc_barrier_distance);

//     StiffnessMatrix contact_hessian = ipc_barrier_weight * ipc_constraints.compute_potential_hessian(collision_mesh, ipc_vertices_deformed, ipc_barrier_distance, false);
    

//     std::vector<Entry> contact_entries = entriesFromSparseMatrix(contact_hessian);
//     entries.insert(entries.end(), contact_entries.begin(), contact_entries.end());

// }

// void FEM3D::updateIPCVertices()
// {
//     deformed = undeformed + u;

//     tbb::parallel_for(0, num_nodes, [&](int i)
//     {
//         ipc_vertices.row(i) = deformed.segment<3>(i * 3);
//     });
//     collision_mesh.build_from_full_mesh(ipc_vertices, ipc_edges, ipc_faces);
// }

// T FEM3D::computeCollisionFreeStepsize(const VectorXT& du)
// {
//     Eigen::MatrixXd current_position = ipc_vertices, 
//         next_step_position = ipc_vertices;

//     tbb::parallel_for(0, num_nodes, [&](int i)
//     {
//         current_position.row(i).segment<3>(0) = undeformed.segment<3>(i * 3) + u.segment<3>(i * 3);
//         next_step_position.row(i).segment<3>(0) = undeformed.segment<3>(i * 3) + u.segment<3>(i * 3) + du.segment<3>(i * 3);
//     });
    
//     T step_size = ipc::compute_collision_free_stepsize(collision_mesh, current_position, 
//             next_step_position, ipc::BroadPhaseMethod::HASH_GRID, 1e-6, 1e7);
//     return step_size;
// }

// void FEM3D::updateBarrierInfo(bool first_step)
// {
    
//     Eigen::MatrixXd ipc_vertices_deformed = ipc_vertices;
    
//     for (int i = 0; i < num_nodes; i++) 
//         ipc_vertices_deformed.row(i) = deformed.segment<3>(i * 3);

//     ipc::CollisionConstraints ipc_constraints;
//     ipc_constraints.set_use_convergent_formulation(true);
//     ipc_constraints.build(collision_mesh, ipc_vertices_deformed, ipc_barrier_distance);
        
//     T current_min_dis = ipc_constraints.compute_minimum_distance(collision_mesh, ipc_vertices_deformed);
//     if (first_step)
//         ipc_min_dis = current_min_dis;
//     else
//     {
//         TV min_corner, max_corner;
//         computeBoundingBox(min_corner, max_corner);
//         T bb_diag = (max_corner - min_corner).norm();
//         ipc::update_barrier_stiffness(ipc_min_dis, current_min_dis, max_barrier_weight, ipc_barrier_weight, bb_diag);
//         ipc_min_dis = current_min_dis;
//         // std::cout << "barrier weight " << barrier_weight << std::endl;
//     }
// }