#include "../include/Scene.h"

static double ROD_A = 2.5e-4;
static double ROD_B = 2.5e-4;

void Scene::buildInterlockingSquareScene(int sub_div)
{
    sim.yarn_map.clear();
        
    clearSimData();

    sim.add_rotation_penalty = false;
    sim.add_pbc_bending = false;
    sim.add_pbc_twisting = false;
    sim.add_pbc = false;

    sim.add_contact_penalty=true;
    sim.new_frame_work = true;
    sim.add_eularian_reg = true;

    sim.ke = 1e-6;

    sim.unit = 0.09;

    auto addCrossingData = [&](int crossing_idx, int rod_idx, int location)
    {
        sim.rod_crossings[crossing_idx]->is_fixed = true;
        sim.rod_crossings[crossing_idx]->rods_involved.push_back(rod_idx);
        sim.rod_crossings[crossing_idx]->on_rod_idx[rod_idx] = location;
        sim.rod_crossings[crossing_idx]->sliding_ranges.push_back(Range::Zero());
    };

    std::vector<Eigen::Triplet<T>> w_entry;
    int full_dof_cnt = 0;
    int node_cnt = 0;
    int rod_cnt = 0;

    std::vector<TV> nodal_positions;

    auto addSquare = [&](const TV& bottom_left, T width)
    {
        
        sim.rod_crossings.push_back(new RodCrossing(node_cnt, std::vector<int>()));
        addPoint(bottom_left, full_dof_cnt, node_cnt);
        nodal_positions.push_back(bottom_left);

        sim.rod_crossings.push_back(new RodCrossing(node_cnt, std::vector<int>()));
        TV bottom_right = bottom_left + TV(width, 0, 0);
        addPoint(bottom_right, full_dof_cnt, node_cnt);
        nodal_positions.push_back(bottom_right);

        sim.rod_crossings.push_back(new RodCrossing(node_cnt, std::vector<int>()));
        TV top_right = bottom_left + TV(width, width, 0);
        addPoint(top_right, full_dof_cnt, node_cnt);
        nodal_positions.push_back(top_right);

        sim.rod_crossings.push_back(new RodCrossing(node_cnt, std::vector<int>()));
        TV top_left = bottom_left + TV(0, width, 0);
        addPoint(top_left, full_dof_cnt, node_cnt);
        nodal_positions.push_back(top_left);
    };

    TV s0 = TV(0, 0, 0) * sim.unit;
    addSquare(s0, 1.0 * sim.unit);

    TV s1 = TV(0.75, 0.75, 0) * sim.unit;
    addSquare(s1, 1.0 * sim.unit);

    TV crossing0 = TV(1.0, 0.75, 0.0) * sim.unit;
    sim.rod_crossings.push_back(new RodCrossing(node_cnt, std::vector<int>()));
    addPoint(crossing0, full_dof_cnt, node_cnt);
    nodal_positions.push_back(crossing0);

    sim.rod_crossings.push_back(new RodCrossing(node_cnt, std::vector<int>()));
    TV crossing1 = TV(0.75, 1.0, 0.0) * sim.unit;
    addPoint(crossing1, full_dof_cnt, node_cnt);
    nodal_positions.push_back(crossing1);

    
    std::vector<TV> passing_points = {nodal_positions[0], nodal_positions[1]};
    std::vector<int> passing_points_id = {0, 1};

    addAStraightRod(passing_points.front(), passing_points.back(), 
        passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);
        
    for (int i = 0; i < passing_points_id.size(); i++) 
        addCrossingData(passing_points_id[i], rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[i]);

    
    passing_points = {nodal_positions[1], nodal_positions[8], nodal_positions[2]};
    passing_points_id = {1, 8, 2};
    addAStraightRod(passing_points.front(), passing_points.back(), passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);
    for (int i = 0; i < passing_points_id.size(); i++) 
        addCrossingData(passing_points_id[i], rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[i]);

    passing_points = {nodal_positions[2], nodal_positions[9], nodal_positions[3]};
    passing_points_id = {2, 9, 3};
    addAStraightRod(passing_points.front(), passing_points.back(), passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);
    for (int i = 0; i < passing_points_id.size(); i++) 
        addCrossingData(passing_points_id[i], rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[i]);

    passing_points = {nodal_positions[3], nodal_positions[0]};
    passing_points_id = {3, 0};
    addAStraightRod(passing_points.front(), passing_points.back(), passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);
    for (int i = 0; i < passing_points_id.size(); i++) 
        addCrossingData(passing_points_id[i], rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[i]);

    passing_points = {nodal_positions[4], nodal_positions[8], nodal_positions[5]};
    passing_points_id = {4, 8, 5};
    addAStraightRod(passing_points.front(), passing_points.back(), passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);
    for (int i = 0; i < passing_points_id.size(); i++) 
        addCrossingData(passing_points_id[i], rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[i]);

    passing_points = {nodal_positions[5], nodal_positions[6]};
    passing_points_id = {5, 6};
    addAStraightRod(passing_points.front(), passing_points.back(), passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);
    for (int i = 0; i < passing_points_id.size(); i++) 
        addCrossingData(passing_points_id[i], rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[i]);

    passing_points = {nodal_positions[6], nodal_positions[7]};
    passing_points_id = {6, 7};
    addAStraightRod(passing_points.front(), passing_points.back(), passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);
    for (int i = 0; i < passing_points_id.size(); i++) 
        addCrossingData(passing_points_id[i], rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[i]);

    passing_points = {nodal_positions[7], nodal_positions[9], nodal_positions[4]};
    passing_points_id = {7, 9, 4};
    addAStraightRod(passing_points.front(), passing_points.back(), passing_points, passing_points_id, sub_div, full_dof_cnt, node_cnt, rod_cnt);
    for (int i = 0; i < passing_points_id.size(); i++) 
        addCrossingData(passing_points_id[i], rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[i]);

    // for (auto& crossing : sim.rod_crossings)
    // {
    //     std::sort( crossing->rods_involved.begin(), crossing->rods_involved.end() );
    //     crossing->rods_involved.erase( std::unique( crossing->rods_involved.begin(), crossing->rods_involved.end() ), crossing->rods_involved.end() );
    // }
    
    for (auto& rod : sim.Rods)
        rod->fixed_by_crossing = std::vector<bool>(rod->dof_node_location.size(), true);

    int dof_cnt = 0;
    markCrossingDoF(w_entry, dof_cnt);

    for (auto& rod : sim.Rods) rod->markDoF(w_entry, dof_cnt);
    
    appendThetaAndJointDoF(w_entry, full_dof_cnt, dof_cnt);
    
    sim.rest_states = deformed_states;
    
    sim.W = StiffnessMatrix(full_dof_cnt, dof_cnt);
    sim.W.setFromTriplets(w_entry.begin(), w_entry.end());
    
    for (auto& rod : sim.Rods)
    {
        rod->fixEndPointEulerian(sim.dirichlet_dof);
        rod->setupBishopFrame();
    }
    
    TV center0;
    T r = 0.05 * sim.unit;
    sim.Rods[0]->x(sim.Rods[0]->indices.front(), center0);
    auto circle0 = [r, center0](const TV& x)->bool
    {
        return (x - center0).norm() < r;
    };

    sim.fixRegion(circle0);

    Offset offset;
    sim.Rods[5]->frontOffsetReduced(offset);
    sim.dirichlet_dof[offset[0]] = 0.3 * sim.unit;
    sim.dirichlet_dof[offset[1]] = 0;
    sim.dirichlet_dof[offset[2]] = 0;

    sim.Rods[5]->backOffsetReduced(offset);
    sim.dirichlet_dof[offset[0]] = 0.3 * sim.unit;
    sim.dirichlet_dof[offset[1]] = 0;
    sim.dirichlet_dof[offset[2]] = 0;

    auto crossing = sim.rod_crossings[8];
    crossing->is_fixed = false;
    crossing->sliding_ranges[1] = Range::Ones();

    crossing = sim.rod_crossings[9];
    crossing->is_fixed = false;
    crossing->sliding_ranges[0] = Range::Ones();

    sim.fixCrossing();

    sim.perturb = VectorXT::Zero(sim.W.cols());
    for (auto& crossing : sim.rod_crossings)
    // for (int i = 0; i < 10; i++)
    {
        // auto crossing = rod_crossings[i];
        Offset off;
        sim.Rods[crossing->rods_involved.front()]->getEntry(crossing->node_idx, off);
        T r = static_cast <T> (rand()) / static_cast <T> (RAND_MAX);
        int z_off = sim.Rods[crossing->rods_involved.front()]->reduced_map[off[3-1]];
        sim.perturb[z_off] += 0.001 * (r - 0.5) * sim.unit;
        // break;
        // sim.perturb[z_off] += 0.001 * r * sim.unit;
        // sim.perturb[z_off] += 0.001 * sim.unit;
        
    }
    sim.boundary_spheres.push_back(std::make_pair(center0, r));   
}

void Scene::clearSimData()
{
    sim.kc = 1e8;
    sim.add_pbc = true;

    if(sim.disable_sliding)
    {
        sim.add_shearing = true;
        sim.add_eularian_reg = false;
        sim.k_pbc = 1e8;
        sim.k_strain = 1e8;
    }
    else
    {
        sim.add_shearing = false;
        sim.add_eularian_reg = true;
        sim.ke = 1e-4;    
        sim.k_yc = 1e8;
    }
    sim.k_pbc = 1e4;
    sim.k_strain = 1e7;
    sim.kr = 1e3;
    sim.yarns.clear();
}

void Scene::markCrossingDoF(std::vector<Eigen::Triplet<T>>& w_entry,
        int& dof_cnt)
{
    for (auto& crossing : sim.rod_crossings)
    {
        int node_idx = crossing->node_idx;
        // std::cout << "node " << node_idx << std::endl;
        std::vector<int> rods_involved = crossing->rods_involved;

        Offset entry_rod0; 
        sim.Rods[rods_involved.front()]->getEntry(node_idx, entry_rod0);

        // push Lagrangian dof first
        for (int d = 0; d < 3; d++)
        {
            
            for (int rod_idx : rods_involved)
            {
                // std::cout << "rods involved " << rod_idx << std::endl;
                // if (node_idx == 21)
                //     std::cout << "rods involved " << rod_idx << std::endl;
                sim.Rods[rod_idx]->reduced_map[entry_rod0[d]] = dof_cnt;
            }    
            w_entry.push_back(Entry(entry_rod0[d], dof_cnt++, 1.0));
        }
        
        // push Eulerian dof for all rods
        for (int rod_idx : rods_involved)
        {
            // std::cout << "3 on rod " <<  rod_idx << std::endl;
            sim.Rods[rod_idx]->getEntry(node_idx, entry_rod0);
            // std::cout << "3 dof on rod " <<  entry_rod0[3] << std::endl;
            sim.Rods[rod_idx]->reduced_map[entry_rod0[3]] = dof_cnt;
            w_entry.push_back(Entry(entry_rod0[3], dof_cnt++, 1.0));
        }
        // std::getchar();
        
    }
}

void Scene::appendThetaAndJointDoF(std::vector<Entry>& w_entry, 
    int& full_dof_cnt, int& dof_cnt)
{
    // for (auto& rod : sim.Rods)
    // {
    //     rod->theta_dof_start_offset = full_dof_cnt;
    //     rod->theta_reduced_dof_start_offset = dof_cnt;
    //     deformed_states.conservativeResize(full_dof_cnt + rod->indices.size() - 1);
    //     for (int i = 0; i < rod->indices.size() - 1; i++)
    //         w_entry.push_back(Entry(full_dof_cnt++, dof_cnt++, 1.0));
    //     deformed_states.template segment(rod->theta_dof_start_offset, 
    //         rod->indices.size() - 1).setZero();
    // }   

    deformed_states.conservativeResize(full_dof_cnt + sim.rod_crossings.size() * 3);
    deformed_states.template segment(full_dof_cnt, sim.rod_crossings.size() * 3).setZero();

    for (auto& crossing : sim.rod_crossings)
    {
        crossing->dof_offset = full_dof_cnt;
        crossing->reduced_dof_offset = dof_cnt;
        for (int d = 0; d < 3; d++)
        {
            w_entry.push_back(Entry(full_dof_cnt++, dof_cnt++, 1.0));
        }
    }

    for (auto& rod : sim.Rods)
    {
        rod->theta_dof_start_offset = full_dof_cnt;
        rod->theta_reduced_dof_start_offset = dof_cnt;
        deformed_states.conservativeResize(full_dof_cnt + rod->indices.size() - 1);
        for (int i = 0; i < rod->indices.size() - 1; i++)
            w_entry.push_back(Entry(full_dof_cnt++, dof_cnt++, 1.0));
        deformed_states.template segment(rod->theta_dof_start_offset, 
            rod->indices.size() - 1).setZero();
    }
}

void Scene::addAStraightRod(const TV& from, const TV& to, 
        const std::vector<TV>& passing_points, 
        const std::vector<int>& passing_points_id, 
        int sub_div, int& full_dof_cnt, int& node_cnt, int& rod_cnt)
{
    
    std::unordered_map<int, Offset> offset_map;

    std::vector<TV> points_on_curve;
    std::vector<int> rod_indices;
    std::vector<int> key_points_location_rod;
    addStraightYarnCrossNPoints(from, to, passing_points, passing_points_id,
                                sub_div, points_on_curve, rod_indices,
                                key_points_location_rod, node_cnt);
                   

    deformed_states.conservativeResize(full_dof_cnt + (points_on_curve.size()) * (3 + 1));

    Rod* rod = new Rod(deformed_states, sim.rest_states, rod_cnt, false, ROD_A, ROD_B);

    for (int i = 0; i < points_on_curve.size(); i++)
    {
        offset_map[node_cnt] = Offset::Zero();
        //push Lagrangian DoF    
        deformed_states.template segment<3>(full_dof_cnt) = points_on_curve[i];
        for (int d = 0; d < 3; d++)
        {
            offset_map[node_cnt][d] = full_dof_cnt++;  
        }
        // push Eulerian DoF
        deformed_states[full_dof_cnt] = (points_on_curve[i] - from).norm() / (to - from).norm();
        offset_map[node_cnt][3] = full_dof_cnt++;
        node_cnt++;
    }
    
    deformed_states.conservativeResize(full_dof_cnt + passing_points.size());

    for (int i = 0; i < passing_points.size(); i++)
    {
        deformed_states[full_dof_cnt] = (passing_points[i] - from).norm() / (to - from).norm();
        offset_map[passing_points_id[i]] = Offset::Zero();
        offset_map[passing_points_id[i]][3] = full_dof_cnt++; 
        Vector<int, 3> offset_dof_lag;
        for (int d = 0; d < 3; d++)
        {
            offset_dof_lag[d] = passing_points_id[i] * 3 + d;
        }
        offset_map[passing_points_id[i]].template segment<3>(0) = offset_dof_lag;
    }
    
    rod->offset_map = offset_map;
    rod->indices = rod_indices;
    Vector<T, 3 + 1> q0, q1;
    rod->frontDoF(q0); rod->backDoF(q1);

    rod->rest_state = new LineCurvature(q0, q1);
    
    rod->dof_node_location = key_points_location_rod;
    
    sim.Rods.push_back(rod);
    rod_cnt++;
}

void Scene::addStraightYarnCrossNPoints(const TV& from, const TV& to,
    const std::vector<TV>& passing_points, 
    const std::vector<int>& passing_points_id, int sub_div,
    std::vector<TV>& sub_points, std::vector<int>& node_idx, 
    std::vector<int>& key_points_location, 
    int start, bool pbc)
{
    
    int cnt = 1;
    if(passing_points.size())
    {
        if ((from - passing_points[0]).norm() < 1e-6 )
        {
            node_idx.push_back(passing_points_id[0]);
            cnt = 0;
        }
        else
        {
            node_idx.push_back(start);
            sub_points.push_back(from);
        }
    }
    else
    {
        node_idx.push_back(start);
        sub_points.push_back(from);
    }
    
    T length_yarn = (to - from).norm();
    TV length_vec = (to - from).normalized();
    
    TV loop_point = from;
    TV loop_left = from;
    for (int i = 0; i < passing_points.size(); i++)
    {
        if ((from - passing_points[i]).norm() < 1e-6 )
        {
            key_points_location.push_back(0);
            continue;
        }
        T fraction = (passing_points[i] - loop_point).norm() / length_yarn;
        int n_sub_nodes = std::ceil(fraction * sub_div);
        T length_sub = (passing_points[i] - loop_point).norm() / T(n_sub_nodes);
        for (int j = 0; j < n_sub_nodes - 1; j++)
        {
            sub_points.push_back(loop_left + length_sub * length_vec);
            loop_left = sub_points.back();
            node_idx.push_back(start + cnt);
            cnt++;
        }
        node_idx.push_back(passing_points_id[i]);
        key_points_location.push_back(cnt + i);
        loop_point = passing_points[i];
        loop_left = passing_points[i];
    }
    if (passing_points.size())
    {
        if ((passing_points.back() - to).norm() < 1e-6)
        {
            
            return;
        }
    }
    T fraction;
    int n_sub_nodes;
    T length_sub;
    if( passing_points.size() )
    {
        fraction = (to - passing_points.back()).norm() / length_yarn;
        n_sub_nodes = std::ceil(fraction * sub_div);
        length_sub = (to - passing_points.back()).norm() / T(n_sub_nodes);
    }
    else
    {
        n_sub_nodes = sub_div + 1;
        length_sub = (to - from).norm() / T(sub_div);
    }
    for (int j = 0; j < n_sub_nodes - 1; j++)
    {
        if (j == 0)
        {
            if(passing_points.size())
            {
                sub_points.push_back(passing_points.back() + length_sub * length_vec);
                loop_left = sub_points.back();
            }
        }
        else
        {
            sub_points.push_back(loop_left + length_sub * length_vec);
            loop_left = sub_points.back();
        }
        if(passing_points.size() == 0 && j == 0)
            continue;
        node_idx.push_back(start + cnt);
        cnt++;
    }
    node_idx.push_back(start + cnt);
    sub_points.push_back(to);
}

void Scene::addPoint(const TV& point, int& full_dof_cnt, int& node_cnt)
{
    deformed_states.conservativeResize(full_dof_cnt + 3);
    deformed_states.template segment<3>(full_dof_cnt) = point;
    full_dof_cnt += 3;
    node_cnt++;
}