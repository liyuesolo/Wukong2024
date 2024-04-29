#include "../include/Scene.h"

static double ROD_A = 2.5e-4;
static double ROD_B = 2.5e-4;

void Scene::buildGridScene(int sub_div)
{
    auto unit_yarn_map = sim.yarn_map;
    sim.yarn_map.clear();
    
    sim.add_rotation_penalty = false;
    sim.add_pbc_bending = false;
    sim.add_pbc_twisting = false;
    sim.add_pbc = false;

    sim.add_contact_penalty=true;
    sim.new_frame_work = true;
    sim.add_eularian_reg = false;

    sim.ke = 1e-4;

    sim.unit = 0.09;
    
    std::vector<Eigen::Triplet<T>> w_entry;
    int full_dof_cnt = 0;
    int node_cnt = 0;

    int n_row = 20, n_col = 20;

    // push crossings first 
    T dy = 1.0 / n_row * sim.unit;
    T dx = 1.0 / n_col * sim.unit;
    
    //num of crossing
    deformed_states.resize(n_col * n_row * 3);
    
    std::unordered_map<int, Offset> crossing_offset_copy;

    auto getXY = [=](int row, int col, T& x, T& y)
    {
        if (row == 0) y = 0.5 * dy;
        else if (row == n_row) y = n_row * dy;
        else y = 0.5 * dy + (row ) * dy;
        if (col == 0) x = 0.5 * dx;
        else if (col == n_col) x = n_col * dx;
        else x = 0.5 * dx + (col ) * dx;
    };


    for (int row = 0; row < n_row; row++)
    {
        for (int col = 0; col < n_col; col++)
        {
            T x, y;
            getXY(row, col, x, y);
            deformed_states.template segment<3>(node_cnt * 3) = TV(x, y, 0);
            
            full_dof_cnt += 3;
            node_cnt ++;       
        }
    }

    int rod_cnt = 0;
    for (int row = 0; row < n_row; row++)
    {
        T x0 = 0.0, x1 = 1.0 * sim.unit;
        T x, y;
        
        std::vector<int> passing_points_id;
        std::vector<TV> passing_points;
        
        for (int col = 0; col < n_col; col++)
        {
            int node_idx = row * n_col + col;
            passing_points_id.push_back(node_idx);
            passing_points.push_back(deformed_states.template segment<3>(node_idx * 3));
        }

        getXY(row, 0, x, y);

        TV from = TV(x0, y, 0);
        TV to = TV(x1, y, 0);
    
        addAStraightRod(from, to, passing_points, passing_points_id, 
            sub_div, full_dof_cnt, node_cnt, rod_cnt);
        
    }
    
    for (int col = 0; col < n_col; col++)
    {
        T y0 = 0.0, y1 = 1.0 * sim.unit;
        T x, y;
        std::vector<int> passing_points_id;
        std::vector<TV> passing_points;
        getXY(0, col, x, y);
        for (int row = 0; row < n_row; row++)
        {
            int node_idx = row * n_col + col;
            passing_points_id.push_back(node_idx);
            passing_points.push_back(deformed_states.template segment<3>(node_idx * 3));
        }
        
        TV from = TV(x, y0, 0);
        TV to = TV(x, y1, 0);

        addAStraightRod(from, to, passing_points, passing_points_id, sub_div, 
                        full_dof_cnt, node_cnt, rod_cnt);
        
    }
    
    for (auto& rod : sim.Rods)
        rod->fixed_by_crossing.resize(rod->dof_node_location.size(), false);

    T dv = 1.0 / n_row;
    T du = 1.0 / n_col;

    int odd_even_cnt = 0;
    for (int row = 0; row < n_row; row++)
    {
        for (int col = 0; col < n_col; col++)
        {
            int node_idx = row * n_col + col;
            RodCrossing* crossing = 
                new RodCrossing(node_idx, {row, n_row + col});

            crossing->on_rod_idx[row] = sim.Rods[row]->dof_node_location[col];
            crossing->on_rod_idx[n_row + col] = sim.Rods[n_row + col]->dof_node_location[row];
            
            // if (odd_even_cnt % 2 == 0)
            //     crossing->is_fixed = true;

            // if (row % 2 == 0)
            //     crossing->is_fixed = true;

            // if (row ==  col)
            //     crossing->is_fixed = true;

            // if (row != col)
            //     crossing->is_fixed = true;

            // crossing->is_fixed = true;

            // if (row == 0 || row == n_row - 1 || col == 0 || col == n_col - 1)                    
                crossing->is_fixed = true;

            sim.Rods[row]->fixed_by_crossing[col] = false;
            sim.Rods[n_row + col]->fixed_by_crossing[row] = false;
            // if (col % 2 == 0)
            {
                // crossing->sliding_ranges.push_back(Range(0, 0));    
                // crossing->sliding_ranges.push_back(Range(1.0/20.0 - 1e3, 1.0/20.0 - 1e3));
                // crossing->sliding_ranges.push_back(Range(1, 1));    
                crossing->sliding_ranges.push_back(Range(1, 1));    
                crossing->sliding_ranges.push_back(Range(0, 0));    
            }
            // else
            // {
            //     crossing->sliding_ranges.push_back(Range(0.02, 0.02));    
            //     crossing->sliding_ranges.push_back(Range(0, 0));
            // }
            
            sim.rod_crossings.push_back(crossing);
            odd_even_cnt++;
        }
    }    

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
    

    T r = 0.05 * sim.unit;
    TV center1, center2;
    sim.getCrossingPosition(0, center1);
    sim.getCrossingPosition(n_row * n_col - 1, center2);

    TV delta1 = TV(-0.05, -0.05, -1e-2) * sim.unit;

    auto circle1 = [r, center1, delta1](const TV& x, TV& delta, Vector<bool, 3>& mask)->bool
    {
        mask = Vector<bool, 3>(true, true, true);
        delta = delta1;
        return (x - center1).norm() < r;
    };

    TV delta2 = TV(0.0, 0.0, 0) * sim.unit;
    auto circle2 = [r, center2, delta2](const TV& x, TV& delta, Vector<bool, 3>& mask)->bool
    {
        mask = Vector<bool, 3>(true, true, true);
        delta = delta2;
        return (x - center2).norm() < r;

    };
    
    TV bottom_left, top_right;
    sim.computeBoundingBox(bottom_left, top_right);
    TV shear_y_left = TV(0.0, 0.5, 0.1) * sim.unit;
    TV shear_y_right = TV(0.0, 0.0, 0) * sim.unit;
    
    T rec_width = 0.1 * sim.unit;

    auto rec1 = [bottom_left, top_right, shear_y_left, rec_width](
        const TV& x, TV& delta, Vector<bool, 3>& mask)->bool
    {
        mask = Vector<bool, 3>(true, true, true);
        delta = shear_y_left;
        T one_third = (top_right[1] - bottom_left[1]) / 3.0;
        if (x[0] < bottom_left[0] + rec_width 
            // &&
            // (x[1] > bottom_left[1] + one_third && x[1] < bottom_left[1] + one_third * 2)
            )
            return true;
        return false;
    };

    auto rec2 = [bottom_left, top_right, shear_y_right, rec_width](
        const TV& x, TV& delta, Vector<bool, 3>& mask)->bool
    {
        mask = Vector<bool, 3>(true, true, true);
        delta = shear_y_right;
        T one_third = (top_right[1] - bottom_left[1]) / 3.0;

        if (x[0] > top_right[0] - rec_width 
            // && 
            // (x[1] > bottom_left[1] + one_third && x[1] < bottom_left[1] + one_third * 2)
            )
            return true;
        return false;
    };

    sim.fixRegionalDisplacement(circle1);
    sim.fixRegionalDisplacement(circle2);

    // sim.fixRegionalDisplacement(rec1);
    // sim.fixRegionalDisplacement(rec2);


    sim.fixCrossing();

    sim.perturb = VectorXT::Zero(sim.W.cols());
    for (auto& crossing : sim.rod_crossings)
    {
        if (crossing->is_fixed)
            continue;
        Offset off;
        sim.Rods[crossing->rods_involved.front()]->getEntry(crossing->node_idx, off);
        T r = static_cast <T> (rand()) / static_cast <T> (RAND_MAX);
        int z_off = sim.Rods[crossing->rods_involved.front()]->reduced_map[off[3-1]];
        sim.perturb[z_off] += 0.001 * (r - 0.5) * sim.unit;
        // sim.perturb[z_off] += 0.001 * r * sim.unit;
    }

    sim.dq = VectorXT::Zero(dof_cnt);
}

void Scene::buildFullScaleSquareScene(int sub_div)
{
    auto unit_yarn_map = sim.yarn_map;
    sim.yarn_map.clear();
    
    clearSimData();

    sim.add_rotation_penalty = true;
    sim.add_pbc_bending = false;
    sim.add_pbc_twisting = false;
    sim.add_pbc = false;

    sim.add_contact_penalty=true;
    sim.new_frame_work = true;
    sim.add_eularian_reg = true;

    sim.ke = 1e-6;

    sim.unit = 0.09;
    sim.visual_R = 0.005;

    std::vector<Eigen::Triplet<T>> w_entry;
    int full_dof_cnt = 0;
    int node_cnt = 0;
    int rod_cnt = 0;

    std::vector<TV> nodal_positions;


    T square_width = 0.012; 
    T overlap = square_width * 0.3;

    auto addCrossingData = [&](int crossing_idx, int rod_idx, int location)
    {
        sim.rod_crossings[crossing_idx]->is_fixed = true;
        sim.rod_crossings[crossing_idx]->rods_involved.push_back(rod_idx);
        sim.rod_crossings[crossing_idx]->on_rod_idx[rod_idx] = location;
        sim.rod_crossings[crossing_idx]->sliding_ranges.push_back(Range::Zero());
    };

    auto addSquare = [&](const TV& bottom_left)
    {
        addCrossingPoint(nodal_positions, bottom_left, full_dof_cnt, node_cnt);
        TV bottom_right = bottom_left + TV(square_width, 0, 0);
        addCrossingPoint(nodal_positions, bottom_right, full_dof_cnt, node_cnt);
        TV top_right = bottom_left + TV(square_width, square_width, 0);
        addCrossingPoint(nodal_positions, top_right, full_dof_cnt, node_cnt);
        TV top_left = bottom_left + TV(0, square_width, 0);
        addCrossingPoint(nodal_positions, top_left, full_dof_cnt, node_cnt);
    };


    auto addInnerSquare = [&](const TV& bottom_left)
    {
        TV v0, v1;

        // add bottom left corner
        addCrossingPoint(nodal_positions, bottom_left, full_dof_cnt, node_cnt);
        v0 = bottom_left + TV(overlap, 0, 0);
        v1 = bottom_left + TV(0, overlap, 0);
        addCrossingPoint(nodal_positions, v0, full_dof_cnt, node_cnt);
        addCrossingPoint(nodal_positions, v1, full_dof_cnt, node_cnt);
        
        // add bottom right corner
        TV bottom_right = bottom_left + TV(square_width, 0, 0);
        addCrossingPoint(nodal_positions, bottom_right, full_dof_cnt, node_cnt);
        v0 = bottom_right - TV(overlap, 0, 0);
        v1 = bottom_right + TV(0, overlap, 0);
        addCrossingPoint(nodal_positions, v0, full_dof_cnt, node_cnt);
        addCrossingPoint(nodal_positions, v1, full_dof_cnt, node_cnt);
        
        TV top_right = bottom_left + TV(square_width, square_width, 0);
        addCrossingPoint(nodal_positions, top_right, full_dof_cnt, node_cnt);
        v0 = top_right - TV(overlap, 0, 0);
        v1 = top_right - TV(0, overlap, 0);
        addCrossingPoint(nodal_positions, v0, full_dof_cnt, node_cnt);
        addCrossingPoint(nodal_positions, v1, full_dof_cnt, node_cnt);
        
        TV top_left = bottom_left + TV(0, square_width, 0);
        addCrossingPoint(nodal_positions, top_left, full_dof_cnt, node_cnt);
        v0 = top_left + TV(overlap, 0, 0);
        v1 = top_left - TV(0, overlap, 0);
        addCrossingPoint(nodal_positions, v0, full_dof_cnt, node_cnt);
        addCrossingPoint(nodal_positions, v1, full_dof_cnt, node_cnt);
        
    };

    auto addRodsForASquare = [&](int bottom_left_node_idx)
    {
        TV bottom_left = nodal_positions[bottom_left_node_idx];
        TV bottom_right = nodal_positions[bottom_left_node_idx + 1];
        TV top_right = nodal_positions[bottom_left_node_idx + 2];
        TV top_left = nodal_positions[bottom_left_node_idx + 3];

        addAStraightRod(bottom_left, bottom_right, 
            {bottom_left, bottom_right}, {bottom_left_node_idx, bottom_left_node_idx + 1},
            sub_div, full_dof_cnt, node_cnt, rod_cnt );
        
        addCrossingData(bottom_left_node_idx, rod_cnt - 1, 0);
        addCrossingData(bottom_left_node_idx + 1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->numSeg());

        addAStraightRod(bottom_right, top_right, 
            {bottom_right, top_right}, {bottom_left_node_idx + 1, bottom_left_node_idx + 2},
            sub_div, full_dof_cnt, node_cnt, rod_cnt );
        addCrossingData(bottom_left_node_idx + 1, rod_cnt - 1, 0);
        addCrossingData(bottom_left_node_idx + 2, rod_cnt - 1, sim.Rods[rod_cnt - 1]->numSeg());
        
        addAStraightRod(top_right, top_left, 
            {top_right, top_left}, {bottom_left_node_idx + 2, bottom_left_node_idx + 3},
            sub_div, full_dof_cnt, node_cnt, rod_cnt );
        addCrossingData(bottom_left_node_idx + 2, rod_cnt - 1, 0);
        addCrossingData(bottom_left_node_idx + 3, rod_cnt - 1, sim.Rods[rod_cnt - 1]->numSeg());

        addAStraightRod(top_left, bottom_left, 
            {top_left, bottom_left}, {bottom_left_node_idx + 3, bottom_left_node_idx},
            sub_div, full_dof_cnt, node_cnt, rod_cnt );
        addCrossingData(bottom_left_node_idx + 3, rod_cnt - 1, 0);
        addCrossingData(bottom_left_node_idx, rod_cnt - 1, sim.Rods[rod_cnt - 1]->numSeg());
    };

    int n_row = 5, n_col = 5;
    
    T length_x = square_width * n_col + (square_width - 2.0 * overlap) * (n_col - 1);
    T length_y = length_x;


    // "longer" rows
    for (int row = 0; row < n_row; row++)
    {
        for (int col = 0; col < n_col; col++)
        {
            T left = col * 2.0 * (square_width - overlap);
            T bottom = row * 2.0 * (square_width - overlap);
            TV bottom_left = TV(left, bottom, 0.0);
            addSquare(bottom_left);
        }    
    }

    for (int row = 0; row < n_row - 1; row++)
    {
        for (int col = 0; col < n_col - 1; col++)
        {
            T left = square_width - overlap + col * 2.0 * (square_width - overlap);
            T bottom = square_width - overlap + row * 2.0 * (square_width - overlap);
            TV bottom_left = TV(left, bottom, 0.0);
            // std::cout << bottom_left.transpose() << std::endl;
            addInnerSquare(bottom_left);
        }    
    }

    // std::cout << nodal_positions.size() << std::endl;

    // add Boundary rods first
    // these are the top row and bottom row

    for (int col = 0; col < n_col; col++)
    {
        int idx0 = (0 * n_col + col) * 4; // 4 is four nodes per square
        int idx1 = ((n_row - 1) * n_col + col) * 4; // 4 is four nodes per square
        TV v0 = nodal_positions[idx0], v1 = nodal_positions[idx1 + 3];
        TV v0_next = nodal_positions[idx0 + 1], v1_next = nodal_positions[idx1 + 2]; 
        addAStraightRod(v0, v0_next, 
            {v0, v0_next}, {idx0, idx0 + 1},
            sub_div, full_dof_cnt, node_cnt, rod_cnt );
        addCrossingData(idx0, rod_cnt - 1, 0);
        addCrossingData(idx0 + 1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->numSeg());
    }

    for (int col = n_col - 1; col > -1; col--)
    {
        
        int idx1 = ((n_row - 1) * n_col + col) * 4; // 4 is four nodes per square
        TV v1 = nodal_positions[idx1 + 2];
        TV v1_next = nodal_positions[idx1 + 3]; 
        
        addAStraightRod(v1, v1_next, 
            {v1, v1_next}, {idx1 + 2, idx1 + 3},
            sub_div, full_dof_cnt, node_cnt, rod_cnt );
        addCrossingData(idx1 + 2, rod_cnt - 1, 0);
        addCrossingData(idx1 + 3, rod_cnt - 1, sim.Rods[rod_cnt - 1]->numSeg());
    }
    
    for (int col = 0; col < n_col; col++)
    {
        auto rod0 = sim.Rods[col];
        auto rod1 = sim.Rods[n_col + n_col - col - 1];
        Offset end0, end1;
        // std::cout << rod0->indices.front() << " " << rod1->indices.back() << std::endl;
        rod0->frontOffset(end0); rod1->backOffset(end1);
        sim.pbc_pairs.push_back(std::make_pair(1, std::make_pair(end0, end1)));
        if (col == 0)
        {
            sim.pbc_pairs_reference[1] = std::make_pair(std::make_pair(end0, end1), 
                std::make_pair(rod0->rod_id, rod1->rod_id));
        }
        rod0->backOffset(end0); rod1->frontOffset(end1);
        // std::cout << rod0->indices.back() << " " << rod1->indices.front() << std::endl;
        sim.pbc_pairs.push_back(std::make_pair(1, std::make_pair(end0, end1)));
    }
    

    // these are the left and right most rows
    for (int row = n_row - 1; row > -1; row--)
    {
        int idx0 = (row * n_col + 0) * 4; // 4 is four nodes per square
        int idx1 = (row * n_col + n_col - 1) * 4; // 4 is four nodes per square
        
        TV v0 = nodal_positions[idx0 + 3];
        TV v0_next = nodal_positions[idx0 + 0];

        addAStraightRod(v0, v0_next, 
            {v0, v0_next}, {idx0 + 3, idx0},
            sub_div, full_dof_cnt, node_cnt, rod_cnt );
        addCrossingData(idx0 + 3, rod_cnt - 1, 0);
        addCrossingData(idx0, rod_cnt - 1, sim.Rods[rod_cnt - 1]->numSeg());
    }

    
    for (int row = 0; row < n_row; row++)
    {
        
        int idx1 = (row * n_col + n_col - 1) * 4; // 4 is four nodes per square
        
        TV v1 = nodal_positions[idx1 + 1];
        TV v1_next = nodal_positions[idx1 + 2]; 

        addAStraightRod(v1, v1_next, 
            {v1, v1_next}, {idx1 + 1, idx1 + 2},
            sub_div, full_dof_cnt, node_cnt, rod_cnt );
        addCrossingData(idx1 + 1, rod_cnt - 1, 0);
        addCrossingData(idx1 + 2, rod_cnt - 1, sim.Rods[rod_cnt - 1]->numSeg());
    }

    for (int row = 0; row < n_row; row++)
    {
        auto rod0 = sim.Rods[row + 2 * n_col];
        auto rod1 = sim.Rods[2 * n_col + 2 * n_row - row - 1];

        Offset end0, end1;
        rod0->frontOffset(end0); rod1->backOffset(end1);
        // std::cout << rod0->indices.front() << " " << rod1->indices.back() << std::endl;
        sim.pbc_pairs.push_back(std::make_pair(0, std::make_pair(end0, end1)));
        if (row == 0)
        {
            sim.pbc_pairs_reference[0] = std::make_pair(std::make_pair(end0, end1), 
                std::make_pair(rod0->rod_id, rod1->rod_id));
        }
        rod0->backOffset(end0); rod1->frontOffset(end1);
        sim.pbc_pairs.push_back(std::make_pair(0, std::make_pair(end0, end1)));
        // std::cout << rod0->indices.back() << " " << rod1->indices.front() << std::endl;
    }
    
    for (int col = 0; col < n_col - 1; col++)
    {
        // vertical ones
        for (int row = 0; row < n_row - 1; row++)    
        {
            int idx0 = (row * n_col + col) * 4 + 1; // 4 is four nodes per square
            int idx1 = (row * n_col + col) * 4 + 2; 

            int idx_middle1 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 1;
            int idx_middle2 = n_row * n_col * 4 + ((row + 1) * (n_col - 1) + col) * 12 + 1;
            // std::cout << idx0 << " " << idx_middle1 << " " << idx1 << std::endl;

            TV v0 = nodal_positions[idx0], v1 = nodal_positions[idx1];
            TV v_middle1 = nodal_positions[idx_middle1];
            TV v_middle2;
            if (row == 0)
            {
                addAStraightRod(v0, v1, 
                    {v0, v_middle1, v1}, {idx0, idx_middle1, idx1},
                    sub_div, full_dof_cnt, node_cnt, rod_cnt );
                addCrossingData(idx0, rod_cnt - 1, 0);
                addCrossingData(idx_middle1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[1]);
                addCrossingData(idx1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->numSeg());
            }
            if (row == n_row - 2)
            {
                idx0 = ((row + 1) * n_col + col) * 4 + 1; // 4 is four nodes per square
                idx1 = ((row + 1) * n_col + col) * 4 + 2; 

                idx_middle1 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 10;
                v0 = nodal_positions[idx0], v1 = nodal_positions[idx1];
                v_middle1 = nodal_positions[idx_middle1];
                addAStraightRod(v0, v1, 
                    {v0, v_middle1, v1}, {idx0, idx_middle1, idx1},
                    sub_div, full_dof_cnt, node_cnt, rod_cnt );
                addCrossingData(idx0, rod_cnt - 1, 0);
                addCrossingData(idx_middle1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[1]);
                addCrossingData(idx1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->numSeg());

            }
            if (row != n_row - 2)
            {
                idx0 = ((row + 1) * n_col + col) * 4 + 1; // 4 is four nodes per square
                idx1 = ((row + 1) * n_col + col) * 4 + 2; 
                
                idx_middle1 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 10;
                idx_middle2 = n_row * n_col * 4 + ((row + 1) * (n_col - 1) + col) * 12 + 1;
                v0 = nodal_positions[idx0], v1 = nodal_positions[idx1];
                v_middle1 = nodal_positions[idx_middle1]; v_middle2 = nodal_positions[idx_middle2];
                
                addAStraightRod(v0, v1, 
                    {v0, v_middle1, v_middle2, v1}, {idx0, idx_middle1, idx_middle2,  idx1},
                    sub_div, full_dof_cnt, node_cnt, rod_cnt );
                addCrossingData(idx0, rod_cnt - 1, 0);
                addCrossingData(idx_middle1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[1]);
                addCrossingData(idx_middle2, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[2]);
                addCrossingData(idx1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->numSeg());
            }
        }
    }

    // add inner rods
    for (int col = 0; col < n_col - 1; col++)
    {
        // vertical ones
        for (int row = 0; row < n_row - 1; row++)
        {
            //these are right column of each square
            int idx0 = (row * n_col + col) * 4 + 1; // 4 is four nodes per square
            int idx1 = (row * n_col + col) * 4 + 2; 

            int idx_middle1 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 1;
            int idx_middle2 = n_row * n_col * 4 + ((row + 1) * (n_col - 1) + col) * 12 + 1;
            // std::cout << idx0 << " " << idx_middle1 << " " << idx1 << std::endl;

            TV v0 = nodal_positions[idx0], v1 = nodal_positions[idx1];
            TV v_middle1 = nodal_positions[idx_middle1];
            TV v_middle2;
            if (row == 0)
            {
                

                idx0 = (row * n_col + col + 1) * 4 + 0; // 4 is four nodes per square
                idx1 = (row * n_col + col + 1) * 4 + 3; 

                idx_middle1 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 4;
                v0 = nodal_positions[idx0], v1 = nodal_positions[idx1];
                v_middle1 = nodal_positions[idx_middle1];
                addAStraightRod(v0, v1, 
                    {v0, v_middle1, v1}, {idx0, idx_middle1, idx1},
                    sub_div, full_dof_cnt, node_cnt, rod_cnt );
                addCrossingData(idx0, rod_cnt - 1, 0);
                addCrossingData(idx_middle1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[1]);
                addCrossingData(idx1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->numSeg());
            }
            if (row == n_row - 2)
            {
                

                idx0 = ((row + 1) * n_col + col + 1) * 4 + 0; // 4 is four nodes per square
                idx1 = ((row + 1) * n_col + col + 1) * 4 + 3; 

                idx_middle1 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 7;
                v0 = nodal_positions[idx0], v1 = nodal_positions[idx1];
                v_middle1 = nodal_positions[idx_middle1];
                addAStraightRod(v0, v1, 
                    {v0, v_middle1, v1}, {idx0, idx_middle1, idx1},
                    sub_div, full_dof_cnt, node_cnt, rod_cnt );
                addCrossingData(idx0, rod_cnt - 1, 0);
                addCrossingData(idx_middle1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[1]);
                addCrossingData(idx1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->numSeg());

            }
            if (row != n_row - 2)
            {
                

                idx0 = ((row + 1) * n_col + col + 1) * 4 + 0; // 4 is four nodes per square
                idx1 = ((row + 1) * n_col + col + 1) * 4 + 3; 

                idx_middle1 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 7;
                idx_middle2 = n_row * n_col * 4 + ((row + 1) * (n_col - 1) + col) * 12 + 4;
                v0 = nodal_positions[idx0], v1 = nodal_positions[idx1];
                v_middle1 = nodal_positions[idx_middle1]; v_middle2 = nodal_positions[idx_middle2];
                addAStraightRod(v0, v1, 
                    {v0, v_middle1, v_middle2, v1}, {idx0, idx_middle1, idx_middle2,  idx1},
                    sub_div, full_dof_cnt, node_cnt, rod_cnt );
                addCrossingData(idx0, rod_cnt - 1, 0);
                addCrossingData(idx_middle1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[1]);
                addCrossingData(idx_middle2, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[2]);
                addCrossingData(idx1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->numSeg());
            }
        }
    }

    
    

    // add inner rods
    for (int row = 0; row < n_row - 1; row++)
    {
        // vertical ones
        for (int col = 0; col < n_col - 1; col++)
        {
            //these are right column of each square
            int idx0 = (row * n_col + col) * 4 + 1; // 4 is four nodes per square
            int idx1 = (row * n_col + col) * 4 + 2; 

            int idx_middle1 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 1;
            int idx_middle2 = n_row * n_col * 4 + ((row + 1) * (n_col - 1) + col) * 12 + 1;
            // std::cout << idx0 << " " << idx_middle1 << " " << idx1 << std::endl;

            TV v0 = nodal_positions[idx0], v1 = nodal_positions[idx1];
            TV v_middle1 = nodal_positions[idx_middle1];
            TV v_middle2;
            

            if (col == 0)
            {
                idx0 = (row * n_col + col) * 4 + 3; // 4 is four nodes per square
                idx1 = (row * n_col + col) * 4 + 2; 

                idx_middle1 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 2;
                v0 = nodal_positions[idx0], v1 = nodal_positions[idx1];
                v_middle1 = nodal_positions[idx_middle1];
                addAStraightRod(v0, v1, 
                    {v0, v_middle1, v1}, {idx0, idx_middle1, idx1},
                    sub_div, full_dof_cnt, node_cnt, rod_cnt );
                addCrossingData(idx0, rod_cnt - 1, 0);
                addCrossingData(idx_middle1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[1]);
                addCrossingData(idx1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->numSeg());
            }
            if (col == n_col - 2)
            {
                idx0 = (row * n_col + col + 1) * 4 + 3; // 4 is four nodes per square
                idx1 = (row * n_col + col + 1) * 4 + 2; 

                idx_middle1 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 5;
                v0 = nodal_positions[idx0], v1 = nodal_positions[idx1];
                v_middle1 = nodal_positions[idx_middle1];
                addAStraightRod(v0, v1, 
                    {v0, v_middle1, v1}, {idx0, idx_middle1, idx1},
                    sub_div, full_dof_cnt, node_cnt, rod_cnt );
                addCrossingData(idx0, rod_cnt - 1, 0);
                addCrossingData(idx_middle1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[1]);
                addCrossingData(idx1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->numSeg());
            }
            if (n_col - 2 != col)
            {
                idx0 = (row * n_col + col + 1) * 4 + 3; // 4 is four nodes per square
                idx1 = (row * n_col + col + 1) * 4 + 2; 
                
                idx_middle1 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 5;
                idx_middle2 = n_row * n_col * 4 + (row * (n_col - 1) + col + 1) * 12 + 2;
                v0 = nodal_positions[idx0], v1 = nodal_positions[idx1];
                v_middle1 = nodal_positions[idx_middle1]; v_middle2 = nodal_positions[idx_middle2];
                
                addAStraightRod(v0, v1, 
                    {v0, v_middle1, v_middle2, v1}, {idx0, idx_middle1, idx_middle2,  idx1},
                    sub_div, full_dof_cnt, node_cnt, rod_cnt );
                addCrossingData(idx0, rod_cnt - 1, 0);
                addCrossingData(idx_middle1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[1]);
                addCrossingData(idx_middle2, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[2]);
                addCrossingData(idx1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->numSeg());
            }
        }
    }




    // add inner rods
    for (int row = 0; row < n_row - 1; row++)
    {
        // vertical ones
        for (int col = 0; col < n_col - 1; col++)
        {
            //these are right column of each square
            int idx0 = (row * n_col + col) * 4 + 1; // 4 is four nodes per square
            int idx1 = (row * n_col + col) * 4 + 2; 

            int idx_middle1 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 1;
            int idx_middle2 = n_row * n_col * 4 + ((row + 1) * (n_col - 1) + col) * 12 + 1;
            // std::cout << idx0 << " " << idx_middle1 << " " << idx1 << std::endl;

            TV v0 = nodal_positions[idx0], v1 = nodal_positions[idx1];
            TV v_middle1 = nodal_positions[idx_middle1];
            TV v_middle2;
            
            if (col == 0)
            {
                

                idx0 = ((row + 1) * n_col + col) * 4 + 0; // 4 is four nodes per square
                idx1 = ((row + 1) * n_col + col) * 4 + 1; 

                idx_middle1 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 11;
                v0 = nodal_positions[idx0], v1 = nodal_positions[idx1];
                v_middle1 = nodal_positions[idx_middle1];
                addAStraightRod(v0, v1, 
                    {v0, v_middle1, v1}, {idx0, idx_middle1, idx1},
                    sub_div, full_dof_cnt, node_cnt, rod_cnt );
                addCrossingData(idx0, rod_cnt - 1, 0);
                addCrossingData(idx_middle1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[1]);
                addCrossingData(idx1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->numSeg());
            }
            if (col == n_col - 2)
            {
                
                idx0 = ((row + 1) * n_col + col + 1) * 4 + 0; // 4 is four nodes per square
                idx1 = ((row + 1) * n_col + col + 1) * 4 + 1; 

                idx_middle1 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 8;
                v0 = nodal_positions[idx0], v1 = nodal_positions[idx1];
                v_middle1 = nodal_positions[idx_middle1];
                addAStraightRod(v0, v1, 
                    {v0, v_middle1, v1}, {idx0, idx_middle1, idx1},
                    sub_div, full_dof_cnt, node_cnt, rod_cnt );
                addCrossingData(idx0, rod_cnt - 1, 0);
                addCrossingData(idx_middle1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[1]);
                addCrossingData(idx1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->numSeg());
            }
            if (n_col - 2 != col)
            {
                

                idx0 = ((row + 1) * n_col + col + 1) * 4 + 0; // 4 is four nodes per square
                idx1 = ((row + 1) * n_col + col + 1) * 4 + 1; 

                idx_middle1 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 8;
                idx_middle2 = n_row * n_col * 4 + (row * (n_col - 1) + col + 1) * 12 + 11;
                v0 = nodal_positions[idx0], v1 = nodal_positions[idx1];
                v_middle1 = nodal_positions[idx_middle1]; v_middle2 = nodal_positions[idx_middle2];
                addAStraightRod(v0, v1, 
                    {v0, v_middle1, v_middle2, v1}, {idx0, idx_middle1, idx_middle2,  idx1},
                    sub_div, full_dof_cnt, node_cnt, rod_cnt );
                addCrossingData(idx0, rod_cnt - 1, 0);
                addCrossingData(idx_middle1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[1]);
                addCrossingData(idx_middle2, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[2]);
                addCrossingData(idx1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->numSeg());
            }
        }
    }



    // std::cout << "inner" << std::endl;
    for (int row = 0; row < n_row - 1; row++)
    {
        for (int col = 0; col < n_col - 1; col++)
        {
            int idx0 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12;
            int idx1 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 3; 
            int middle1 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 1; 
            int middle2 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 4; 

            
            TV v0 = nodal_positions[idx0], v1 = nodal_positions[idx1];
            TV v_middle1 = nodal_positions[middle1], v_middle2 = nodal_positions[middle2];
            addAStraightRod(v0, v1, 
                    {v0, v_middle1, v_middle2, v1}, {idx0, middle1, middle2, idx1},
                    sub_div, full_dof_cnt, node_cnt, rod_cnt );
            addCrossingData(idx0, rod_cnt - 1, 0);
            addCrossingData(middle1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[1]);
            addCrossingData(middle2, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[2]);
            addCrossingData(idx1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->numSeg());

            idx0 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 3;
            idx1 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 6; 
            middle1 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 5; 
            middle2 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 8; 

            v0 = nodal_positions[idx0], v1 = nodal_positions[idx1];
            v_middle1 = nodal_positions[middle1], v_middle2 = nodal_positions[middle2];
            addAStraightRod(v0, v1, 
                    {v0, v_middle1, v_middle2, v1}, {idx0, middle1, middle2, idx1},
                    sub_div, full_dof_cnt, node_cnt, rod_cnt );
            addCrossingData(idx0, rod_cnt - 1, 0);
            addCrossingData(middle1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[1]);
            addCrossingData(middle2, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[2]);
            addCrossingData(idx1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->numSeg());

            idx0 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 6;
            idx1 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 9; 
            middle1 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 7; 
            middle2 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 10; 

            v0 = nodal_positions[idx0], v1 = nodal_positions[idx1];
            v_middle1 = nodal_positions[middle1], v_middle2 = nodal_positions[middle2];
            addAStraightRod(v0, v1, 
                    {v0, v_middle1, v_middle2, v1}, {idx0, middle1, middle2, idx1},
                    sub_div, full_dof_cnt, node_cnt, rod_cnt );
            addCrossingData(idx0, rod_cnt - 1, 0);
            addCrossingData(middle1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[1]);
            addCrossingData(middle2, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[2]);
            addCrossingData(idx1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->numSeg());
            
            idx0 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 9;
            idx1 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 0; 
            middle1 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 11; 
            middle2 = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12 + 2; 

            v0 = nodal_positions[idx0], v1 = nodal_positions[idx1];
            v_middle1 = nodal_positions[middle1], v_middle2 = nodal_positions[middle2];
            addAStraightRod(v0, v1, 
                    {v0, v_middle1, v_middle2, v1}, {idx0, middle1, middle2, idx1},
                    sub_div, full_dof_cnt, node_cnt, rod_cnt );
            addCrossingData(idx0, rod_cnt - 1, 0);
            addCrossingData(middle1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[1]);
            addCrossingData(middle2, rod_cnt - 1, sim.Rods[rod_cnt - 1]->dof_node_location[2]);
            addCrossingData(idx1, rod_cnt - 1, sim.Rods[rod_cnt - 1]->numSeg());
        }
    }

    for (auto& rod : sim.Rods)
        rod->fixed_by_crossing = std::vector<bool>(rod->dof_node_location.size(), true);

    // for (int row = 0; row < n_row - 1; row++)
    // {
    //     for (int col = 0; col < n_col - 1; col++)
    //     {
    //         int base = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12;
    //         for (int corner = 0; corner < 4; corner++)
    //         {
                
    //             auto crossing = sim.rod_crossings[base + corner * 3 + 1];
    //             crossing->is_fixed = false;
    //             crossing->sliding_ranges[1] = Range::Ones();
    //             sim.Rods[crossing->rods_involved[1]]->fixed_by_crossing[1] = false;
    //             sim.Rods[crossing->rods_involved[1]]->fixed_by_crossing[2] = false;

    //             sim.Rods[crossing->rods_involved[0]]->fixed_by_crossing[1] = false;

    //             crossing = sim.rod_crossings[base + corner * 3 + 2];
    //             crossing->is_fixed = false;
    //             crossing->sliding_ranges[0] = Range::Ones();
    //             sim.Rods[crossing->rods_involved[0]]->fixed_by_crossing[1] = false;

    //             sim.Rods[crossing->rods_involved[1]]->fixed_by_crossing[1] = false;
    //             sim.Rods[crossing->rods_involved[1]]->fixed_by_crossing[2] = false;

    //         }
    //     }
    // }

    

    

    int dof_cnt = 0;
    markCrossingDoF(w_entry, dof_cnt);
    // std::cout << "mark dof" << std::endl;
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
    
    if (sim.add_pbc)
        sim.Rods[0]->fixPointLagrangianByID(0, TV::Zero(), Mask::Ones(), sim.dirichlet_dof);

    TV bottom_left, top_right;
    sim.computeBoundingBox(bottom_left, top_right);

    TV shear_x_right = TV(0.1, 0.0, 0.0) * length_x;
    TV shear_x_left = TV(0.0, 0.0, 0) * sim.unit;


    T rec_width = 0.015 * sim.unit;

    auto rec1 = [bottom_left, top_right, shear_x_left, rec_width](
        const TV& x, TV& delta, Vector<bool, 3>& mask)->bool
    {
        mask = Vector<bool, 3>(true, true, true);
        delta = shear_x_left;
        if (x[0] < bottom_left[0] + rec_width 
            )
            return true;
        return false;
    };

    auto rec2 = [bottom_left, top_right, shear_x_right, rec_width](
        const TV& x, TV& delta, Vector<bool, 3>& mask)->bool
    {
        mask = Vector<bool, 3>(true, true, true);
        delta = shear_x_right;

        if (x[0] > top_right[0] - rec_width)
            return true;
        return false;
    };

    if (!sim.add_pbc)
    {
        sim.fixRegionalDisplacement(rec2);
        sim.fixRegionalDisplacement(rec1);
    }

    sim.fixCrossing();
    // std::cout << "fix crossing" << std::endl;
    sim.perturb = VectorXT::Zero(sim.W.cols());
    for (auto& crossing : sim.rod_crossings)
    {
        Offset off;
        sim.Rods[crossing->rods_involved.front()]->getEntry(crossing->node_idx, off);
        T r = static_cast <T> (rand()) / static_cast <T> (RAND_MAX);
        int z_off = sim.Rods[crossing->rods_involved.front()]->reduced_map[off[3-1]];
        sim.perturb[z_off] += 0.001 * (r - 0.5) * sim.unit;
        
    }
    sim.dq = VectorXT::Zero(dof_cnt);
}

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

    sim.dq = VectorXT::Zero(dof_cnt);
}

void Scene::buildStraightRodScene(int sub_div)
{
    auto unit_yarn_map = sim.yarn_map;
    sim.yarn_map.clear();
    sim.add_rotation_penalty = false;
    sim.add_pbc_bending = false;
    sim.new_frame_work = true;

    clearSimData();

    std::vector<Eigen::Triplet<T>> w_entry;
    int full_dof_cnt = 0;
    int node_cnt = 0;
    int rod_cnt = 0;
    
    TV from = TV(0, 0.5, 0) * sim.unit;
    TV to = TV(2, 0.5001, 0.001) * sim.unit;

    std::vector<int> passing_points_id;
    std::vector<TV> passing_points;

    addAStraightRod(from, to, passing_points, passing_points_id, 
            sub_div, full_dof_cnt, node_cnt, rod_cnt);
    
    int dof_cnt;
    for (auto& rod : sim.Rods) rod->markDoF(w_entry, dof_cnt);
    
    appendThetaAndJointDoF(w_entry, full_dof_cnt, dof_cnt);
    
    sim.rest_states = deformed_states;

    sim.W = StiffnessMatrix(full_dof_cnt, dof_cnt);
    sim.W.setFromTriplets(w_entry.begin(), w_entry.end());


    int cnt = 0;
    for (auto& rod : sim.Rods)
    {
        // std::cout << rod->kt << std::endl;
        // rod->kt  =0 ;
        rod->fixEndPointEulerian(sim.dirichlet_dof);
        sim.dirichlet_dof[rod->theta_reduced_dof_start_offset] = 0;
        // sim.dirichlet_dof[rod->theta_reduced_dof_start_offset + rod->indices.size()-1] = 0;
        rod->setupBishopFrame();
        Offset end0, end1;
        rod->frontOffset(end0); rod->backOffset(end1);
        sim.pbc_pairs.push_back(std::make_pair(0, std::make_pair(end0, end1)));
    }

    Offset end0, end1;
    sim.Rods[0]->frontOffset(end0); sim.Rods[0]->backOffset(end1);
    sim.pbc_pairs_reference[0] = std::make_pair(std::make_pair(end0, end1), std::make_pair(0, 0));

    sim.Rods[0]->fixPointLagrangian(0, TV::Zero(), sim.dirichlet_dof);
    sim.Rods[0]->fixPointLagrangian(1, TV::Zero(), sim.dirichlet_dof);

    sim.dq = VectorXT::Zero(dof_cnt);
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

void Scene::addCrossingPoint(std::vector<TV>& existing_nodes, 
        const TV& point, int& full_dof_cnt, int& node_cnt)
{
    sim.rod_crossings.push_back(new RodCrossing(node_cnt, std::vector<int>())); 
    deformed_states.conservativeResize(full_dof_cnt + 3);
    deformed_states.template segment<3>(full_dof_cnt) = point;
    existing_nodes.push_back(point);
    full_dof_cnt += 3;
    node_cnt++;
}