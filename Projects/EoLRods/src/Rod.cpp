#include "../include/Rod.h"

void Rod::setupBishopFrame()
{
    // rest_tangents.resize(indices.size() - 1);
    reference_frame_us.resize(indices.size() - 1);
    reference_twist.resize(indices.size() - 1);
    prev_tangents.resize(indices.size() - 1);
    reference_twist.setZero();
    reference_angles = full_states.segment(theta_dof_start_offset, indices.size() - 1);
    TV prev_tangent;
    iterateSegments([&](int node_i, int node_j, int segment){
        TV xi, xj;
        x(node_i, xi); x(node_j, xj);
        TV tangent = (xj - xi).normalized();
        // rest_tangents[segment] = tangent;
        TV u, v;
        if (segment == 0)
        {
            if(colinear(tangent, TV(0, 0, 1)))
                v = tangent.cross(TV(0, 1, 0));
            else
                v = tangent.cross(TV(0, 0, 1));
            u = v.cross(tangent);
            
        }
        else
        {
            u = parallelTransportOrthonormalVector(
                reference_frame_us[segment - 1], 
                prev_tangent, 
                tangent);
        }
        prev_tangent = tangent;
        prev_tangents[segment] = tangent;
        reference_frame_us[segment] = u;
        if (segment)
            reference_twist[segment] = computeReferenceTwist(tangent, prev_tangent, segment);
    });    
    
    rest_tangents = prev_tangents;
    rest_normals = reference_frame_us;
}


void Rod::rotateReferenceFrameToLastNewtonStepAndComputeReferenceTwsit()
{
    TV prev_tangent;
    iterateSegments([&](int node_i, int node_j, int segment){
        TV xi, xj;
        x(node_i, xi); x(node_j, xj);
        TV tangent = (xj - xi).normalized();
        TV tangent_last_time_step = prev_tangents[segment];
        reference_frame_us[segment] = parallelTransportOrthonormalVector(
                                            reference_frame_us[segment], 
                                            tangent_last_time_step,
                                            tangent);
        if (segment)
            reference_twist[segment] = computeReferenceTwist(tangent, prev_tangent, segment);
            // std::cout << computeReferenceTwist(tangent, prev_tangent, segment) << std::endl;
        prev_tangent = tangent;
        prev_tangents[segment] = tangent;
            // reference_twist[segment] = computeReferenceTwist(tangent, prev_tangent, segment);
    });
}




void Rod::curvatureBinormal(const TV& t1, const TV& t2, TV& kb)
{
    T denominator = 1. + t1.dot(t2);
    if (denominator <= 0. || denominator < std::numeric_limits<T>::epsilon()) 
    {
        if (denominator <= 0.) 
            denominator = 1. + t1.normalized().dot(t2.normalized());
        

        if (denominator <= 0.) 
        {
            std::cerr << "CurvatureBinormals::compute() denominator == "
                        << denominator << " t1 = " << t1
                        << " t2 = " << t2 << std::endl;

            kb = TV::Constant(
                std::numeric_limits<T>::infinity());  // Should not be
                                                            // accepted.
        } 
        else 
        {
            TV normal;
            if(colinear(t1, TV(0, 0, 1)))
                normal = t1.cross(TV(0, 1, 0));
            else
                normal = t1.cross(TV(0, 0, 1));
            kb = 4. * std::tan(.5 * std::acos(denominator - 1.)) * normal;
        }
    } 
    else 
    {
        kb = 2.0 * t1.cross(t2) / denominator;
    }
    
}


T Rod::computeReferenceTwist(const TV& tangent, const TV& prev_tangent, int rod_idx)
{
    TV ut = parallelTransportOrthonormalVector(reference_frame_us[rod_idx-1], prev_tangent, tangent);
    TV u = reference_frame_us[rod_idx];
    // rotate by current value of reference twist
    T before_twist = reference_twist[rod_idx];
    rotateAxisAngle(ut, tangent, before_twist);
    T after_twist = before_twist + signedAngle(ut, u, tangent);
    // std::cout << after_twist << std::endl;
    return after_twist;
}

// the following codes assume all dof node has been added/marked

void Rod::markDoF(std::vector<Entry>& w_entry, int& dof_cnt)
{
    
    // std::cout << "[Rod" << rod_id << "]" << std::endl;
    int loop_id = 0;
    if(dof_node_location.size())
    {
        if (dof_node_location.front() != 0)
        {
            // add first node
            Offset offset = offset_map[indices.front()];
            // std::cout << "node " <<  indices.front() << " added all dof " << std::endl;
            for(int d = 0; d < 3 + 1; d++)
            {
                reduced_map[offset[d]] = dof_cnt;
                w_entry.push_back(Eigen::Triplet<T>(offset[d], dof_cnt++, 1.0));
            }
        }
    }
    else
    {
        Offset offset = offset_map[indices.front()];
        // std::cout << "node " <<  indices.front() << " added all dof " << std::endl;
        for(int d = 0; d < 3 + 1; d++)
        {
            reduced_map[offset[d]] = dof_cnt;
            w_entry.push_back(Eigen::Triplet<T>(offset[d], dof_cnt++, 1.0));
        }
    }
    
    // loop over each segment
    for (int i = 0; i < dof_node_location.size(); i++)
    {   
        // compute weight value for in-between nodes
        // using linear interpolation
        // std::cout << "left node " << indices[loop_id] << " right node " << indices[dof_node_location[i]] << std::endl;
        
        Offset offset_left_node = offset_map[indices[loop_id]];
        Offset offset_right_node = offset_map[indices[dof_node_location[i]]];
        T ui = full_states[offset_left_node[3]];
        T uj = full_states[offset_right_node[3]];

        for (int j = loop_id + 1; j < dof_node_location[i]; j++)
        {
            int current_global_idx = indices[j];
            // std::cout << "current node " << current_global_idx << std::endl;
            //push Lagrangian DoF first
            Offset offset = offset_map[current_global_idx];
            for(int d = 0; d < 3; d++)
            {
                reduced_map[offset[d]] = dof_cnt;
                // std::cout << "add lagrangian entry to " << offset[d] << " " << dof_cnt << std::endl;
                w_entry.push_back(Eigen::Triplet<T>(offset[d], dof_cnt++, 1.0));
            }
            // compute Eulerian weight
            T u = full_states[offset[3]];
            T alpha = (u - ui) / (uj - ui);
            // std::cout << "alpha " << alpha << std::endl;
            // std::cout << "add eulerian entry to " << offset[3] << " " << reduced_map[offset_left_node[3]] << std::endl;
            w_entry.push_back(Eigen::Triplet<T>(offset[3], reduced_map[offset_left_node[3]], 1.0 - alpha));
            w_entry.push_back(Eigen::Triplet<T>(offset[3], reduced_map[offset_right_node[3]], alpha));
        }
        
        loop_id = dof_node_location[i];
        
    }
    // last segment
    if (loop_id != indices.size() - 1 && !closed)
    {
        //last node is not a crossing node
        Offset offset = offset_map[indices.back()];
        // std::cout << "node " <<  indices.back() << " added all dof " << std::endl;
        for(int d = 0; d < 3 + 1; d++)
        {
            reduced_map[offset[d]] = dof_cnt;
            w_entry.push_back(Eigen::Triplet<T>(offset[d], dof_cnt++, 1.0));
        }
    }
    
    for (int j = loop_id + 1; j < indices.size() - 1; j++)
    {
        // std::cout << "rod id " << rod_id << std::endl;
        // std::cout << j << " " <<  indices.size() - 1 << std::endl;

        // std::cout << "left node " << indices[loop_id] << " right node " << indices.back() << std::endl;

        Offset offset_left_node = offset_map[indices[loop_id]];
        Offset offset_right_node = offset_map[indices.back()];

        T ui = full_states[offset_left_node[3]];
        T uj = full_states[offset_right_node[3]];
        int current_global_idx = indices[j];
        // std::cout << "current node " << current_global_idx << std::endl;
        //push Lagrangian DoF first
        Offset offset = offset_map[current_global_idx];
        for(int d = 0; d < 3; d++)
        {
            reduced_map[offset[d]] = dof_cnt;
            w_entry.push_back(Eigen::Triplet<T>(offset[d], dof_cnt++, 1.0));
        }
        // compute Eulerian weight
        T u = full_states[offset[3]];
        T alpha = (u - ui) / (uj - ui);
        
        if (closed) alpha *= -1;
        
        // std::cout << "alpha " << alpha << std::endl;
        // std::cout << "add eulerian entry to " << offset[3] << " " << reduced_map[offset_left_node[3]] << std::endl;
        w_entry.push_back(Eigen::Triplet<T>(offset[3], reduced_map[offset_left_node[3]], 1.0 - alpha));
        w_entry.push_back(Eigen::Triplet<T>(offset[3], reduced_map[offset_right_node[3]], alpha));
    }
}


bool Rod::isFixedNodeForPrinting(int node_idx, int rod_idx)
{
    std::vector<int> dof_node_idx;
    for (int location : dof_node_location)
        dof_node_idx.push_back(indices[location]);
    
    // if (dof_node_location.size())
    // {
    //     if (rod_idx < dof_node_location.front())
    //         return true;
    //     if (rod_idx >= dof_node_location.back())
    //         return true;
    // }

    
    
    auto iter = std::find(dof_node_idx.begin(), dof_node_idx.end(), node_idx);
    if (iter == dof_node_idx.end())
    {
        return true;
        int left_node = 0;
        int right_node = indices.size() - 1;
        for (int i = 0; i < dof_node_location.size(); i++)
        {
            if (dof_node_location[i] <= rod_idx)
            {
                left_node = i;
            }
            if ( dof_node_location[i] >= rod_idx)
            {
                right_node = i;
                break;
            }
        }

        // for (auto x : fixed_by_crossing)
        //         std::cout << x << " ";
        //     std::cout << std::endl;
            

        // for (int location : dof_node_location)
        //     std::cout << fixed_by_crossing[location] << " " << location << std::endl;
        // std::cout << "dd" << std::endl;
        // std::cout << left_node << " " << right_node << " " << rod_idx << std::endl;
        // std::cout << fixed_by_crossing[left_node] << " " << fixed_by_crossing[right_node] << std::endl;
        // std::getchar();
        
        // if((fixed_by_crossing[left_node] || left_node == 0) && (fixed_by_crossing[right_node] || right_node == indices.size() - 1))
        //     return true;
        // else
        //     return false;
        bool fix_left = fixed_by_crossing[left_node] || left_node == 0;
        bool fix_right = fixed_by_crossing[right_node] || right_node == indices.size() - 1;
        if (!fix_left && !fix_right)
            return false;
        return true;
    }
    int location = std::distance(dof_node_idx.begin(), iter);
    if (fixed_by_crossing[location])
        return true;
    return false;
}
