#include "../include/IMLSContact.h"
#include "../include/VecMatDef.h"

VectorXa IMLSContact::divisionGradient(AScalar& fx, AScalar& gx, VectorXa& dfxdx, VectorXa& dgxdx)
{
  return (gx*dfxdx-fx*dgxdx)/(pow(gx,2));
}


void IMLSContact::IMLS_local_gradient_to_global_gradient(
    VectorXa& local_grad,
    Eigen::VectorXi ids,
    int dim1,
    Eigen::VectorXd& grad)
{
    for (int i = 0; i < ids.size(); i++) 
    {
        if(ids(i) < 0) continue;
        grad.segment(dim * ids(i), dim) += local_grad.segment(dim * i, dim);
    }
}


void IMLSContact:: IMLS_local_gradient_to_global_gradient_and_multiply(
    VectorXa& local_grad,
    Eigen::VectorXi ids,
    int dim1,
    AScalar scale,
    VectorXa& global_vector,
    std::vector<Eigen::Triplet<double>>& triplets)
{
    for (int i = 0; i < ids.size(); i++) {
        if(ids(i) < 0) continue;
        for (int j = 0; j < global_vector.size(); j++) {
            if(global_vector(j) == 0) continue;
            for (int k = 0; k < dim; k++) {
                triplets.emplace_back(
                    dim * ids(i) + k, j,
                    scale*local_grad(dim * i + k)*global_vector(j));
            }
        }
    }

    for (int i = 0; i < global_vector.size(); i++) {
        if(global_vector(i) == 0) continue;
        for (int j = 0; j < ids.size(); j++) {
            if(ids(j) < 0) continue;
            for (int k = 0; k < dim; k++) {
                triplets.emplace_back(
                    i,  dim * ids(j) + k,
                     scale*local_grad(dim*j+k)*global_vector(i));
            }
        }
    }
}


void IMLSContact::IMLS_local_hessian_to_global_triplets(
    VectorXa& local_hessian,
    Eigen::VectorXi ids,
    int dim1,
    std::vector<Eigen::Triplet<double>>& triplets)
{
    for (int i = 0; i < ids.size(); i++) {
        if(ids(i) < 0) continue;
        for (int j = 0; j < ids.size(); j++) {
            if(ids(j) < 0) continue;
            for (int k = 0; k < dim; k++) {
                for (int l = 0; l < dim; l++) {
                    triplets.emplace_back(
                        dim * ids(i) + k, dim * ids(j) + l,
                        local_hessian(dim*ids.size()*(dim * i + k)+dim * j + l));
                }
            }
        }
    }
}


void IMLSContact::IMLS_local_hessian_matrix_to_global_triplets(
    MatrixXa& local_hessian,
    Eigen::VectorXi ids,
    int dim1,
    std::vector<Eigen::Triplet<double>>& triplets)
{
    for (int i = 0; i < ids.size(); i++) {
        if(ids(i) < 0) continue;
        for (int j = 0; j < ids.size(); j++) {
            if(ids(j) < 0) continue;
            for (int k = 0; k < dim; k++) {
                for (int l = 0; l < dim; l++) {
                    triplets.emplace_back(
                        dim * ids(i) + k, dim * ids(j) + l,
                        local_hessian(dim*i+k,dim*j+l));
                }
            }
        }
    }
}



void IMLSContact::IMLS_vector_muliplication_to_triplets(
    VectorXa& v1,
    VectorXa& v2,
    AScalar scale,
    std::vector<Eigen::Triplet<double>>& triplets)
{
    assert(v1.size() == v2.size());
    for(int i=0; i<v1.size(); ++i)
    {
        if(v1(i) == 0) continue;
        for(int j=0; j<v2.size(); ++j)
        {
            if(v2(j) == 0) continue;
            triplets.emplace_back(
            i, j, scale*v1(i)*v2(j));
        }
    }
}


void IMLSContact::BuildAcceleratorTreeSC()
{
    accelerator.clear();

    for(int i=0; i<num_nodes; ++i)
    {
        if(is_surface_vertex[i] == 0) continue;
        Eigen::VectorXd point = deformed.segment<3>(3*i);
        accelerator.insert(std::make_tuple(Point_d(point[0], point[1], point[2]), i));
    }
}


void IMLSContact::BuildAcceleratorTreeSCTest()
{
    accelerator.clear();

    for(int i=0; i<4225; ++i)
    {
        if(is_surface_vertex[i] == 0) continue;
        Eigen::VectorXd point = deformed.segment<3>(3*i);
        accelerator.insert(std::make_tuple(Point_d(point[0], point[1], point[2]), i));
    }
}


void IMLSContact::divisionHessian(AScalar& fx, AScalar& gx, VectorXa& dfxdx, VectorXa& dgxdx, std::vector<Eigen::Triplet<double>>& d2fxdx2, std::vector<Eigen::Triplet<double>>& d2gxdx2, std::vector<Eigen::Triplet<double>>& out)
{
   // h1 = d2f/dx2
   StiffnessMatrix h1(dim*num_nodes,dim*num_nodes);
   h1.setFromTriplets(d2fxdx2.begin(), d2fxdx2.end());

   // h2 = d2g/dx2
   StiffnessMatrix h2(dim*num_nodes,dim*num_nodes);
   h2.setFromTriplets(d2gxdx2.begin(), d2gxdx2.end());

   //h3 = dfdx*dgdx^T
   std::vector<Eigen::Triplet<double>> h3_t;
   IMLS_vector_muliplication_to_triplets(dfxdx,dgxdx,1.,h3_t);
   StiffnessMatrix h3(dim*num_nodes,dim*num_nodes);
   h3.setFromTriplets(h3_t.begin(), h3_t.end());

   //h4 = dgdx*dgdx^T
   std::vector<Eigen::Triplet<double>> h4_t;
   IMLS_vector_muliplication_to_triplets(dgxdx,dgxdx,1.,h4_t);
   StiffnessMatrix h4(dim*num_nodes,dim*num_nodes);
   h4.setFromTriplets(h4_t.begin(), h4_t.end());

   //h6 = dgdx*dfdx^T
   std::vector<Eigen::Triplet<double>> h6_t;
   IMLS_vector_muliplication_to_triplets(dgxdx,dfxdx,1.,h6_t);
   StiffnessMatrix h6(dim*num_nodes,dim*num_nodes);
   h6.setFromTriplets(h6_t.begin(), h6_t.end());

   //h5 = 3(f/g)^2 1/g h1 -2/g^2 h3+2f/g^3 h4 -f/g^2 h2
   StiffnessMatrix h5(dim*num_nodes,dim*num_nodes);
   h5 =1./gx*h1 - 1./pow(gx,2)*h3 -1./pow(gx,2)*h6 + 2*fx/pow(gx,3)*h4-fx/pow(gx,2)*h2;
   std::vector<Entry> h5e = entriesFromSparseMatrix(h5.block(0, 0, num_nodes * dim , num_nodes * dim));
   out.insert(out.end(), h5e.begin(), h5e.end());
}

void IMLSContact::addFastRIMLSSCTestEnergy(T& energy)
{
   BuildAcceleratorTreeSCTest();
   current_indices.resize(num_nodes);
   dist_info.clear();

   std::vector<std::vector<double>> tbb_current_indices(num_nodes);
   tbb::parallel_for(
   tbb::blocked_range<size_t>(size_t(0), num_nodes),
   [&](const tbb::blocked_range<size_t>& r) {
       for (size_t i = r.begin(); i < r.end(); i++) {
           auto& local_current_index = tbb_current_indices[i];
           local_current_index.clear();
           //std::cout<<i<<std::endl;
           //if(is_surface_vertex[i] != 0)
           if(i>=4225)
           {
               Vector3a xi = deformed.segment<3>(3*i);

               Fuzzy_sphere fs(Point_d(xi[0], xi[1], xi[2]), radius+radius*0.2, radius*0.2);

               std::vector<std::tuple<Point_d, int>> points_query;
               accelerator.search(std::back_inserter(points_query), fs);

               for(int k=0; k<points_query.size(); ++k)
               {
                   int index = std::get<1>(points_query[k]);
                   //if((geodist_close_matrix.coeff(i,index) == 0))
                       local_current_index.push_back(index);
               }
           }
       }
   });

   std::vector<double> temp_energies(num_nodes,0);
   tbb::enumerable_thread_specific<double> tbb_temp_energies(0);
   std::vector<std::pair<int,double>> tbb_dist_info(num_nodes);

   tbb::parallel_for(
   tbb::blocked_range<size_t>(size_t(0), num_nodes),
   [&](const tbb::blocked_range<size_t>& r) {
       auto& local_energy_temp = tbb_temp_energies.local();
       for (size_t i = r.begin(); i < r.end(); i++) {
           auto& local_current_index = tbb_current_indices[i];
           auto& local_dist_info = tbb_dist_info[i];
           local_dist_info = (std::pair<int,double>(i,1e8));
           //if(is_surface_vertex[i] != 0)
           if(i>=4225)
           {
               AScalar fx = 0;
               AScalar gx = 0;
               VectorXa N(3); N.setZero();
               VectorXa GF(3); GF.setZero();
               VectorXa GW(3); GW.setZero();
               AScalar fx2 = 0;
               AScalar gx2 = 0;

               Vector3a xi = deformed.segment<3>(3*i);

               for(int k=0; k<local_current_index.size(); ++k)
               {
                   Vector3a ck = deformed.segment<3>(local_current_index[k]*3);
                   int valence = nodes_adjacency[local_current_index[k]].size();
                   if(valence > 11) std::cout<<"WTF!!!!"<<std::endl;
                   VectorXa vs;

                   vs.resize(33);
                   for(int a=0; a<11; ++a)
                   {
                       vs.segment<3>(a*3) = ck;
                   }


                   for(int a=0; a<nodes_adjacency[local_current_index[k]].size(); ++a)
                   {
                       vs.segment<3>(a*3) = deformed.segment<3>(nodes_adjacency[local_current_index[k]][a]*3);
                   }

                   if((ck-xi).norm()<=radius)
                   {
                       fx += fx_func11(xi, vs, ck, radius);
                       gx += gx_func11(xi, ck, radius);
                       GF(0) += sumGF1_func(xi, vs, ck, radius); GF(1) += sumGF2_func(xi, vs, ck, radius); GF(2) += sumGF3_func(xi, vs, ck, radius);
                       GW(0) += sumGW1_func(xi, vs, ck, radius); GW(1) += sumGW2_func(xi, vs, ck, radius); GW(2) += sumGW3_func(xi, vs, ck, radius);
                       N(0) += sumN1_func(xi, vs, ck, radius); N(1) += sumN2_func(xi, vs, ck, radius); N(2) += sumN3_func(xi, vs, ck, radius);
                   }
               }

               if(abs(gx)<1e-9)
                   continue;

               AScalar f = fx/gx;
               VectorXa grad_f = (GF-f*GW+N)/gx;

               for(int k=0; k<local_current_index.size(); ++k)
               {
                   Vector3a ck = deformed.segment<3>(local_current_index[k]*3);
                   int valence = nodes_adjacency[local_current_index[k]].size();
                   VectorXa vs;

                   vs.resize(33);
                   for(int a=0; a<11; ++a)
                   {
                       vs.segment<3>(a*3) = ck;
                   }


                   for(int a=0; a<nodes_adjacency[local_current_index[k]].size(); ++a)
                   {
                       vs.segment<3>(a*3) = deformed.segment<3>(nodes_adjacency[local_current_index[k]][a]*3);
                   }

                   if((ck-xi).norm()<=radius)
                   {
                       fx2 += fx_func_2(xi, vs, ck, radius,f,grad_f,sigma_r,sigma_n);
                       gx2 += gx_func_2(xi, vs, ck, radius,f,grad_f,sigma_r,sigma_n);
                   }
               }

               AScalar dist = fx2/gx2;

               if(abs(gx2)>1e-9)
               {
                   //std::cout<<i<<" "<<fx<<" "<<gx<<" "<<dist<<std::endl;
                   if(BARRIER_ENERGY)
                   {
                       if(dist > 0)
                       {
                           if(dist <= barrier_distance)
                               local_energy_temp += -IMLS_param*pow((dist-barrier_distance),2)*log(dist/barrier_distance);
                       }
                       // else
                       //     std::cout<<"ERROR!!! Negative Distance!"<<std::endl;
                   }
                   else if(!BARRIER_ENERGY && dist<0)
                       local_energy_temp += -IMLS_param*pow((fx/gx),3);
                   local_dist_info = (std::pair<int,double>(i,fx/gx));
                   // close_slave_nodes.push_back(it->second);
               }

           }
       }
   });

   for(auto& local_energy_temp: tbb_temp_energies)
   {
       energy += local_energy_temp;
   }
   int i=0;
   for(auto& local_dist_info: tbb_dist_info)
   {
       // if(local_dist_info.first > 0)
       if(i>=4225)
           dist_info.push_back(local_dist_info);
       i++;
   }
}

void IMLSContact::addFastRIMLSSCEnergy(T& energy)
{
   BuildAcceleratorTreeSC();
   current_indices.resize(num_nodes);
   dist_info.clear();

   std::vector<std::vector<double>> tbb_current_indices(num_nodes);
   tbb::parallel_for(
   tbb::blocked_range<size_t>(size_t(0), num_nodes),
   [&](const tbb::blocked_range<size_t>& r) {
       for (size_t i = r.begin(); i < r.end(); i++) {
           auto& local_current_index = tbb_current_indices[i];
           local_current_index.clear();
           //std::cout<<i<<std::endl;
           if(is_surface_vertex[i] != 0)
           {
               Vector3a xi = deformed.segment<3>(3*i);

               Fuzzy_sphere fs(Point_d(xi[0], xi[1], xi[2]), radius+radius*0.2, radius*0.2);

               std::vector<std::tuple<Point_d, int>> points_query;
               accelerator.search(std::back_inserter(points_query), fs);

               for(int k=0; k<points_query.size(); ++k)
               {
                   int index = std::get<1>(points_query[k]);
                   if((geodist_close_matrix.coeff(i,index) == 0))
                       local_current_index.push_back(index);
               }
           }
       }
   });

   std::vector<double> temp_energies(num_nodes,0);
   tbb::enumerable_thread_specific<double> tbb_temp_energies(0);
   //tbb::enumerable_thread_specific<std::vector<std::pair<int,double>>> tbb_dist_info(0);
   std::vector<std::pair<int,double>> tbb_dist_info(num_nodes, std::pair<int,double>(-1,0));

   tbb::parallel_for(
   tbb::blocked_range<size_t>(size_t(0), num_nodes),
   [&](const tbb::blocked_range<size_t>& r) {
       auto& local_energy_temp = tbb_temp_energies.local();
       for (size_t i = r.begin(); i < r.end(); i++) {
           auto& local_current_index = tbb_current_indices[i];
           auto& local_dist_info = tbb_dist_info[i];
           //local_dist_info = (std::pair<int,double>(i,1e8));
           if(is_surface_vertex[i] != 0)
           {
            // std::cout<<i<<std::endl;
               AScalar fx = 0;
               AScalar gx = 0;
               VectorXa N(3); N.setZero();
               VectorXa GF(3); GF.setZero();
               VectorXa GW(3); GW.setZero();
               AScalar fx2 = 0;
               AScalar gx2 = 0;

               Vector3a xi = deformed.segment<3>(3*i);

               for(int k=0; k<local_current_index.size(); ++k)
               {
                   Vector3a ck = deformed.segment<3>(local_current_index[k]*3);
                   int valence = nodes_adjacency[local_current_index[k]].size();
                   if(valence > 11) std::cout<<"WTF!!!!"<<std::endl;
                   VectorXa vs;

                   vs.resize(33);
                   for(int a=0; a<11; ++a)
                   {
                       vs.segment<3>(a*3) = ck;
                   }


                   for(int a=0; a<nodes_adjacency[local_current_index[k]].size(); ++a)
                   {
                       vs.segment<3>(a*3) = deformed.segment<3>(nodes_adjacency[local_current_index[k]][a]*3);
                   }

                   if((ck-xi).norm()<=radius)
                   {
                       fx += fx_func11(xi, vs, ck, radius);
                       gx += gx_func11(xi, ck, radius);
                       GF(0) += sumGF1_func(xi, vs, ck, radius); GF(1) += sumGF2_func(xi, vs, ck, radius); GF(2) += sumGF3_func(xi, vs, ck, radius);
                       GW(0) += sumGW1_func(xi, vs, ck, radius); GW(1) += sumGW2_func(xi, vs, ck, radius); GW(2) += sumGW3_func(xi, vs, ck, radius);
                       N(0) += sumN1_func(xi, vs, ck, radius); N(1) += sumN2_func(xi, vs, ck, radius); N(2) += sumN3_func(xi, vs, ck, radius);
                   }
               }

               if(abs(gx)<1e-9)
                   continue;

               AScalar f = fx/gx;
               VectorXa grad_f = (GF-f*GW+N)/gx;

               for(int k=0; k<local_current_index.size(); ++k)
               {
                   Vector3a ck = deformed.segment<3>(local_current_index[k]*3);
                   int valence = nodes_adjacency[local_current_index[k]].size();
                   VectorXa vs;

                   vs.resize(33);
                   for(int a=0; a<11; ++a)
                   {
                       vs.segment<3>(a*3) = ck;
                   }


                   for(int a=0; a<nodes_adjacency[local_current_index[k]].size(); ++a)
                   {
                       vs.segment<3>(a*3) = deformed.segment<3>(nodes_adjacency[local_current_index[k]][a]*3);
                   }

                   if((ck-xi).norm()<=radius)
                   {
                       fx2 += fx_func_2(xi, vs, ck, radius,f,grad_f,sigma_r,sigma_n);
                       gx2 += gx_func_2(xi, vs, ck, radius,f,grad_f,sigma_r,sigma_n);
                   }
               }

               AScalar dist = pow(fx2/gx2,2);

               if(abs(gx2)>1e-9)
               {
                   //std::cout<<i<<" "<<fx<<" "<<gx<<" "<<dist<<std::endl;
                   if(BARRIER_ENERGY)
                   {
                       if(dist > 0)
                       {
                           if(dist <= barrier_distance)
                               local_energy_temp += -IMLS_param*pow((dist-barrier_distance),2)*log(dist/barrier_distance);
                       }
                       // else
                       //     std::cout<<"ERROR!!! Negative Distance!"<<std::endl;
                   }
                   else if(!BARRIER_ENERGY && dist<0)
                       local_energy_temp += -IMLS_param*pow((fx/gx),3);
                   //local_energy_temp += dist;
                   //local_energy_temp += grad_f(2);
                   local_dist_info = (std::pair<int,double>(i,fx2/gx2));
                   // close_slave_nodes.push_back(it->second);
               }

           }
       }
   });

   for(auto& local_energy_temp: tbb_temp_energies)
   {
       energy += local_energy_temp;
   }
   for(auto& local_dist_info: tbb_dist_info)
   {
       if(local_dist_info.first > 0)
       {
            // std::cout<<local_dist_info.first<<" "<<local_dist_info.second<<std::endl;
            dist_info.push_back(local_dist_info);
       }
           
        
   }
}

void IMLSContact::addFastRIMLSSCForceEntries(VectorXT& residual)
{
   BuildAcceleratorTreeSC();
    current_indices.resize(num_nodes);
   dist_info.clear();

   std::vector<std::vector<double>> tbb_current_indices(num_nodes);
   tbb::parallel_for(
   tbb::blocked_range<size_t>(size_t(0), num_nodes),
   [&](const tbb::blocked_range<size_t>& r) {
       for (size_t i = r.begin(); i < r.end(); i++) {
           auto& local_current_index = tbb_current_indices[i];
           local_current_index.clear();
           //std::cout<<i<<std::endl;
           if(is_surface_vertex[i] != 0)
           {
               Vector3a xi = deformed.segment<3>(3*i);

               Fuzzy_sphere fs(Point_d(xi[0], xi[1], xi[2]), radius+radius*0.2, radius*0.2);

               std::vector<std::tuple<Point_d, int>> points_query;
               accelerator.search(std::back_inserter(points_query), fs);

               for(int k=0; k<points_query.size(); ++k)
               {
                   int index = std::get<1>(points_query[k]);
                   if((geodist_close_matrix.coeff(i,index) == 0))
                       local_current_index.push_back(index);
               }
           }
       }
   });

   std::vector<std::vector<std::pair<int,double>>> vector_temps(num_nodes);
   tbb::enumerable_thread_specific<std::vector<std::pair<int,double>>> tbb_vector_temps;

   tbb::parallel_for(
   tbb::blocked_range<size_t>(size_t(0), num_nodes),
   [&](const tbb::blocked_range<size_t>& r) {
       auto& local_vector_temp = tbb_vector_temps.local();
       for (size_t i = r.begin(); i < r.end(); i++) {
           auto& local_current_index = tbb_current_indices[i];
           if(is_surface_vertex[i] != 0)
           {
               VectorXa gradient(num_nodes*dim);
               gradient.setZero();
               VectorXa sum_dfdx(dim*num_nodes);
               sum_dfdx.setZero();
               VectorXa sum_dgdx(dim*num_nodes);
               sum_dgdx.setZero();
               VectorXa sum_dGF1dx(dim*num_nodes);
               sum_dGF1dx.setZero();
               VectorXa sum_dGF2dx(dim*num_nodes);
               sum_dGF2dx.setZero();
               VectorXa sum_dGF3dx(dim*num_nodes);
               sum_dGF3dx.setZero();
               VectorXa sum_dGW1dx(dim*num_nodes);
               sum_dGW1dx.setZero();
               VectorXa sum_dGW2dx(dim*num_nodes);
               sum_dGW2dx.setZero();
               VectorXa sum_dGW3dx(dim*num_nodes);
               sum_dGW3dx.setZero();
               VectorXa sum_dN1dx(dim*num_nodes);
               sum_dN1dx.setZero();
               VectorXa sum_dN2dx(dim*num_nodes);
               sum_dN2dx.setZero();
               VectorXa sum_dN3dx(dim*num_nodes);
               sum_dN3dx.setZero();
               VectorXa sum_df2dx(dim*num_nodes);
               sum_df2dx.setZero();
               VectorXa sum_dg2dx(dim*num_nodes);
               sum_dg2dx.setZero();

               //std::cout<<i<<" out of "<<slave_nodes_3d[0].size()<<std::endl;
               AScalar fx = 0;
               AScalar gx = 0;
               VectorXa N(3); N.setZero();
               VectorXa GF(3); GF.setZero();
               VectorXa GW(3); GW.setZero();
               AScalar fx2 = 0;
               AScalar gx2 = 0;
               Vector3a xi = deformed.segment<3>(3*i);

               VectorXa ele_dfdx;
               VectorXa ele_dgdx;
               VectorXa ele_df2dx;
               VectorXa ele_dg2dx;
               VectorXa ele_dGF1dx,ele_dGF2dx,ele_dGF3dx;
               VectorXa ele_dGW1dx,ele_dGW2dx,ele_dGW3dx;
               VectorXa ele_dN1dx,ele_dN2dx,ele_dN3dx;
               //std::cout<<i<<" out of "<<slave_nodes_3d[0].size()<<std::endl;

               for(int k=0; k<local_current_index.size(); ++k)
               {
                   Vector3a ck = deformed.segment<3>(local_current_index[k]*3);
                   //f(i == 1567) std::cout<<k<<" k out of "<<local_current_index.size()<<std::endl;
                   int valence = nodes_adjacency[local_current_index[k]].size();
                   Eigen::VectorXi valence_indices;
                   VectorXa vs;

                   ele_dfdx.resize(39);
                   ele_dgdx.resize(39);
                   vs.resize(33);
                   valence_indices.resize(13);
                   for(int a=0; a<11; ++a)
                   {
                       vs.segment<3>(a*3) = ck;
                       valence_indices(a+2)=(local_current_index[k]);
                   }


                   valence_indices(0)=(i);
                   valence_indices(1)=(local_current_index[k]);


                   for(int a=0; a<nodes_adjacency[local_current_index[k]].size(); ++a)
                   {
                       vs.segment<3>(a*3) = deformed.segment<3>(nodes_adjacency[local_current_index[k]][a]*3);
                       valence_indices(a+2) = nodes_adjacency[local_current_index[k]][a];
                   }

                   if((ck-xi).norm()<=radius) {
                       ele_dfdx = dfdx_func11(xi, vs, ck, radius);
                       ele_dgdx = dgdx_func11(xi, ck, radius);
                       ele_dGF1dx = dsumGF1dx_func(xi, vs, ck, radius); ele_dGF2dx = dsumGF2dx_func(xi, vs, ck, radius); ele_dGF3dx = dsumGF3dx_func(xi, vs, ck, radius);
                       ele_dGW1dx = dsumGW1dx_func(xi, vs, ck, radius); ele_dGW2dx = dsumGW2dx_func(xi, vs, ck, radius); ele_dGW3dx = dsumGW3dx_func(xi, vs, ck, radius);
                       ele_dN1dx = dsumN1dx_func(xi, vs, ck, radius); ele_dN2dx = dsumN2dx_func(xi, vs, ck, radius); ele_dN3dx = dsumN3dx_func(xi, vs, ck, radius);
                       fx += fx_func11(xi, vs, ck, radius);
                       gx += gx_func11(xi, ck, radius);
                       GF(0) += sumGF1_func(xi, vs, ck, radius); GF(1) += sumGF2_func(xi, vs, ck, radius); GF(2) += sumGF3_func(xi, vs, ck, radius);
                       GW(0) += sumGW1_func(xi, vs, ck, radius); GW(1) += sumGW2_func(xi, vs, ck, radius); GW(2) += sumGW3_func(xi, vs, ck, radius);
                       N(0) += sumN1_func(xi, vs, ck, radius); N(1) += sumN2_func(xi, vs, ck, radius); N(2) += sumN3_func(xi, vs, ck, radius);
                       IMLS_local_gradient_to_global_gradient(ele_dfdx,valence_indices,dim,sum_dfdx);
                       IMLS_local_gradient_to_global_gradient(ele_dgdx,valence_indices,dim,sum_dgdx);
                       IMLS_local_gradient_to_global_gradient(ele_dGF1dx,valence_indices,dim,sum_dGF1dx);
                       IMLS_local_gradient_to_global_gradient(ele_dGF2dx,valence_indices,dim,sum_dGF2dx);
                       IMLS_local_gradient_to_global_gradient(ele_dGF3dx,valence_indices,dim,sum_dGF3dx);
                       IMLS_local_gradient_to_global_gradient(ele_dGW1dx,valence_indices,dim,sum_dGW1dx);
                       IMLS_local_gradient_to_global_gradient(ele_dGW2dx,valence_indices,dim,sum_dGW2dx);
                       IMLS_local_gradient_to_global_gradient(ele_dGW3dx,valence_indices,dim,sum_dGW3dx);
                       IMLS_local_gradient_to_global_gradient(ele_dN1dx,valence_indices,dim,sum_dN1dx);
                       IMLS_local_gradient_to_global_gradient(ele_dN2dx,valence_indices,dim,sum_dN2dx);
                       IMLS_local_gradient_to_global_gradient(ele_dN3dx,valence_indices,dim,sum_dN3dx);

                   }
               }
               if(abs(gx)<1e-9)
                   continue;

               AScalar f = fx/gx;
               VectorXa grad_f = (GF-f*GW+N)/gx;

               VectorXa dfdx = (gx*sum_dfdx-fx*sum_dgdx)/(pow(gx,2));
               // A = GF/gx
               VectorXa dA1dx = (gx*sum_dGF1dx-GF(0)*sum_dgdx)/(pow(gx,2)); VectorXa dA2dx = (gx*sum_dGF2dx-GF(1)*sum_dgdx)/(pow(gx,2)); VectorXa dA3dx = (gx*sum_dGF3dx-GF(2)*sum_dgdx)/(pow(gx,2));
               // B = GW/gx
               VectorXa dB1dx = (gx*sum_dGW1dx-GW(0)*sum_dgdx)/(pow(gx,2)); VectorXa dB2dx = (gx*sum_dGW2dx-GW(1)*sum_dgdx)/(pow(gx,2)); VectorXa dB3dx = (gx*sum_dGW3dx-GW(2)*sum_dgdx)/(pow(gx,2));
               // C = N/gx
               VectorXa dC1dx = (gx*sum_dN1dx-N(0)*sum_dgdx)/(pow(gx,2)); VectorXa dC2dx = (gx*sum_dN2dx-N(1)*sum_dgdx)/(pow(gx,2)); VectorXa dC3dx = (gx*sum_dN3dx-N(2)*sum_dgdx)/(pow(gx,2));

               VectorXa dgrad_f1dx = dA1dx - dfdx*(GW(0)/gx) - f*dB1dx + dC1dx;
               VectorXa dgrad_f2dx = dA2dx - dfdx*(GW(1)/gx) - f*dB2dx + dC2dx;
               VectorXa dgrad_f3dx = dA3dx - dfdx*(GW(2)/gx) - f*dB3dx + dC3dx;

               for(int k=0; k<local_current_index.size(); ++k)
               {
                   Vector3a ck = deformed.segment<3>(local_current_index[k]*3);
                   int valence = nodes_adjacency[local_current_index[k]].size();
                   Eigen::VectorXi valence_indices;
                   VectorXa vs;

                   ele_df2dx.resize(43); ele_dg2dx.resize(43);
                   vs.resize(33);
                   valence_indices.resize(13);
                   for(int a=0; a<11; ++a)
                   {
                       vs.segment<3>(a*3) = ck;
                       valence_indices(a+2)=(local_current_index[k]);
                   }


                   valence_indices(0)=(i);
                   valence_indices(1)=(local_current_index[k]);


                   for(int a=0; a<nodes_adjacency[local_current_index[k]].size(); ++a)
                   {
                       vs.segment<3>(a*3) = deformed.segment<3>(nodes_adjacency[local_current_index[k]][a]*3);
                       valence_indices(a+2) = nodes_adjacency[local_current_index[k]][a];
                   }

                   AScalar df2df = 0; AScalar dg2df = 0;
                   Vector3a df2dgrad_f; Vector3a dg2dgrad_f;
                   if((ck-xi).norm()<=radius)
                   {
                       // df2dx contains df2/df df2/dgrad_f1 ...
                       VectorXa temp1 = dfdx_func_2(xi, vs, ck, radius,f,grad_f,sigma_r,sigma_n);
                       VectorXa temp2 = dgdx_func_2(xi, vs, ck, radius,f,grad_f,sigma_r,sigma_n);
                       df2df = temp1(39);  dg2df = temp2(39);
                       df2dgrad_f = temp1.segment<3>(40); dg2dgrad_f = temp2.segment<3>(40);
                       ele_df2dx = temp1.segment<39>(0); ele_dg2dx = temp2.segment<39>(0);
                       fx2 += fx_func_2(xi, vs, ck, radius,f,grad_f,sigma_r,sigma_n);
                       gx2 += gx_func_2(xi, vs, ck, radius,f,grad_f,sigma_r,sigma_n);
                       IMLS_local_gradient_to_global_gradient(ele_df2dx,valence_indices,dim,sum_df2dx);
                       IMLS_local_gradient_to_global_gradient(ele_dg2dx,valence_indices,dim,sum_dg2dx);
                       sum_df2dx += (df2df * dfdx + df2dgrad_f(0) * dgrad_f1dx + df2dgrad_f(1) * dgrad_f2dx + df2dgrad_f(2) * dgrad_f3dx);
                       sum_dg2dx += (dg2df * dfdx + dg2dgrad_f(0) * dgrad_f1dx + dg2dgrad_f(1) * dgrad_f2dx + dg2dgrad_f(2) * dgrad_f3dx);
                   }
               }

               AScalar dist = pow(fx2/gx2,2);

               if(abs(gx2) > 1e-9)
               {
                   gradient.setZero();
                   if(BARRIER_ENERGY)
                   {
                       if(dist > 0)
                       {
                           if(dist<=barrier_distance)
                               gradient = 2*(fx2/gx2)*IMLS_param*((barrier_distance-dist)*(2*log(dist/barrier_distance)-barrier_distance/dist+1))*(gx2*sum_df2dx-fx2*sum_dg2dx)/(pow(gx2,2));
                       }
                       else
                           std::cout<<"ERROR!!! Negative Distance!"<<std::endl;
                   }
                   else if(!BARRIER_ENERGY && dist<0)
                       gradient = -3*IMLS_param*pow(fx/gx,2)*(gx*sum_dfdx-fx*sum_dgdx)/(pow(gx,2));


                //    gradient = (gx2*sum_df2dx-fx2*sum_dg2dx)/(pow(gx2,2));
                   for(int j=0; j<dim*num_nodes; ++j)
                   {
                       if(gradient(j) != 0) local_vector_temp.push_back(std::pair<int, double>(j,gradient(j)));
                   }
               }
           }
       }
   });

   for (const auto& local_vector_temp : tbb_vector_temps)
   {
       for(int j=0; j<local_vector_temp.size(); ++j)
       {
           residual(local_vector_temp[j].first) -= local_vector_temp[j].second;
       }
   }
}

void IMLSContact::addFastRIMLSSCHessianEntries(std::vector<Entry>& entries)
{
   BuildAcceleratorTreeSC();
   current_indices.resize(num_nodes);
   close_slave_nodes.clear();

   std::vector<std::vector<double>> tbb_current_indices(num_nodes);
   tbb::parallel_for(
   tbb::blocked_range<size_t>(size_t(0), num_nodes),
   [&](const tbb::blocked_range<size_t>& r) {
       for (size_t i = r.begin(); i < r.end(); i++) {
           auto& local_current_index = tbb_current_indices[i];
           local_current_index.clear();
           //std::cout<<i<<std::endl;
           if(is_surface_vertex[i] != 0)
           {
               Vector3a xi = deformed.segment<3>(3*i);

               Fuzzy_sphere fs(Point_d(xi[0], xi[1], xi[2]), radius+radius*0.2, radius*0.2);

               std::vector<std::tuple<Point_d, int>> points_query;
               accelerator.search(std::back_inserter(points_query), fs);

               for(int k=0; k<points_query.size(); ++k)
               {
                   int index = std::get<1>(points_query[k]);
                   if((geodist_close_matrix.coeff(i,index) == 0))
                       local_current_index.push_back(index);
               }
           }
       }
   });


   std::vector<std::vector<Entry>> temp_triplets(num_nodes);
   tbb::enumerable_thread_specific<std::vector<Eigen::Triplet<double>>> tbb_temp_triplets;
   for (auto& local_temp_triplets : tbb_temp_triplets)
   {
       local_temp_triplets.reserve(100000);
   }

//    tbb::task_scheduler_init init(4);

   tbb::parallel_for(
   tbb::blocked_range<size_t>(size_t(0), num_nodes),
   [&](const tbb::blocked_range<size_t>& r) {
       auto& local_temp_triplets = tbb_temp_triplets.local();
       for (size_t i = r.begin(); i < r.end(); i++) {
           //for (int i = 0; i < num_nodes; i++) {
           auto& local_current_index = tbb_current_indices[i];
           if(is_surface_vertex[i] != 0)
           {
               VectorXa sum_dfdx(dim*num_nodes);
               sum_dfdx.setZero();
               VectorXa sum_dgdx(dim*num_nodes);
               sum_dgdx.setZero();
               std::vector<Entry> sum_d2fdx2;
               std::vector<Entry> sum_d2gdx2;

               //std::cout<<i<<" out of "<<slave_nodes_3d[0].size()<<std::endl;
               AScalar fx = 0;
               AScalar gx = 0;
               Vector3a GW; GW.setZero();
               Vector3a GF; GF.setZero();
               Vector3a N; N.setZero();

               std::vector<VectorXa> sum_dGWdx(3,VectorXa(dim*num_nodes));
               std::vector<VectorXa> sum_dGFdx(3,VectorXa(dim*num_nodes));
               std::vector<VectorXa> sum_dNdx(3,VectorXa(dim*num_nodes));
               for(int l=0; l<3; l++)
               {
                   sum_dGWdx[l].setZero();
                   sum_dGFdx[l].setZero();
                   sum_dNdx[l].setZero();
               }
               std::vector<std::vector<Entry>> sum_d2GWdx2(3,std::vector<Entry>());
               std::vector<std::vector<Entry>> sum_d2GFdx2(3,std::vector<Entry>());
               std::vector<std::vector<Entry>> sum_d2Ndx2(3,std::vector<Entry>());


               Vector3a xi = deformed.segment<3>(3*i);

               VectorXa ele_dfdx;
               VectorXa ele_dgdx;
               VectorXa ele_d2fdx;
               VectorXa ele_d2gdx;
               std::vector<VectorXa> ele_dgwdx(3); std::vector<VectorXa> ele_dgfdx(3); std::vector<VectorXa> ele_dndx(3);
               std::vector<VectorXa> ele_d2gwdx2(3); std::vector<VectorXa> ele_d2gfdx2(3); std::vector<VectorXa> ele_d2ndx2(3);
               //std::cout<<i<<" out of "<<slave_nodes_3d[0].size()<<std::endl;

               Eigen::MatrixXd valence_indices;
               valence_indices.resize(local_current_index.size(),8);

               for(int k=0; k<local_current_index.size(); ++k)
               {
                   Vector3a ck = deformed.segment<3>(local_current_index[k]*3);
                   //f(i == 1567) std::cout<<k<<" k out of "<<local_current_index.size()<<std::endl;
                   int valence = nodes_adjacency[local_current_index[k]].size();
                   Eigen::VectorXi valence_indices;
                   VectorXa vs;

                   ele_dfdx.resize(39);
                   ele_dgdx.resize(39);
                   ele_d2fdx.resize(39*39);
                   ele_d2gdx.resize(39*39);
                   vs.resize(33);
                   valence_indices.resize(13);
                   for(int a=0; a<11; ++a)
                   {
                       vs.segment<3>(a*3) = ck;
                       valence_indices(a+2)=(local_current_index[k]);
                   }


                   valence_indices(0)=i;
                   valence_indices(1)=(local_current_index[k]);


                   for(int a=0; a<nodes_adjacency[local_current_index[k]].size(); ++a)
                   {
                       vs.segment<3>(a*3) = deformed.segment<3>(nodes_adjacency[local_current_index[k]][a]*3);
                       valence_indices(a+2) = nodes_adjacency[local_current_index[k]][a];
                   }


                   if((ck-xi).norm()<=radius) {

                       ele_d2fdx = (d2fdx_func11(xi, vs, ck, radius));
                       ele_d2gdx = (d2gdx_func11(xi, ck, radius));
                       ele_dfdx = dfdx_func11(xi, vs, ck, radius);
                       ele_dgdx = dgdx_func11(xi, ck, radius);
                       fx += fx_func11(xi, vs, ck, radius);
                       gx += gx_func11(xi, ck, radius);


                       IMLS_local_gradient_to_global_gradient(ele_dfdx,valence_indices,dim,sum_dfdx);
                       IMLS_local_gradient_to_global_gradient(ele_dgdx,valence_indices,dim,sum_dgdx);
                       IMLS_local_hessian_to_global_triplets(ele_d2fdx,valence_indices,dim,sum_d2fdx2);
                       IMLS_local_hessian_to_global_triplets(ele_d2gdx,valence_indices,dim,sum_d2gdx2);

                       GW(0) += sumGW1_func(xi, vs, ck, radius); ele_dgwdx[0] = dsumGW1dx_func(xi, vs, ck, radius); ele_d2gwdx2[0] = d2sumGW1dx2_func(xi, vs, ck, radius);
                       GW(1) += sumGW2_func(xi, vs, ck, radius); ele_dgwdx[1] = dsumGW2dx_func(xi, vs, ck, radius); ele_d2gwdx2[1] = d2sumGW2dx2_func(xi, vs, ck, radius);
                       GW(2) += sumGW3_func(xi, vs, ck, radius); ele_dgwdx[2] = dsumGW3dx_func(xi, vs, ck, radius); ele_d2gwdx2[2] = d2sumGW3dx2_func(xi, vs, ck, radius);
                       IMLS_local_gradient_to_global_gradient(ele_dgwdx[0],valence_indices,dim,sum_dGWdx[0]);
                       IMLS_local_hessian_to_global_triplets(ele_d2gwdx2[0],valence_indices,dim,sum_d2GWdx2[0]);
                       IMLS_local_gradient_to_global_gradient(ele_dgwdx[1],valence_indices,dim,sum_dGWdx[1]);
                       IMLS_local_hessian_to_global_triplets(ele_d2gwdx2[1],valence_indices,dim,sum_d2GWdx2[1]);
                       IMLS_local_gradient_to_global_gradient(ele_dgwdx[2],valence_indices,dim,sum_dGWdx[2]);
                       IMLS_local_hessian_to_global_triplets(ele_d2gwdx2[2],valence_indices,dim,sum_d2GWdx2[2]);

                       GF(0) += sumGF1_func(xi, vs, ck, radius); ele_dgfdx[0] = dsumGF1dx_func(xi, vs, ck, radius); ele_d2gfdx2[0] = d2sumGF1dx2_func(xi, vs, ck, radius);
                       GF(1) += sumGF2_func(xi, vs, ck, radius); ele_dgfdx[1] = dsumGF2dx_func(xi, vs, ck, radius); ele_d2gfdx2[1] = d2sumGF2dx2_func(xi, vs, ck, radius);
                       GF(2) += sumGF3_func(xi, vs, ck, radius); ele_dgfdx[2] = dsumGF3dx_func(xi, vs, ck, radius); ele_d2gfdx2[2] = d2sumGF3dx2_func(xi, vs, ck, radius);
                       IMLS_local_gradient_to_global_gradient(ele_dgfdx[0],valence_indices,dim,sum_dGFdx[0]);
                       IMLS_local_hessian_to_global_triplets(ele_d2gfdx2[0],valence_indices,dim,sum_d2GFdx2[0]);
                       IMLS_local_gradient_to_global_gradient(ele_dgfdx[1],valence_indices,dim,sum_dGFdx[1]);
                       IMLS_local_hessian_to_global_triplets(ele_d2gfdx2[1],valence_indices,dim,sum_d2GFdx2[1]);
                       IMLS_local_gradient_to_global_gradient(ele_dgfdx[2],valence_indices,dim,sum_dGFdx[2]);
                       IMLS_local_hessian_to_global_triplets(ele_d2gfdx2[2],valence_indices,dim,sum_d2GFdx2[2]);

                       N(0) += sumN1_func(xi, vs, ck, radius); ele_dndx[0] = dsumN1dx_func(xi, vs, ck, radius); ele_d2ndx2[0] = d2sumN1dx2_func(xi, vs, ck, radius);
                       N(1) += sumN2_func(xi, vs, ck, radius); ele_dndx[1] = dsumN2dx_func(xi, vs, ck, radius); ele_d2ndx2[1] = d2sumN2dx2_func(xi, vs, ck, radius);
                       N(2) += sumN3_func(xi, vs, ck, radius); ele_dndx[2] = dsumN3dx_func(xi, vs, ck, radius); ele_d2ndx2[2] = d2sumN3dx2_func(xi, vs, ck, radius);
                       IMLS_local_gradient_to_global_gradient(ele_dndx[0],valence_indices,dim,sum_dNdx[0]);
                       IMLS_local_hessian_to_global_triplets(ele_d2ndx2[0],valence_indices,dim,sum_d2Ndx2[0]);
                       IMLS_local_gradient_to_global_gradient(ele_dndx[1],valence_indices,dim,sum_dNdx[1]);
                       IMLS_local_hessian_to_global_triplets(ele_d2ndx2[1],valence_indices,dim,sum_d2Ndx2[1]);
                       IMLS_local_gradient_to_global_gradient(ele_dndx[2],valence_indices,dim,sum_dNdx[2]);
                       IMLS_local_hessian_to_global_triplets(ele_d2ndx2[2],valence_indices,dim,sum_d2Ndx2[2]);
                   }
               }


               if(abs(gx)<1e-9) continue;
               AScalar f = fx/gx;
               VectorXa grad_f = (GF-f*GW+N)/gx;

               VectorXa dfdx = divisionGradient(fx,gx,sum_dfdx,sum_dgdx);
               std::vector<VectorXa> dAdx(3),dBdx(3),dCdx(3), dgrad_fdx(3);

               // A = GF/gx; B = GW/gx; C = N/gx
               for(int l=0; l<3; l++)
               {
                   dAdx[l] = divisionGradient(GF(l),gx,sum_dGFdx[l],sum_dgdx);
                   dBdx[l] = divisionGradient(GW(l),gx,sum_dGWdx[l],sum_dgdx);
                   dCdx[l] = divisionGradient(N(l),gx,sum_dNdx[l],sum_dgdx);
                   dgrad_fdx[l]  = dAdx[l] - dfdx*(GW(l)/gx) - f*dBdx[l] + dCdx[l];
               }

               std::vector<Eigen::Triplet<double>> d2fdx2;
               std::vector<std::vector<Eigen::Triplet<double>>> d2Adx2(3,std::vector<Eigen::Triplet<double>>());
               std::vector<std::vector<Eigen::Triplet<double>>> d2Bdx2(3,std::vector<Eigen::Triplet<double>>());
               std::vector<std::vector<Eigen::Triplet<double>>> d2Cdx2(3,std::vector<Eigen::Triplet<double>>());
               std::vector<Eigen::SparseMatrix<double>> d2grad_fdx2(3,Eigen::SparseMatrix<double>(num_nodes*dim,num_nodes*dim));

               divisionHessian(fx,gx,sum_dfdx,sum_dgdx,sum_d2fdx2,sum_d2gdx2,d2fdx2);

               for(int l=0; l<3; l++)
               {
                   divisionHessian(GF(l),gx,sum_dGFdx[l],sum_dgdx,sum_d2GFdx2[l],sum_d2gdx2,d2Adx2[l]);
                   divisionHessian(GW(l),gx,sum_dGWdx[l],sum_dgdx,sum_d2GWdx2[l],sum_d2gdx2,d2Bdx2[l]);
                   divisionHessian(N(l),gx,sum_dNdx[l],sum_dgdx,sum_d2Ndx2[l],sum_d2gdx2,d2Cdx2[l]);

                   Eigen::SparseMatrix<double> h1(dim*num_nodes,dim*num_nodes), h2(dim*num_nodes,dim*num_nodes), h3(dim*num_nodes,dim*num_nodes), h4(dim*num_nodes,dim*num_nodes), h5(dim*num_nodes,dim*num_nodes), h6(dim*num_nodes,dim*num_nodes);
                   h1.setFromTriplets(d2Adx2[l].begin(),d2Adx2[l].end());
                   h2.setFromTriplets(d2Bdx2[l].begin(),d2Bdx2[l].end());
                   h3.setFromTriplets(d2Cdx2[l].begin(),d2Cdx2[l].end());
                   h4.setFromTriplets(d2fdx2.begin(),d2fdx2.end());

                   std::vector<Eigen::Triplet<double>> h5_t;
                   IMLS_vector_muliplication_to_triplets(dfdx,dBdx[l],1.,h5_t);
                   h5.setFromTriplets(h5_t.begin(), h5_t.end());

                   std::vector<Eigen::Triplet<double>> h6_t;
                   IMLS_vector_muliplication_to_triplets(dBdx[l],dfdx,1.,h6_t);
                   h6.setFromTriplets(h6_t.begin(), h6_t.end());

                   Eigen::SparseMatrix<double> h7(dim*num_nodes,dim*num_nodes);
                   h7 = h1 - h4*(GW(l)/gx) - h5 - h6 - f*h2 + h3;

                   d2grad_fdx2[l] = h7;
               }
               //int cnt = 0;

               // for(int j=0; j<dfdx.size(); ++j)
               // {
               //     if(dgrad_fdx[0](j) !=0 ) cnt++;
               // }
               // std::cout<<"Nonzeros: "<<cnt<<" out of "<<dfdx.size()<<std::endl;

               AScalar fx2 = 0;
               AScalar gx2 = 0;

               VectorXa ele_df2dx;
               VectorXa ele_dg2dx;
               MatrixXa ele_d2f2dx;
               MatrixXa ele_d2g2dx;
               MatrixXa ele_d2f2duvwz2;
               MatrixXa ele_d2g2duvwz2;
               std::vector<VectorXa> ele_d2f2dxuvwz(4);
               std::vector<VectorXa> ele_d2g2dxuvwz(4);

               VectorXa sum_df2dx(dim*num_nodes);
               sum_df2dx.setZero();
               VectorXa sum_dg2dx(dim*num_nodes);
               sum_dg2dx.setZero();
               std::vector<Entry> sum_d2f2dx2;
               std::vector<Entry> sum_d2g2dx2;

               std::vector<VectorXa> duvwzdx = {dfdx,dgrad_fdx[0],dgrad_fdx[1],dgrad_fdx[2]};
               std::vector<Eigen::SparseMatrix<double>> d2uvwzdx2(4,Eigen::SparseMatrix<double>(num_nodes*dim,num_nodes*dim));
               d2uvwzdx2[0].setFromTriplets(d2fdx2.begin(),d2fdx2.end());
               d2uvwzdx2[1] = d2grad_fdx2[0];
               d2uvwzdx2[2] = d2grad_fdx2[1];
               d2uvwzdx2[3] = d2grad_fdx2[2];
               Vector4a smf_coeff; smf_coeff.setZero();
               Vector4a smg_coeff; smg_coeff.setZero();

               MatrixXa vvtf_coeff(4,4); vvtf_coeff.setZero();
               MatrixXa vvtg_coeff(4,4); vvtg_coeff.setZero();
               std::vector<VectorXa> vvtf(4,VectorXa(num_nodes*dim));
               std::vector<VectorXa> vvtg(4,VectorXa(num_nodes*dim));

               for(int a=0; a<4; ++a)
               {
                   vvtf[a].setZero();
                   vvtg[a].setZero();
               }

               for(int k=0; k<local_current_index.size(); ++k)
               {
                   Vector3a ck = deformed.segment<3>(local_current_index[k]*3);
                   int valence = nodes_adjacency[local_current_index[k]].size();
                   Eigen::VectorXi valence_indices;
                   VectorXa vs;

                   vs.resize(33);
                   valence_indices.resize(13);
                   for(int a=0; a<11; ++a)
                   {
                       vs.segment<3>(a*3) = ck;
                       valence_indices(a+2)=(local_current_index[k]);
                   }


                   valence_indices(0)=i;
                   valence_indices(1)=(local_current_index[k]);


                   for(int a=0; a<nodes_adjacency[local_current_index[k]].size(); ++a)
                   {
                       vs.segment<3>(a*3) = deformed.segment<3>(nodes_adjacency[local_current_index[k]][a]*3);
                       valence_indices(a+2) = nodes_adjacency[local_current_index[k]][a];
                   }


                   if((ck-xi).norm()<=radius) {
                       Vector4a df2duvwz; Vector4a dg2duvwz;

                       fx2 += fx_func_2(xi, vs, ck, radius,f,grad_f,sigma_r,sigma_n);
                       gx2 += gx_func_2(xi, vs, ck, radius,f,grad_f,sigma_r,sigma_n);

                       VectorXa temp1 = dfdx_func_2(xi, vs, ck, radius,f,grad_f,sigma_r,sigma_n);
                       VectorXa temp2 = dgdx_func_2(xi, vs, ck, radius,f,grad_f,sigma_r,sigma_n);
                       df2duvwz = temp1.segment<4>(39); dg2duvwz = temp2.segment<4>(39);
                       ele_df2dx = temp1.segment<39>(0); ele_dg2dx = temp2.segment<39>(0);

                       IMLS_local_gradient_to_global_gradient(ele_df2dx,valence_indices,dim,sum_df2dx);
                       IMLS_local_gradient_to_global_gradient(ele_dg2dx,valence_indices,dim,sum_dg2dx);

                       for(int l=0; l<4; ++l)
                       {
                           sum_df2dx += (df2duvwz(l) * duvwzdx[l]);
                           sum_dg2dx += (dg2duvwz(l) * duvwzdx[l]);
                           smf_coeff(l) += df2duvwz(l);
                           smg_coeff(l) += dg2duvwz(l);
                       }
                    
                       MatrixXa temp3 = d2fdx2_func_2(xi, vs, ck, radius,f,grad_f,sigma_r,sigma_n);
                       MatrixXa temp4 = d2gdx2_func_2(xi, vs, ck, radius,f,grad_f,sigma_r,sigma_n);

                       ele_d2f2dx = temp3.block(0,0,39,39);
                       ele_d2g2dx = temp4.block(0,0,39,39);
                       for(int l=0; l<4; ++l)
                       {
                           ele_d2f2dxuvwz[l] = temp3.block(0,39+l,39,1);
                           ele_d2g2dxuvwz[l] = temp4.block(0,39+l,39,1);
                       }
                       ele_d2f2duvwz2 = temp3.block(39,39,4,4);
                       ele_d2g2duvwz2 = temp4.block(39,39,4,4);


                       IMLS_local_hessian_matrix_to_global_triplets(ele_d2f2dx,valence_indices,dim,sum_d2f2dx2);
                       IMLS_local_hessian_matrix_to_global_triplets(ele_d2g2dx,valence_indices,dim,sum_d2g2dx2);

                       for(int a=0; a<4; ++a)
                       {
                           VectorXa temp1(num_nodes*dim); temp1.setZero();
                           VectorXa temp2(num_nodes*dim); temp2.setZero();
                           IMLS_local_gradient_to_global_gradient(ele_d2f2dxuvwz[a],valence_indices,dim,temp1);
                           IMLS_local_gradient_to_global_gradient(ele_d2g2dxuvwz[a],valence_indices,dim,temp2);
                           vvtf[a] += temp1;
                           vvtg[a] += temp2;

                           for(int b=0; b<4; b++)
                           {
                               vvtf_coeff(a,b) += ele_d2f2duvwz2(a,b);
                               vvtg_coeff(a,b) += ele_d2g2duvwz2(a,b);
                           }
                       }
                   }
               }

                if(abs(gx2) < 1e-9) continue;

                AScalar dist = pow(fx2/gx2,2);

                if(BARRIER_ENERGY && dist <= 0)
                {
                    std::cout<<"NEGATIVE DISTANCE!!!"<<std::endl;
                    continue;
                } 
                if(BARRIER_ENERGY && dist > barrier_distance)
                {
                    continue;
                } 
                if(!BARRIER_ENERGY && dist > 0)
                {
                    continue;
                }


                Eigen::SparseMatrix<double> additional_smf(num_nodes*dim, num_nodes*dim); additional_smf.setZero();
                Eigen::SparseMatrix<double> additional_smg(num_nodes*dim, num_nodes*dim); additional_smg.setZero();
                for(int a=0; a<4; ++a)
                {
                    additional_smf += smf_coeff(a)*d2uvwzdx2[a];
                    additional_smg += smg_coeff(a)*d2uvwzdx2[a];

                    IMLS_vector_muliplication_to_triplets(vvtf[a],duvwzdx[a],1.,sum_d2f2dx2);
                    IMLS_vector_muliplication_to_triplets(duvwzdx[a],vvtf[a],1.,sum_d2f2dx2);

                    IMLS_vector_muliplication_to_triplets(vvtg[a],duvwzdx[a],1.,sum_d2g2dx2);
                    IMLS_vector_muliplication_to_triplets(duvwzdx[a],vvtg[a],1.,sum_d2g2dx2);

                    for(int b=0; b<4; b++)
                    {
                        IMLS_vector_muliplication_to_triplets(duvwzdx[a],duvwzdx[b],vvtf_coeff(a,b),sum_d2f2dx2);
                        IMLS_vector_muliplication_to_triplets(duvwzdx[a],duvwzdx[b],vvtg_coeff(a,b),sum_d2g2dx2);
                    }
                }

                Eigen::SparseMatrix<double> sum_d2f2dx2m(num_nodes*dim, num_nodes*dim); sum_d2f2dx2m.setZero();
                Eigen::SparseMatrix<double> sum_d2g2dx2m(num_nodes*dim, num_nodes*dim); sum_d2g2dx2m.setZero();
                sum_d2f2dx2m.setFromTriplets(sum_d2f2dx2.begin(), sum_d2f2dx2.end());
                sum_d2g2dx2m.setFromTriplets(sum_d2g2dx2.begin(), sum_d2g2dx2.end());
                sum_d2f2dx2m += additional_smf;
                sum_d2g2dx2m += additional_smg;

                sum_d2f2dx2.clear(); sum_d2f2dx2 = entriesFromSparseMatrix(sum_d2f2dx2m.block(0, 0, num_nodes * dim , num_nodes * dim));
                sum_d2g2dx2.clear(); sum_d2g2dx2 = entriesFromSparseMatrix(sum_d2g2dx2m.block(0, 0, num_nodes * dim , num_nodes * dim));

                VectorXa ddist_dx = divisionGradient(fx2,gx2,sum_df2dx,sum_dg2dx);
                VectorXa ddist_dx1 = ddist_dx*(2*fx2/gx2);

                std::vector<Eigen::Triplet<double>> d2f2dx2;

                divisionHessian(fx2,gx2,sum_df2dx,sum_dg2dx,sum_d2f2dx2,sum_d2g2dx2,d2f2dx2);
               

                std::vector<Eigen::Triplet<double>> Hessian_t;
                if(!BARRIER_ENERGY)
                    IMLS_vector_muliplication_to_triplets(ddist_dx,ddist_dx,-IMLS_param*6*(dist),Hessian_t);
                else
                {
                    IMLS_vector_muliplication_to_triplets(ddist_dx1,ddist_dx1,IMLS_param*((barrier_distance/dist+2)*(barrier_distance/dist)-2*log(dist/barrier_distance)-3),Hessian_t);
                    IMLS_vector_muliplication_to_triplets(ddist_dx,ddist_dx,2*((barrier_distance-dist)*(2*log(dist/barrier_distance)-barrier_distance/dist+1)),Hessian_t);
                }
                    
                std::vector<Entry> local;
                local.insert(local.end(), Hessian_t.begin(), Hessian_t.end());

                StiffnessMatrix h5(dim*num_nodes,dim*num_nodes);
                h5.setFromTriplets(d2f2dx2.begin(),d2f2dx2.end());
                if(!BARRIER_ENERGY)
                    h5 *= -IMLS_param*3*pow(dist,2);
                else
                    h5 *= IMLS_param*((barrier_distance-dist)*(2*log(dist/barrier_distance)-barrier_distance/dist+1))*(2*fx2/gx2);
                std::vector<Entry> h5e = entriesFromSparseMatrix(h5.block(0, 0, num_nodes * dim , num_nodes * dim));


                local.insert(local.end(),h5e.begin(),h5e.end());

                h5.setFromTriplets(local.begin(),local.end());
                std::vector<Entry> h5f = entriesFromSparseMatrix(h5.block(0, 0, num_nodes * dim , num_nodes * dim));
                local_temp_triplets.insert(local_temp_triplets.end(),h5f.begin(),h5f.end());
           }
       }
   });

   for (const auto& local_temp_triplets : tbb_temp_triplets)
   {
       entries.insert(entries.end(),local_temp_triplets.begin(),local_temp_triplets.end());
   }
}