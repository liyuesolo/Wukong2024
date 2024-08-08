#include "../include/IMLSContact.h"

void IMLSContact::checkTotalGradient(bool perturb)
{   
    run_diff_test = true;
    VectorXT du(num_nodes * 3);
    du.setRandom();
    du *= 1.0 / du.norm();
    du *= 0.001;
    if (perturb)
        u += du;

    std::cout << "======================== CHECK GRADIENT ========================" << std::endl;
    int n_dof = num_nodes * 3;
    T epsilon = 1e-6;
    VectorXT gradient(n_dof);
    gradient.setZero();
    computeResidual(gradient);
    gradient *= -1.0;
    
    VectorXT gradient_FD(n_dof);
    gradient_FD.setZero();

    int cnt = 0;
    for(int dof_i = 0; dof_i < n_dof; dof_i++)
    {
        u(dof_i) += epsilon;
        // std::cout << W * dq << std::endl;
        T E0 = computeTotalEnergy();
        
        u(dof_i) -= 2.0 * epsilon;
        T E1 = computeTotalEnergy();
        u(dof_i) += epsilon;
        // std::cout << "E1 " << E1 << " E0 " << E0 << std::endl;
        gradient_FD(dof_i) = (E0 - E1) / (2.0 *epsilon);
        if( gradient_FD(dof_i) == 0 && gradient(dof_i) == 0)
            continue;
        if (std::abs( gradient_FD(dof_i) - gradient(dof_i)) < 1e-3 * std::abs(gradient(dof_i)))
            continue;
        std::cout << " dof " << dof_i << " " << gradient_FD(dof_i) << " " << gradient(dof_i) << std::endl;
        std::getchar();
        cnt++;   
    }
    
    run_diff_test = false;
}

void IMLSContact::checkTotalGradientScale(bool perturb)
{
    run_diff_test = true;
    add_cubic_plane = false;
    std::cout << "======================== CHECK GRADIENT Scale ========================" << std::endl;
    T epsilon = 1e-5;
    VectorXT du(num_nodes * 3);
    du.setRandom();
    du *= 1.0 / du.norm();
    du *= 0.001;
    if (perturb)
        u += du;
    
    int n_dof = num_nodes * 3;

    VectorXT gradient(n_dof);
    gradient.setZero();
    computeResidual(gradient);
    gradient *= -1.0;

    T E0 = computeTotalEnergy();
    VectorXT dx(n_dof);
    dx.setRandom();
    dx *= 1.0 / dx.norm();
    dx *= 0.001;
    T previous = 0.0;
    VectorXT u_current = u;
    for (int i = 0; i < 10; i++)
    {
        u = u_current + dx;
        T E1 = computeTotalEnergy();
        T dE = E1 - E0;
        
        dE -= gradient.dot(dx);
        // std::cout << "dE " << dE << std::endl;
        if (i > 0)
        {
            std::cout << (previous/dE) << std::endl;
        }
        previous = dE;
        dx *= 0.5;
    }
    run_diff_test = false;
}

void IMLSContact::checkTotalHessianScale(bool perturb)
{
    run_diff_test = true;
    std::cout << "===================== check Hessian 2nd Scale =====================" << std::endl;

    VectorXT du(num_nodes * 3);
    du.setRandom();
    du *= 1.0 / du.norm();
    du *= 0.001;
    if (perturb)
        u += du;
    
    int n_dof = num_nodes * 3;

    StiffnessMatrix A;
    
    buildSystemMatrix(A);

    std::cout << "build matrix" << std::endl;
    VectorXT f0(n_dof);
    f0.setZero();
    computeResidual(f0);
    f0 *= -1;
    
    VectorXT dx(n_dof);
    dx.setRandom();
    dx *= 1.0 / dx.norm();
    for(int i = 0; i < n_dof; i++) dx[i] += 0.5;
    dx *= 0.001;
    T previous = 0.0;
    VectorXT u_current = u;
    for (int i = 0; i < 10; i++)
    {
        VectorXT f1(n_dof);
        f1.setZero();
        u = u_current+dx;
        computeResidual(f1);
        f1 *= -1;
        T df_norm = (f0 + (A * dx) - f1).norm();
        // std::cout << "df_norm " << df_norm << std::endl;
        if (i > 0)
        {
            std::cout << (previous/df_norm) << std::endl;
        }
        previous = df_norm;
        dx *= 0.5;
        u = u_current;
    }
    run_diff_test = false;
}