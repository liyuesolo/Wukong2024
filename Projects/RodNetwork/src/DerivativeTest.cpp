#include "../include/RodNetwork.h"

void RodNetwork::testGradientFD()
{
    run_diff_test = true;
    std::cout
        << "======================== CHECK GRADIENT ========================"
        << std::endl;
    T epsilon = 1e-7;
    int n_dof = deformed_states.rows();

    VectorXT gradient(n_dof);
    gradient.setZero();
    VectorXT dx(n_dof);
    dx.setRandom();
    dx *= 1.0 / dx.norm();
    dx *= 0.001;
    dq += dx;
    VectorXT dq0 = dq;
    computeResidual(gradient);
    gradient *= -1;

    for (int i = 0; i < n_dof; i++)
    {
        dq[i] = dq0[i] + epsilon;
        T E0 = computeTotalEnergy();
        dq[i] = dq0[i] - 2.0 * epsilon;
        T E1 = computeTotalEnergy();
        dq[i] = dq0[i] + epsilon;
        T fd = (E0 - E1) / (2.0 * epsilon);
        if (std::abs(fd - gradient[i]) > 1e-3 * std::abs(gradient[i]))
        {
            std::cout << " dof " << i << " fd " << fd << " symbolic "
                      << gradient[i] << std::endl;
            std::getchar();
        }
    }

    run_diff_test = false;
}

void RodNetwork::testGradient2ndOrderTerm()
{
    run_diff_test = true;
    std::cout
        << "======================== CHECK GRADIENT ========================"
        << std::endl;
    T epsilon = 1e-6;
    int n_dof = deformed_states.rows();

    VectorXT gradient(n_dof);
    gradient.setZero();
    VectorXT dx(n_dof);
    dx.setRandom();
    dx *= 1.0 / dx.norm();
    dx *= 0.01;
    dq += dx;
    VectorXT dq0 = dq;
    T E0 = computeTotalEnergy();
    computeResidual(gradient);

    gradient *= -1;
    T previous = 0.0;
    for (int i = 0; i < 10; i++)
    {
        dq = dq0 + dx;
        T E1 = computeTotalEnergy();
        T dE = E1 - E0;

        dE -= gradient.dot(dx);
        // std::cout << "dE " << dE << std::endl;
        if (i > 0)
        {
            std::cout << (previous / dE) << std::endl;
        }
        previous = dE;
        dx *= 0.5;
    }
    run_diff_test = false;
}

void RodNetwork::testHessian2ndOrderTerm()
{
    run_diff_test = true;
    std::cout
        << "======================== CHECK Hessian ========================"
        << std::endl;
    T epsilon = 1e-6;
    int n_dof = deformed_states.rows();

    VectorXT f0(n_dof);
    f0.setZero();
    VectorXT dx(n_dof);
    dx.setRandom();
    dx *= 1.0 / dx.norm();
    dx *= 0.001;
    dq += dx;
    VectorXT dq0 = dq;

    StiffnessMatrix A;
    buildSystemDoFMatrix(A);
    computeResidual(f0);
    f0 *= -1;

    T previous = 0.0;
    for (int i = 0; i < 10; i++)
    {
        dq = dq0 + dx;
        VectorXT f1(n_dof);
        f1.setZero();
        computeResidual(f1);
        f1 *= -1;
        T df_norm = (f0 + (A * dx) - f1).norm();
        // std::cout << "df_norm " << df_norm << std::endl;
        if (i > 0)
        {
            std::cout << (previous / df_norm) << std::endl;
        }
        previous = df_norm;
        dx *= 0.5;
    }

    run_diff_test = false;
}