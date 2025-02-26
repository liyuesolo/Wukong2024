#include "../include/Autodiff.h"

void Autodiff::generateDiscreteShellGravitionalEnergy()
{
    
    std::vector<A> x(9);
    Independent(x);

    // add_argument_expr(function, undeformedDensity);
	// add_argument_expr(function, thickness);
	// add_argument_vector(function, gravity, 3);
	// add_argument_vector(function, x1, 3);
	// add_argument_vector(function, x2, 3);
	// add_argument_vector(function, x3, 3);
	// add_argument_vector(function, x1Undef, 3);
	// add_argument_vector(function, x2Undef, 3);
	// add_argument_vector(function, x3Undef, 3);
}


   // redefine cpp ad scalar
    using A  = CppAD::AD<CppAD::cg::CG<T>>;
    using AV = Vector<A, 3>;
    using AV2 = Vector<A, 2>;
    using AM = Matrix<A, 3, 3>;
    using AM2 = Matrix<A, 2, 2>;

    bool use_CppAD = false;
	std::vector<ADCGWrapper::Autodiff::Function> autodiff_energies;

// VectorXT xp(9 + 11);
        // xp.segment<3>(0) = x0; xp.segment<3>(3) = x1;xp.segment<3>(6) = x2;
        // xp[9] = nu; xp[10] = k_s; 
        // xp.segment<3>(11) = X0; xp.segment<3>(14) = X1; xp.segment<3>(17) = X2;
        // energy += autodiff_energies[1](xp)[0];

// VectorXT xp(9 + 11);
        // xp.segment<3>(0) = x0; xp.segment<3>(3) = x1;xp.segment<3>(6) = x2;
        // xp[9] = nu; xp[10] = k_s; 
        // xp.segment<3>(11) = X0; xp.segment<3>(14) = X1; xp.segment<3>(17) = X2;
        // Vector<T, 9> dedx; dedx.setZero();
        // StiffnessMatrix jacobian_sparse = autodiff_energies[1].Jacobian(xp);
        // for (int k = 0; k < jacobian_sparse.outerSize(); ++k)
        //     for (StiffnessMatrix::InnerIterator it(jacobian_sparse, k); it; ++it)
        //         dedx(it.col()) = it.value();

// VectorXT xp(9 + 14);
        // xp.segment<3>(0) = x0; xp.segment<3>(3) = x1;xp.segment<3>(6) = x2;
        // xp[9] = density; xp[10] = thickness; 
        // xp.segment<3>(11) = gravity;
        // xp.segment<3>(14) = X0; xp.segment<3>(17) = X1; xp.segment<3>(20) = X2;
        // energy += autodiff_energies[0](xp)[0];

// VectorXT xp(9 + 14);
        // xp.segment<3>(0) = x0; xp.segment<3>(3) = x1;xp.segment<3>(6) = x2;
        // xp[9] = density; xp[10] = thickness; 
        // xp.segment<3>(11) = gravity;
        // xp.segment<3>(14) = X0; xp.segment<3>(17) = X1; xp.segment<3>(20) = X2;
        // Vector<T, 9> dedx = autodiff_energies[0].Jacobian(xp).toDense().reshaped(9, 1);

// VectorXT xp(9 + 14);
        // xp.segment<3>(0) = x0; xp.segment<3>(3) = x1;xp.segment<3>(6) = x2;
        // xp[9] = density; xp[10] = thickness; 
        // xp.segment<3>(11) = gravity;
        // xp.segment<3>(14) = X0; xp.segment<3>(17) = X1; xp.segment<3>(20) = X2;
        // Matrix<T, 9, 9> hessian = autodiff_energies[0].Hessian(xp).toDense().reshaped(9, 9);
        // addHessianEntry<3, 3>(entries, {indices[0], indices[1], indices[2]}, hessian);

VectorXT xp(9 + 11);
            xp.segment<3>(0) = x0; xp.segment<3>(3) = x1;xp.segment<3>(6) = x2;
            xp[9] = nu; xp[10] = k_s; 
            xp.segment<3>(11) = X0; xp.segment<3>(14) = X1; xp.segment<3>(17) = X2;
            Matrix<T, 9, 9> hessian; hessian.setZero();
            StiffnessMatrix hessian_sparse = autodiff_energies[1].Hessian(xp);
            for (int k = 0; k < hessian_sparse.outerSize(); ++k)
                for (StiffnessMatrix::InnerIterator it(hessian_sparse, k); it; ++it)
                    hessian(it.row(), it.col()) = it.value();

void DiscreteShell::generateAutodiffRuntime()
{
    
    using namespace ADCGWrapper;

    auto normalized = [](const VectorXad &x)  -> VectorXad
    {
        A z = x.squaredNorm();
        return x / CppAD::CondExpGt(z, A{0}, sqrt(z), A{1.0});
    };

    // variable are those we take derivatives of
    // parameters are constants

    // index_t variable_size = 9_idx;
    // index_t parameter_size = 14_idx;
    // xp => x[x0, x1, x2], p[undeformedDensity, thickness, gravity, X0, X1, X2]
    const Autodiff::Function::Blueprint gravitional_energy_blueprint{
        [&](const VectorXad& xp, VectorXad& y) {
            y = VectorXr::Zero(1);
            VectorXad xUndef = xp.segment(9_idx + 5, 9);
            auto undeformedDensity = xp[9_idx + 0];
            auto thickness = xp[9_idx + 1];
            auto x = xp.segment<9_idx>(0);
            auto gravity = xp.segment<3>(9_idx + 2);
            auto undefArea = 0.5 * (xUndef.segment<3>(3) - xUndef.segment<3>(0)).cross(xUndef.segment<3>(6) - xUndef.segment<3>(0)).norm();

            y[0] += (undefArea * thickness) * undeformedDensity * gravity.dot(1.0/3.0 * 
                x.segment<3>(0) + 1.0/3.0 * x.segment<3>(3) + 1.0/3.0 * x.segment<3>(6)
            );
        },
        9_idx,
        14_idx,
        "computeGravitionalEnergyDerivatives"sv,
        EnabledDerivatives::JACOBIAN | EnabledDerivatives::HESSIAN};
    
    autodiff_energies.emplace_back(Autodiff::MakeFunction(gravitional_energy_blueprint, true));

    auto barycentricJacobian = [&](const AV &a, const AV &b, const AV &c)
    {
        AV v0 = b - a, v1 = c - a;
        //AV v2 = p - a;
        A d00 = v0.dot(v0); // (b - a).dot(b - a) => d d00 / d p = 0
        A d01 = v0.dot(v1); // (b - a).dot(c - a) => d d01 / d p = 0
        A d11 = v1.dot(v1); // (c - a).dot(c - a) => d d11 / dp = 0
        //A d20 = v2.dot(v0); // (p - a).dot(b - a) => d d20 / dp = (b - a)^T
        //A d21 = v2.dot(v1); //  ( p - a).dot(c - a) = > d21 / dp = (c - a)^T
        A denom = d00 * d11 - d01 * d01; // => d00 * d11 constant in => drops out, d01 constant in p => derivative is 0
        //v = (d11 * d20 - d01 * d21) / denom;
        //w = (d00 * d21 - d01 * d20) / denom;
        //u = 1.0f - v - w;
        //AV dvdp = (d11 * dd20 / dp - d01 * d d21 / dp) / denom;
        AV dvdp = (d11 * (b - a) - d01 * (c - a)) / denom;
        AV dwdp = (d00 * (c - a) - d01 * (b - a)) / denom;
        Matrix<A, 2, 3> result;
        result.row(0) = dvdp.transpose();
        result.row(1) = dwdp.transpose();
        return result;
    };

    auto computeDeformationGradient = [&](const AV& x1Undef, const AV& x2Undef, const AV& x3Undef,
        const AV& x1, const AV& x2, const AV& x3)
    {
        auto tUndef = normalized(x2Undef - x1Undef) ;
        auto e2Undef = (x3Undef - x1Undef);
        auto qUndef = normalized(e2Undef - tUndef * e2Undef.dot(tUndef));

        AM x;
         x << x1, x2, x3;
        // x.row(0) = x1;
        // x.row(1) = x2;
        // x.row(2) = x3;

        //N(b) = [1 - b1 - b2, b1, b2]
        Matrix<A, 3, 2> dNdb;
        dNdb.row(0) = AV2(-1.0, -1.0);
        dNdb.row(1) = AV2(1.0, 0.0);
        dNdb.row(2) = AV2(0.0, 1.0);
        

        Matrix<A, 2, 3> dBdX = barycentricJacobian(x1Undef, x2Undef, x3Undef);
        Matrix<A, 3, 2> dXdXStar;
        dXdXStar << tUndef, qUndef;
        Matrix<A, 3, 2> def_grad = x * dNdb * dBdX * dXdXStar;
        return def_grad;
    };

    auto computeGreenStrain = [&](const Matrix<A, 3, 2> &defGrad)
    {
        AM2 greenStrain = 0.5 * (defGrad.transpose() * defGrad - AM2::Identity());
        return greenStrain;
    };

    // variable_size = 9_idx;
    // parameter_size = 11_idx;
    // xp => x[x0, x1, x2], p[nu, ks, X0, X1, X2]
    const Autodiff::Function::Blueprint shell_stretching_energy_blueprint{
        [&](const VectorXad& xp, VectorXad& y) 
        {
            y = VectorXr::Zero(1);

            VectorXad xUndef = xp.segment(9_idx + 2, 9);
            AV x1Undef = xUndef.segment<3>(0);
            AV x2Undef = xUndef.segment<3>(3);
            AV x3Undef = xUndef.segment<3>(6);

            A nu = xp[9_idx + 0];
            A ks = xp[9_idx + 1];
            AV x1 = xp.segment<3>(0);
            AV x2 = xp.segment<3>(3);
            AV x3 = xp.segment<3>(6);
            
            A undefArea = 0.5 * (x2Undef - x1Undef).cross(x3Undef - x1Undef).norm();

            Matrix<A, 3, 2> defGrad = computeDeformationGradient(x1Undef, x2Undef, x3Undef, x1, x2, x3);
            AM2 greenStrain = computeGreenStrain(defGrad);

            A trE = greenStrain.trace();
            A trESq = greenStrain.squaredNorm();
            y[0] = ks * undefArea * ((1.0 - nu) * trESq + nu * (trE * trE));
        },
        9_idx,
        11_idx,
        "computeShellStretchingEnergyDerivatives"sv,
        EnabledDerivatives::JACOBIAN | EnabledDerivatives::HESSIAN};
    
    autodiff_energies.emplace_back(Autodiff::MakeFunction(shell_stretching_energy_blueprint, true));

    // index_t variable_size = 12_idx;
    // index_t parameter_size = 13_idx;
    // const Autodiff::Function::Blueprint shell_bending_energy_blueprint{
    //     [&](const VectorXad& xp, VectorXad& y) {
    //         y = VectorXr::Zero(1);
    //         VectorXad xUndef = xp.segment(variable_size + 5, 9);
    //         auto undeformedDensity = xp[variable_size + 0];
    //         auto thickness = xp[variable_size + 1];
    //         auto x = xp.segment<variable_size>(0);
    //         auto gravity = xp.segment<3>(variable_size + 2);
    //         auto undefArea = 0.5 * (xUndef.segment<3>(3) - xUndef.segment<3>(0)).cross(xUndef.segment<3>(6) - xUndef.segment<3>(0)).norm();

    //         y[0] += (undefArea * thickness) * undeformedDensity * gravity.dot(1.0/3.0 * 
    //             x.segment<3>(0) + 1.0/3.0 * x.segment<3>(3) + 1.0/3.0 * x.segment<3>(6)
    //         );
    //     },
    //     variable_size,
    //     parameter_size,
    //     "computeShellBendingEnergyDerivatives"sv,
    //     EnabledDerivatives::JACOBIAN | EnabledDerivatives::HESSIAN};
    
    // autodiff_energies.emplace_back(Autodiff::MakeFunction(shell_bending_energy_blueprint, true));
}