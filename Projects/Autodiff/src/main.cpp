#include "../include/Autodiff.h"
#include "../include/App.h"
#include "adcg_wrapper/autodiff/function.hpp"
#include "../autodiff/CST3DShell.h"
#include "../include/Timer.h"

int main()
{
	using namespace ADCGWrapper;


    // add_argument_expr(function, undeformedDensity);
	// add_argument_expr(function, thickness);
	// add_argument_vector(function, gravity, 3);
	// add_argument_vector(function, x1, 3);
	// add_argument_vector(function, x2, 3);
	// add_argument_vector(function, x3, 3);
	// add_argument_vector(function, x1Undef, 3);
	// add_argument_vector(function, x2Undef, 3);
	// add_argument_vector(function, x3Undef, 3);

	// Expr third = Expr(1.0) / Expr(3.0);

	// Expr undefArea = Expr(0.5) * (x2Undef - x1Undef).cross(x3Undef - x1Undef).norm();
	// Expr energy = (undefArea * thickness) * undeformedDensity * gravity.dot(third * x1 + third * x2 + third * x3);


    /***************   Generate derivatives for functions of variable maps. */
    const index_t xSize = 9_idx;
    const index_t pSize = 16_idx;
    const Autodiff::Function::Blueprint functionBlueprint{
        [&](const VectorXad& xp, VectorXad& y) {
            y = VectorXr::Zero(1);
            VectorXad xUndef = xp.segment(xSize + 5, 9);
            auto undeformedDensity = xp[xSize + 0];
            auto thickness = xp[xSize + 1];
            auto x = xp.segment<xSize>(0);
            auto gravity = xp.segment<3>(xSize + 2);
            auto undefArea = 0.5 * (xUndef.segment<3>(3) - xUndef.segment<3>(0)).cross(xUndef.segment<3>(6) - xUndef.segment<3>(0)).norm();

            y[0] += (undefArea * thickness) * undeformedDensity * gravity.dot(1.0/3.0 * 
                x.segment<3>(0) + 1.0/3.0 * x.segment<3>(3) + 1.0/3.0 * x.segment<3>(6)
            );
        },
        xSize,
        pSize,
        "function_example"sv,
        EnabledDerivatives::JACOBIAN | EnabledDerivatives::HESSIAN};
    const auto function = Autodiff::MakeFunction(functionBlueprint, true);

    VectorXr xp                         = VectorXr::Random(xSize + pSize);
    const VectorXr value                = function(xp);
    
    const SparseMatrix<T> jacobian = function.Jacobian(xp);
    START_TIMING(grad_cppad)
    // for (int i = 0; i < 10000; i++)
    tbb::parallel_for(0, 10000, [&](int i)
    {
        const SparseMatrix<T> jac = function.Jacobian(xp);
        const SparseMatrix<T> hessian  = function.Hessian(xp);
    });
    FINISH_TIMING_PRINT(grad_cppad)

    std::cout << "Jonas AD" << std::endl;
    Eigen::Matrix<T, 9, 1, 0, 9, 1> dedx;
    Eigen::Matrix<T, 9, 9, 0, 9, 9> d2edx2;
    START_TIMING(grad_jonas)
    // for (int i = 0; i < 10000; i++)
    tbb::parallel_for(0, 10000, [&](int i)
    {
        compute3DCSTGravitationalEnergyGradient(
            xp[xSize], xp[xSize + 1], xp.segment<3>(xSize + 2), xp.segment<3>(0),
            xp.segment<3>(3), xp.segment<3>(6), 
            xp.segment<3>(xSize + 5),
            xp.segment<3>(xSize + 8),
            xp.segment<3>(xSize + 11), dedx
        );
        compute3DCSTGravitationalEnergyHessian(
            xp[xSize], xp[xSize + 1], xp.segment<3>(xSize + 2), xp.segment<3>(0),
            xp.segment<3>(3), xp.segment<3>(6), 
            xp.segment<3>(xSize + 5),
            xp.segment<3>(xSize + 8),
            xp.segment<3>(xSize + 11), d2edx2
        );
    });
    FINISH_TIMING_PRINT(grad_jonas)
    std::cout << dedx.transpose() << std::endl;
    std::cout << "Cpp AD" << std::endl;
    std::cout << jacobian.toDense() << std::endl;

	return 0;
}