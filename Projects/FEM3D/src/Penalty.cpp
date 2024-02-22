#include "../include/FEM3D.h"

T height = 0.2;

T FEM3D::addCubicPlaneEnergy(T w)
{
    T energy = 0.0;
    for (int i = 0; i < num_nodes; i++)
    {
        TV xi = deformed.segment<3>(i * 3);
        if (xi[1] > height)
            energy += w * std::pow(xi[1] - height, 3);   
    }
    return energy;
}

void FEM3D::addCubicPlaneForceEntry(VectorXT& residual, T w)
{
    for (int i = 0; i < num_nodes; i++)
    {
        TV xi = deformed.segment<3>(i * 3);
        if (xi[1] > height)
        {
            residual[i * 3 + 1] += -w * 3.0 * std::pow(xi[1] - height, 2);
        }
    }
}

void FEM3D::addCubicPlaneHessianEntry(std::vector<Entry>& entries, T w)
{
    for (int i = 0; i < num_nodes; i++)
    {
        TV xi = deformed.segment<3>(i * 3);
        if (xi[1] > height)
        {
            entries.emplace_back(i * 3 + 1, i * 3 + 1, w * 6.0 * (xi[1] - height));
        }
    }
}