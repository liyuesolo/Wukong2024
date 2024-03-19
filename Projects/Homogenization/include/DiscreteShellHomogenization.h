#ifndef DISCRETE_SHELL_HOMOGENIZATION_H
#define DISCRETE_SHELL_HOMOGENIZATION_H

#include "../../DiscreteShell/include/DiscreteShell.h"

class DiscreteShellHomogenization : public DiscreteShell
{
public:
    
public:
    void imposeUniaxialBending(T angle, T curvature) {}

    DiscreteShellHomogenization() {}
    ~DiscreteShellHomogenization() {}
};

#endif