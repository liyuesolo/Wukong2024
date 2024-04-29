#include "../include/Homogenization.h"
#include "../include/DiscreteShellHomogenization.h"

template<class NativeScaleModel>
void Homogenization<NativeScaleModel>::imposeUniaxialBending(T angle, T curvature)
{
    native_scale_model.imposeUniaxialBending(angle, curvature);
}

template class Homogenization<DiscreteShellHomogenization>;