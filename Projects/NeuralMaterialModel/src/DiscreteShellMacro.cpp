#include "../include/NeuralMaterialModel.h"
#include "../include/DiscreteShellMacro.h"

template<class NativeScaleModel>
void DiscreteShellMacro<NativeScaleModel>::addShellInplaneEnergy(T& energy)
{
    
}

template<class NativeScaleModel>
void DiscreteShellMacro<NativeScaleModel>::addShellInplaneForceEntries(VectorXT& residual)
{

}

template<class NativeScaleModel>
void DiscreteShellMacro<NativeScaleModel>::addShellInplaneHessianEntries(std::vector<Entry>& entries)
{

}

template class DiscreteShellMacro<NeuralMaterialModel>;