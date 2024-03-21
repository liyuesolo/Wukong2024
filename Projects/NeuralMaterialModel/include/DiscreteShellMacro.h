#ifndef DISCRETE_SHELL_MACRO_H
#define DISCRETE_SHELL_MACRO_H

#include "../../DiscreteShell/include/DiscreteShell.h"


// using discrete elastic shell as macro model
// and the native scale model is from the template simulation
template<class NativeScaleModel>
class DiscreteShellMacro : public DiscreteShell
{
public:
    NativeScaleModel& native_scale_model;

public:
    void addShellInplaneEnergy(T energy);
    void addShellInplaneForceEntry(const VectorXT& residual);

    DiscreteShellMacro(NativeScaleModel& _native_scale_model) 
        : native_scale_model(_native_scale_model) {}
    ~DiscreteShellMacro() {}
};

#endif