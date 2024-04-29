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
    virtual void addShellInplaneEnergy(T& energy) override;
    virtual void addShellInplaneForceEntries(VectorXT& residual) override;
    virtual void addShellInplaneHessianEntries(std::vector<Entry>& entries) override;

    DiscreteShellMacro(NativeScaleModel& _native_scale_model) 
        : native_scale_model(_native_scale_model) {}
    ~DiscreteShellMacro() {}
};

#endif