#ifndef __SOP_TETEMBEDDEDDEFORM_H__
#define __SOP_TETEMBEDDEDDEFORM_H__

#include <SOP/SOP_Node.h>

namespace TetEmbeddedDeform
{
    class DeformationTracker;
}

/// SOP for interpolating the deformation of a tetrahedral mesh onto a set of 
/// embedded points.
class SOP_TetEmbeddedDeform : public SOP_Node
{
public:
    SOP_TetEmbeddedDeform(OP_Network *net, const char *name, OP_Operator *op);
    ~SOP_TetEmbeddedDeform() override;

    static OP_Node *myConstructor(
        OP_Network *net,
        const char *name,
        OP_Operator *entry);

    static PRM_Template myTemplateList[];

protected:
    OP_ERROR cookMySop(OP_Context &context) override;

private:
    TetEmbeddedDeform::DeformationTracker *_tracker;  
};

#endif
