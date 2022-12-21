#include "SOP_TetEmbeddedDeform.h"
#include "DeformationTracker.h"

#include <GA/GA_Handle.h>
#include <GU/GU_Detail.h>
#include <OP/OP_AutoLockInputs.h>
#include <OP/OP_Operator.h>
#include <OP/OP_OperatorTable.h>
#include <PRM/PRM_Include.h>
#include <PRM/PRM_Range.h>
#include <UT/UT_DSOVersion.h>

#include <iostream>

using namespace TetEmbeddedDeform;

void
newSopOperator(OP_OperatorTable *table)
{
    const char *inputLabels[] = {
        "Mesh to deform",
        "Rest tet mesh",
        "Deformed tet mesh",
        nullptr
    };

    table->addOperator(new OP_Operator(
        "tet_embedded_deform",
        "Tet Embedded Deform",
        SOP_TetEmbeddedDeform::myConstructor,
        SOP_TetEmbeddedDeform::myTemplateList,
        3,
        3,
        0,
        0,
        inputLabels));
}

static PRM_Name modeName("mode", "Deformation Mode");
static PRM_Name modeChoices[] = {
    PRM_Name("linear", "Linear"),
    PRM_Name("phong", "Phong"),
    PRM_Name(0)
};
static PRM_ChoiceList modeMenu(PRM_CHOICELIST_SINGLE, modeChoices);

PRM_Template
SOP_TetEmbeddedDeform::myTemplateList[] = {
    PRM_Template(PRM_STRING, 1, &modeName, 0, &modeMenu),
    PRM_Template()
};

OP_Node *
SOP_TetEmbeddedDeform::myConstructor(
    OP_Network *net,
    const char *name, OP_Operator *op)
{
    return new SOP_TetEmbeddedDeform(net, name, op);
}

SOP_TetEmbeddedDeform::SOP_TetEmbeddedDeform(
    OP_Network *net,
    const char *name,
    OP_Operator *op) :
    SOP_Node(net, name, op)
{
    mySopFlags.setManagesDataIDs(true);
}

SOP_TetEmbeddedDeform::~SOP_TetEmbeddedDeform()
{
    if (_tracker) {
        delete _tracker;
        _tracker = nullptr;
    }
}

OP_ERROR
SOP_TetEmbeddedDeform::cookMySop(OP_Context &context)
{
    OP_AutoLockInputs inputs(this);
    if (inputs.lock(context) >= UT_ERROR_ABORT) {
        return error();
    }

    OP_Node::flags().setTimeDep(true);

    if (!_tracker) {
        _tracker = new DeformationTracker();
    }

    const GU_Detail *detail = inputGeo(0, context);
    const GU_Detail *restTetMeshDetail = inputGeo(1, context);
    const GU_Detail *deformedTetMeshDetail = inputGeo(2, context);

    _tracker->update(detail, restTetMeshDetail, deformedTetMeshDetail);

    UT_String modeString;
    evalString(modeString, "mode", 0, context.getTime()); 
    bool doPhong = modeString == "phong";

    duplicatePointSource(0, context);
    GA_RWHandleV3 pointsAttrHandle(gdp->getP());

    for (GA_Size i = 0; i < gdp->getNumPoints(); ++i) {
        const GA_Offset offset = gdp->pointOffset(i);
        UT_Vector3 point = pointsAttrHandle.get(offset);
        if (doPhong) {
            point = _tracker->applyPhongDeformation(point, i);
        } else {
            point = _tracker->applyLinearDeformation(point, i);
        }
        pointsAttrHandle.set(offset, point);
    }

    pointsAttrHandle.bumpDataId();
    
    return error();
}


