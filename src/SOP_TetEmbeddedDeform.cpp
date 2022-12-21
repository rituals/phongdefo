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

static PRM_Name radiusName("searchradius", "Search radius");
static PRM_Default radiusDefault(0.2f);
static PRM_Range radiusRange(PRM_RANGE_UI, 0.1f, PRM_RANGE_UI, 10.0f);

static PRM_Name enablePhongName(
    "enablephongdeformation", "Phong Deformation");

PRM_Template
SOP_TetEmbeddedDeform::myTemplateList[] = {
    PRM_Template(PRM_FLT_J, 1, &radiusName, &radiusDefault, 0, &radiusRange),
    PRM_Template(PRM_TOGGLE, 1, &enablePhongName),
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

    //float now = context.getTime(); 
    //const bool phong = evalInt("enablephongdeformation", 0, now); 

    duplicatePointSource(0, context);
    GA_RWHandleV3 pointsAttrHandle(gdp->getP());

    for (GA_Size i = 0; i < gdp->getNumPoints(); ++i) {
        const GA_Offset offset = gdp->pointOffset(i);
        UT_Vector3 point = pointsAttrHandle.get(offset);
        //if (phong) {
            //_tracker->applyPhongDeformation(&point, i);
        //} else {
        //    _tracker->applyLinearDeformation(&point, i);
        //}
        point = _tracker->applyPhongDeformation(point, i);
        pointsAttrHandle.set(offset, point);
    }

    pointsAttrHandle.bumpDataId();
    
    return error();
}


