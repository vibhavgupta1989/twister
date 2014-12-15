#include "MyWorld.h"
#include "dart/utils/Paths.h"
#include "dart/utils/SkelParser.h"
#include "dart/dynamics/BodyNode.h"
#include "dart/dynamics/Joint.h"
#include "dart/dynamics/Marker.h"
#include <iostream>

using namespace Eigen;
using namespace dart::dynamics;

MyWorld::MyWorld() {
    // Load a skeleton from file

    mSkel = dart::utils::SkelParser::readSkeleton(DART_DATA_PATH"skel/human.skel");
    numConstraints=15;
    mJ = MatrixXd::Zero(3*numConstraints, mSkel->getNumDofs());
    mC = VectorXd::Zero(3*numConstraints);
    for(int i=0;i<numConstraints;i++)
    {    
        mConstrainedMarker[i] = false;
        mTarget[i].x() = 0.0;
        mTarget[i].y() = 0.0;
        mTarget[i].z() = 0.0;
    }
    for (int i = 0; i < 15; i++){
        BodyNode *node = mSkel->getMarker(i)->getBodyNode();
        bodyNode[i] = node;
    }
}

MyWorld::~MyWorld() {
    delete mSkel;
}

void MyWorld::solve(int activeMarker) {

    bool noMarkerSet = true;
    int numIter = 300;
    double alpha = 0.01;
    int nDof = mSkel->getNumDofs();
    VectorXd gradients(nDof);
    VectorXd newPose(nDof);

    for (int i = 0; i < numConstraints; i++){
        if(mConstrainedMarker[i]==true){
            noMarkerSet = false;
            break;
        }
    }

    if(noMarkerSet==true){
        return;
    }
    
    for (int i = 0; i < numIter; i++) {
        gradients = updateGradients();
        newPose = mSkel->getPositions() - alpha * gradients;
        mSkel->setPositions(newPose); 
        mSkel->computeForwardKinematics(true, false, false); // DART updates all the transformations based on newPose
        mJ = MatrixXd::Zero(3*numConstraints, mSkel->getNumDofs());
    }
}

VectorXd MyWorld::updateGradients() {

    VectorXd gradients;
    VectorXd v;
    Vector3d diff;
    Vector4d offset;
    Vector4d jCol;
    int nDofs;
    Matrix4d worldToParent;
    Matrix4d parentToJoint;
    Matrix4d jointToChild;
    Matrix4d dR;
    Matrix4d R, R1, R2;
    int colIndex;

    for (int i=0; i < numConstraints; i++){
        if(mConstrainedMarker[i]==true){

            diff = mSkel->getMarker(i)->getWorldPosition() - mTarget[i];

            mC[3*i] = diff.x();
            mC[3*i+1] = diff.y();
            mC[3*i+2] = diff.z();

            offset << mSkel->getMarker(i)->getLocalPosition(), 1;
            BodyNode *node = mSkel->getMarker(i)->getBodyNode();
            while(node){
                Joint *joint = node->getParentJoint();
                if(joint){
                    nDofs = joint->getNumDofs();

                    if(node->getParentBodyNode()){
                        worldToParent = node->getParentBodyNode()->getTransform().matrix();
                        parentToJoint = joint->getTransformFromParentBodyNode().matrix();
                    }
                    jointToChild =joint->getTransformFromChildBodyNode().inverse().matrix();

                    if(nDofs==1){ // joint has 1 degree of freedom
                        dR = joint->getTransformDerivative(0);
                        if(node->getParentBodyNode()){
                            jCol = worldToParent * parentToJoint * dR * jointToChild * offset;
                        }
                        else{
                            jCol = dR * jointToChild * offset;
                        }
                        colIndex = joint->getIndexInSkeleton(0);
                        mJ.col(colIndex)[3*i] = jCol.x();
                        mJ.col(colIndex)[3*i+1] = jCol.y();
                        mJ.col(colIndex)[3*i+2] = jCol.z();
                        // reset offset when moving towards root
                        if(node->getParentBodyNode()){
                            offset = parentToJoint * joint->getTransform(0) * jointToChild * offset;
                        }
                        else{
                            offset = joint->getTransform(0) * jointToChild * offset;
                        }
                    }
                    else if(nDofs==2){ // joint has 2 degrees of freedom
                        dR = joint->getTransformDerivative(0); // Doesn't need .matrix() because it returns a Matrix4d instead of Isometry3d
                        R = joint->getTransform(1).matrix();
                        if(node->getParentBodyNode()){
                            jCol = worldToParent * parentToJoint * dR * R * jointToChild * offset;
                        }else{
                            jCol = dR * R * jointToChild * offset;
                        }
                        colIndex = joint->getIndexInSkeleton(0);
                        mJ.col(colIndex)[3*i] = jCol.x();
                        mJ.col(colIndex)[3*i+1] = jCol.y();
                        mJ.col(colIndex)[3*i+2] = jCol.z();

                        dR = joint->getTransformDerivative(1);
                        R = joint->getTransform(0).matrix();
                        if(node->getParentBodyNode()){
                            jCol = worldToParent * parentToJoint * R * dR * jointToChild * offset;
                        }
                        else{
                            jCol = R * dR * jointToChild * offset;
                        }
                        colIndex = joint->getIndexInSkeleton(1);
                        mJ.col(colIndex)[3*i] = jCol.x();
                        mJ.col(colIndex)[3*i+1] = jCol.y();
                        mJ.col(colIndex)[3*i+2] = jCol.z();
                        // reset offset when moving towards root
                        if(node->getParentBodyNode()){
                            offset = parentToJoint * joint->getTransform(0) * joint->getTransform(1) * jointToChild * offset;
                        }
                        else{
                            offset = joint->getTransform(0) * joint->getTransform(1) * jointToChild * offset;
                        }
                    }
                    else{ // joint has 3 degrees of freedom
                        dR = joint->getTransformDerivative(0); // Doesn't need .matrix() because it returns a Matrix4d instead of Isometry3d
                        R1 = joint->getTransform(1).matrix();
                        R2 = joint->getTransform(2).matrix();
                        if(node->getParentBodyNode()){
                            jCol = worldToParent * parentToJoint * dR * R1 * R2 * jointToChild * offset;
                        }
                        else{
                            jCol = dR * R1 * R2 * jointToChild * offset;
                        }
                        colIndex = joint->getIndexInSkeleton(0);
                        mJ.col(colIndex)[3*i] = jCol.x();
                        mJ.col(colIndex)[3*i+1] = jCol.y();
                        mJ.col(colIndex)[3*i+2] = jCol.z();

                        dR = joint->getTransformDerivative(1);
                        R1 = joint->getTransform(0).matrix();
                        R2 = joint->getTransform(2).matrix();
                        if(node->getParentBodyNode()){
                            jCol = worldToParent * parentToJoint * R1 * dR * R2 * jointToChild * offset;
                        }
                        else{
                            jCol = R1 * dR * R2 * jointToChild * offset;
                        }
                        colIndex = joint->getIndexInSkeleton(1);
                        mJ.col(colIndex)[3*i] = jCol.x();
                        mJ.col(colIndex)[3*i+1] = jCol.y();
                        mJ.col(colIndex)[3*i+2] = jCol.z();

                        dR = joint->getTransformDerivative(2);
                        R1 = joint->getTransform(0).matrix();
                        R2 = joint->getTransform(1).matrix();
                        if(node->getParentBodyNode()){
                            jCol = worldToParent * parentToJoint * R1 * R2 * dR * jointToChild * offset;
                        }
                        else{
                            jCol = R1 * R2 * dR * jointToChild * offset;
                        }
                        colIndex = joint->getIndexInSkeleton(2);
                        mJ.col(colIndex)[3*i] = jCol.x();
                        mJ.col(colIndex)[3*i+1] = jCol.y();
                        mJ.col(colIndex)[3*i+2] = jCol.z();

                        // reset offset when moving towards root
                        if(node->getParentBodyNode()){
                            offset = parentToJoint * joint->getTransform(0) * joint->getTransform(1) * joint->getTransform(2) * jointToChild * offset;
                        }
                        else{
                            offset = joint->getTransform(0) * joint->getTransform(1) * joint->getTransform(2) * jointToChild * offset;
                        }
                    }
                }
                node = node->getParentBodyNode();
            }
        }
    }
    gradients = 2.0 * mJ.transpose() * mC;
    return gradients;
}

void MyWorld::createConstraint(int _index) {
    std::cout << "index = " << _index << std::endl;
    mTarget[_index] = mSkel->getMarker(_index)->getWorldPosition();
    mConstrainedMarker[_index] = true;
}

void MyWorld::modifyConstraint(int _index, Vector3d _deltaP) {
    if (mConstrainedMarker[_index] == true){
        mTarget[_index] += _deltaP;
    }
}

void MyWorld::removeConstraint(int _index) {
    mConstrainedMarker[_index] = false;
}


