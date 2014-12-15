#ifndef _MYWORLD_
#define _MYWORLD_

#include <vector>
#include "dart/dynamics/Skeleton.h"


class MyWorld {
 public:
    MyWorld();
    virtual ~MyWorld();
    dart::dynamics::Skeleton* getSkel() {
        return mSkel;
    }

    void solve(int activeMarker);
    void createConstraint(int _index);
    void modifyConstraint(int _index, Eigen::Vector3d _deltaP);
    void removeConstraint(int _index);

 protected:
    Eigen::VectorXd updateGradients();

    dart::dynamics::Skeleton *mSkel;
    Eigen::VectorXd mC;      // Hard-coding that numConstraints=15
    Eigen::MatrixXd mJ;
    int numConstraints;
    Eigen::Vector3d mTarget[15]; // The target location of the constriant
    bool mConstrainedMarker[15]; // The index of the constrained marker

    void *bodyNode[15];
};

#endif
