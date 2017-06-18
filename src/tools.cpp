#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {

	VectorXd rmse(4);
	rmse << 0,0,0,0;

	if (estimations.size() == 0) {
	    cout << "Estimations size 0" << endl;
	    return rmse;
	}

	if (estimations.size() != ground_truth.size()) {
	    cout << "GT & estimation size mismatch" << endl;
	    return rmse;
	}

	std::vector<VectorXd> vecs;

	//accumulate squared residuals
	for(int i=0; i < estimations.size(); ++i){
    VectorXd err = estimations[i] - ground_truth[i];
		vecs.push_back((err.array() * err.array()).matrix());
	}

	//calculate the mean

	for (int i = 0; i < vecs.size(); ++i) {
	    rmse += vecs[i];
	}

	rmse /= vecs.size();

	//calculate the squared root

	rmse = rmse.array().sqrt().matrix();
	return rmse;
}
