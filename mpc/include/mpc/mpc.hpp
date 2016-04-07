/*
 * mpc.hpp
 *
 *  Created on: Apr 7, 2016
 *      Author: filip
 */

#ifndef MPC_HPP_
#define MPC_HPP_


#include <Eigen/Dense>

/*********************************************************************
 ***  MPC class definition
 ********************************************************************/


class MPC
{

public:

	typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrix;
	typedef Eigen::Matrix<double, Eigen::Dynamic, 1> vector;



	/*****************************************************************
	 ***  Class functions
	 ****************************************************************/

	//predictModelState(), getCost(), setModelState(), getModelState(), setControlHorizon(), getControlHorizon() itd.

	MPC();

	~MPC();

	void setModelState(vector state);

	double getCost();

	vector getModelState();

	void setControlHorizon(int value);


	/*********************************************************************
	 ***  Class variables
	 ********************************************************************/

	vector state_;

	/*** Control horizon ***/
	int control_horizon_;



};



MPC::MPC()
{
	state_ = Eigen::VectorXd::Zero(5);
	control_horizon_ = 5;

}

MPC::~MPC()
{

}

void MPC::setModelState(vector state)
{
	state_ = state;
}

#endif /* MPC_HPP_ */

