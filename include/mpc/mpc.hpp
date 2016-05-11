/*
 * mpc.hpp
 *
 *  Created on: Apr 7, 2016
 *      Author: filip
 */

#ifndef MPC_HPP_
#define MPC_HPP_

#include <symbolicc++.h>
#include <nlopt.h>
//#include <nlopt.hpp>
#include <Eigen/Dense>
//#include <cmath>
# define pi 3.14159265358979323846


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

	void setParameters(int Km, float K, float Tm, float Ts, float d, float alpha);

	void setReference(Symbolic ref);

	void setModel();

	void setDiscreteModel();

	void setModelState(Symbolic x);

	void setHorizon(int Hw, int Hp, int Hu);

	void setMatrixDimension();

	void setControlMatrices(Symbolic Q, Symbolic R);

	void setInitialState(double *u);

	void setNlopt();

	void setBounds(double *lb, double *ub);

	void getTotalMatrices();

	void getPredictModelState();

	void getReferenceWindow();

	void setMinObjective();

	double myvfunc(unsigned dd, const double *u, double *grad, void *my_func_data);

	void substState(const double *u);

	void getPredictOutput();

	void getOutputWindow();

	void substCost(const double *u);

	void getControlWindow(const double *u);

	double getCost();

	Symbolic getControlValue(const double *u);

	Symbolic getModelState(Symbolic u_real);

	Symbolic getModelOutput(Symbolic x_real);


	void test();

	//void setModelState(vector state);
	//void setModelState();

	//double getCost();

	//vector getModelState();

	//void setControlHorizon(int value);

	/*********************************************************************
	 ***  Class variables
	 ********************************************************************/
	vector state_;

	/*** Control horizon ***/
	int control_horizon_;


	/*** Parameters ***/
	int Km_;
	float K_, Tm_, Ts_, d_ , alpha_;

	int i_, j_, n_, m_;

	Symbolic ref_;
	Symbolic A_c_, B_c_, I_;
	Symbolic A_d_, B_d_, C_d_;

	Symbolic x_;

	int Hw_, Hp_, Hu_;
	Symbolic Q_, R_;

	double u_[12];

	Symbolic u11_, u12_, u13_;
	Symbolic u21_, u22_, u23_;
	Symbolic u31_, u32_, u33_;
	Symbolic u41_, u42_, u43_;
	Symbolic ulaz_;
	int vel_;

	double lb_[12], ub_[12];

    nlopt_opt opt_;

    double minf_;

    int k_, ii_, kk_, br_;

    Symbolic A_d_uk1_;
    Symbolic A_d_uk2_;

    Symbolic xp_;
    Symbolic x_uk_;

    //Symbolic r_p_;
    Symbolic r_k_;

    Symbolic nelin_;
    Symbolic z_;
    Symbolic z_k_;

    Symbolic pom1_, pom2_;

    Symbolic u_k_;

    int km_;
    Symbolic V_;
    double vv_;

    //Symbolic u_real_ ;

    int iter_;

};



MPC::MPC()
{
	////////////////////////////////////////

	state_ = Eigen::VectorXd::Zero(5);
	control_horizon_ = 5;
	////////////////////////////////////////

	i_ = 0, j_ = 0, n_ = 0, m_ = 0;
	Km_ = 0, K_ = 0, Tm_ = 0, Ts_ = 0, d_ = 0, alpha_ = 0;

    A_c_ = Symbolic ("A_c_", 4, 4);
    B_c_ = Symbolic("B_c_", 4, 4);
    I_ = Symbolic("I_", 4, 4);
    A_d_ = Symbolic ("A_d_", 4, 4);
    B_d_ = Symbolic("B_d_", 4, 4);
    C_d_ = Symbolic("C_d_", 3, 4);

    x_ = Symbolic("x_", 4, 1);

	Hw_ = 0, Hp_ = 0, Hu_ = 0;

	//Q_ = Symbolic ("Q_", 3*(Hp_-Hw_+1), 3*(Hp_-Hw_+1));
	//R_ = Symbolic ("R_", 4*Hu_, 4*Hu_);
    //x_uk_ = Symbolic ("x_uk_", 4, Hp_);


	u11_ = Symbolic ("u11_"), u12_ = Symbolic ("u12_"), u13_ = Symbolic ("u13_");
	u21_ = Symbolic ("u21_"), u22_ = Symbolic ("u22_"), u23_ = Symbolic ("u23_");
	u31_ = Symbolic ("u31_"), u32_ = Symbolic ("u32_"), u33_ = Symbolic ("u33_");
	u41_ = Symbolic ("u41_"), u42_ = Symbolic ("u42_"), u43_ = Symbolic ("u43_");
    ulaz_ = ((u11_, u12_, u13_),
             (u21_, u22_, u23_),
             (u31_, u32_, u33_),
             (u41_, u42_, u43_));
    vel_ = ulaz_.columns();

    minf_ = 0;
    k_ = 0, ii_ = 0, kk_ = 0, br_ = 0, km_ = 0, vv_ = 0;


    xp_ = Symbolic ("xp_", 4, 1);

    //r_p_ = Symbolic ("r_p_", Hp_-Hw_+1, 3);
    //r_k_ = Symbolic ("r_k_", 3*(Hp_-Hw_+1), 1);

    z_ = Symbolic ("z_", 3, 1);

    //u_real_ = Symbolic ("u_real_", 4, 1);

    iter_ = 1;

}

MPC::~MPC()
{


}

void MPC::setParameters(int Km, float K, float Tm, float Ts, float d, float alpha){

	Km_ = Km;
	K_ = K;
	Tm_ = Tm;
	Ts_ = Ts;
	d_ = d;
	alpha_ = alpha;

}

void MPC::setReference(Symbolic ref){

	ref_ = ref;

}

void MPC::setModel(){

    for (i_ = 0; i_ < 4; i_++){
            for (j_ = 0; j_ < 4; j_++){
                    if (i_ == j_){
                            A_c_(i_,j_) = -1/Tm_;
                            B_c_(i_,j_) = Km_/Tm_;
                            I_(i_,j_) = 1;
                    }
                    else {
                            A_c_(i_,j_) = B_c_(i_,j_) = 0;
                            I_(i_,j_) = 0;
                    }
            }
    }

}

void MPC::setDiscreteModel(){

    A_d_ = I_ + Ts_ * A_c_;
    B_d_ = Ts_ * B_c_;
    C_d_(0,0) = C_d_(0,1) = cos(alpha_);
    C_d_(0,2) = C_d_(0,3) = -cos(alpha_);
    C_d_(1,0) = C_d_(1,2) = sin(alpha_);
    C_d_(1,1) = C_d_(1,3) = -sin(alpha_);
    C_d_(2,0) = C_d_(2,3) = d_;
    C_d_(2,1) = C_d_(2,2) = -d_;

}

void MPC::setModelState(Symbolic x){

	x_ = x;

}

void MPC::setHorizon(int Hw, int Hp, int Hu){

	Hw_ = Hw;
	Hp_ = Hp;
	Hu_ = Hu;

}

void MPC::setMatrixDimension(){

	Q_ = Symbolic ("Q_", 3*(Hp_-Hw_+1), 3*(Hp_-Hw_+1));
	R_ = Symbolic ("R_", 4*Hu_, 4*Hu_);
	x_uk_ = Symbolic ("x_uk_", 4, Hp_);
    //r_p_ = Symbolic ("r_p_", Hp_-Hw_+1, 3);
    r_k_ = Symbolic ("r_k_", 3*(Hp_-Hw_+1), 1);
    nelin_ = Symbolic ("nelin_", 4, (Hp_-Hw_+1));
    z_k_ = Symbolic ("z_k_", 3*(Hp_-Hw_+1), 1);
    u_k_ = Symbolic ("u_k_", 4*Hu_, 1);

}

void MPC::setControlMatrices(Symbolic Q, Symbolic R){

	Q_ = Q;
	R_ = R;

}

void MPC::setInitialState(double *u){

	for (i_ = 0; i_ < 12 ; i_++) {
		u_[i_] = *(u + i_);
	}

}

void MPC::setNlopt(){

    opt_ = nlopt_create(NLOPT_LN_COBYLA, 12); // algorithm and dimensionality //

    nlopt_set_xtol_rel(opt_, 1e-3);
    nlopt_set_ftol_rel(opt_, 1e-3);

}

void MPC::setBounds(double *lb, double *ub){

	for (i_ = 0; i_ < 12 ; i_++) {
		lb_[i_] = *(lb + i_);
		ub_[i_] = *(ub + i_);
	}


    nlopt_set_lower_bounds(opt_, lb_);
    nlopt_set_upper_bounds(opt_, ub_);


}

void MPC::getTotalMatrices(){

    A_d_uk1_ = I_;
    A_d_uk2_ = I_;
    for (ii_ = 0; ii_ < i_+1; ii_++){
            A_d_uk1_ = A_d_ * A_d_uk1_;
            if(ii_ >= i_) break;
            A_d_uk2_ = A_d_uk2_ + A_d_uk1_;
    }


}

void MPC::getPredictModelState(){

    xp_ = A_d_uk1_ * x_ + A_d_uk2_ * B_d_ * ulaz_.column(kk_);

	x_uk_(0,i_) = xp_(0);
	x_uk_(1,i_) = xp_(1);
	x_uk_(2,i_) = xp_(2);
	x_uk_(3,i_) = xp_(3);

}

void MPC::getReferenceWindow(){

    //vektor dimenzija 1x3*(Hp-Hw+1)
	br_ = 0;
    for (j_ = Hw_; j_ <= Hp_; j_++){
            for (i_ = 0; i_ < 3; i_++){
                    r_k_(br_,0) = ref_(i_,j_+k_);
                    br_++;
            }
    }

}

void MPC::setMinObjective(){


}


double MPC::myvfunc(unsigned dd, const double *u, double *grad, void *my_func_data){

	double ret = 8.0;
	return ret;
}


void MPC::substState(const double *u){

    x_uk_ = x_uk_.subst(u11_ == u[0*vel_+0]);
    x_uk_ = x_uk_.subst(u12_ == u[0*vel_+1]);
    x_uk_ = x_uk_.subst(u13_ == u[0*vel_+2]);
    x_uk_ = x_uk_.subst(u21_ == u[1*vel_+0]);
    x_uk_ = x_uk_.subst(u22_ == u[1*vel_+1]);
    x_uk_ = x_uk_.subst(u23_ == u[1*vel_+2]);
    x_uk_ = x_uk_.subst(u31_ == u[2*vel_+0]);
    x_uk_ = x_uk_.subst(u32_ == u[2*vel_+1]);
    x_uk_ = x_uk_.subst(u33_ == u[2*vel_+2]);
    x_uk_ = x_uk_.subst(u41_ == u[3*vel_+0]);
    x_uk_ = x_uk_.subst(u42_ == u[3*vel_+1]);
	x_uk_ = x_uk_.subst(u43_ == u[3*vel_+2]);

}


void MPC::getPredictOutput(){

	/*
	for (m_ = 0; m_ < Hp_; m_++){
	    for (n_ = 0; n_ < 4; n_++){

	    	dd_ = (double)(x_uk_(n_,m_));
	    	nelin_(n_,m_) = K_ * (abs(dd_)) * dd_;

	    	}

    	z_ = C_d_ * nelin_.column(m_);

    	z_uk_(0,m_) = z_(0);
    	z_uk_(1,m_) = z_(1);
    	z_uk_(2,m_) = z_(2);
	}
	*/

	br_ = 0;
	for (m_ = Hw_; m_ <= Hp_; m_++){
	    for (n_ = 0; n_ < 4; n_++){

	    	nelin_(n_,m_-Hw_) = K_ * (abs((double)(x_uk_(n_,m_-1)))) * ((double)(x_uk_(n_,m_-1)));

	    	}

    	z_ = C_d_ * nelin_.column(m_-Hw_);

    	/*
    	z_uk_(0,m_-Hw_) = z_(0);
    	z_uk_(1,m_-Hw_) = z_(1);
    	z_uk_(2,m_-Hw_) = z_(2);
    	*/

    	z_k_(br_,0) = z_(0);
    	z_k_(br_+1,0) = z_(1);
    	z_k_(br_+2,0) = z_(2);
    	br_ = br_ + 3;
	}

}

void MPC::getOutputWindow(){


}

void MPC::substCost(const double *u){

    pom1_ = ((z_k_.transpose())-(r_k_.transpose()))*Q_*(z_k_-r_k_);
    pom1_ = pom1_.subst(u11_ == u[0*vel_+0]);
    pom1_ = pom1_.subst(u12_ == u[0*vel_+1]);
    pom1_ = pom1_.subst(u13_ == u[0*vel_+2]);
    pom1_ = pom1_.subst(u21_ == u[1*vel_+0]);
    pom1_ = pom1_.subst(u22_ == u[1*vel_+1]);
    pom1_ = pom1_.subst(u23_ == u[1*vel_+2]);
    pom1_ = pom1_.subst(u31_ == u[2*vel_+0]);
    pom1_ = pom1_.subst(u32_ == u[2*vel_+1]);
    pom1_ = pom1_.subst(u33_ == u[2*vel_+2]);
    pom1_ = pom1_.subst(u41_ == u[3*vel_+0]);
    pom1_ = pom1_.subst(u42_ == u[3*vel_+1]);
	pom1_ = pom1_.subst(u43_ == u[3*vel_+2]);

}

void MPC::getControlWindow(const double *u){

    br_ = 0;
    for (m_ = 0; m_ < Hu_; m_++){
    	for (n_ = 0; n_ < 4; n_++){
    		if (m_ > 2) km_ = 2;
    		else km_ = m_;

    		u_k_(br_,0) = u[n_*vel_+km_];
            br_++;
        }
    }

}

double MPC::getCost(){

    pom2_ = (u_k_.transpose())*R_*(u_k_);

    return ((double)(pom1_+pom2_));

}

Symbolic MPC::getControlValue(const double *u){

    Symbolic u_real_("u_real_", 4, 1);

    u_real_(0,0) = u[0*vel_+0];
    u_real_(1,0) = u[1*vel_+0];
    u_real_(2,0) = u[2*vel_+0];
    u_real_(3,0) = u[3*vel_+0];

    return u_real_;

}

Symbolic MPC::getModelState(Symbolic u_real){

    x_ = A_d_ * x_ + B_d_ * u_real;

	return x_;
}

Symbolic MPC::getModelOutput(Symbolic x_real){

	Symbolic z_("z_", 3, 1);

	nelin_(0,0) = K_ * (abs((double)(x_real(0)))) * ((double)(x_real(0)));
	nelin_(1,0) = K_ * (abs((double)(x_real(1)))) * ((double)(x_real(1)));
	nelin_(2,0) = K_ * (abs((double)(x_real(2)))) * ((double)(x_real(2)));
	nelin_(3,0) = K_ * (abs((double)(x_real(3)))) * ((double)(x_real(3)));

	z_ = C_d_ * nelin_.column(0);

	return z_;
}


void MPC::test(){

	cout << "\n test \n";
}

////////////////////////////////////////

//void MPC::setModelState(vector state)


////////////////////////////////////////

#endif /* MPC_HPP_ */

