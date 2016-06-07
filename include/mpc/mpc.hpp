/*
 * mpc.hpp
 *
 *  Created on: Apr 7, 2016
 *      Author: filip
 */

#ifndef MPC_HPP_
#define MPC_HPP_


//#include <nlopt.h>
#include <Eigen/Dense>
//#include <cmath>
#include <nlopt.hpp>
#include <symbolicc++.h>

# define pi 3.14159265358979323846


/*********************************************************************
 ***  MPC class definition
 ********************************************************************/


class MPC
{

public:

	typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrix;
	typedef Eigen::Matrix<double, Eigen::Dynamic, 1> vector;

	enum {state_num = 12};



	/*****************************************************************
	 ***  Class functions
	 ****************************************************************/

	MPC();

	~MPC();

	void setParameters(int Km, float K, float Tm, float Ts, float d, float alpha);

	void setReference(matrix ref);

	void setModel();

	void setDiscreteModel();

	void getSymbolicMatrices();

	void setModelState(vector x);

	void setHorizon(int Hw, int Hp, int Hu);

	void setMatrixDimension();

	void setControlMatrices(matrix Q, matrix R);

	void setInitialState(double *u);

	void setNlopt();

	void setBounds(double *lb, double *ub);

	void getSymbolicState();

	void getPredictModelState();

	void getReferenceWindow();

	double myfunc(const std::vector<double> &u, std::vector<double> &grad, void *my_func_data);

	Symbolic substState(const double *u, Symbolic x_uk_pom);

	void getPredictOutput(Symbolic x_uk_pom);

	void substCost();

	void getControlWindow(const double *u);

	double getCost();

	vector getControlValue(const double *u);

	vector getModelState(vector u_real);

	vector getModelOutput(vector x_real);


	/*********************************************************************
	 ***  Class variables
	 ********************************************************************/

	/*** Parameters ***/
	int Km_;
	float K_, Tm_, Ts_, d_ , alpha_;

	/*** Counters ***/
	int i_, j_, n_, m_;

	/*** Matrices ***/
	matrix ref_;				//reference matrix
	matrix A_c_, B_c_, I_;		//continuous model
	matrix A_d_, B_d_, C_d_;	//discrete model

	/*** State vector ***/
	vector x_;
	Symbolic x_sym_;

	/*** Horizon, control matrices ***/
	int Hw_, Hp_, Hu_;
	matrix Q_, R_;

	/*** Control value ***/

	std::vector<double> u_;
	//double u_[12];
	Symbolic u11_, u12_, u13_;
	Symbolic u21_, u22_, u23_;
	Symbolic u31_, u32_, u33_;
	Symbolic u41_, u42_, u43_;
	Symbolic ulaz_;
	int vel_;					//number of columns

	/*** Nlopt ***/
	//double lb_[12], ub_[12];	//lower and upper_bounds
	std::vector<double> lb_;
	std::vector<double> ub_;
    //nlopt_opt opt_;
    double minf_;				//minimum

	/*** Additional counters ***/
    int k_, ii_, kk_, br_, km_;

	/*** Total matrices ***/
    Symbolic A_d_sym_;		//A_d^i	,	i = 1, ... , Hp
    Symbolic B_d_sym_;		//(sum(A_d^i))*B_d  ,  i = 0, ... , Hp -1

	/*** State vector ***/
    Symbolic xp_;				//predict state [4, 1]
    Symbolic x_uk_;				//total predict matrix [4, Hp]

	/*** Reference vector ***/
    vector r_k_;				// dimensions:  [3*(Hp_-Hw_+1), 1]

	/*** Nonlinear matrix ***/
    matrix nelin_;				// dimensions:  [4, Hp_-Hw_+1]

	/*** Output vector ***/
    vector zp_;					//predict vector; [3, 1]
    vector z_k_;				//total predict vector; [3*(Hp_-Hw_+1), 1]
    //vector z_;					//real predict vector [3, 1]

    vector pom1_, pom2_;

	/*** Control value ***/
    vector u_k_;				//total predict vector; [4*Hu, 1]
    //vector u_real_;				//real predict vector [4, 1]

    int iter_;					//the number of iterations

    nlopt::opt opt_;
};



MPC::MPC()
{

	/*** Parameters ***/
	Km_ = 0, K_ = 0, Tm_ = 0, Ts_ = 0, d_ = 0, alpha_ = 0;

	/*** Counters ***/
	i_ = 0, j_ = 0, n_ = 0, m_ = 0;

	/*** Matrices ***/
    A_c_ = Eigen::MatrixXd::Zero(4,4);
    B_c_ = Eigen::MatrixXd::Zero(4,4);
    I_ = Eigen::MatrixXd::Zero(4,4);
    A_d_ = Eigen::MatrixXd::Zero(4,4);
    B_d_ = Eigen::MatrixXd::Zero(4,4);
    C_d_ = Eigen::MatrixXd::Zero(3,4);

	/*** State vector ***/
	x_ = Eigen::VectorXd::Zero(4);
    x_sym_ = Symbolic ("x_sym_", 4, 1);

	/*** Horizon ***/
	Hw_ = 0, Hp_ = 0, Hu_ = 0;

	/*** Control value ***/

    u_ = std::vector<double> (state_num);
	u11_ = Symbolic ("u11_"), u12_ = Symbolic ("u12_"), u13_ = Symbolic ("u13_");
	u21_ = Symbolic ("u21_"), u22_ = Symbolic ("u22_"), u23_ = Symbolic ("u23_");
	u31_ = Symbolic ("u31_"), u32_ = Symbolic ("u32_"), u33_ = Symbolic ("u33_");
	u41_ = Symbolic ("u41_"), u42_ = Symbolic ("u42_"), u43_ = Symbolic ("u43_");
    ulaz_ = ((u11_, u12_, u13_),
             (u21_, u22_, u23_),
             (u31_, u32_, u33_),
             (u41_, u42_, u43_));
    vel_ = ulaz_.columns();

	/*** Nlopt ***/
    lb_ = std::vector<double> (state_num);
    ub_ = std::vector<double> (state_num);

    minf_ = 0;

	/*** Additional counters ***/
    k_ = 0, ii_ = 0, kk_ = 0, br_ = 0, km_ = 0;

	/*** Total matrices ***/
    A_d_sym_ = Symbolic ("A_d_sym_", 4, 4);
    B_d_sym_ = Symbolic ("A_d_sym_", 4, 4);

	/*** State vector ***/
    xp_ = Symbolic ("xp_", 4, 1);

	/*** Output vector ***/
	zp_ = Eigen::VectorXd::Zero(3);
	//z_ = Eigen::VectorXd::Zero(3);

	/*** Control value ***/
	//u_real_ = Eigen::VectorXd::Zero(4);

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

void MPC::setReference(matrix ref){

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

void MPC::getSymbolicMatrices(){


	for (i_ = 0; i_ < 4; i_++ ){
		for (j_ = 0; j_ < 4; j_++){
			A_d_sym_(i_,j_) = A_d_(i_,j_);
			B_d_sym_(i_,j_) = B_d_(i_,j_);
		}
	}
}

void MPC::setModelState(vector x){

	x_ = x;

}

void MPC::setHorizon(int Hw, int Hp, int Hu){

	Hw_ = Hw;
	Hp_ = Hp;
	Hu_ = Hu;

}

void MPC::setMatrixDimension(){

    Q_ = Eigen::MatrixXd::Zero(3*(Hp_-Hw_+1), 3*(Hp_-Hw_+1));
    R_ = Eigen::MatrixXd::Zero(4*Hu_, 4*Hu_);
	x_uk_ = Symbolic ("x_uk_", 4, Hp_);
	r_k_ = Eigen::VectorXd::Zero(3*(Hp_-Hw_+1));
    nelin_ = Eigen::MatrixXd::Zero(4, Hp_-Hw_+1);
	z_k_ = Eigen::VectorXd::Zero(3*(Hp_-Hw_+1));
    u_k_ = Eigen::VectorXd::Zero(4*Hu_);


}

void MPC::setControlMatrices(matrix Q, matrix R){

	Q_ = Q;
	R_ = R;

}

void MPC::setInitialState(double *u){

	for (i_ = 0; i_ < state_num ; i_++) {
		u_[i_] = *(u + i_);
	}

}

void MPC::setNlopt(){

	/*
    opt_ = nlopt_create(NLOPT_LN_COBYLA, 12); // algorithm and dimensionality //

    nlopt_set_xtol_rel(opt_, 1e-3);
    nlopt_set_ftol_rel(opt_, 1e-3);
	*/

	opt_ = nlopt::opt(nlopt::LN_COBYLA, state_num);			// algorithm and dimensionality //

	opt_.set_xtol_rel(1e-3);
	opt_.set_ftol_rel(1e-3);


}

void MPC::setBounds(double *lb, double *ub){

	for (i_ = 0; i_ < state_num ; i_++) {
		lb_[i_] = *(lb + i_);
		ub_[i_] = *(ub + i_);
	}

	opt_.set_lower_bounds(lb_);
	opt_.set_upper_bounds(ub_);

    //nlopt_set_lower_bounds(opt_, lb_);
    //nlopt_set_upper_bounds(opt_, ub_);


}

void MPC::getSymbolicState(){

	xp_(0) = x_(0);
	xp_(1) = x_(1);
	xp_(2) = x_(2);
	xp_(3) = x_(3);


}

void MPC::getPredictModelState(){



    xp_ = A_d_sym_ * xp_ + B_d_sym_ * ulaz_.column(kk_);

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


double MPC::myfunc(const std::vector<double> &u, std::vector<double> &grad, void *my_func_data){
	//double myfunc(unsigned dd, const double *u, double *grad, void *my_func_data){

		/********** Define x_uk_ matrix **********/

		Symbolic x_uk_pom;
		//x_uk_pom = Symbolic ("x_uk_pom", 4, controller.Hp_);
		x_uk_pom = x_uk_;

		/********** Substitute control variable in state matrix **********/

		x_uk_pom = substState(&u[0], x_uk_pom);

		/*
	    cout <<"Dimensions x_uk: "<< controller.x_uk_.rows() << " x " << controller.x_uk_.columns() << "\n";
	    cout << "x_uk_pom = " << x_uk_pom << "\n";
		*/

		/********** Get predictive output (Hw to Hp) **********/

		getPredictOutput(x_uk_pom);

		/*
	    cout <<"Dimensions nelin: "<< controller.nelin_.rows() << " x " << controller.nelin_.cols() << "\n";
	    cout << "nelin = \n" << controller.nelin_ << "\n";

	    cout <<"Dimensions zp: "<< controller.zp_.rows() << " x " << controller.zp_.cols() << "\n";
	    cout << "zp = \n" << controller.zp_ << "\n";

	    cout <<"Dimensions z_k: "<< controller.z_k_.rows() << " x " << controller.z_k_.cols() << "\n";
	    cout << "z_k = \n" << controller.z_k_ << "\n";
		*/

		/********** Substitute control variable in cost **********/

		substCost();

		/********** Get control window **********/

		getControlWindow(&u[0]);

		/*
	    cout <<"Dimensions u_k: "<< controller.u_k_.rows() << " x " << controller.u_k_.cols() << "\n";
	    cout << "u_k = \n" << controller.u_k_ << "\n";
		*/

		/*****************************************/

	    cout << "Iteracija br.: " << iter_ << "\n";
	    iter_++;

		/********** Get cost and return value **********/

		return getCost();
}


Symbolic MPC::substState(const double *u, Symbolic x_uk_pom){

	 Equations rules = (u11_ == u[0*vel_+0], u12_ == u[0*vel_+1], u13_ == u[0*vel_+2], u21_ == u[1*vel_+0], u22_ == u[1*vel_+1],u23_ == u[1*vel_+2], u31_ == u[2*vel_+0], u32_ == u[2*vel_+1],u33_ == u[2*vel_+2], u41_ == u[3*vel_+0], u42_ == u[3*vel_+1], u43_ == u[3*vel_+2]);

	 x_uk_pom = x_uk_pom.subst_all(rules);

	return x_uk_pom;

}


void MPC::getPredictOutput(Symbolic x_uk_pom){

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

	    	nelin_(n_,m_-Hw_) = K_ * (abs((double)(x_uk_pom(n_,m_-1)))) * ((double)(x_uk_pom(n_,m_-1)));

	    	}

    	zp_ = C_d_ * nelin_.col(m_-Hw_);

    	/*
    	z_uk_(0,m_-Hw_) = z_(0);
    	z_uk_(1,m_-Hw_) = z_(1);
    	z_uk_(2,m_-Hw_) = z_(2);
    	*/

    	z_k_(br_,0) = zp_(0);
    	z_k_(br_+1,0) = zp_(1);
    	z_k_(br_+2,0) = zp_(2);
    	br_ = br_ + 3;
	}

}


void MPC::substCost(){

    pom1_ = ((z_k_.transpose())-(r_k_.transpose()))*Q_*(z_k_-r_k_);

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

    return ((double)(pom1_(0)+pom2_(0)));

}

Eigen::VectorXd MPC::getControlValue(const double *u){

	vector u_real;
	u_real = Eigen::VectorXd::Zero(4);

    u_real(0) = u[0*vel_+0];
    u_real(1) = u[1*vel_+0];
    u_real(2) = u[2*vel_+0];
    u_real(3) = u[3*vel_+0];

    return u_real;

}

Eigen::VectorXd MPC::getModelState(vector u_real){

    x_ = A_d_ * x_ + B_d_ * u_real;

	return x_;
}

Eigen::VectorXd MPC::getModelOutput(vector x_real){

	nelin_(0,0) = K_ * (abs((double)(x_real(0)))) * ((double)(x_real(0)));
	nelin_(1,0) = K_ * (abs((double)(x_real(1)))) * ((double)(x_real(1)));
	nelin_(2,0) = K_ * (abs((double)(x_real(2)))) * ((double)(x_real(2)));
	nelin_(3,0) = K_ * (abs((double)(x_real(3)))) * ((double)(x_real(3)));

	vector z_;
	z_ = C_d_ * nelin_.col(0);

	return z_;
}


////////////////////////////////////////


////////////////////////////////////////

#endif /* MPC_HPP_ */
