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

	void setHorizon(int Hw, int Hp, int Hu);

	void setMatrixDimension();

	void setControlMatrices(matrix Q, matrix R);

	void setInitialState(double *tau);

	void setNlopt();

	void setBounds(double *lb, double *ub);

	void getPredictModelState();

	void getReferenceWindow();

	double myfunc(const std::vector<double> &tau, std::vector<double> &grad, void *my_func_data);

	Symbolic substState(const double *tau, Symbolic z_uk_pom);

	void getPredictOutput(Symbolic z_uk_pom);

	void substCost();

	void getControlWindow(const double *u);

	double getCost();

	vector getControlValue(const double *tau);

	vector getModelOutput(vector tau_real);

	void setVelocityParameters(double beta_uu, double beta_vv, double beta_rr, double a_u, double a_v, double a_r);

	void setInitialVelocity(vector velocity);

	void setInitialPosition(vector position);

	void getAdditionalVelocity();

	void getAdditionalPosition();

	void getPredictVelocity();

	void getPredictPosition();

	vector getRealVelocity(vector z_real);

	vector getRealPosition();


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

	/*** Horizon, control matrices ***/
	int Hw_, Hp_, Hu_;
	matrix Q_, R_;

	/*** Control value ***/

	std::vector<double> tau_;
	//double u_[12];
	Symbolic tau11_, tau12_, tau13_;
	Symbolic tau21_, tau22_, tau23_;
	Symbolic tau31_, tau32_, tau33_;
	Symbolic u41_, u42_, u43_;
	Symbolic tau_ulaz_;
	int vel_;					//number of columns

	/*** Nlopt ***/
	//double lb_[12], ub_[12];	//lower and upper_bounds
	std::vector<double> lb_;
	std::vector<double> ub_;
    //nlopt_opt opt_;
    double minf_;				//minimum

	/*** Additional counters ***/
    int k_, ii_, kk_, br_, km_;

	/*** Symbolic matrices ***/
    Symbolic A_d_sym_;
    Symbolic B_d_sym_;
    Symbolic C_d_sym_;

	/*** State vector ***/
    Symbolic z_p_;				//predict state [4, 1]
    Symbolic z_uk_;				//total predict matrix [4, Hp]

	/*** Reference vector ***/
    vector r_k_;				// dimensions:  [3*(Hp_-Hw_+1), 1]


	/*** Output vector ***/
    vector z_k_;				//total predict vector; [3*(Hp_-Hw_+1), 1]
    //vector z_;					//real predict vector [3, 1]

    vector pom1_, pom2_;

	/*** Control value ***/
    vector tau_k_;				//total predict vector; [4*Hu, 1]
    //vector tau_real_;				//real predict vector [4, 1]

    int iter_;					//the number of iterations

    nlopt::opt opt_;

    /*** Velocity parameters ***/
	double beta_uu_, beta_vv_, beta_rr_ , alpha_u_, alpha_v_, alpha_r_;

    /*** Velocity  ***/
	vector velocity_;
	vector velocity_p1_;
	vector velocity_p2_;

    /*** Position  ***/
	vector position_;
	vector position_p1_;
	vector position_p2_;

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

	/*** Horizon ***/
	Hw_ = 0, Hp_ = 0, Hu_ = 0;

	/*** Control value ***/

    tau_ = std::vector<double> (state_num);
	tau11_ = Symbolic ("tau11_"), tau12_ = Symbolic ("tau12_"), tau13_ = Symbolic ("tau13_");
	tau21_ = Symbolic ("tau21_"), tau22_ = Symbolic ("tau22_"), tau23_ = Symbolic ("tau23_");
	tau31_ = Symbolic ("tau31_"), tau32_ = Symbolic ("tau32_"), tau33_ = Symbolic ("tau33_");
	u41_ = Symbolic ("u41_"), u42_ = Symbolic ("u42_"), u43_ = Symbolic ("u43_");
    tau_ulaz_ = ((tau11_, tau12_, tau13_),
                 (tau21_, tau22_, tau23_),
                 (tau31_, tau32_, tau33_),
                 (u41_, u42_, u43_));
    vel_ = tau_ulaz_.columns();


	/*** Nlopt ***/
    lb_ = std::vector<double> (state_num);
    ub_ = std::vector<double> (state_num);

    minf_ = 0;

	/*** Additional counters ***/
    k_ = 0, ii_ = 0, kk_ = 0, br_ = 0, km_ = 0;

	/*** Total matrices ***/
    A_d_sym_ = Symbolic ("A_d_sym_", 4, 4);
    B_d_sym_ = Symbolic ("A_d_sym_", 4, 4);
    C_d_sym_ = Symbolic ("C_d_sym_", 3, 4);

	/*** Output vector ***/
    z_p_ = Symbolic ("z_p_", 3, 1);

	//zp_ = Eigen::VectorXd::Zero(3);
	//z_ = Eigen::VectorXd::Zero(3);

	/*** Control value ***/
	//tau_real_ = Eigen::VectorXd::Zero(4);

    iter_ = 1;

    /*** Velocity parameters ***/
	beta_uu_ = 0, beta_vv_ = 0, beta_rr_ = 0;
	alpha_u_ = 0, alpha_v_ = 0, alpha_r_ = 0;

    /*** Velocity  ***/
	velocity_ = Eigen::VectorXd::Zero(3);
	velocity_p1_ = Eigen::VectorXd::Zero(3);
	velocity_p2_ = Eigen::VectorXd::Zero(3);

    /*** Position  ***/
	position_ = Eigen::VectorXd::Zero(3);
	position_p1_ = Eigen::VectorXd::Zero(3);
	position_p2_ = Eigen::VectorXd::Zero(3);

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
			if (i_ != 3) C_d_sym_(i_,j_) = C_d_(i_,j_);
		}
	}
}

void MPC::setHorizon(int Hw, int Hp, int Hu){

	Hw_ = Hw;
	Hp_ = Hp;
	Hu_ = Hu;

}

void MPC::setMatrixDimension(){

    Q_ = Eigen::MatrixXd::Zero(3*(Hp_-Hw_+1), 3*(Hp_-Hw_+1));
    R_ = Eigen::MatrixXd::Zero(4*Hu_, 4*Hu_);
	z_uk_ = Symbolic ("z_uk_", 3, Hp_);
	r_k_ = Eigen::VectorXd::Zero(3*(Hp_-Hw_+1));
	z_k_ = Eigen::VectorXd::Zero(3*(Hp_-Hw_+1));
    tau_k_ = Eigen::VectorXd::Zero(4*Hu_);

}

void MPC::setControlMatrices(matrix Q, matrix R){

	Q_ = Q;
	R_ = R;

}

void MPC::setInitialState(double *tau){

	for (i_ = 0; i_ < 12 ; i_++) {
		tau_[i_] = *(tau + i_);
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

void MPC::getPredictModelState(){

    z_p_ = C_d_sym_ * tau_ulaz_.column(kk_);

	z_uk_(0,i_) = z_p_(0);
	z_uk_(1,i_) = z_p_(1);
	z_uk_(2,i_) = z_p_(2);

}

void MPC::getReferenceWindow(){

    //vektor dimenzija 1x3*(Hp-Hw+1)
	br_ = 0;
    for (j_ = Hw_; j_ <= Hp_; j_++){
            for (i_ = 0; i_ < 3; i_++){
                    r_k_(br_,0) = ref_(i_,j_); //umjesto k_ ide 1
                    br_++;
            }
    }

}


double MPC::myfunc(const std::vector<double> &tau, std::vector<double> &grad, void *my_func_data){
	//double myfunc(unsigned dd, const double *u, double *grad, void *my_func_data){

		/********** Define x_uk_ matrix **********/

		Symbolic z_uk_pom;
		z_uk_pom = z_uk_;

		/********** Substitute control variable in state matrix **********/

		z_uk_pom = substState(&tau[0], z_uk_pom);

		/*
	    cout <<"Dimensions z_uk: "<< controller.z_uk_.rows() << " x " << controller.z_uk_.columns() << "\n";
	    cout << "z_uk_pom = " << z_uk_pom << "\n";
		 */

		/********** Get predictive output (Hw to Hp) **********/

		getPredictOutput(z_uk_pom);

		/*
	    cout <<"Dimensions z_p: "<< controller.z_p_.rows() << " x " << controller.z_p_.cols() << "\n";
	    cout << "z_p = \n" << controller.zp_ << "\n";

	    cout <<"Dimensions z_k: "<< controller.z_k_.rows() << " x " << controller.z_k_.cols() << "\n";
	    cout << "z_k = \n" << controller.z_k_ << "\n";
		*/

		/********** Substitute control variable in cost **********/

		substCost();

		/********** Get control window **********/

		getControlWindow(&tau[0]);

		/*
	    cout <<"Dimensions tau_k: "<< controller.tau_k_.rows() << " x " << controller.tau_k_.cols() << "\n";
	    cout << "tau_k = \n" << controller.tau_k_ << "\n";
		*/

		/*****************************************/

	    cout << "Iteracija br.: " << iter_ << "\n";
	    iter_++;

		/********** Get cost and return value **********/

		return getCost();
}


Symbolic MPC::substState(const double *tau, Symbolic z_uk_pom){


	 Equations rules = (tau11_ == tau[0*vel_+0], tau12_ == tau[0*vel_+1], tau13_ == tau[0*vel_+2], tau21_ == tau[1*vel_+0], tau22_ == tau[1*vel_+1],tau23_ == tau[1*vel_+2], tau31_ == tau[2*vel_+0], tau32_ == tau[2*vel_+1],tau33_ == tau[2*vel_+2], u41_ == tau[3*vel_+0], u42_ == tau[3*vel_+1], u43_ == tau[3*vel_+2]);
	 z_uk_pom = z_uk_pom.subst_all(rules);


	return z_uk_pom;

}


void MPC::getPredictOutput(Symbolic z_uk_pom){

	getAdditionalVelocity();
    getAdditionalPosition();

	br_ = 0;
	for (m_ = 0; m_ < Hp_; m_++){


    	z_p_ = z_uk_pom.column(m_);


    	getPredictVelocity();
    	getPredictPosition();

    	if((Hw_ - 1) <= m_){

    		z_k_(br_,0) = position_p2_(0);
    		z_k_(br_ + 1,0) = position_p2_(1);
    		z_k_(br_ + 2,0) = position_p2_(2);
    		br_ = br_ + 3;

    	}

    	velocity_p1_ = velocity_p2_;
    	position_p1_ = position_p2_;

	}
	//cout << "z_k_ = \n" << z_k_ << "\n";
	//exit(2);

}



void MPC::substCost(){

    pom1_ = ((z_k_.transpose())-(r_k_.transpose()))*Q_*(z_k_-r_k_);

}

void MPC::getControlWindow(const double *tau){

    br_ = 0;
    for (m_ = 0; m_ < Hu_; m_++){
    	for (n_ = 0; n_ < 4; n_++){
    		if (m_ > 2) km_ = 2;
    		else km_ = m_;

    		tau_k_(br_,0) = tau[n_*vel_+km_];
            br_++;
        }
    }

}


double MPC::getCost(){

    pom2_ = (tau_k_.transpose())*R_*(tau_k_);

    return ((double)(pom1_(0)+pom2_(0)));

}

Eigen::VectorXd MPC::getControlValue(const double *tau){

	vector tau_real;
	tau_real = Eigen::VectorXd::Zero(4);

    tau_real(0) = tau[0*vel_+0];
    tau_real(1) = tau[1*vel_+0];
    tau_real(2) = tau[2*vel_+0];
    tau_real(3) = tau[3*vel_+0];

    return tau_real;

}

Eigen::VectorXd MPC::getModelOutput(vector tau_real){

	vector z_real;

	z_real = C_d_ * tau_real;

	return z_real;
}

void MPC::setVelocityParameters(double beta_uu, double beta_vv, double beta_rr, double alpha_u, double alpha_v, double alpha_r){

	beta_uu_ = beta_uu;
	beta_vv_ = beta_vv;
	beta_rr_ = beta_rr;
	alpha_u_ =  alpha_u;
	alpha_v_ = alpha_v;
	alpha_r_ = alpha_r;

}

void MPC::setInitialVelocity(vector velocity){

	velocity_ = velocity;

}

void MPC::setInitialPosition(vector position){

	position_ = position;

}

void MPC::getAdditionalVelocity(){

	velocity_p1_ = velocity_;

}

void MPC::getAdditionalPosition(){

	position_p1_ = position_;

}

void MPC::getPredictVelocity(){

	velocity_p2_(0) = velocity_p1_(0) - beta_uu_ / alpha_u_ * (abs((double)(velocity_p1_(0)))) * ((double)(velocity_p1_(0))) * Ts_ + 1 / alpha_u_ * z_p_(0) * Ts_;
	velocity_p2_(1) = velocity_p1_(1) - beta_vv_ / alpha_v_ * (abs((double)(velocity_p1_(1)))) * ((double)(velocity_p1_(1))) * Ts_ + 1 / alpha_v_ * z_p_(1) * Ts_;
	velocity_p2_(2) = velocity_p1_(2) - beta_rr_ / alpha_r_ * (abs((double)(velocity_p1_(2)))) * ((double)(velocity_p1_(2))) * Ts_ + 1 / alpha_r_ * z_p_(2) * Ts_;


}

void MPC::getPredictPosition(){


	position_p2_(0) = position_p1_(0) + velocity_p1_(0) * cos((double)(position_p1_(2))) * Ts_ - velocity_p1_(1) * sin((double)(position_p1_(2))) * Ts_;
	position_p2_(1) = position_p1_(1) + velocity_p1_(0) * sin((double)(position_p1_(2))) * Ts_ + velocity_p1_(1) * cos((double)(position_p1_(2))) * Ts_;
	position_p2_(2) = position_p1_(2) + velocity_p1_(2) * Ts_;

}

Eigen::VectorXd MPC::getRealVelocity(vector z_real){

	velocity_(0) = velocity_p1_(0) - beta_uu_ / alpha_u_ * (abs((double)(velocity_p1_(0)))) * ((double)(velocity_p1_(0))) * Ts_ + 1 / alpha_u_ * z_real(0) * Ts_;
	velocity_(1) = velocity_p1_(1) - beta_vv_ / alpha_v_ * (abs((double)(velocity_p1_(1)))) * ((double)(velocity_p1_(1))) * Ts_ + 1 / alpha_v_ * z_real(1) * Ts_;
	velocity_(2) = velocity_p1_(2) - beta_rr_ / alpha_r_ * (abs((double)(velocity_p1_(2)))) * ((double)(velocity_p1_(2))) * Ts_ + 1 / alpha_r_ * z_real(2) * Ts_;

	return velocity_;
}

Eigen::VectorXd MPC::getRealPosition(){


	position_(0) = position_p1_(0) + velocity_p1_(0) * cos((double)(position_p1_(2))) * Ts_ - velocity_p1_(1) * sin((double)(position_p1_(2))) * Ts_;
	position_(1) = position_p1_(1) + velocity_p1_(0) * sin((double)(position_p1_(2))) * Ts_ + velocity_p1_(1) * cos((double)(position_p1_(2))) * Ts_;
	position_(2) = position_p1_(2) + velocity_p1_(2) * Ts_;

	return position_;
}


////////////////////////////////////////


////////////////////////////////////////

#endif /* MPC_HPP_ */
