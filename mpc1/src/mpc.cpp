#include <mpc/mpc.hpp>
//#include <ros/ros.h>


/*********************************************************************
 ***  Main function
 ********************************************************************/

MPC controller;

double myfunc(const std::vector<double> &u, std::vector<double> &grad, void *my_func_data){
//double myfunc(unsigned dd, const double *u, double *grad, void *my_func_data){

	/********** Define x_uk_ matrix **********/

	Symbolic x_uk_pom;
	//x_uk_pom = Symbolic ("x_uk_pom", 4, controller.Hp_);
	x_uk_pom = controller.x_uk_;

	/********** Substitute control variable in state matrix **********/

	x_uk_pom = controller.substState(&u[0], x_uk_pom);

	/*
    cout <<"Dimensions x_uk: "<< controller.x_uk_.rows() << " x " << controller.x_uk_.columns() << "\n";
    cout << "x_uk_pom = " << x_uk_pom << "\n";
	*/

	/********** Get predictive output (Hw to Hp) **********/

	controller.getPredictOutput(x_uk_pom);

	/*
    cout <<"Dimensions nelin: "<< controller.nelin_.rows() << " x " << controller.nelin_.cols() << "\n";
    cout << "nelin = \n" << controller.nelin_ << "\n";

    cout <<"Dimensions zp: "<< controller.zp_.rows() << " x " << controller.zp_.cols() << "\n";
    cout << "zp = \n" << controller.zp_ << "\n";

    cout <<"Dimensions z_k: "<< controller.z_k_.rows() << " x " << controller.z_k_.cols() << "\n";
    cout << "z_k = \n" << controller.z_k_ << "\n";
	*/

	/********** Substitute control variable in cost **********/

	controller.substCost();

	/********** Get control window **********/

	controller.getControlWindow(&u[0]);

	/*
    cout <<"Dimensions u_k: "<< controller.u_k_.rows() << " x " << controller.u_k_.cols() << "\n";
    cout << "u_k = \n" << controller.u_k_ << "\n";
	*/

	/*****************************************/

    cout << "Iteracija br.: " << controller.iter_ << "\n";
    controller.iter_++;

	/********** Get cost and return value **********/

	return controller.getCost();
}


int main(int argc, char** argv)
{

    int Km = 1000;
    float K = 1e-5, Tm = 50e-3, Ts = 0.1, d = 0.5, alpha = pi/4;
	int i = 0, j = 0;

	typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrix;
	typedef Eigen::Matrix<double, Eigen::Dynamic, 1> vector;


	//ros::init(argc, argv, "mpc");
	//MPC controller;
	//ros::spin();


	/********** Set parameters **********/

	controller.setParameters(Km, K, Tm, Ts, d, alpha);

	/*
	cout << "Km = " << controller.Km_ << "\n";
	cout << "K = " << controller.K_ << "\n";
	cout << "Tm = " << controller.Tm_ << "\n";
	cout << "Ts = " << controller.Ts_ << "\n";
	cout << "d = " << controller.d_ << "\n";
	cout << "alpha = " << controller.alpha_ << "\n";
	*/

	/********** Set reference **********/

	int kor = 50;
    matrix ref;
    ref = Eigen::MatrixXd::Zero(3, kor);


    for (i = 0; i < 3; i++){
            for (j = 0; j < kor; j++){
                    if(i == 0) {
                    	//if (j < 24) ref(i,j) =  8 * sin(pi*j/1/kor);
                    	//else ref(i,j) = 8;
                    	ref(i,j) = 5;           		//referent value
                    }
            }
    }

	controller.setReference(ref);

	/*
	cout << "ref = \n" << controller.ref_ << "\n";
	cout << "Dimensions ref: " << controller.ref_.rows() << " x " << controller.ref_.cols() << "\n";
	*/

	/********** Set model **********/

	controller.setModel();

	/*
	cout << "A_c = \n" << controller.A_c_ << "\n";
	cout << "B_c = \n" << controller.B_c_ << "\n";
	cout << "I = \n" << controller.I_ << "\n";
	*/

	/********** Set discrete model **********/

	controller.setDiscreteModel();

	/*
	cout << "A_d = \n" << controller.A_d_ << "\n";
	cout << "B_d = \n" << controller.B_d_ << "\n";
	cout << "C_d = \n" << controller.C_d_ << "\n";
	*/

	/********** Set model state **********/

	vector x;
	x = Eigen::VectorXd::Zero(4);

	x(0) = 400;
	x(1) = 400;
	x(2) = -400;
	x(3) = -400;


	controller.setModelState(x);

	//cout << "x = \n" << controller.x_ << "\n";

	/********** Set horizon **********/

	int Hw, Hp, Hu;

	Hw = 1;
	Hp = 4;
	Hu = 4;

	controller.setHorizon(Hw, Hp, Hu);

	/*
	cout << "Hw = " << controller.Hw_ << "\n";
	cout << "Hp = " << controller.Hp_ << "\n";
	cout << "Hu = " << controller.Hu_ << "\n";
	*/

	/********** Set matrix dimension **********/

	controller.setMatrixDimension();

	/********** Set control matrices **********/

	matrix Q, R;
    Q = Eigen::MatrixXd::Zero(3*(Hp-Hw+1), 3*(Hp-Hw+1));
    R = Eigen::MatrixXd::Zero(4*Hu, 4*Hu);

    for (i = 0; i < 3*(Hp-Hw+1); i++){
            for (j = 0; j < 3*(Hp-Hw+1); j++){
                    if (i == j) Q(i,j) = 10;		//value of Q(i,i)
            }
    }

    for (i = 0; i < 4*Hu; i++){
            for (j = 0; j < 4*Hu; j++){
                    if (i == j) R(i,j) = 1;			//value of R(i,i)
            }
    }

    controller.setControlMatrices(Q, R);

    /*
	cout << "Dimensions Q: " << controller.Q_.rows() << " x " << controller.Q_.cols() << "\n";
    cout << "Dimensions R: " << controller.R_.rows() << " x " << controller.R_.cols() << "\n";

    cout << "Q = \n" << controller.Q_ << "\n";
    cout << "R = \n" << controller.R_ << "\n";
     */

	/********** Set initial state **********/

    ////// {u11, u12, u13, u21, u22, u23, u31, u32, u33, u41, u42, u43} //////

    //double u[12] = {0};
    //double u[12] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1};
    double u[12] = {0.4, 0.4, 0.4, 0.4, 0.4, 0.4, -0.4, -0.4, -0.4, -0.4, -0.4, -0.4};

    controller.setInitialState(&u[0]);

    /*
    cout << "u0 = \n";
    for (i = 0; i < 4; i++){
            for (j = 0; j < 3; j++){
                    cout << controller.u_[i*3+j] << "  ";
            }
            cout << "\n";
    }
	*/

	/********** Set NLopt **********/

    controller.setNlopt();

	/********** Set Bounds **********/

    // lower bounds //
    double lb[12] = { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};

    // upper bounds //
    double ub[12] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    controller.setBounds(&lb[0], &ub[0]);

    /*
    cout << "lb = \n";
    for (i = 0; i < 4; i++){
            for (j = 0; j < 3; j++){
                    cout << controller.lb_[i*3+j] << "  ";
            }
            cout << "\n";
    }

    cout << "ub = \n";
    for (i = 0; i < 4; i++){
            for (j = 0; j < 3; j++){
                    cout << controller.ub_[i*3+j] << "  ";
            }
            cout << "\n";
    }
	*/

	/********************************/

    int counter;

    counter = kor - Hp;

    controller.k_ = 0;



    for (controller.k_ = 0; controller.k_ < counter; controller.k_++){

    	for (controller.i_ = 0; controller.i_ < Hp; controller.i_++){//x(k+1|k) to x(k+Hp|k)


    		if(controller.i_  >= 3) controller.kk_ = 2;
    		else controller.kk_  = controller.i_ ;

    		controller.getTotalMatrices();

    		/*
            cout << "A_d_uk1 = \n" << controller.A_d_uk1_ << "\n";
            cout << "A_d_uk2 = \n" << controller.A_d_uk2_ << "\n";
			*/

       		controller.getPredictModelState();

       		//cout << "xp = " << controller.xp_ << "\n";

    	}

    	/*
        cout <<"Dimensions x_uk: "<< controller.x_uk_.rows() << " x " << controller.x_uk_.columns() << "\n";
        cout << "x_uk = " << controller.x_uk_ << "\n";
		*/

    	/********** Get Reference Window **********/

    	controller.getReferenceWindow();

    	/*
        cout <<"Dimensions r_k: "<< controller.r_k_.rows() << " x " << controller.r_k_.cols() << "\n";
        cout << "r_k = \n" << controller.r_k_ << "\n";
		*/

    	/********** Function **********/

        cout << "\n" << "Korak br.: " << controller.k_+1 << "\n";

        //nlopt_set_min_objective(controller.opt_, myfunc, NULL);
        //nlopt_set_min_objective(controller.opt_, controller.myfunc, NULL);

        controller.opt_.set_min_objective(myfunc, NULL);

        /*
        if (nlopt_optimize(controller.opt_, &controller.u_[0], &controller.minf_) < 0) {
            cout <<"Nlopt failed!" << "\n";
        }
        else {
            cout <<"Found minimum: " << controller.minf_ << "\n";
        }
		*/

        if ( controller.opt_.optimize(controller.u_, controller.minf_) < 0) {
            cout <<"Nlopt failed!" << "\n";
        }
        else {
            cout <<"Found minimum: " << controller.minf_ << "\n";
        }

        controller.iter_ = 1;

     	/********** Get control value **********/

        vector u_real;
        //u_real = Eigen::VectorXd::Zero(4);

        u_real = controller.getControlValue(&controller.u_[0]);

        cout << "u_real = \n" << u_real << "\n";

      	/********** Get model state **********/

        vector x_real;
        //x_real = Eigen::VectorXd::Zero(4);

        x_real = controller.getModelState(u_real);

        cout << "x = \n" << x_real << "\n";

        /********** Get model output **********/

        vector z_real;
        //z = Eigen::VectorXd::Zero(3);

        z_real = controller.getModelOutput(controller.x_);

        cout << "z = \n" << z_real << "\n";


        /********** Initial control value for next step **********/

     	for (int row = 0; row < 4; row++){
     		for (int col = 0; col < 3; col++){

     			if(col < 2) controller.u_[row*controller.vel_+col] = controller.u_[row*controller.vel_+col+1];
     			//cout << controller.u_[row*controller.vel_+col] << "  ";

     			}
     	    //cout << "\n";
     	}


     	//exit(4);

    }//main


	/********************************/


	cout << "\n!!STOP()!!\n";
	return 0;

}
