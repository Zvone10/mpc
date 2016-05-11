#include <mpc/mpc.hpp>
//#include <ros/ros.h>
//#include <symbolicc++.h>


/*********************************************************************
 ***  Main function
 ********************************************************************/

MPC controller;

double myvfunc(unsigned dd, const double *u, double *grad, void *my_func_data){

	/********** Substitute control variable in state matrix **********/

	controller.substState(&u[0]);

	 /*
    cout <<"Dimensions x_uk: "<< controller.x_uk_.rows() << " x " << controller.x_uk_.columns() << "\n";
    cout << "x_uk = " << controller.x_uk_ << "\n";
	 */

	/********** Get predictive output (Hw to Hp) **********/

	controller.getPredictOutput();

	/*
    cout <<"Dimensions nelin: "<< controller.nelin_.rows() << " x " << controller.nelin_.columns() << "\n";
    cout << "nelin = " << controller.nelin_ << "\n";

    cout <<"Dimensions z: "<< controller.z_.rows() << " x " << controller.z_.columns() << "\n";
    cout << "z = " << controller.z_ << "\n";

    cout <<"Dimensions z_k: "<< controller.z_k_.rows() << " x " << controller.z_k_.columns() << "\n";
    cout << "z_k = " << controller.z_k_ << "\n";
	*/

	/********** Substitute control variable in cost **********/

	controller.substCost(&u[0]);

	/*
    cout <<"Dimensions pom1: "<< controller.pom1_.rows() << " x " << controller.pom1_.columns() << "\n";
    cout << "pom1 = " << controller.pom1_ << "\n";
	*/

	/********** Get control window **********/

	controller.getControlWindow(&u[0]);

	/*
    cout <<"Dimensions u_k: "<< controller.u_k_.rows() << " x " << controller.u_k_.columns() << "\n";
    cout << "u_k = " << controller.u_k_ << "\n";
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

	int kor = 30;
    Symbolic ref("ref", 3, kor);

    for (i = 0; i < 3; i++){
            for (j = 0; j < kor; j++){
                    if(i == 0) {
                    	if (j < 16) ref(i,j) = 8 * sin(pi*j/1/kor);
                    	else ref(i,j) = 8;
                    }
                    else ref(i,j) = 0;
            }
    }

	controller.setReference(ref);

	//cout << controller.ref_;

	/********** Set model **********/

	controller.setModel();

	/*
	cout << "A_c = " << controller.A_c_ << "\n";
	cout << "B_c = " << controller.B_c_ << "\n";
	cout << "I = " << controller.I_ << "\n";
	*/

	/********** Set discrete model **********/

	controller.setDiscreteModel();

	/*
	cout << "A_d = " << controller.A_d_ << "\n";
	cout << "B_d = " << controller.B_d_ << "\n";
	cout << "C_d = " << controller.C_d_ << "\n";
	*/

	/********** Set model state **********/

	Symbolic x("x", 4, 1);

	x(0,0) = 600;
	x(1,0) = 500;
	x(2,0) = -300;
	x(3,0) = -400;

	controller.setModelState(x);

	//cout << "x = " << controller.x_ << "\n";

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

	Symbolic Q("Q", 3*(Hp-Hw+1), 3*(Hp-Hw+1));
	Symbolic R("R", 4*Hu, 4*Hu);

    for (i = 0; i < 3*(Hp-Hw+1); i++){
            for (j = 0; j < 3*(Hp-Hw+1); j++){
                    if (i == j) Q(i,j) = 1;
                    else Q(i,j) = 0;
            }
    }

    for (i = 0; i < 4*Hu; i++){
            for (j = 0; j < 4*Hu; j++){
                    if (i == j) R(i,j) = 10;
                    else R(i,j) = 0;
            }
    }

    controller.setControlMatrices(Q, R);

    /*
    cout << "Q = " << controller.Q_ << "\n";
    cout << "R = " << controller.R_ << "\n";

    cout << "Dimensions Q: " << controller.Q_.rows() << " x " << controller.Q_.columns() << "\n";
    cout << "Dimensions R: " << controller.R_.rows() << " x " << controller.R_.columns() << "\n";
	*/

	/********** Set initial state **********/

    // {u11, u12, u13, u21, u22, u23, u31, u32, u33, u41, u42, u43} //
    double u[12] = {0};
    //double u[12] = {11, 12, 13, 21, 22, 23, 31, 32, 33, 41, 42, 43};

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

    int k, counter;

    counter = kor - Hp;

    controller.k_ = 0;



    for (controller.k_ = 0; controller.k_ < counter; controller.k_++){

    	 for (controller.i_ = 0; controller.i_ < Hp; controller.i_++){//x(k+1|k) do x(k+Hp|k)

    		 if(controller.i_  >= 3) controller.kk_ = 2;
    		 else controller.kk_  = controller.i_ ;

    		 controller.getTotalMatrices();

    		 /*
             cout << "A_d_uk1 =" << controller.A_d_uk1_ << "\n";
             cout << "A_d_uk2 =" << controller.A_d_uk2_ << "\n";
             */

       		 controller.getPredictModelState();

     		 /*
             cout <<"Dimensions x_uk: "<< controller.x_uk_.rows() << " x " << controller.x_uk_.columns() << "\n";
             cout << "x_uk = " << controller.x_uk_ << "\n";
             cout <<"Dimensions xp: "<< controller.xp_.rows() << " x " << controller.xp_.columns() << "\n";
             cout << "xp = " << controller.xp_ << "\n";
			 */

    	 }

    	 controller.getReferenceWindow();

    	 /*
         cout <<"Dimensions r_k: "<< controller.r_k_.rows() << " x " << controller.r_k_.columns() << "\n";
         cout << "r_k = " << controller.r_k_ << "\n";
		 */


         cout << "\n" << "Korak br.: " << controller.k_+1 << "\n";

    	 /********** Function **********/

         nlopt_set_min_objective(controller.opt_, myvfunc, NULL);


         if (nlopt_optimize(controller.opt_, &controller.u_[0], &controller.minf_) < 0) {
             cout <<"Nlopt failed!" << "\n";
         }
         else {
             cout <<"Found minimum: " << controller.minf_ << "\n";
         }



         controller.iter_ = 1;
    	 /********************************/



     	/********** Get real control value **********/

         Symbolic u_real("u_real", 4, 1);

         u_real = controller.getControlValue(&controller.u_[0]);

         //cout << "u_real = " << u_real << "\n";

      	/********** Get real model state **********/

         x = controller.getModelState(u_real);

         //cout << "x = " << x << "\n";

       	/********** Get real model state **********/

         Symbolic z("z", 3, 1);

         z = controller.getModelOutput(x);

         //cout << "z = " << z << "\n";



    }//main counter







	/********************************/


	cout << "\n!!STOP(-971.967)!!\n";
	return 0;

}

