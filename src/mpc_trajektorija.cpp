#include <mpc/mpc_trajektorija.hpp>
#include <ros/ros.h>
#include <boost/bind.hpp>
#include <fstream>
#include<cmath>

#include <auv_msgs/NavSts.h>
#include <auv_msgs/BodyForceReq.h>
#include <std_msgs/Float32MultiArray.h>


ofstream myfile_ref, myfile_z, myfile_velocity, myfile_position, myfile_tau;


class MPCNode
{
public:

	enum {kor = 50};

	typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrix;
	typedef Eigen::Matrix<double, Eigen::Dynamic, 1> vector;

	MPCNode();

	~MPCNode();

	void initialize();

	void step(matrix ref);

	matrix get_reference();

	void onReference(const auv_msgs::NavSts::ConstPtr& data);

	void onStateHat(const auv_msgs::NavSts::ConstPtr& data);

	static double wrap(const std::vector<double> &x, std::vector<double> &grad, void *data)
	{
		MPC *obj = static_cast<MPC *>(data);
	    return obj->myfunc(x,grad,data);
	}

	MPC controller;

	ros::Subscriber sub_reference, sub_state;

	ros::Publisher pub_tau_i, pub_tau, pub_vel, pub_pos;

	matrix reference;

	auv_msgs::NavSts state;

};

MPCNode::MPCNode():
		reference(Eigen::MatrixXd::Zero(3, kor))
{

	ros::NodeHandle nh;

	sub_state = nh.subscribe<auv_msgs::NavSts>("position",1,&MPCNode::onStateHat,this);
	sub_reference = nh.subscribe<auv_msgs::NavSts>("stateRef",1,&MPCNode::onReference,this);

	pub_tau_i = nh.advertise<std_msgs::Float32MultiArray>("tau_i",1);
	pub_tau = nh.advertise<auv_msgs::BodyForceReq>("tau",1);
	pub_vel = nh.advertise<auv_msgs::NavSts>("velocity",1);
	pub_pos = nh.advertise<auv_msgs::NavSts>("position" ,1);

}

MPCNode::~MPCNode()
{

}

void MPCNode::initialize()
{

		myfile_z.open ("z_plot.txt");
		myfile_ref.open ("ref_plot.txt");
		myfile_velocity.open ("velocity_plot.txt");
		myfile_position.open ("position_plot.txt");
		myfile_tau.open ("tau_plot.txt");


	    int Km = 1000;
	    float K = 1.13137e-5, Tm = 50e-3, Ts = 0.1, d = 0.5, alpha = pi/4;
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

		/********** Set symbolic A_d and B_d matrices **********/

		controller.getSymbolicMatrices();

		/*
		cout << "A_d_sym = \n" << controller.A_d_sym_ << "\n";
		cout << "B_d_sym = \n" << controller.B_d_sym_ << "\n";
		cout << "C_d_sym = \n" << controller.C_d_sym_ << "\n";
		*/


		/********** Set horizon **********/

		int Hw, Hp, Hu;

		Hw = 1;
		Hp = 4;
		Hu = 3;

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
	                    if (i == j) Q(i,j) = 150;		//value of Q(i,i)
	            }
	    }

	    for (i = 0; i < 4*Hu; i++){
	            for (j = 0; j < 4*Hu; j++){
	                    if (i == j) R(i,j) = 0.01;			//value of R(i,i)
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

	    double tau[12] = {0};
	    //double tau[12] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1};
	    //double tau[12] = {0.4, 0.4, 0.4, 0.4, 0.4, 0.4, -0.4, -0.4, -0.4, -0.4, -0.4, -0.4};
	    //double tau[12] = {0.4, 0.4, 0.4, 0.1, 0.1, 0.1, -0.4, -0.4, -0.4, -0.4, -0.4, -0.4};
	    //double tau[12] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5};

	    controller.setInitialState(&tau[0]);

	    /*
	    cout << "u0 = \n";
	    for (i = 0; i < 4; i++){
	            for (j = 0; j < 3; j++){
	                    cout << controller.tau_[i*3+j] << "  ";
	            }
	            cout << "\n";
	    }
		*/

		/********** Set NLopt **********/

	    controller.setNlopt();

		/********** Set Bounds **********/

	    // lower bounds //
	    //double lb[12] = { -1, -1, -1, -0.3, -0.3, -0.3, -1, -1, -1, -1, -1, -1};
	    double lb[12] = { -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15};
	    //double lb[12] = { -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10};
	    //double lb[12] = { -10, -10, -10, -10, -10, -10, -10, -10, -10, -1, -1, -1};

	    // upper bounds //
	    //double ub[12] = { 1, 1, 1, 0.3, 0.3, 0.3, 1, 1, 1, 1, 1, 1};
	    double ub[12] = {15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15};
	    //double ub[12] = { 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10};
	    //double ub[12] = { 10, 10, 10, 10, 10, 10, 10, 10, 10, 1, 1, 1};

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


		/********** Set Velocity Parameters **********/

		double beta_uu, beta_vv, beta_rr , alpha_u, alpha_v, alpha_r;

		beta_uu = 1;
		beta_vv = 1;
		beta_rr = 1;
		alpha_u = 1;
		alpha_v = 1;
		alpha_r = 1;

		controller.setVelocityParameters(beta_uu, beta_vv, beta_rr, alpha_u, alpha_v, alpha_r);

		/********** Set Initial Velocity **********/

		vector velocity;
		velocity = Eigen::VectorXd::Zero(3);

		velocity(0) = 0;
		velocity(1) = 0;
		velocity(2) = 0;

		controller.setInitialVelocity(velocity);

		//cout << "velocity = \n" << controller.velocity_ << "\n";

		/********** Set Initial Position **********/

		vector position;
		position = Eigen::VectorXd::Zero(3);

		position(0) = 0;
		position(1) = 0;
		position(2) = 0;

		controller.setInitialPosition(position);

		//cout << "position = \n" << controller.position_ << "\n";

		/********************************/
		/********************************/


}

void MPCNode::step(matrix ref)
{

    int Hp = controller.Hp_;


	controller.setReference(ref);

	/*
	cout << "ref = \n" << controller.ref_ << "\n";
	cout << "Dimensions ref: " << controller.ref_.rows() << " x " << controller.ref_.cols() << "\n";
	*/

    	for (controller.i_ = 0; controller.i_ < Hp; controller.i_++){//x(k+1|k) to x(k+Hp|k)


    		if(controller.i_  >= 3) controller.kk_ = 2;
    		else controller.kk_  = controller.i_ ;


       		controller.getPredictModelState();

       		//cout << "z_p = " << controller.z_p_ << "\n";
    	}

    	/*
        cout <<"Dimensions z_uk: "<< controller.z_uk_.rows() << " x " << controller.z_uk_.columns() << "\n";
        cout << "z_uk = " << controller.z_uk_ << "\n";
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

        //controller.opt_.set_min_objective(boost::bind(&MPCNode::munc,this,_1,_2,_3), NULL);
        //controller.opt_.set_min_objective(myfunc, NULL);
        //opt.set_min_objective(MyFunction::wrap, &some_MyFunction)
        controller.opt_.set_min_objective(wrap, this);



        /*
        if (nlopt_optimize(controller.opt_, &controller.u_[0], &controller.minf_) < 0) {
            cout <<"Nlopt failed!" << "\n";
        }
        else {
            cout <<"Found minimum: " << controller.minf_ << "\n";
        }
		*/

        if ( controller.opt_.optimize(controller.tau_, controller.minf_) < 0) {
            cout <<"Nlopt failed!" << "\n";
        }
        else {
            cout <<"Found minimum: " << controller.minf_ << "\n";
        }

        controller.iter_ = 1;


    	controller.getAdditionalVelocity();
    	controller.getAdditionalPosition();



    	/********** Get control value **********/

        vector tau_real;

        tau_real = controller.getControlValue(&controller.tau_[0]);

        cout << "tau_real = \n" << tau_real << "\n";


        /// OVDJE SE IZRACUNAVA X,Y,N ZA JASNIJI PRIKAZ

        matrix B(Eigen::Matrix<double,3,4>::Zero());
        matrix Binv(Eigen::Matrix<double,4,3>::Zero());


        /*
        B << cos(M_PI/4),cos(M_PI/4),-cos(M_PI/4),-cos(M_PI/4),
        		sin(M_PI/4),-sin(M_PI/4),sin(M_PI/4),-sin(M_PI/4),
				0.5,-0.5,-0.5,0.5;
		*/

        /*
        B(0,0) = B(0,1) = cos(controller.alpha_);
        B(0,2) = B(0,3) = -cos(controller.alpha_);
        B(1,0) = B(1,2) = sin(controller.alpha_);
        B(1,1) = B(1,3) = -sin(controller.alpha_);
        B(2,0) = B(2,3) = controller.d_;
        B(2,1) = B(2,2) = -controller.d_;
		*/

        //Binv = B.transpose()*(B*B.transpose()).inverse();

        //vector tau_1 = Binv*tau_real;
        /*
        vector tau_1=(controller.C_d_*tau_real);

        auv_msgs::BodyForceReq tau_pub_1;

        tau_pub_1.header.stamp = ros::Time::now();
        tau_pub_1.wrench.force.x = tau_1[0];
        tau_pub_1.wrench.force.y = tau_1[1];
        tau_pub_1.wrench.force.z = 0;
		tau_pub_1.wrench.torque.z = tau_1[2];
		tau_pub_1.wrench.torque.x = 0;
		tau_pub_1.wrench.torque.y = 0;

        std_msgs::Float32MultiArray tau_pub;
        tau_pub.data.push_back(tau_real[0]);
        tau_pub.data.push_back(tau_real[1]);
        tau_pub.data.push_back(tau_real[2]);
        tau_pub.data.push_back(tau_real[3]);

        pub_tau_i.publish(tau_pub);
        pub_tau.publish(tau_pub_1);
	*/
        /// KRAJ //////////////////////////////////////////


        std_msgs::Float32MultiArray tau_pub;

        tau_pub.data.push_back(tau_real[0]);
        tau_pub.data.push_back(tau_real[1]);
        tau_pub.data.push_back(tau_real[2]);
        tau_pub.data.push_back(tau_real[3]);

        pub_tau_i.publish(tau_pub);

        /********** Get model output **********/

        vector z_real;

        z_real = controller.getModelOutput(tau_real);

        cout << "z_real = \n" << z_real << "\n";

        auv_msgs::BodyForceReq tau_pub_1;

        tau_pub_1.header.stamp = ros::Time::now();
        tau_pub_1.wrench.force.x = z_real[0];
        tau_pub_1.wrench.force.y = z_real[1];
        tau_pub_1.wrench.force.z = 0;
  		tau_pub_1.wrench.torque.z = z_real[2];
  		tau_pub_1.wrench.torque.x = 0;
  		tau_pub_1.wrench.torque.y = 0;

        pub_tau.publish(tau_pub_1);

        /********** Get real velocity **********/

     	vector velocity_real;

        velocity_real = controller.getRealVelocity(z_real);

        cout << "velocity = \n" << velocity_real << "\n";

        auv_msgs::NavSts vel_pub;

        vel_pub.header.stamp = ros::Time::now();
        vel_pub.body_velocity.x = velocity_real[0];
        vel_pub.body_velocity.y= velocity_real[1];
  		vel_pub.body_velocity.z = 0;
        vel_pub.orientation_rate.roll = 0;
        vel_pub.orientation_rate.pitch = 0;
        vel_pub.orientation_rate.yaw = velocity_real[2];

        pub_vel.publish(vel_pub);

        /********** Get real position **********/

        vector position_real;

        position_real = controller.getRealPosition();

        cout << "position = \n" << position_real << "\n";

        auv_msgs::NavSts pos_pub;

        pos_pub.header.stamp = ros::Time::now();
        pos_pub.position.north = position_real[0];
        pos_pub.position.east = position_real[1];
        pos_pub.position.depth = 0;
   		pos_pub.orientation.roll = 0;
   		pos_pub.orientation.pitch  = 0;
   		pos_pub.orientation.yaw  = position_real[2];

        pub_pos.publish(pos_pub);

        /********** Initial control value for next step **********/

     	for (int row = 0; row < 4; row++){
     		for (int col = 0; col < 3; col++){

     			if(col < 2) controller.tau_[row*controller.vel_+col] = controller.tau_[row*controller.vel_+col+1];
     			//cout << controller.tau_[row*controller.vel_+col] << "  ";

     			}
     	    //cout << "\n";
     	}



     	//exit(4);

     	myfile_z << z_real.transpose() << " ;\n";
     	//myfile_ref << ref.transpose().row(controller.k_ + controller.Hw_)<< " ;\n";
     	myfile_ref << reference.block<1,kor>(0,0) << " ";
     	myfile_ref << reference.block<1,kor>(1,0) << " ";
     	myfile_ref << reference.block<1,kor>(2,0) << " ;\n";
     	myfile_velocity << velocity_real.transpose() << " ;\n";
     	myfile_position << position_real.transpose() << " ;\n";
     	myfile_tau << tau_real.transpose() << " ;\n";

    //main


	/********************************/


}

MPCNode::matrix MPCNode::get_reference()
{
	/********** Set reference **********/

/*	int kor = 50;
	int i = 0, j = 0;
    matrix ref;
    ref = Eigen::MatrixXd::Zero(3, kor);


    double rr = 0;

    for (i = 0; i < 3; i++){
           for (j = 0; j < kor; j++){
                    if(i == 0) {
                    	//if (j < 24) ref(i,j) =  1 * sin(pi*j/1/kor);
                    	//else ref(i,j) = 1;
                    	//ref(i,j) = 1 * sin(pi*j/2/kor) + 0;           		//referent value
                    	ref(i,j) = rr;
                    	rr =  rr + ((6.5217e-3)*j);
                    }
            }
    }

    return ref;
*/	  return reference;
}
void MPCNode::onReference(const auv_msgs::NavSts::ConstPtr& data)
{
/*
	reference.block<1,kor>(0,0) = data->body_velocity.x*Eigen::MatrixXd::Ones(1, kor);
	reference.block<1,kor>(1,0) = data->body_velocity.y*Eigen::MatrixXd::Ones(1, kor);
	reference.block<1,kor>(2,0) = data->body_velocity.z*Eigen::MatrixXd::Ones(1, kor);
*/

	reference.block<1,kor>(0,0) = data->position.north*Eigen::MatrixXd::Ones(1, kor);
	reference.block<1,kor>(1,0) = data->position.east*Eigen::MatrixXd::Ones(1, kor);
	reference.block<1,kor>(2,0) = data->orientation.yaw*Eigen::MatrixXd::Ones(1, kor);

}

void MPCNode::onStateHat(const auv_msgs::NavSts::ConstPtr& data)
{
	state = *data;
}

/*********************************************************************
 ***  Main function
 ********************************************************************/

int main(int argc, char** argv)
{
	ros::init(argc, argv, "mpc_node");
	ros::NodeHandle nh;
	MPCNode mpc_node;

	/*** Time sample (in seconds) ***/
	double Ts(1.0);
	ros::Rate rate(1/Ts);

	mpc_node.initialize();


	// PRIVREMENO RJESENJE
    int counter;
    //int kor = 50; ///////// Ovo je privremeno!!!!!
    int Hp = mpc_node.controller.Hp_;

    counter = MPCNode::kor - Hp;

    mpc_node.controller.k_ = 0;


    //for (controller.k_ = 0; controller.k_ < counter; mpc_node.controller.k_++){

	//while (ros::ok() && mpc_node.controller.k_ < counter )
	while (ros::ok())
	{

		mpc_node.step(mpc_node.get_reference());

		mpc_node.controller.k_++;

		if(mpc_node.controller.k_ == 46){
			   myfile_z.close();
			    myfile_ref.close();
			    myfile_velocity.close();
			    myfile_position.close();
			    myfile_tau.close();
			    cout << "\n SAVED \n";
		}
		ros::spinOnce();
		rate.sleep();
	}





	return 0;

}
