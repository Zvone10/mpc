#include <mpc/mpc_trajektorija.hpp>
#include <ros/ros.h>
#include <boost/bind.hpp>
#include <fstream>
#include<cmath>

#include <auv_msgs/NavSts.h>
#include <auv_msgs/BodyForceReq.h>
#include <std_msgs/Float32MultiArray.h>




class Relej
{
public:

	Relej();

	~Relej();


	void onStateHat(const auv_msgs::NavSts::ConstPtr& data);
	void onTauSlow(const auv_msgs::BodyForceReq::ConstPtr& data);


	ros::Subscriber  sub_state,sub_tau_slow;
	ros::Publisher  pub_tau;
	auv_msgs::NavSts state;
	auv_msgs::BodyForceReq tau_out_;

};

Relej::Relej()
{

	ros::NodeHandle nh;

	sub_state = nh.subscribe<auv_msgs::NavSts>("/cres/position",1,&Relej::onStateHat,this);
	sub_tau_slow = nh.subscribe<auv_msgs::BodyForceReq>("/cres/tauOutSlow",1,&Relej::onTauSlow,this);

	pub_tau = nh.advertise<auv_msgs::BodyForceReq>("/cres/tauOut",1);

}

Relej::~Relej()
{

}

void Relej::onTauSlow(const auv_msgs::BodyForceReq::ConstPtr& data)
{
	//ROS_ERROR("DEBUG");
	tau_out_ = *data;
}


void Relej::onStateHat(const auv_msgs::NavSts::ConstPtr& data)
{
	//ROS_ERROR("DEBUG");
	state = *data;
	tau_out_.header.stamp = ros::Time::now();
	pub_tau.publish(tau_out_);
}

/*********************************************************************
 ***  Main function
 ********************************************************************/

int main(int argc, char** argv)
{
	ros::init(argc, argv, "relej_node");
	ros::NodeHandle nh;
	Relej relej_node;

	ros::spin();

	/*** Time sample (in seconds) ***/
	double Ts(1.0);
	ros::Rate rate(1/Ts);

    //ros::MultiThreadedSpinner spinner(2);
	while (ros::ok())
	{

		ros::spinOnce();
		//spinner.spinOnce();
		rate.sleep();
	}





	return 0;

}
