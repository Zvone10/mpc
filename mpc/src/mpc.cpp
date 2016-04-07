#include <mpc/mpc.hpp>
#include <ros/ros.h>

/*********************************************************************
 ***  Main function
 ********************************************************************/

int main(int argc, char** argv)
{
	ros::init(argc, argv, "mpc");
	MPC controller;
	ros::spin();
	return 0;
}





