# CreepyCrawler-Quadruped-Robot
Quadruped robots have been of interest to the robotics community as they have high maneuverability and they can be designed to be very robust to the environment they are placed in. The use of quadrupeds as service robots is being realized by the community. Understanding the dynamics of a quadruped robot is a complex task as there are many factors that need to be taken into consideration while designing the robot. This project intends to develop an eight degree of freedom (DOF) sprawling-type quadruped robot. Control of high DOF robots can be highly challenging and computationally heavy. This project intends to work on bringing down the complexity of the robot control by decoupling the dynamics of the robot into four separate dynamic systems, i.e one for each leg, allowing us to achieve a decentralized control of the robot. Gazebo simulation of the robot is achieved with a few gait implementations using a decentralized control strategy.

# Creepy_Crawler

* Author: Kunal Nandanwar  <kgnandanwar@wpi.edu>
* License: GNU General Public License, version 3 (GPL-3.0)

Example robots and code for interfacing Gazebo with ROS

## Quick Start
One_Leg:

        Rviz:

            roslaunch creepy_crawler one_leg_rviz.launch

        Gazebo:

            roslaunch creepy_crawler one_leg_world.launch

        ROS Control:

            roslaunch creepy_crawler one_leg_control.launch

        Example of Moving Joints:

            rostopic pub /creepy_crawler/joint1_1_effort_controller/command std_msgs/Float64 "data: -0.9"

Creepy_Crawler:

        Rviz:

            roslaunch creepy_crawler creepy_crawler_rviz.launch

        Gazebo:

            roslaunch creepy_crawler creepy_crawler_world.launch

        ROS Control:

            roslaunch creepy_crawler creepy_crawler_effort_control.launch

        Example of Moving Joints:

            rostopic pub /creepy_crawler/joint1_1_effort_controller/command std_msgs/Float64 "data: -0.9"


Creepy_Spot:
        Gazebo:

            roslaunch creepy_crawler creepy_spot_world.launch

        ROS Control:

            roslaunch creepy_crawler creepy_spot_effort_control.launch

        Example of Moving Joints:

            rostopic pub /creepy_spots/joint1_1_effort_controller/command std_msgs/Float64 "data: -0.9"
            
            
src folder contains all the neccassary matlab codes for the project:

	Trajectory folder contains trajectpry generation
	controls contains code for feedback linearization control algo and lagrangian dynamics model
	dynamics_RNE contains RNE based dynamics model and matlab implementation of one leg movement.
	

