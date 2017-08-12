// -- BEGIN LICENSE BLOCK -----------------------------------------------------
// -- END LICENSE BLOCK -------------------------------------------------------

//-----------------------------------------------------------------------------
/*!\file    cartesian_controller_base.h
 *
 * \author  Stefan Scherzinger <scherzin@fzi.de>
 * \date    2017/07/27
 *
 */
//-----------------------------------------------------------------------------

#ifndef CARTESIAN_CONTROLLER_BASE_H_INCLUDED
#define CARTESIAN_CONTROLLER_BASE_H_INCLUDED

// ROS
#include <ros/node_handle.h>
#include <trajectory_msgs/JointTrajectoryPoint.h>
#include <geometry_msgs/WrenchStamped.h>

// ros_controls
#include <controller_interface/controller.h>
#include <hardware_interface/joint_command_interface.h>

// KDL
#include <kdl/treefksolverpos_recursive.hpp>

// Project
#include <cartesian_controller_base/ForwardDynamicsSolver.h>
#include <cartesian_controller_base/SpatialPIDController.h>
#include <cartesian_controller_base/Utility.h>

// Other
#include <vector>
#include <string>

namespace cartesian_controller_base
{

template <class HardwareInterface>
class CartesianControllerBase : public controller_interface::Controller<HardwareInterface>
{
  public:
    CartesianControllerBase();
    virtual ~CartesianControllerBase<HardwareInterface>(){};

    virtual bool init(HardwareInterface* hw, ros::NodeHandle& nh);

    virtual void starting(const ros::Time& time);

  protected:
    void writeJointControlCmds();

    void computeJointControlCmds(const ctrl::Vector6D& error, const ros::Duration& period);

    ctrl::Vector6D displayInBaseLink(const geometry_msgs::WrenchStamped& wrench, const std::string& from);

    std::vector<hardware_interface::JointHandle>      m_joint_handles;
    std::vector<std::string>                          m_joint_names;
    ForwardDynamicsSolver                             m_forward_dynamics_solver;
    boost::shared_ptr<KDL::TreeFkSolverPos_recursive> m_forward_kinematics_solver;
    trajectory_msgs::JointTrajectoryPoint             m_simulated_joint_motion;
    SpatialPIDController                              m_spatial_controller;
    ctrl::Vector6D                                    m_cartesian_input;
    std::string                                       m_robot_base_link;
    std::string                                       m_end_effector_link;
};

}

#include <cartesian_controller_base/cartesian_controller_base.hpp>

#endif