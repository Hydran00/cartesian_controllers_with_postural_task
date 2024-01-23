////////////////////////////////////////////////////////////////////////////////
// Copyright 2019 FZI Research Center for Information Technology
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived from this
// software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
////////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
/*!\file    ForwardDynamicsSolver.cpp
 *
 * \author  Stefan Scherzinger <scherzin@fzi.de>
 * \date    2020/03/24
 *
 */
//-----------------------------------------------------------------------------

#include <algorithm>
#include <cartesian_controller_base/ForwardDynamicsSolver.h>
#include <kdl/framevel.hpp>
#include <kdl/jntarrayvel.hpp>
#include <map>
#include <pluginlib/class_list_macros.hpp>
#include <sstream>

/**
 * \class cartesian_controller_base::ForwardDynamicsSolver
 *
 * Users may explicitly specify it with \a "forward_dynamics" as \a ik_solver
 * in their controllers.yaml configuration file for each controller:
 *
 * \code{.yaml}
 * <name_of_your_controller>:
 *   ros__parameters:
 *     ik_solver: "forward_dynamics"
 *     ...
 *
 *     solver:
 *         ...
 *         forward_dynamics:
 *             link_mass: 0.5
 * \endcode
 *
 */
PLUGINLIB_EXPORT_CLASS(cartesian_controller_base::ForwardDynamicsSolver, cartesian_controller_base::IKSolver)

namespace cartesian_controller_base
{

  ForwardDynamicsSolver::ForwardDynamicsSolver()
  {
  }

  ForwardDynamicsSolver::~ForwardDynamicsSolver() {}

  trajectory_msgs::msg::JointTrajectoryPoint ForwardDynamicsSolver::getJointControlCmds(
      rclcpp::Duration period,
      const ctrl::Vector6D &net_force)
  {

    // Compute joint space inertia matrix with actualized link masses
    buildGenericModel();
    m_jnt_space_inertia_solver->JntToMass(m_current_positions, m_jnt_space_inertia);

    // Compute joint jacobian
    m_jnt_jacobian_solver->JntToJac(m_current_positions, m_jnt_jacobian);

    //////////Compute postural task///////////////////////////////////////////
    Eigen::MatrixXd J = m_jnt_jacobian.data;

    USING_NAMESPACE_QPOASES
    /* Setup data of first QP. */
 
// [         J1_1^2 + J2_1^2 + J3_1^2, J1_1*J1_2 + J2_1*J2_2 + J3_1*J3_2, J1_1*J1_3 + J2_1*J2_3 + J3_1*J3_3, J1_1*J1_4 + J2_1*J2_4 + J3_1*J3_4, J1_1*J1_5 + J2_1*J2_5 + J3_1*J3_5, J1_1*J1_6 + J2_1*J2_6 + J3_1*J3_6]
// [J1_1*J1_2 + J2_1*J2_2 + J3_1*J3_2,          J1_2^2 + J2_2^2 + J3_2^2, J1_2*J1_3 + J2_2*J2_3 + J3_2*J3_3, J1_2*J1_4 + J2_2*J2_4 + J3_2*J3_4, J1_2*J1_5 + J2_2*J2_5 + J3_2*J3_5, J1_2*J1_6 + J2_2*J2_6 + J3_2*J3_6]
// [J1_1*J1_3 + J2_1*J2_3 + J3_1*J3_3, J1_2*J1_3 + J2_2*J2_3 + J3_2*J3_3,          J1_3^2 + J2_3^2 + J3_3^2, J1_3*J1_4 + J2_3*J2_4 + J3_3*J3_4, J1_3*J1_5 + J2_3*J2_5 + J3_3*J3_5, J1_3*J1_6 + J2_3*J2_6 + J3_3*J3_6]
// [J1_1*J1_4 + J2_1*J2_4 + J3_1*J3_4, J1_2*J1_4 + J2_2*J2_4 + J3_2*J3_4, J1_3*J1_4 + J2_3*J2_4 + J3_3*J3_4,          J1_4^2 + J2_4^2 + J3_4^2, J1_4*J1_5 + J2_4*J2_5 + J3_4*J3_5, J1_4*J1_6 + J2_4*J2_6 + J3_4*J3_6]
// [J1_1*J1_5 + J2_1*J2_5 + J3_1*J3_5, J1_2*J1_5 + J2_2*J2_5 + J3_2*J3_5, J1_3*J1_5 + J2_3*J2_5 + J3_3*J3_5, J1_4*J1_5 + J2_4*J2_5 + J3_4*J3_5,          J1_5^2 + J2_5^2 + J3_5^2, J1_5*J1_6 + J2_5*J2_6 + J3_5*J3_6]
// [J1_1*J1_6 + J2_1*J2_6 + J3_1*J3_6, J1_2*J1_6 + J2_2*J2_6 + J3_2*J3_6, J1_3*J1_6 + J2_3*J2_6 + J3_3*J3_6, J1_4*J1_6 + J2_4*J2_6 + J3_4*J3_6, J1_5*J1_6 + J2_5*J2_6 + J3_5*J3_6,          J1_6^2 + J2_6^2 + J3_6^2]
    real_t H[6 * 6] = {
      pow(J(0,0),2) + pow(J(1,0),2) + pow(J(2,0),2), J(0,0)*J(0,1) + J(1,0)*J(1,1) + J(2,0)*J(2,1), J(0,0)*J(0,2) + J(1,0)*J(1,2) + J(2,0)*J(2,2), J(0,0)*J(0,3) + J(1,0)*J(1,3) + J(2,0)*J(2,3), J(0,0)*J(0,4) + J(1,0)*J(1,4) + J(2,0)*J(2,4), J(0,0)*J(0,5) + J(1,0)*J(1,5) + J(2,0)*J(2,5),
      J(0,0)*J(0,1) + J(1,0)*J(1,1) + J(2,0)*J(2,1), pow(J(0,1),2) + pow(J(1,1),2) + pow(J(2,1),2), J(0,1)*J(0,2) + J(1,1)*J(1,2) + J(2,1)*J(2,2), J(0,1)*J(0,3) + J(1,1)*J(1,3) + J(2,1)*J(2,3), J(0,1)*J(0,4) + J(1,1)*J(1,4) + J(2,1)*J(2,4), J(0,1)*J(0,5) + J(1,1)*J(1,5) + J(2,1)*J(2,5),
      J(0,0)*J(0,2) + J(1,0)*J(1,2) + J(2,0)*J(2,2), J(0,1)*J(0,2) + J(1,1)*J(1,2) + J(2,1)*J(2,2), pow(J(0,2),2) + pow(J(1,2),2) + pow(J(2,2),2), J(0,2)*J(0,3) + J(1,2)*J(1,3) + J(2,2)*J(2,3), J(0,2)*J(0,4) + J(1,2)*J(1,4) + J(2,2)*J(2,4), J(0,2)*J(0,5) + J(1,2)*J(1,5) + J(2,2)*J(2,5),
      J(0,0)*J(0,3) + J(1,0)*J(1,3) + J(2,0)*J(2,3), J(0,1)*J(0,3) + J(1,1)*J(1,3) + J(2,1)*J(2,3), J(0,2)*J(0,3) + J(1,2)*J(1,3) + J(2,2)*J(2,3), pow(J(0,3),2) + pow(J(1,3),2) + pow(J(2,3),2), J(0,3)*J(0,4) + J(1,3)*J(1,4) + J(2,3)*J(2,4), J(0,3)*J(0,5) + J(1,3)*J(1,5) + J(2,3)*J(2,5),
      J(0,0)*J(0,4) + J(1,0)*J(1,4) + J(2,0)*J(2,4), J(0,1)*J(0,4) + J(1,1)*J(1,4) + J(2,1)*J(2,4), J(0,2)*J(0,4) + J(1,2)*J(1,4) + J(2,2)*J(2,4), J(0,3)*J(0,4) + J(1,3)*J(1,4) + J(2,3)*J(2,4), pow(J(0,4),2) + pow(J(1,4),2) + pow(J(2,4),2), J(0,4)*J(0,5) + J(1,4)*J(1,5) + J(2,4)*J(2,5),
      J(0,0)*J(0,5) + J(1,0)*J(1,5) + J(2,0)*J(2,5), J(0,1)*J(0,5) + J(1,1)*J(1,5) + J(2,1)*J(2,5), J(0,2)*J(0,5) + J(1,2)*J(1,5) + J(2,2)*J(2,5), J(0,3)*J(0,5) + J(1,3)*J(1,5) + J(2,3)*J(2,5), J(0,4)*J(0,5) + J(1,4)*J(1,5) + J(2,4)*J(2,5), pow(J(0,5),2) + pow(J(1,5),2) + pow(J(2,5),2)
    };
    Eigen::MatrixXd ddq = m_current_accelerations;
    Eigen::MatrixXd ddx = 
    // J1_1*(J1_1*ddq_0 - ddx_x + J1_2*ddq_1 + J1_3*ddq_2 + J1_4*ddq_3 + J1_5*ddq_4 + J1_6*ddq_5 + dJ1_1*dq_0 + dJ1_2*dq_1 + dJ1_3*dq_2 + dJ1_4*dq_3 + dJ1_5*dq_4 + dJ1_6*dq_5) + J2_1*(J2_1*ddq_0 - ddx_y + J2_2*ddq_1 + J2_3*ddq_2 + J2_4*ddq_3 + J2_5*ddq_4 + J2_6*ddq_5 + dJ2_1*dq_0 + dJ2_2*dq_1 + dJ2_3*dq_2 + dJ2_4*dq_3 + dJ2_5*dq_4 + dJ2_6*dq_5) + J3_1*(J3_1*ddq_0 - ddx_z + J3_2*ddq_1 + J3_3*ddq_2 + J3_4*ddq_3 + J3_5*ddq_4 + J3_6*ddq_5 + dJ3_1*dq_0 + dJ3_2*dq_1 + dJ3_3*dq_2 + dJ3_4*dq_3 + dJ3_5*dq_4 + dJ3_6*dq_5)
    // J1_2*(J1_1*ddq_0 - ddx_x + J1_2*ddq_1 + J1_3*ddq_2 + J1_4*ddq_3 + J1_5*ddq_4 + J1_6*ddq_5 + dJ1_1*dq_0 + dJ1_2*dq_1 + dJ1_3*dq_2 + dJ1_4*dq_3 + dJ1_5*dq_4 + dJ1_6*dq_5) + J2_2*(J2_1*ddq_0 - ddx_y + J2_2*ddq_1 + J2_3*ddq_2 + J2_4*ddq_3 + J2_5*ddq_4 + J2_6*ddq_5 + dJ2_1*dq_0 + dJ2_2*dq_1 + dJ2_3*dq_2 + dJ2_4*dq_3 + dJ2_5*dq_4 + dJ2_6*dq_5) + J3_2*(J3_1*ddq_0 - ddx_z + J3_2*ddq_1 + J3_3*ddq_2 + J3_4*ddq_3 + J3_5*ddq_4 + J3_6*ddq_5 + dJ3_1*dq_0 + dJ3_2*dq_1 + dJ3_3*dq_2 + dJ3_4*dq_3 + dJ3_5*dq_4 + dJ3_6*dq_5)
    // J1_3*(J1_1*ddq_0 - ddx_x + J1_2*ddq_1 + J1_3*ddq_2 + J1_4*ddq_3 + J1_5*ddq_4 + J1_6*ddq_5 + dJ1_1*dq_0 + dJ1_2*dq_1 + dJ1_3*dq_2 + dJ1_4*dq_3 + dJ1_5*dq_4 + dJ1_6*dq_5) + J2_3*(J2_1*ddq_0 - ddx_y + J2_2*ddq_1 + J2_3*ddq_2 + J2_4*ddq_3 + J2_5*ddq_4 + J2_6*ddq_5 + dJ2_1*dq_0 + dJ2_2*dq_1 + dJ2_3*dq_2 + dJ2_4*dq_3 + dJ2_5*dq_4 + dJ2_6*dq_5) + J3_3*(J3_1*ddq_0 - ddx_z + J3_2*ddq_1 + J3_3*ddq_2 + J3_4*ddq_3 + J3_5*ddq_4 + J3_6*ddq_5 + dJ3_1*dq_0 + dJ3_2*dq_1 + dJ3_3*dq_2 + dJ3_4*dq_3 + dJ3_5*dq_4 + dJ3_6*dq_5)
    // J1_4*(J1_1*ddq_0 - ddx_x + J1_2*ddq_1 + J1_3*ddq_2 + J1_4*ddq_3 + J1_5*ddq_4 + J1_6*ddq_5 + dJ1_1*dq_0 + dJ1_2*dq_1 + dJ1_3*dq_2 + dJ1_4*dq_3 + dJ1_5*dq_4 + dJ1_6*dq_5) + J2_4*(J2_1*ddq_0 - ddx_y + J2_2*ddq_1 + J2_3*ddq_2 + J2_4*ddq_3 + J2_5*ddq_4 + J2_6*ddq_5 + dJ2_1*dq_0 + dJ2_2*dq_1 + dJ2_3*dq_2 + dJ2_4*dq_3 + dJ2_5*dq_4 + dJ2_6*dq_5) + J3_4*(J3_1*ddq_0 - ddx_z + J3_2*ddq_1 + J3_3*ddq_2 + J3_4*ddq_3 + J3_5*ddq_4 + J3_6*ddq_5 + dJ3_1*dq_0 + dJ3_2*dq_1 + dJ3_3*dq_2 + dJ3_4*dq_3 + dJ3_5*dq_4 + dJ3_6*dq_5)
    // J1_5*(J1_1*ddq_0 - ddx_x + J1_2*ddq_1 + J1_3*ddq_2 + J1_4*ddq_3 + J1_5*ddq_4 + J1_6*ddq_5 + dJ1_1*dq_0 + dJ1_2*dq_1 + dJ1_3*dq_2 + dJ1_4*dq_3 + dJ1_5*dq_4 + dJ1_6*dq_5) + J2_5*(J2_1*ddq_0 - ddx_y + J2_2*ddq_1 + J2_3*ddq_2 + J2_4*ddq_3 + J2_5*ddq_4 + J2_6*ddq_5 + dJ2_1*dq_0 + dJ2_2*dq_1 + dJ2_3*dq_2 + dJ2_4*dq_3 + dJ2_5*dq_4 + dJ2_6*dq_5) + J3_5*(J3_1*ddq_0 - ddx_z + J3_2*ddq_1 + J3_3*ddq_2 + J3_4*ddq_3 + J3_5*ddq_4 + J3_6*ddq_5 + dJ3_1*dq_0 + dJ3_2*dq_1 + dJ3_3*dq_2 + dJ3_4*dq_3 + dJ3_5*dq_4 + dJ3_6*dq_5)
    // J1_6*(J1_1*ddq_0 - ddx_x + J1_2*ddq_1 + J1_3*ddq_2 + J1_4*ddq_3 + J1_5*ddq_4 + J1_6*ddq_5 + dJ1_1*dq_0 + dJ1_2*dq_1 + dJ1_3*dq_2 + dJ1_4*dq_3 + dJ1_5*dq_4 + dJ1_6*dq_5) + J2_6*(J2_1*ddq_0 - ddx_y + J2_2*ddq_1 + J2_3*ddq_2 + J2_4*ddq_3 + J2_5*ddq_4 + J2_6*ddq_5 + dJ2_1*dq_0 + dJ2_2*dq_1 + dJ2_3*dq_2 + dJ2_4*dq_3 + dJ2_5*dq_4 + dJ2_6*dq_5) + J3_6*(J3_1*ddq_0 - ddx_z + J3_2*ddq_1 + J3_3*ddq_2 + J3_4*ddq_3 + J3_5*ddq_4 + J3_6*ddq_5 + dJ3_1*dq_0 + dJ3_2*dq_1 + dJ3_3*dq_2 + dJ3_4*dq_3 + dJ3_5*dq_4 + dJ3_6*dq_5)

    real_t g[6]{
      J(0,0)*
    }
    Eigen::VectorXd ddq_min = Eigen::VectorXd::Zero(6);
    ddq_min << -0.1, -0.1, -0.1, -0.1, -0.1, -0.1;
    Eigen::VectorXd ddq_max = Eigen::VectorXd::Zero(6);
    ddq_max << 0.1, 0.1, 0.1, 0.1, 0.1, 0.1;
    real_t lb_ddq[6] = {ddq_min(0), ddq_min(1), ddq_min(2), ddq_min(3), ddq_min(4), ddq_min(5)};
    real_t ub_ddq[6] = {ddq_max(0), ddq_max(1), ddq_max(2), ddq_max(3), ddq_max(4), ddq_max(5)};
    // Imposing dynamics M ddq = J^T f -> ddq = M^-1 J^T f
    Eigen::MatrixXd dyn_constraint = m_jnt_space_inertia.data.inverse() * m_jnt_jacobian.data.transpose() * net_force;
    real_t lb_A[6] = {dyn_constraint[0], dyn_constraint[1], dyn_constraint[2], dyn_constraint[3], dyn_constraint[4], dyn_constraint[5]};
    real_t ub_A[6] = {dyn_constraint[0], dyn_constraint[1], dyn_constraint[2], dyn_constraint[3], dyn_constraint[4], dyn_constraint[5]};

    /* Setting up QProblem object. */
    QProblem min_problem(6, 3);

    /* Solve QP. */
    int nWSR = 10;
    min_problem.init(H, g, A, lb, ub, lbA, ubA, nWSR);

    // real_t xOpt[3];
    // min_problem.getPrimalSolution(xOpt);

    // ImpedanceBase::m_cartesian_stiffness(0, 0) = xOpt(0);
    // ImpedanceBase::m_cartesian_stiffness(1, 1) = xOpt[1];
    // ImpedanceBase::m_cartesian_stiffness(2, 2) = xOpt[2];


    //////////////////////////////////////////////////////////////////////////

    // Computes the joint accelerations according to: \f$ \ddot{q} = H^{-1} ( J^T f) \f$
    // m_current_accelerations.data = m_jnt_space_inertia.data.inverse() * m_jnt_jacobian.data.transpose() * net_force + ddq0;
    real_t xOpt[6];
    min_problem.getPrimalSolution( xOpt );
    m_current_accelerations.data << xOpt[0], xOpt[1], xOpt[2], xOpt[3], xOpt[4], xOpt[5];

    // Numerical time integration with the Euler forward method
    m_current_positions.data = m_last_positions.data + m_last_velocities.data * period.seconds();
    m_current_velocities.data = m_last_velocities.data + m_current_accelerations.data * period.seconds();
    m_current_velocities.data *= 0.9; // 10 % global damping against unwanted null space motion.
                                      // Will cause exponential slow-down without input.
    // Make sure positions stay in allowed margins
    applyJointLimits();

    // Apply results
    trajectory_msgs::msg::JointTrajectoryPoint control_cmd;
    for (int i = 0; i < m_number_joints; ++i)
    {
      control_cmd.positions.push_back(m_current_positions(i));
      control_cmd.velocities.push_back(m_current_velocities(i));

      // Accelerations should be left empty. Those values will be interpreted
      // by most hardware joint drivers as max. tolerated values. As a
      // consequence, the robot will move very slowly.
    }
    control_cmd.time_from_start = period; // valid for this duration

    // Update for the next cycle
    m_last_positions = m_current_positions;
    m_last_velocities = m_current_velocities;

    // m_last_delta_q0 = delta_q0;

    return control_cmd;
  }

#if defined CARTESIAN_CONTROLLERS_HUMBLE || defined CARTESIAN_CONTROLLERS_IRON
  bool ForwardDynamicsSolver::init(std::shared_ptr<rclcpp_lifecycle::LifecycleNode> nh,
#else
  bool ForwardDynamicsSolver::init(std::shared_ptr<rclcpp::Node> nh,
#endif
                                   const KDL::Chain &chain,
                                   const KDL::JntArray &upper_pos_limits,
                                   const KDL::JntArray &lower_pos_limits)
  {
    IKSolver::init(nh, chain, upper_pos_limits, lower_pos_limits);

    if (!buildGenericModel())
    {
      RCLCPP_ERROR(nh->get_logger(), "Something went wrong in setting up the internal model.");
      return false;
    }

    // Forward dynamics
    m_jnt_jacobian_solver.reset(new KDL::ChainJntToJacSolver(m_chain));
    m_jnt_space_inertia_solver.reset(new KDL::ChainDynParam(m_chain, KDL::Vector::Zero()));
    m_jnt_jacobian.resize(m_number_joints);
    m_jnt_space_inertia.resize(m_number_joints);
    m_postural_joints.resize(m_number_joints);

    m_postural_joints = Eigen::VectorXd::Zero(m_number_joints);
    m_postural_conf = Eigen::VectorXd::Zero(m_number_joints);
    m_postural_joints << 0.0, 0.0, 1.0, 0.0, 0.0, 0.0;
    m_postural_conf << 0, 0, 0.0, 0, 0, 0;
    m_last_delta_q0 = Eigen::VectorXd::Zero(m_number_joints);
    m_last_delta_q0 << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    std::cout << "m_postural_joints: \n"
              << m_last_delta_q0 << std::endl;

    // Set the initial value if provided at runtime, else use default value.
    m_min = nh->declare_parameter<double>(m_params + "/link_mass", 0.1);

    RCLCPP_INFO(nh->get_logger(), "Forward dynamics solver initialized");
    RCLCPP_INFO(nh->get_logger(), "Forward dynamics solver has control over %i joints", m_number_joints);

    return true;
  }

  bool ForwardDynamicsSolver::buildGenericModel()
  {
    // Set all masses and inertias to minimal (yet stable) values.
    double ip_min = 0.000001;
    for (size_t i = 0; i < m_chain.segments.size(); ++i)
    {
      // Fixed joint segment
      if (m_chain.segments[i].getJoint().getType() == KDL::Joint::None)
      {
        m_chain.segments[i].setInertia(
            KDL::RigidBodyInertia::Zero());
      }
      else // relatively moving segment
      {
        m_chain.segments[i].setInertia(
            KDL::RigidBodyInertia(
                m_min,               // mass
                KDL::Vector::Zero(), // center of gravity
                KDL::RotationalInertia(
                    ip_min, // ixx
                    ip_min, // iyy
                    ip_min  // izz
                    // ixy, ixy, iyz default to 0.0
                    )));
      }
    }

    // Only give the last segment a generic mass and inertia.
    // See https://arxiv.org/pdf/1908.06252.pdf for a motivation for this setting.
    double m = 1;
    double ip = 1;
    m_chain.segments[m_chain.segments.size() - 1].setInertia(
        KDL::RigidBodyInertia(
            m,
            KDL::Vector::Zero(),
            KDL::RotationalInertia(ip, ip, ip)));

    return true;
  }

} // namespace
