#!/usr/bin/env python

# Copyright (C) 2017 Udacity Inc.
#
# This file is part of Robotic Arm: Pick and Place project for Udacity
# Robotics nano-degree program
#
# All Rights Reserved.

# Author: Harsh Pandya

# import modules
import rospy
import tf
from kuka_arm.srv import *
from trajectory_msgs.msg import JointTrajectory, JointTrajectoryPoint
from geometry_msgs.msg import Pose
from mpmath import *
from sympy import *


class FK:
    """Pre-compute forward kinematics parameters"""

    def __init__(self):
        # Create symbols
        # link offset
        d1, d2, d3, d4, d5, d6, d7 = symbols('d1:8')
        # link length
        a0, a1, a2, a3, a4, a5, a6 = symbols('a0:7')
        # Twist angles
        alpha0, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6 = symbols('alpha0:7')
        # Joint angle symbols
        self.q1, self.q2, self.q3, q4, q5, q6, q7 = symbols('q1:8')

        # Create Modified DH parameters
        # DH params
        DH_Table = {alpha0: 0, a0: 0, d1: 0.75, q1: q1,
                    alpha1: -pi / 2., a1: 0.35, d2: 0, q2: q2 - pi / 2.,
                    alpha2: 0, a2: 1.25, d3: 0, q3: q3,
                    alpha3: -pi / 2., a3: -0.054, d4: 1.5, q4: q4,
                    alpha4: pi / 2., a4: 0, d5: 0, q5: q5,
                    alpha5: -pi / 2., a5: 0, d6: 0, q6: q6,
                    alpha6: 0, a6: 0, d7: 0.303, q7: 0,
                    }

        # Create individual transformation matrices
        self.t0_1 = tf_matrix(alpha0, a0, d1, q1).subs(DH_Table)
        self.t1_2 = tf_matrix(alpha1, a1, d2, q2).subs(DH_Table)
        self.t2_3 = tf_matrix(alpha2, a2, d3, q3).subs(DH_Table)
        t3_4 = tf_matrix(alpha3, a3, d4, q4).subs(DH_Table)
        t4_5 = tf_matrix(alpha4, a4, d5, q5).subs(DH_Table)
        t5_6 = tf_matrix(alpha5, a5, d6, q6).subs(DH_Table)
        t6_ee = tf_matrix(alpha6, a6, d7, q7).subs(DH_Table)

        t0_ee = t0_1 * t1_2 * t2_3 * t3_4 * t4_5 * t5_6 * t6_ee

        # Extract rotation matrices from the transformation matrices
        r, p, y = symbols('r p y')

        # Roll
        rot_x = Matrix([[1, 0, 0],
                        [0, cos(r), -sin(r)],
                        [0, sin(r), cos(r)]])
        # Pitch
        rot_y = Matrix([[cos(p), 0, sin(p)],
                        [0, 1, 0],
                        [-sin(p), 0, cos(p)]])
        # Yaw
        rot_z = Matrix([[cos(y), -sin(y), 0],
                        [sin(y), cos(y), 0],
                        [0, 0, 1]])

        rot_ee = rot_z * rot_y * rot_x
        # Error correction - DH parameters to URDF
        rot_error = rot_z.subs(y, radians(180)) * rot_y.subs(p, radians(-90))
        self.corrected_ee_rot_matrix = corrected_ee_rot_matrix * rot_error


# Define Modified DH Transformation matrix
def tf_matrix(alpha, a, d, q):
    """Return a Denavit-Hartenberg transformation matrix based on DH params"""

    TF = Matrix([[cos(q), -sin(q), 0, a],
                 [sin(q) * cos(alpha), cos(q) * cos(alpha), -sin(alpha), -sin(alpha) * d],
                 [sin(q) * sin(alpha), cos(q) * sin(alpha), cos(alpha), cos(alpha) * d],
                 [0, 0, 0, 1]])
    return TF


print("Starting FK ..")
_forward_kinematics = FK()  # type: FK
print("FK initialization completed")


def handle_calculate_IK(req):
    rospy.loginfo("Received %s eef-poses from the plan" % len(req.poses))
    if len(req.poses) < 1:
        print "No valid poses received"
        return -1
    else:
        # Initialize service response
        joint_trajectory_list = []
        for x in xrange(0, len(req.poses)):
            # IK code starts here
            joint_trajectory_point = JointTrajectoryPoint()

            # Extract end-effector position and orientation from request
            # px,py,pz = end-effector position
            # roll, pitch, yaw = end-effector orientation
            px = req.poses[x].position.x
            py = req.poses[x].position.y
            pz = req.poses[x].position.z

            (roll, pitch, yaw) = tf.transformations.euler_from_quaternion(
                [req.poses[x].orientation.x, req.poses[x].orientation.y,
                 req.poses[x].orientation.z, req.poses[x].orientation.w])

            ### Your IK code here
            # Compensate for rotation discrepancy between DH parameters and Gazebo
            rot_ee = _forward_kinematics.corrected_ee_rot_matrix.subs({'r': roll, 'p': pitch, 'y': yaw})

            # Calculate wrist center position from EE position.
            # 0_r_WC/0 = 0_r_EE/0 - d * 0-6R * ([0][0][1])
            EE = Matrix([[px], [py], [pz]])
            WC = EE - (0.303) * rot_ee[:, 2]

            # Calculate joint angles using Geometric IK method
            theta1 = atan2(WC[1], WC[0])
            # SSS triangle for theta2, theta3
            side_a = 1.501
            side_b = sqrt(pow((sqrt(WC[0] * WC[0] + WC[1] * WC[1]) - 0.35), 2) + pow((WC[2] - 0.75), 2))
            side_c = 1.25

            angle_a = acos((side_b * side_b + side_c * side_c - side_a * side_a) / (2 * side_b * side_c))
            angle_b = acos((side_a * side_a + side_c * side_c - side_b * side_b) / (2 * side_a * side_c))
            angle_c = acos((side_a * side_a + side_b * side_b - side_c * side_c) / (2 * side_a * side_b))

            theta2 = pi / 2. - angle_a - atan2(WC[2] - 0.75, sqrt(WC[0] * WC[0] + WC[1] * WC[1]) - 0.35)
            ## 0.036 accounts for sag in link4 of -0.054m
            theta3 = pi / 2. - (angle_b + 0.036)

            r0_3 = _forward_kinematics.t0_1[0:3, 0:3] \
                   * _forward_kinematics.t1_2[0:3, 0:3] \
                   * _forward_kinematics.t2_3[0:3, 0:3]
            r0_3 = r0_3.evalf(
                subs={_forward_kinematics.q1: theta1, _forward_kinematics.q2: theta2, _forward_kinematics.q3: theta3})
            r3_6 = r0_3.inv("LU") * rot_ee

            theta4 = atan2(r3_6[2, 2], -r3_6[0, 2])
            theta5 = atan2(sqrt(r3_6[0, 2] * r3_6[0, 2] + r3_6[2, 2] * r3_6[2, 2]), r3_6[1, 2])
            theta6 = atan2(-r3_6[1, 1], r3_6[1, 0])

            ###

            # Populate response for the IK request
            # In the next line replace theta1,theta2...,theta6 by your joint angle variables
            joint_trajectory_point.positions = [theta1, theta2, theta3, theta4, theta5, theta6]
            joint_trajectory_list.append(joint_trajectory_point)

        rospy.loginfo("length of Joint Trajectory List: %s" % len(joint_trajectory_list))
        return CalculateIKResponse(joint_trajectory_list)


def IK_server():
    # initialize node and declare calculate_ik service

    rospy.init_node('IK_server', log_level=rospy.DEBUG)
    s = rospy.Service('calculate_ik', CalculateIK, handle_calculate_IK)
    print "Ready to receive an IK request"
    rospy.spin()

    if __name__ == "__main__":
        IK_server()
