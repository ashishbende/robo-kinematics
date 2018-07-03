#!/usr/bin/env python

# Copyright (C) 2017 Electric Movement Inc.
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


### Define functions for homogeneous transform matrices of rotations about x, y, and z given specific angle.
def rot_x(q):
	R_x = Matrix([[1, 0, 0],
	              [0, cos(q), -sin(q)],
	              [0, sin(q), cos(q)]])

	return R_x


def rot_y(q):
	R_y = Matrix([[cos(q), 0, sin(q)],
	              [0, 1, 0],
	              [-sin(q), 0, cos(q)]])
	return R_y


def rot_z(q):
	R_z = Matrix([[cos(q), -sin(q), 0],
	              [sin(q), cos(q), 0],
	              [0, 0, 1]])
	return R_z


### Creates Homogeneous Transform Matrix
def homogeneous_transform(q, d, a, alpha):
	T = Matrix([[cos(q), -sin(q), 0, a],
	            [sin(q) * cos(alpha), cos(q) * cos(alpha), -sin(alpha), -sin(alpha) * d],
	            [sin(q) * sin(alpha), cos(q) * sin(alpha), cos(alpha), cos(alpha) * d],
	            [0, 0, 0, 1]])
	return T


def handle_calculate_IK(req):
	rospy.loginfo("Received %s eef-poses from the plan" % len(req.poses))
	if len(req.poses) < 1:
		print "No valid poses received"
		return -1
	else:
		### Your FK code here
		# Create symbols
		print 'Creating symbols'
		# Joint angle (variable for revolute)
		q1, q2, q3, q4, q5, q6, q7 = symbols('q1:8')
		# Link offset (variable for prismatic)
		d1, d2, d3, d4, d5, d6, d7 = symbols('d1:8')
		# Link length
		a0, a1, a2, a3, a4, a5, a6 = symbols('a0:7')
		# Twist angles
		alpha0, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6 = symbols('alpha0:7')

		# Create Modified DH parameters
		DH = {alpha0: 0, a0: 0, d1: 0.75, q1: q1,
		      alpha1: -pi / 2, a1: 0.35, d2: 0, q2: q2 - pi / 2,
		      alpha2: 0, a2: 2, d3: 0, q3: q3,
		      alpha3: -pi / 2, a3: 0.06, d4: 1.5, q4: q4,
		      alpha4: pi / 2, a4: 0, d5: 0, q5: q5,
		      alpha5: -pi / 2, a5: 0, d6: 0, q6: q6,
		      alpha6: 0, a6: 0, d7: 0.453, q7: 0}

		# Define Modified DH Transformation matrix
		# Create individual transformation matrices
		print 'Creating HTs'
		T0_1 = homogeneous_transform(q1, d1, a0, alpha0).subs(DH)
		T1_2 = homogeneous_transform(q2, d2, a1, alpha1).subs(DH)
		T2_3 = homogeneous_transform(q3, d3, a2, alpha2).subs(DH)
		T3_4 = homogeneous_transform(q4, d4, a3, alpha3).subs(DH)
		T4_5 = homogeneous_transform(q5, d5, a4, alpha4).subs(DH)
		T5_6 = homogeneous_transform(q6, d6, a5, alpha5).subs(DH)
		T6_G = homogeneous_transform(q7, d7, a6, alpha6).subs(DH)

		# Extract rotation matrices from the transformation matrices
		# Generalized homogeneous transform
		T0_G = T0_1 * T1_2 * T2_3 * T3_4 * T4_5 * T5_6 * T6_G

		# Correcting difference in orientation between DH convention and URDF file
		r, p, y = symbols('r p y')
		R_x = rot_x(r)
		R_y = rot_y(p)
		R_z = rot_z(y)

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
			print 'Processing pose: ' + str(x)
			R_corr = R_z.subs(y, pi) * R_y.subs(p, -pi / 2)
			Rrpy = R_z * R_y * R_x * R_corr
			Rrpy = Rrpy.subs({'r': roll, 'p': pitch, 'y': yaw})

			# Calculate joint angles using Geometric IK method
			# Calculate the wrist center position first
			print 'Calculating wrist center'
			nx = Rrpy[0, 2]
			ny = Rrpy[1, 2]
			nz = Rrpy[2, 2]

			# d7 = 0.453
			wx = px - 0.453 * nx
			wy = py - 0.453 * ny
			wz = pz - 0.453 * nz

			print 'Calculating position angles'
			# Theta 1 from above is a simple atan2
			theta1 = atan2(wy, wx)
			# Calculate Radius from above
			r = sqrt(wx ** 2 + wy ** 2) - 0.35  # a1: 0.35

			# Calculating Theta 2 and 3 using cosine law. A, B and C are sides of the triangle
			A = 1.501
			B = sqrt(r ** 2 + (wz - 0.75) ** 2)  # d1: 0.75
			C = 1.25

			# a corresponds to angle apha
			a = acos((B ** 2 + C ** 2 - A ** 2) / (2 * B * C))
			theta2 = pi / 2 - a - atan2(wz - 0.75, r)  # d1: 0.75

			# b corresponds to angle beta
			b = acos((A ** 2 + C ** 2 - B ** 2) / (2 * A * C))
			theta3 = pi / 2 - (b + 0.036)

			print 'Calculating orientation angles'
			# Calculating Euler angles from orientation
			R0_3 = T0_1[0:3, 0:3] * T1_2[0:3, 0:3] * T2_3[0:3, 0:3]
			R0_3 = R0_3.evalf(subs={'q1': theta1, 'q2': theta2, 'q3': theta3})
			# R3_6 = R0_3.inv("LU")*Rrpy # using transpose instead:
			R3_6 = R0_3.T * Rrpy

			theta5 = atan2(sqrt(R3_6[0, 2] ** 2 + R3_6[2, 2] ** 2), R3_6[1, 2])
			# Choosing between multiple possible solutions:
			if sin(theta5) < 0:
				theta4 = atan2(-R3_6[2, 2], R3_6[0, 2])
				theta6 = atan2(R3_6[1, 1], -R3_6[1, 0])
			else:
				theta4 = atan2(R3_6[2, 2], -R3_6[0, 2])
				theta6 = atan2(-R3_6[1, 1], R3_6[1, 0])

			print 'Appending calculated values'
			# Populate response for the IK request
			# In the next line replace theta1,theta2...,theta6 by your joint angle variables
			joint_trajectory_point.positions = [theta1, theta2, theta3, theta4, theta5, theta6]
			joint_trajectory_list.append(joint_trajectory_point)

		rospy.loginfo("length of Joint Trajectory List: %s" % len(joint_trajectory_list))
		return CalculateIKResponse(joint_trajectory_list)


def IK_server():
	# initialize node and declare calculate_ik service
	rospy.init_node('IK_server')
	s = rospy.Service('calculate_ik', CalculateIK, handle_calculate_IK)
	print "Ready to receive an IK request"
	rospy.spin()


if __name__ == "__main__":
	IK_server()
