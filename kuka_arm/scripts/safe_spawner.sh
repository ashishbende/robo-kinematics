#!/bin/bash
# This script safely launches ros nodes with buffer time to allow param server population
mate-terminal --tab -e "roslaunch kuka_arm target_description.launch" &
sleep 3 &&
mate-terminal --tab -e "roslaunch kuka_arm cafe.launch" &
sleep 3 &&
mate-terminal --tab -e "roslaunch kuka_arm spawn_target.launch" &
sleep 5 &&
#x-terminal-emulator -e roslaunch kuka_arm inverse_kinematics.launch
roslaunch kuka_arm inverse_kinematics.launch
