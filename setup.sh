#!/bin/bash
# File: setup.sh
# Author: Hyunbin Kim (khb7840@gmail.com)

# Variables
DOCKER_ID="khb7840"
PREFIX="introbioinfo"
EXERCISE_NAME="exercise04"
EXERCISE_DIR=$PWD

# Run docker
docker run -it --mount type=bind,src=$PWD,target=/home/$EXERCISE_NAME/$EXERCISE_NAME \
--name $EXERCISE_NAME "$DOCKER_ID/$PREFIX-$EXERCISE_NAME"
