#!/bin/bash
# File: restart.sh
# Author: Hyunbin Kim (khb7840@gmail.com)

# Variables
DOCKER_ID="khb7840"
PREFIX="introbioinfo"
EXERCISE_NAME="exercise04"
EXERCISE_DIR=$PWD

# Run docker
docker restart $EXERCISE_NAME
docker exec -it $EXERCISE_NAME /bin/bash
