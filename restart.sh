#!/bin/bash
# File: restart.sh
# Author: Seongeun Kim (eunbelivable@snu.ac.kr)

# Variables
DOCKER_ID="seamustard52"
PREFIX="introbioinfo"
EXERCISE_NAME="exercise06"
EXERCISE_DIR=$PWD

# Run docker
docker restart $EXERCISE_NAME
docker exec -it $EXERCISE_NAME /bin/bash