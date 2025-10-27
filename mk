#!/bin/bash

function check_okay {
	if [ $? -ne 0 ]
	then
		echo
		echo "FAILED"
		echo
		exit 1
	fi
}


cd make; make TAG=$1 USE_MPI=$2; check_okay

