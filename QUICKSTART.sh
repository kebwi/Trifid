#!/bin/bash

# Assuming you have successfully built Trifid, this bash script will automatically generate a series of demonstrative images that you can use to quickly get a feel for Trifid's capabilities.  Run it by simply typing "./QUICKSTART.sh" on the command line from the project directory.

# This script can be run on any 512x512 ppm image by passing the filename to the script (% ./QUICKSTART.sh inputImageFilename).  If no filename is provided, the included sample image 'gazebo.ppm' will be used by default.

# Look for an argument to this script and use it as the input filename.  If no argument is found, use 'gazebo.ppm'.
if [ -n "$1" ]  # If command-line argument present...
then
	input="$1"
else
	input="images/gazebo.ppm"
fi

# If the Trifid program is found and is an executable file, generate a series of noisy images on the input image.
if [ -e ./Trifid ]
	then
	if [ -x ./Trifid ]
		then
		echo "Generating a series of simulated images from $input..."
		./Trifid -N noiseParameters/1_day -b $input
		./Trifid -N noiseParameters/2_duskHandheld -b $input
		./Trifid -N noiseParameters/3_eveningHandheld -b $input
		./Trifid -N noiseParameters/5_fullMoonHandheld -b $input
		./Trifid -N noiseParameters/7_quarterMoonHandheld -b $input
		echo "...Done generating a series of simulated images from $input."
	else echo "./Trifid is not executable.  Try typing 'chmod 755 Trifid'."
	fi
else echo "./Trifid was not found.  Have you built it yet?  Try typing 'make'."
fi

# After this script has run, you might consider investigating the PSNR and Lab Delta E statistics of the results.  See the README and TUTORIAL for directions on how to do this.
