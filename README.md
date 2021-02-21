# ffi_dia
The master branch for the updated Vanderbilt DIA pipeline.

This code is an updated version of the difference imaging pipeline released in 
Oelkers & Stassun 2018. 

Please note, in order ot run the code successfully, users will need to update 
various components of config.py. These include:

Directories:

    1. WORKING_DIRECTORY - this is the location of the code on your local machine.
    2. DATA_DIRECTORY - this is the location of the TESS data on your local machine.

Machine:

    1. Update the machine to be the specific machine you are using. For stassunlab
    users, this is the stassunlab machine.
    
Skip Steps:

    1. There are various steps to the process. If you would like to skip a step,
    change the step skip flag from N to Y. 
    
Sector:
    
    1. TESS has various numerical codes for sectors, cameras, and ccds. Please make
    sure they are the correct before running the code.

difference.py

    1. Please note, you may need to update line 134 in difference.py in order to point
    to the correct location of the cfitsio libraries.
    
I do not recommend changing any other configuration flags.  