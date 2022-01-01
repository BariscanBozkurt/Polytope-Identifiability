# Polytope-Identifiability

Official Codes for "On Identifiable Polytope Characterization for Polytopic Matrix Factorization"

# Identifiable Polytope Characterization

# Implementation Environment Details

-> Python Version: Python 3.7.11 (default, Jul 27 2021, 14:32:16)

-> Required Python Packages: Specified in requirements.txt file.

-> Platform Info : "OS: Linux (x86_64-pc-linux-gnu) CPU: Intel(R) Xeon(R) Gold 6248 CPU @ 2.50GHz"

# SageMath Environment Usage
######################################################################

ALTERNATIVE 1 - USING WINDOWS SUBSYSTEM FOR LINUX (wsl) IN WINDOWS 10

STEP1 : Download Ubuntu 18.04 LTS from Microsoft Store. We will call this app as wsl from now on.

STEP2 : Follow the steps from https://gist.github.com/kauffmanes/5e74916617f9993bc3479f401dfec7da
	to install Anaconda in wsl.

STEP3 : Follow the steps from https://doc.sagemath.org/html/en/installation/conda.html to create
	an Anaconda environment with SageMath installed. Moreover, you may need to install the 
	python libraries submitted in requirements.txt

STEP4: Activate environment using one of the following command depending on your environment
	
	"source activate <your_env_name>"
	"conda activate <your_env_name>"


######################################################################

ALTERNATIVE 2 - USE CoCalc (Collaborative Calculation and Data Science)

STEP1 : Go to https://cocalc.com/projects?session=default

STPE2 : Create a new project with an extension of a Jupyter Notebook

STEP3 : Select the kernel as SageMath 9.2

STEP4 : They will give a Jupyter Notebook with SageMath kernel which has 30 min connection.
	You can copy paste the required library and functions from the submitted codes. Then, you can play with those codes in a web-basd cloud computing
	service provided.


