'''
Written by Justin Cook
Initially created on 9 October 2019
Released on 11 October 2019
Purpose:
- This is a helper program to use with the Denavit-Hartenberg convention that:
-- Calculates the  transformation matrix for each rotation/translation
-- Calculates the combined transformation matrix, displaying each step in finding it
- Performs all matrix calculations with any combination of numerical parameters and joint variable parameters by using Sympy
- Outputs the calculated matrices in a LaTeX amsmath bmatrix compatible format
- Allows easy file based input where the input needs only consist of (# of frames) lines of the DH parameters (d,a,theta,alpha)
Requirements:
- Python 3
- Sympy (available through pip3)
'''

from math import *
import sympy as sym
import sys

def info():
	print("""
Info on kinematics_DH_solver
-------------------------------
General Info
--------------------------
- This method called with "python kinematics_DH_solver.py info"
- For a short version of a help message (requirements, arguments, input, output), use "python kinematics_DH_solver.py help"
- For a detailed version of a help message, use "python kinematics_DH_solver.py help_long"

Intended Use
--------------------------
- The solver is designed to not be used piecewise (i.e. used through main) through an input file
- The input file is to be a list of several lines with each line describing one set of comma-separated joint parameters
- The parameters can be specified either numerically (i.e. 10) or symbolically (i.e. theta)
- By default, the solver is designed to output matrices in a LaTeX amsmath compatible way
- For this reason, the solver automatically prepends "\\" for symbolic expressions

Possible Use
--------------------------
- The allocation and printing functions can be used independently
- To use the allocation functions then just simply pass the 4 parameters or an array of the 4 parameters
- allocate_AiM is a convince function that just calls allocate_Ai by splitting the input array into the 4 parameters
- The printing functions simply split a 1d or 2d matrix into a LaTeX friendly and so can be used for any Sympy Matrix
	""")

def help():
	print("""
Help for using kinematics_DH_solver
-----------------------------------
- Python3 and Sympy are required
- Call with "python kinematics_DH_solver.py <optional input file> <optional output file (input file must be specified)>"
- If no arguments are given, then the program will ask for the file to be typed
- Input file specified as # of frames lines of d,a,theta,alpha
- Output file will consist of A matrices (# = # of joints), T matrices (# = # of A-1)
- Matrices are outputted in a way designed for that they can be directly copied into a LaTeX amsmath bmatrix environment
	""")

def help_long():
	print("""
Detailed Help for using kinematics_DH_solver
--------------------------------------------
General Info
--------------------------
- For information on the program, use "python kinematics_DH_solver.py info"
- For a short version of a help message (requirements, arguments, input, output), use "python kinematics_DH_solver.py help"
- This method called with "python kinematics_DH_solver.py help_long"
- Python3 and Sympy are required
- Call with "python kinematics_DH_solver.py <optional input file> <optional output file (input file must be specified)>"
- If no arguments are given, then the program will ask for the file to be typed
- Input file specified as # of frames lines of d,a,theta,alpha
- Output file will consist of A matrices (# = # of joints), T matrices (# = # of A-1) (Note: file will be overwritten)
- Matrices are outputted in a way designed for that they can be directly copied into a LaTeX amsmath bmatrix environment

Parameter Specification
-----------------------
- Each variable can be specified numerically or symbolically
- If a parameter has a numerical value then just enter the value (will be converted to a float)
- If a parameter has a symbolic value then enter the variable name without the joint number
- Symbolic variables will be converted to display properly in LaTeX (i.e. theta will become \\theta_n)
- Since the program has no way to know if an asterisk should be used, these must be specified (i.e. theta^*)
- Angles should be specified in degrees

Examples
-----------------------
- Example call: "python kinematics_DH_solver.py DATA OUT" will use input file of DATA and output (and overwrite) file OUT
- TODO
	""")

def allocate_AiM(matrix):
	return allocate_Ai(matrix[0], matrix[1], matrix[2], matrix[3])

def allocate_Ai(d, a, theta, alpha):
	try:
		ctheta = round(cos(radians(theta)), 4)
		stheta = round(sin(radians(theta)), 4)
	except Exception:
		ctheta = sym.Symbol("C" + str(theta))
		stheta = sym.Symbol("S" + str(theta))
	try:
		calpha = round(cos(radians(alpha)), 4)
		salpha = round(sin(radians(alpha)), 4)
	except Exception:
		calpha = sym.Symbol("C" + str(alpha))
		salpha = sym.Symbol("S" + str(alpha))
	Ai = sym.Matrix(
		[[ctheta, -stheta*calpha, stheta*salpha, a*ctheta],
		[stheta, ctheta*calpha, -ctheta*salpha, a*stheta],
		[0, salpha, calpha, d],
		[0, 0, 0, 1]]
	)
	return Ai

def matrix1d_print(matrix, sep):
	for i in range(matrix.shape[0]):
		print(matrix[i],sep,end='')

def matrix2d_print(matrix, sep="& ", fname=None):
	matrix_rows = matrix.shape[0]
	matrix_cols = matrix.shape[1]
	for i in range(matrix_rows):
		for j in range(matrix_cols-1):
			print((matrix.row(i))[j],sep,end='',file=fname)
		if i == (matrix_rows - 1):
			print((matrix.row(matrix_rows-1))[matrix_cols-1], "",file=fname)
		else:
			print((matrix.row(i))[matrix_cols-1], "\\\\",file=fname)

def main():
	args = sys.argv[1:]
	# Check if user is getting help or info
	if len(args) == 1:
		if args[0] == "help":
			help()
			return
		elif args[0] == "help_long":
			help_long()
			return
		elif args[0] == "info":
			info()
			return
	# Get file and then create a matrix of the parameters, adding joint identifier to each symbolic variable
	try:
		in_file = open(args[0], "r")
	except Exception:
		in_file = open(input("Please enter the input filepath: "), "r")
	parameters = []
	line_number = 1
	for line in in_file:
		values = line.rstrip().split(",")
		for i in range(len(values)):
			try:
				values[i] = float(values[i])
			except Exception:
				if values[i] == "theta" or values[i] == "alpha":
					values[i] = sym.Symbol("\\" + values[i] + "_" + str(line_number))
				else:
					values[i] = sym.Symbol(values[i] + "_" + str(line_number))
		parameters.append(values)
		line_number += 1
	in_file.close()

	# Calculate the A matrices and then the Transformation matrices, outputs to a file if specified
	try:
		out_file = open(args[1], "w")
	except Exception:
		out_file = None
	As=[]
	for Avalues in parameters:
		A = allocate_AiM(Avalues)
		As.append(A)
		matrix2d_print(A, fname=out_file)
	T=As[0]
	for i in range(1,len(As)):
		T=T*As[i]
		matrix2d_print(T, fname=out_file)
	if out_file != None:
		out_file.close()

main()

