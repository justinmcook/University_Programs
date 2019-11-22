'''
Written by Justin Cook
Initially created on 9 October 2019
Released on 11 October 2019
Purpose:

Requirements:
- Python 3
- Sympy (available through pip3)
'''

from math import *
import sympy as sym
import sys

def info():
	print("""
Info on jacobian_solver
-------------------------------
General Info
--------------------------
- This method called with "python jacobian_solver.py info"
- For a short version of a help message (requirements, arguments, input, output), use "python jacobian_solver.py help"
- For a detailed version of a help message, use "python jacobian_solver.py help_long"

Intended Use
--------------------------


Possible Use
--------------------------

	""")

def help():
	print("""
Help for using jacobian_solver
-----------------------------------
- Python3 and Sympy are required
	""")

def help_long():
	print("""
Detailed Help for using jacobian_solver
--------------------------------------------
General Info
--------------------------
- For information on the program, use "python jacobian_solver.py info"
- For a short version of a help message (requirements, arguments, input, output), use "python jacobian_solver.py help"
- This method called with "python jacobian_solver.py help_long"
- Python3 and Sympy are required

Parameter Specification
-----------------------


Examples
-----------------------

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

