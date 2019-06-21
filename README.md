# 2_5DGreenFunctionRadar
This repository contains all the files to calculated the Green function
for the 2.5 D geometry of radar wave modeling and to reproduced all the results
presented in the paper:

Accurate 2.5D frequency domain radar waves modeling using 
weighted-averaging diﬀerence operators 

Brief description of files:

1- main.m 
			The main file. Most of the other functions are called from this main file
			
2- GenerateModel.m
			A useful function to generate three simple parameter models to test the code
			The numerical solution can be calculated for the homogeneous model and the 3 
			layer model and compare with the analytical solution.
			
3- InitializeForwardModeling.m
			A function that will take the input modeling setup of the user and translate it	
			for the function EM_NumericalGreenFunction(). This InitializeForwardModeling()
			function can plot the initial parameter model and the numerical grid that is used
			to solve for the Green's functions. Position of sources and receptors are also 
			indicated on the plot. 

4- EM_NumericalGreenFunction.m
			This function solves equation 6 and 7. The sparse matrix is constructed 
			with RadarSparseMatrix (standard coefficients) or 
			RadarSparseMatrixOpt (optimal coefficients)
			
5a RadarSparseMatrix.m
			This function calculates the elements of the sparse matrix associated with eq 6
			(standard coefficients formulation)
			
5b RadarSparseMatrix_mex.mexw64
			MEX file associated to RadarSparseMatrix.m 

5c RadarSparseMatrix.prj
			Project coder to build the mex file from RadarSparseMatrix.m 
			
5d Call_RadarSparseMatrix,m
			Dummy function to call RadarSparseMatrix.m when building the mex file
			(useful only when building the mex file the first time)
			
6a RadarSparseMatrixOpt
			This function calculates the elements of the sparse matrix associated with eq 6
			(optimal coefficients formulation)
			
6b RadarSparseMatrixOpt.mexw64
			MEX file associated to RadarSparseMatrixOpt.m 
			
6c RadarSparseMatrixOpt.prj
			Project coder to build the mex file from RadarSparseMatrixOpt.m 
			
6d Call_RadarSparseMatrixOpt.m
			Dummy function to call RadarSparseMatrixOpt.m when building the mex file
			(useful only when building the mex file the first time)
			
7- AnalyticalSolutionReceptorsHomogeneous.m
			Analytical solution for the homogeneous parameters
			
8- AnalyticalSolutionReceptorsThreeLayers.m
			Analytical solution for the 3 layers parameters
			
9- GreenFunctionErrorReceptors.m
			A plot to compare analytical solution to numerical solution
			
10- GreenFunctionErrorReceptorsComp.m
			A plot to compare the solution obtained with the standard 
			coefficients and the optimal coefficients
