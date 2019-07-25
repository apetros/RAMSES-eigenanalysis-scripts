% EIGS PSCC 2016
%
% Files
%   analyze_results - Analyzes the results of the eigenvalue/eigenvector computation
%   eigenvals_eig   - Computing egenvalues with the use of the state matrix and Matlab's eig function
%   eigenvals_eigs  - Computing egenvalues with the use of the descriptor matrix and Matlab's Eigs function (ARPACK Arnoldi)
%   eigs_solver     - Implements y = A\x, needed for the Arnoldi method (eigs) using different solvers
%   init            - Reads the descriptor matrix in raw format and creates the internal data structure
%   loop_analysis   - Creates an interactive selection list to analyze the output of the eigenvalue computations
%   ssa             - Receives the descriptor form of an eigenvalue problem and calculates the eigenvalues/eigenvectors with different methods
