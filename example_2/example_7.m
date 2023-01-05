DREAMPar.d = 30;                          % Dimension of the problem
DREAMPar.N = 10;                         % Number of Markov chains
DREAMPar.T = 40000;                       % Number of generations
DREAMPar.lik = 2;                       % Model output is simulation: Gaussian likelihood function
Par_info.prior = 'latin';                 % Latin hypercube sampling
Par_info.boundhandling = 'reflect';         % Explicit boundary handling
A=[0.7 0.7];B=[1.2 1.2]; 
Par_info.min =repmat(A,1,15);
Par_info.max =repmat(B,1,15);  
%%RUN PROGRAM
Func_name = 'objfunc';
%% Run the DREAM algorithm
[chain,output,fx] = DREAM(Func_name,DREAMPar,Par_info);
