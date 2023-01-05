DREAMPar.d = 15;                          % Dimension of the problem
DREAMPar.N = 10;                         % Number of Markov chains
DREAMPar.T = 15000;                       % Number of generations
DREAMPar.lik = 2;                       % Model output is simulation: Gaussian likelihood function
Par_info.prior = 'latin';                 % Latin hypercube sampling
Par_info.boundhandling = 'reflect';         % Explicit boundary handling
% A1=[1 1];B1=[3 3];A2=[2 2];B2=[4 4]; 
% Par_info.min =[repmat(A1,1,5) repmat(A2,1,5)];
% Par_info.max =[repmat(B1,1,5) repmat(B2,1,5)];  

A1=[0.5 0.5];B1=[1.5 1.5];
Par_info.min =[repmat(A1,1,7) 0.5];
Par_info.max =[repmat(B1,1,7) 1.5]; 
%%RUN PROGRAM
Func_name = 'objfunc';
%% Run the DREAM algorithm
[chain,output,fx] = DREAM(Func_name,DREAMPar,Par_info);
