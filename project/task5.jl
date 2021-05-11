# Load the JuMP related packages and several solvers
# ECOS - for solving SOCP problems
# Juniper & Ipopt - for solving MI-NLP problems
using JuMP, Juniper, ECOS, Ipopt
# Load the data/file processing related packages
using CSV, Glob, DataFrames, Statistics
# Load the Plot package for illustrating the solution
using Plots
# Load the custom functions for benchmarking  
include("./reusablefunc.jl");
# Change "subgroup1" to other names according to the subgroup you are assigned.
path_subgroup = "./ftec_project_subgroup2/" 
files = glob( "*_train.csv", path_subgroup );
dfs = DataFrame.( CSV.File.( files ) );
T = 800; n = length(dfs);
stocks_retur = zeros(T,n);
# calculate r_i and Sigma
bar_R = [ mean( stocks_retur[:,i] ) for i in 1:length(dfs) ];
Sigma = [ mean( (stocks_retur[:,i].-bar_R[i]).*(stocks_retur[:,j].-bar_R[j]) ) for i=1:n, j=1:n ]; 
# your code here
B = 20;
one_vec = ones(20);
Rd = 1.01 * sum(bar_R);
r0 = bar_R'*Sigma^-1*bar_R;
r1 = one_vec'*Sigma^-1*bar_R;
r2 = one_vec'*Sigma^-1*one_vec;
portfolio_opt = Sigma^-1 *((r0*B-r1*Rd)*one_vec+(r2*Rd-r1*B)*bar_R)/(r0*r2-r1^2);
plot( names_stocks, portfolio_opt, labels = "portfolio", xticks = :all )
# you may adjust the scale factor "1000" to scale up/down the expected return to make it comparable with 
# the value of the portfolio (*for improved visualization only*).
plot!( names_stocks, 1000*bar_R , labels = "(scaled) expected return") 
plot!( names_stocks, 1000*[Sigma[i,i] for i in 1:n], labels = "(scaled) variance" )
