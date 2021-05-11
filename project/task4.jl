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
for i = 1:n
    # compute the realized return R_i(t)
    stocks_retur[:,i] = (dfs[i].close-dfs[i].open) ./ dfs[i].open;
end
names_stocks = [ dfs[i].Name[1] for i in 1:n ];
stock_id = 1; # You may change this stock_id to different numbers from 1 to 20
plot( dfs[stock_id].date, stocks_retur[:,stock_id] , label = dfs[stock_id].Name[1], title = dfs[stock_id].Name[1]*"'s return over time" )
stock_id = 8; # You may change this stock_id to different numbers from 1 to 20
plot( dfs[stock_id].date, stocks_retur[:,stock_id] , label = dfs[stock_id].Name[1], title = dfs[stock_id].Name[1]*"'s return over time" )
stock_id = 15; # You may change this stock_id to different numbers from 1 to 20
plot( dfs[stock_id].date, stocks_retur[:,stock_id] , label = dfs[stock_id].Name[1], title = dfs[stock_id].Name[1]*"'s return over time" )
stock_id = 20; # You may change this stock_id to different numbers from 1 to 20
plot( dfs[stock_id].date, stocks_retur[:,stock_id] , label = dfs[stock_id].Name[1], title = dfs[stock_id].Name[1]*"'s return over time" )