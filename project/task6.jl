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
# the following code specifies the constants as described in the problem
M = 20; B = 20; c = 2; w = 1*ones(n); Rd = 1.01*sum( w.*bar_R ); 

# the following code setup the JuMP model with the right solver
nl_solver = optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0)
optimizer = Juniper.Optimizer
model = Model(optimizer_with_attributes(optimizer, "nl_solver"=>nl_solver, "atol"=>1e-10));
# your code here
n = 20;
@variable(model, k[1:n]);
@variable(model, y[1:n], Bin);
@variable(model, x_mip[1:n]);
@constraint(model, k .== sqrt(Sigma)*(x_mip+w));
@constraint(model, sum((x_mip+w) .* bar_R) >= Rd);
@constraint(model, sum(x_mip+y.*c)<=B);
@constraint(model, x_mip .>= -M.*y);
@constraint(model, x_mip .<= M.*y);
@constraint(model, x_mip .>= -w);
@NLobjective(model, Min, sum(k[i]^2 for i in 1:n));
optimize!(model);
# specify the problem parameters
M = 20; B = 20; a = 2; w = ones(n); Rd = 1.01*sum( w.*bar_R ); 
# specify the JuMP model with ECOS as the optimizer
m_socp = Model( ECOS.Optimizer );
# your code here
using LinearAlgebra
n = 20;
a1 = zeros(n+1,n+1);
b1 = [Rd;zeros(n)];
c1 = [0;bar_R];
d1 = bar_R'*w;
a2 = [zeros(n) sqrt(a).*Array{Float64}(I,n,n)];
b2 = 1/(2*sqrt(a)).*ones(n);
d2 = sqrt(B+n/(4*a));
a3 = [zeros(n) sqrt(Sigma)];
b3 = sqrt(Sigma)*w;
c3 = [1;zeros(n)];
@variable(m_socp,x_socp[1:n+1]);
@constraint(m_socp,[c1'*x_socp+d1;a1*x_socp+b1] in SecondOrderCone());
@constraint(m_socp,[d2;a2*x_socp+b2] in SecondOrderCone());
@constraint(m_socp,[c3'*x_socp;a3*x_socp+b3] in SecondOrderCone());
@constraint(m_socp,[i=2:n+1],-M <= x_socp[i] <= M);
@constraint(m_socp,[i=2:n+1],x_socp[i] >= -w[i-1]);
@objective(m_socp, Min, x_socp[1]);
optimize!(m_socp)
plot( names_stocks, portfolio_opt, labels = "unconstrained portfolio" )
plot!( names_stocks, JuMP.value.(x_mip + w), labels = "MI-SOCP")
plot!( names_stocks, JuMP.value.(x_socp[2:n+1] + w), labels = "SOCP")
plot!( names_stocks, 1000*bar_R, labels ="(Scaled) Expected Return")
plot!( names_stocks, 1000*[Sigma[i,i] for i in 1:n], labels ="(Scaled) Variance")
sharpe_IP = sharpe_ratio( path_subgroup, JuMP.value.(x_mip), ones(n), 2 );
sharpe_SOCP = sharpe_ratio( path_subgroup, JuMP.value.(x_socp[2:n+1]), ones(n) , 2 );
sharpe_Opt = sharpe_ratio( path_subgroup, portfolio_opt, zeros(n), 2 );