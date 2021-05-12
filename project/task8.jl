# using XLXS
# Load the JuMP related packages and several solvers
# ECOS - for solving SOCP problems
# Juniper & Ipopt - for solving MI-NLP problems
using JuMP, Juniper, ECOS, Ipopt
# Load the data/file processing related packages
using CSV, Glob, DataFrames, Statistics
# Load the Plot package for illustrating the solution
using Plots
# Load the custom functions for benchmarking  
include("./func.jl");
# load the full data set!
files = glob( "*_train.csv", "./ftec_project_files/");
dfs = DataFrame.( CSV.File.( files ) );
T = 800; n = length(dfs);
stocks_retur_full = zeros(T,n);
for i = 1:n
    # compute the realized return R_i(t)
    stocks_retur_full[:,i] = (dfs[i].close-dfs[i].open) ./ dfs[i].open;
end
names_stocks_full = [ dfs[i].Name[1] for i in 1:n ];
# calculate r_i and Sigma
bar_R_full = [ mean( stocks_retur_full[:,i] ) for i in 1:length(dfs) ];
Sigma_full = [ mean( (stocks_retur_full[:,i].-bar_R_full[i]).*(stocks_retur_full[:,j].-bar_R_full[j]) ) for i=1:n, j=1:n ]; 
# your code/functions here
B = 20;
Rd = 1.01*sum(bar_R_full);
function calc_gamma(delta)
    delta = 0.01;
    L = a/delta + max(diag(Sigma_full));
    gamma = 1/L;
    return gamma;
end

function gd(x,w,a,lambda,upsilon,delta)
    @assert length(x)==n;
    grad = Sigma_full*x+Sigma_full*w - lambda.*bar_R_full+upsilon.*ones(n);
    for i=1:n
        if abs(x[i])> delta
            grad[i] += sign(x[i])*upsilon*a;
        else
            grad[i] += upsilon*(a/delta)*x[i];
        end
    end
    return grad;
end

function obj_v(x,w,a,lambda,upsilon,delta)
    @assert length(x)==n;
    obj = 1/2*((x+w)'*Sigma_full*(x+w))+lambda*(Rd-bar_R_full'*x-bar_R_full'*w)+(upsilon*sum(x)-upsilon*B);
    for i=1:n
        if abs(x[i])>delta
            obj += upsilon*(a*abs(x[i])-a*delta/2);
        else
            obj += upsilon*(a/(2*delta)*x[i]^2);
        end
    end
    return obj;
end
# set the parameters as specified by the problem
M = 20; w = 1*ones(n); a = 1; delta = 0.01;

# calculate l and u - your code here (should be a simple formula)
u = M.*ones(n);
l = -w;
outfile = "output.txt";
# initialize the algorithm
for lambda = 1:10:70000
    for upsilon = 1:1000
        x_custom = zeros(n); 
        store_obj = []
        push!(store_obj, obj_v(x_custom,w,a,lambda,upsilon,delta) ) # replace ".." with the function you wrote for computing the objective val.

        for iteration_no = 1 : 50000 # feel free to adjust the number of iterations run here.
            # your code here
            gamma = 1/(iteration_no+1);
            x_custom -= gamma* gd(x_custom,w,a,lambda,upsilon,delta);
            for i=1:n
                if x_custom[i] > u[i]
                    x_custom[i] = u[i];
                elseif x_custom[i] < l[i]
                    x_custom[i] = l[i];
                end
            end
            push!(store_obj, obj_v(x_custom,w,a,lambda,upsilon,delta) ) # replace ".." with the function you wrote for computing the objective val.
        end
        x_custom_pp = copy( x_custom )
        x_custom_pp[ abs.(x_custom_pp) .< delta ] .= 0;
        sharpe_PPGD = sharpe_ratio( "./ftec_project_files/", x_custom_pp, ones(n), 2 )
        open(outfile, "w") do f
            for da in sharpe_PPGD
              println(f, da)
            end
        end
    end
end