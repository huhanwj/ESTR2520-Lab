using CSV, Glob, DataFrames
using Statistics

function sharpe_ratio( path_to_test, x, w, c )
    files = glob( "*_test.csv", path_to_test );
    dfs = DataFrame.( CSV.File.( files ) ); 
    T = 459; n = length(dfs);
    stocks_retur = zeros(T,n);
    for i = 1:n
        # compute the realized return R_i(t)
        stocks_retur[:,i] = (dfs[i].close-dfs[i].open) ./ dfs[i].open;
    end
    names_stocks = [ dfs[i].Name[1] for i in 1:n ];
    # calculate r_i and Sigma
    bar_R = [ mean( stocks_retur[:,i] ) for i in 1:length(dfs) ];
    Sigma = [ mean( (stocks_retur[:,i].-bar_R[i]).*(stocks_retur[:,j].-bar_R[j]) ) for i=1:n, j=1:n ]; 
    # calculate Sharpe Ratio
    portfolio = w + x;
    if sum( portfolio ) < 1e-10
        sharpe_ratio = 0;
    else
        sharpe_ratio = sum( bar_R.* portfolio ) / sqrt(portfolio'*Sigma*portfolio);
    end
    cost = sum( x + w );
    tx_cost = c*sum( abs.(x) .> 1e-5 );
    print("Sharpe Ratio = ", sharpe_ratio, ", Return = ", sum(bar_R.*portfolio), ", Tx Cost = ", tx_cost, ", Portfo Value = ", cost );
    k=[sharpe_ratio;sum(bar_R.*portfolio);tx_cost];
    return k
end