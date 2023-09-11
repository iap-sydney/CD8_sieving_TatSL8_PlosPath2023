%Steffen Docken
%25-11-21
%Function to calculate Log-Likelihood of restricted Poisson model based on 
%parameters of model and model definition vector

function logL = loglikelihood_Poisson_restricted(params, model_def, ...
    react_det, preART_prop, avg_pred)

if (model_def(1) == 1) %if r_j != 0
    r_j = params(1);
    params = params(2:end);
else
    r_j = 0;
end

if (model_def(2) == 1) %if m_j != 0
    m_j = params(1);
    params = params(2:end);
else
    m_j = 0;
end

if (model_def(3) == 1) %if b_j != 0
    b_j = params(1);
else
    b_j = 0;
end

logL = loglikelihood_Poisson_full(react_det, preART_prop, r_j, m_j, b_j,...
    avg_pred);