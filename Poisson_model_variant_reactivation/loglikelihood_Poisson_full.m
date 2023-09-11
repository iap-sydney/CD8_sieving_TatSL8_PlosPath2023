%Steffen Docken
%8-11-21
%Code to calculate log-likelihood of reactivation data for full Poisson
%model. Based on code from Yuhuang Wu.

% probability of reactivation is modeled as 1 - exp(-r_j - m_j*f_k - b_j*q_k,j))

%r_j is constant term for probability of reactivation

%m_j is coefficient of average size of variant on day 14 in all animals

%b_j is coefficient of predictor within animal

%preART_prop is vector of proportion in VL before ART)

%avg_pred is vector of average size of variant on day 14 in all animals

%react_det is the vector that contains whether a variant was detected in
%reactivation or not.

function logL = loglikelihood_Poisson_full(react_det, preART_prop, ...
    r_j, m_j, b_j, avg_pred)

num_SL8_var = length(react_det); %number of SL8 variants to consider

% log likelihood function
logL = 0;
for ii = 1:num_SL8_var
    p_reactivate = react_Poisson_Model_full_211118(preART_prop(ii),...
        r_j, m_j, b_j, avg_pred(ii));
    loglike_var = react_det(ii)*log(p_reactivate) +...
        (1-react_det(ii))*log(1-p_reactivate);
    logL = logL + loglike_var;
end