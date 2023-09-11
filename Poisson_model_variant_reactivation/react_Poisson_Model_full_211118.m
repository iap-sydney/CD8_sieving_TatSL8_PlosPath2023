%Steffen Docken
%10-11-21
%Code to calculate probability of reactivation given model parameters and
%reactivation predictors (react_pred and avg_pred)

function p_react = react_Poisson_Model_full_211118(preART_prop,...
    r_j, m_j, b_j, avg_pred)

p_react = 1 - exp(-r_j - m_j*avg_pred - b_j*preART_prop);