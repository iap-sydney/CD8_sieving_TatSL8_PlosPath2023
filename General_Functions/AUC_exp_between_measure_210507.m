%Steffen Docken
%7-5-21
%Code to calculate the integral of viral load (or some other data) when 
%we assume the variable changes exponentially between time points

function AUC_exp = AUC_exp_between_measure_210507(x_data, y_data)

num_datapoints = length(x_data);

AUC_exp_vec = zeros(num_datapoints-1, 1); %will hold the integral for each
%interval between measurements

for ii = 1:num_datapoints-1
    x_old = x_data(ii);
    x_new = x_data(ii+1);
    y_old = y_data(ii);
    y_new = y_data(ii+1);
    
    if (y_new == y_old) %no change in y_data
        AUC_exp_vec(ii) = y_old*(x_new - x_old);
    else
        
        AUC_exp_vec(ii) = (y_new - y_old)*(x_new - x_old)/(log(y_new/y_old));
    end
    
end

AUC_exp = sum(AUC_exp_vec);