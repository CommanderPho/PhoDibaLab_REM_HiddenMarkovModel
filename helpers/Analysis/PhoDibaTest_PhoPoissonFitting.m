

temp.curr_data = across_experiment_results{1, 1}.active_processing.spikes.time{1};  


[lambdahat,lambdaci] = poissfit(temp.curr_data);
