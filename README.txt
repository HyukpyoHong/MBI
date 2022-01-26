See the section "Computational package of MBI" in Supplementary Text G in the paper for the detailed instruction to use the code.

- List of R code functions
1. 'Queueing_functions.R': the set of all source functions needed to run the MBI main functions below. This is automatically called by the line 'source("Queueing_functions.R")' when the MBI main functions are executed.
The next three functions estimate parameters under the assumption that observed data is generated from the delayed birth-death process.
2. 'MBI_single_cell_main.R': the main function to estimate parameters from multiple time traces obtained from a single cell. In other words, all the time traces are assumed to have the same parameters.
3. 'MBI_multiple_cells_main.R': the main function to estimate parameters from multiple time traces obtained from multiple cells (i.e., each time trace is measured from each of single cells in a heterogeneous population).
4. 'MBI_multiple_cells_same_n_main.R': the main function to estimate parameters from multiple time traces obtained from multiple cells. In this function, the number of rate-limiting steps (i.e., the shape parameter of the gamma delay distribution) are assumed to be the same among the single cells.
The next three functions estimate parameters under the assumption that observed data is generated from the delayed birth-death process including negative feedback regulation with a fixed delay.
5. 'MBI_feedback_single_cell_main.R': the main function to estimate parameters from multiple time traces obtained from a single cell. In other words, all the time traces are assumed to have the same parameters.
6. 'MBI_feedback_multiple_cells_main.R': the main function to estimate parameters from multiple time traces obtained from multiple cells (i.e., each time trace is measured from each of single cells in a heterogeneous population).
7. 'MBI_feedback_multiple_cells_same_n_main.R': the main function to estimate parameters from multiple time traces obtained from multiple cells. In this function, the number of rate-limiting steps (i.e., the shape parameter of the gamma delay distribution) are assumed to be the same among the single cells.

- Input data example
1. 'input_data_example.csv': the example file for input data from the delayed birth-death process. The first column contains the time points at which the data were measured. From the second to the last column, the measured time traces are contained.
2. 'input_data_feedback_example.csv': the example file for input data from the delayed birth-death process including negative feedback regulation with a fixed delay. The first column contains the time points at which the data were measured. From the second to the last column, the measured time traces are contained.

- Output data example
1.  'post_samples_MBI_single_cell_example.csv': the example file for output data of the main function 'MBI_single_cell_main.R'. The four columns contain posterior samples of lambda_b, lambda_d, n, and t_r, respectively. 
2.  'mean_std_post_samples_MBI_single_cell_example.csv': the example file for output data of the main function 'MBI_single_cell_main.R'. The four columns contain posterior means and standard deviations of lambda_b, lambda_d, n, and t_r, respectively.
