# New SCP Method
This is a fast trajectory planning method via a new sequential convex programming based on first-order method framework. Each SCP subproblem is handled by applying an augmented Lagrangian approach to incorporate inequality constraints into the cost function. The resulting equality-constrained problem is then solved via an iterative scheme that avoids costly matrix inversions.

Some suggestions for selecting the values of important hyperparameters

rho : 0.5 - 0.8               
eta = 1.6 - 2            
selfsigma = 1e-3            
rho_update_inteval : 50 / 25            
adaptive_rho_tolerance : 5            
omega = 0.3
