# This function is a simple yet reasonably robust implementation of a
# primal-dual interior-point solver for convex programs with convex
# inequality constraints (it does not handle equality constraints).
# It will compute the solution to the following optimization problem:
#
#     minimize    f(x)
#     subject to  c(x) < 0,
#
# where f(x) is a convex objective and c(x) is a vector-valued
# function with outputs that are convex in x. There are many
# optimization problems that can be framed in this form (e.g., see the
# book "Convex Optimization" by Boyd & Vandenberghe). The code is
# mostly based on:
#
#    Armand, Gilbert & Jan-Jegou. A feasible BFGS interior point
#    algorithm for solving convex minimization problems. SIAM Journal
#    on Optimization 11, 199-222.
#
# However, to understand what is going on you will need to read up on
# interior-point methods for constrained optimization. A good starting
# point is the "Convex Optmization" book.
#
# Input argument x0 is the initial point for the solver. Input "tol"
# is the tolerance of the convergence criterion; it determines when
# the solver should stop. Input "maxiter" is the maximum number of
# iterations.
#
# The inputs "obj", "grad", "constr" and "jac" are functions:
# 
#   obj takes input x, the vector of optimization variables, and
#   returns the value of the objective function f(x) at x.
#
#   grad takes the same input x and returns a list with
#   two elements: the gradient and n x n Hessian of the objective at x.
#
#   constr returns the value of the constraint function c(x) at x.
#
#   jac takes two inputs: the primal varaibles x and the dual
#   variables z. The return value is a list with two elements: the m x
#   n Jacobian matrix (containing the first-order partial derivatives
#   of the inequality constraint functions), and W is the n x n
#   Hessian of the Lagrangian (minus the Hessian of the objective),
#   equal to
#
#     W = z[1]*W1 + z[2]*W2 + ... + z[m]*Wm,
#
#   where Wi is the Hessian of the ith constraint.
#
# If you set "verbose" to true, then at each iteration the solver will
# output the following information (from left to right): 1. the
# iteration number, 2. the value of the objective, 3. the barrier
# parameter mu, 4. the centering parameter sigma, 4. the residuals of
# the perturbed Karush-Kuhn-Tucker system (rx, rc), 5. the step size,
# and the number of iterations in the line search before we found a
# suitable descent step.
#
# Note that the interior-point solver may not work very well if your
# problem is very poorly scaled (i.e. the Hessian of the objective or
# the Hessian of one of the constraint functions is poorly
# conditioned). It is up to you to make sure you look at the
# conditioning of your problem.
ipsolver <- function (x, obj, grad, constr, jac, tol, maxiter, verbose) {

  # Some algorithm parameters.
  eps       <- 1e-8   # A number close to zero.
  sigmamax  <- 0.5    # Maximum centering parameter.
  etamax    <- 0.25   # Maximum forcing number.
  mumin     <- 1e-9   # Minimum barrier parameter.
  alphamax  <- 0.995  # Maximum step size.
  alphamin  <- 1e-6   # Minimum step size.
  beta      <- 0.75   # Granularity of backtracking search.
  tau       <- 0.01   # Decrease we will accept in line search.

  # INITIALIZATION
  # --------------
  # Get the number of primal variables (n), the number of constraints (m),
  # and the total number of primal-dual optimization variables (nv).
  # Initialize the Lagrange multipliers. Also, initialize the second-order
  # information.
  c  <- constraints(x)
  n  <- length(x)
  m  <- length(c)
  nv <- n + m
  z  <- rep(1,m)
  B  <- diag(rep(1,n))
  
  if verbose
    fprintf('  i f(x)       lg(mu) sigma   ||rx||  ||rc||  alpha   #ls\n');
  end
  
  % Repeat while the convergence criterion has not been satisfied, and
  % we haven't reached the maximum number of iterations.
  alpha = 0;
  ls    = 0;
  for iter = 1:maxiter

    % COMPUTE OBJECTIVE, GRADIENT, CONSTRAINTS, ETC.  
    % Compute the response of the objective function, the gradient of the
    % objective, the response of the inequality constraints, the Jacobian of
    % the inequality constraints, the Hessian of the Lagrangian (minus the
    % Hessian of the objective) and, optionally, the Hessian of the
    % objective.
    f     = objective(x);
    c     = constraints(x);
    [J W] = jacobian(x,z);
    if strcmp(descentdir,'newton')
      [g B] = gradient(x);
    else
      g = gradient(x);
    end
    
    % Compute the responses of the unperturbed Karush-Kuhn-Tucker
    % optimality conditions.
    rx = g + J'*z;  % Dual residual.
    rc = c.*z;      % Complementarity.
    r0 = [rx; rc]; 
    
    % Set some parameters that affect convergence of the primal-dual
    % interior-point method.
    eta        = min(etamax,norm(r0)/nv);
    sigma      = min(sigmamax,sqrt(norm(r0)/nv));
    dualitygap = -c'*z;
    mu         = max(mumin,sigma*dualitygap/m);
    
    % Print the status of the algorithm.
    if verbose
      fprintf('%3d %+0.3e  %+5.2f %0.1e %0.1e %0.1e %0.1e %3d\n',...
	      iter,f,log10(mu),sigma,norm(rx),norm(rc),alpha,ls);
    end

    % CONVERGENCE CHECK.
    % If the norm of the responses is less than the specified tolerance,
    % we are done. 
    if norm(r0)/nv < tolerance
      break
    end
    
    % SOLUTION TO PERTURBED KKT SYSTEM.
    % Compute the search direction of x and z.
    S  = diag(sparse(z./(c-eps)));
    gb = g - mu*J'*(1./(c-eps));
    px = (B + W - J'*S*J) \ (-gb);
    pz = -(z + mu./(c-eps) + S*J*px);
    
    % BACKTRACKING LINE SEARCH.
    % To ensure global convergence, execute backtracking line search to
    % determine the step length. First, we have to find the largest step
    % size which ensures that z remains feasible. Next, we perform
    % backtracking line search.
    alpha = alphamax;
    is    = find(z + pz < 0);
    if length(is)
      alpha = alphamax * min(1,min(z(is) ./ -pz(is)));
    end
    
    % Compute the response of the merit function and the directional
    % gradient at the current point and search direction.
    psi  = merit(x,z,f,c,mu,eps);
    dpsi = gradmerit(x,z,px,pz,g,c,J,mu,eps);
    ls   = 0;
    while true

      % Compute the candidate point, the constraints, and the response of
      % the objective function and merit function at the candidate point.
      ls     = ls + 1;
      xnew   = x + alpha * px;
      znew   = z + alpha * pz;
      f      = objective(xnew);
      c      = constraints(xnew);
      psinew = merit(xnew,znew,f,c,mu,eps);
      
      % Stop backtracking search if we've found a candidate point that
      % sufficiently decreases the merit function and satisfies all the
      % constraints.
      if sum(c > 0) == 0 & psinew < psi + tau*eta*alpha*dpsi
	x     = xnew;
	z     = znew;
	gprev = g;
	break
      end
      
      % The candidate point does not meet our criteria, so decrease the step
      % size for 0 < beta < 1.
      alpha = alpha * beta;
      if alpha < alphamin
	error('Step size too small');
      end
    end
  end
  
% ------------------------------------------------------------------
% Compute the response of the merit function at (x,z).
function psi = merit (x, z, f, c, mu, eps)
  psi = f - c'*z - mu*sum(log(c.^2.*z+eps));
  
% ------------------------------------------------------------------
% Compute the directional derivative of the merit function at (x,z).
function dpsi = gradmerit (x, z, px, pz, g, c, J, mu, eps)
  dpsi = px'*(g - J'*z - 2*mu*J'*(1./(c-eps))) - pz'*(c + mu./(z+eps));

% ------------------------------------------------------------------
% Update the quasi-Newton approximation using Broyden-Fletcher-
% Goldfarb-Shanno (BFGS) formula.
function B = bfgsupdate (B, s, y)  
  if y'*s < 0
    error('dot(y,s) > 0 is not satisfied');
  end
  x = B*s;
  B = B - x*x'/(x'*s) + y*y'/(y'*s);
  
