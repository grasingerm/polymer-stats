 _ALGO_HELP_MSG_ = """algorithm (ISRES | PRAXIS | SBPLX | SLSQP | LBFGS | AUGLAG | Newton | Anneal)
  
  Separate algorithms by commas to run them in succession; e.g. LBGFS,Newton

  Algorithms:
  ISRES   -   Improved stochastic ranking evolution strategy; 
              nonlinear rained global optimization;
              note: requires arbitrary bounds on unknowns
  PRAXIS  -   PRincipal AXIS; gradient-free local optimization
  SBPLX   -   Subplex; variant of Nelder-Mead that is more efficient and robust;
              gradient-free local optimization
  SLSQP   -   Sequential quadratic programming for nonlineraly rained,
              gradient-based optimizationa
  LBFGS   -   low-storage BFGS; gradient-based unrained optimization
  AUGLAG  -   Augmented Lagrangian; rained optimization with LBFGS
  Newton  -   My implementation of the Newton-Raphson method
  Anneal  -   My implementation of simulated annealing

  Special syntax: try algorithm1,algorithm2 then backup_algorithms
  This first tries an algorithm that may be less robust but faster and
  if it fails then an algorithm that is more robust and slower.

  Note: in general, I've found SBPLX to be incredibly robust (for this problem).
  If you have no clue at an initial guess, my advice is to start with some
  simulated anneeling, hit with the SBPLX, and then finish it off with a little
  of Newton's method. I've tested this with intentionally pathological initial
  guesses and often converge to the solution anyway :)
""";
