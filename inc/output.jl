include("free_energy.jl");

function standard_solution_output(xs::Vector)
  println();
  println("A = $(free_energy(xs))");
  println("lambda = $(xs[2])");
  println("alpha = $(xs[3])");
  println("C = $(xs[1])");
  println("|lambda| = $(sqrt(xs[2]^2 + xs[3]^2))");
  μ = polarization(xs);
  println("mu1 = $(μ[1])");
  println("mu2 = $(μ[2])");
  println("mu3 = $(μ[3])");
end

function residuals_output(rs::Real)
  println("L2 = $rs");
end

function residuals_output(rs::Vector)
  println("L2 = $(norm(rs, 2))");
end
