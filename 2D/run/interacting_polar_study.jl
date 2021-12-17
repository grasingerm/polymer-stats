println("Hello world!");

@everywhere using Printf;

@show workdir = if length(ARGS) > 0
  ARGS[1];
else
  ".";
end

@everywhere fmt(x) = @sprintf("%07d", round(Int, 1e3*x));
@everywhere fmt_int(x) = @sprintf("%03d", x);

@everywhere function prefix(case)
  "E0-$(fmt(case[:E0]))_mu-$(fmt(case[:mu]))_kT-$(fmt(case[:kT]))_Fz-$(fmt(case[:Fz]))_Fx-$(fmt(case[:Fx]))_n-$(fmt(case[:n]))_b-$(fmt(case[:b]))";
end

cases = Any[];
E0s = [0.0; 1e-1; 1.0; 10.0];
μs = [1e-6; 1e-3; 1e-2; 1e-1; 0.5; 1.0; 2.0; 10.0]; 
kTs = [1.0];
Fzs = sort(vcat(-[0.25; 0.5; 0.75; 1.0; 2.0; 3.0; 4.0; 5.0; 10.0; 20.0],
                [0.0; 0.25; 0.5; 0.75; 1.0; 2.0; 3.0; 4.0; 5.0; 10.0; 20.0]));
Fxs = [0.0; 0.25; 0.5; 0.75; 1.0; 2.0; 3.0; 4.0; 5.0; 10.0; 20.0];
ns = Int[100; 200];
bs = [0.5; 1.0; 2.0];
for b in bs, n in ns, Fx in Fxs, Fz in Fzs, kT in kTs, E0 in E0s, μ in μs
    push!(cases, Dict(:E0 => E0, :mu => μ,
                      :kT => kT, :Fz => Fz,
                      :Fx => Fx, :n => n, :b => b));
end

@info "total number of cases to run: $(length(cases))";

pmap(case -> begin;
  outfile = joinpath(workdir, "$(prefix(case)).out");
  if !isfile(outfile)
    println("Running case: $case.");
    command = `julia -O 3 mcmc_eap_chain.jl --chain-type polar --energy-type interacting -b $(case[:b]) --E0 $(case[:E0]) --mu $(case[:mu]) --kT $(case[:kT]) --Fz $(case[:Fz]) --Fx $(case[:Fx]) -n $(case[:n]) --num-steps 1000000 -v 2 --prefix $(joinpath(workdir, prefix(case))) --numeric-type big`;
    output = read(command, String);
    write(outfile, output); 
  else
    println("Case: $case has already been run.");
  end
end, cases);
