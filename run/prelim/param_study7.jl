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
  "E0-$(fmt(case[:E0]))_K1-$(fmt(case[:K1]))_K2-$(fmt(case[:K2]))_kT-$(fmt(case[:kT]))_Fz-$(fmt(case[:Fz]))_Fx-$(fmt(case[:Fx]))_n-$(fmt(case[:n]))_b-$(fmt(case[:b]))";
end

cases = Any[];
E0s = [1.0; 1.25; 1.5; 2.0; 2.5; 5.0; 7.5; 10.0];
Ks = zip([1.0; 0.0], [0.0; 1.0]);
kTs = [1.0];
Fzs = vcat(collect(range(0.0, 1.0; step=0.1)), range(2.0, 4.0; step=1.0), range(5.0, 20.0; step=2.5), 40.0);
Fxs = [0.0];
ns = Int[100];
bs = [1.0];
for b in bs, n in ns, Fx in Fxs, Fz in Fzs, kT in kTs, Kpair in Ks, 
    E0 in E0s
    push!(cases, Dict(:E0 => E0, :K1 => Kpair[1], :K2 => Kpair[2], 
                      :kT => kT, :Fz => Fz,
                      :Fx => Fx, :n => n, :b => b));
end

@info "total number of cases to run: $(length(cases))";

pmap(case -> begin;
  outfile = joinpath(workdir, "$(prefix(case)).out");
  if !isfile(outfile)
    println("Running case: $case.");
    command = `julia -O 3 mcmc_eap_chain.jl --chain-type dielectric --energy-type noninteracting --do-flips -b $(case[:b]) --E0 $(case[:E0]) --K1 $(case[:K1]) --K2 $(case[:K2]) --kT $(case[:kT]) --Fz $(case[:Fz]) --Fx $(case[:Fx]) -n $(case[:n]) --num-steps 1500000 -v 2 --prefix $(joinpath(workdir, prefix(case))) --umbrella-sampling --numeric-type big`;
    output = read(command, String);
    write(outfile, output); 
  else
    println("Case: $case has already been run.");
  end
end, cases);
