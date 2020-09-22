include(joinpath(@__DIR__, "..", "eap_chain.jl"));

const stepout = 500;
outpath(pargs) = "$(pargs["prefix"])_trajectory.csv";

[
  (chain::EAPChain, step::Int, startup::Bool, cleanup::Bool, pargs) -> begin;
    global outfile;
    if startup
      outfile = open(outpath(pargs), "a");
    end
  end,

  (chain::EAPChain, step::Int, startup::Bool, cleanup::Bool, pargs) -> begin;
    global outfile;
    if cleanup
      close(outfile);
    end
  end,

  (chain::EAPChain, step::Int, startup::Bool, cleanup::Bool, pargs) -> begin;
    global outfile;
    if step % stepout == 0
      writedlm(outfile, hcat(step, chain.r, chain_μ(chain), chain.U), ',');
    end
  end

]
