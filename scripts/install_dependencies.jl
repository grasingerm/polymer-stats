import Pkg;

@show depends = ([
                  "ArgParse";
                  "Distributions";
                  "Logging";
                  "DecFP";
                  "Quadmath";
                  "NLsolve";
                  "Plots";
                  "Optim";
                  "NLopt";
                  "Cubature";
                 ]);

foreach(Pkg.add, depends);

Pkg.update();
