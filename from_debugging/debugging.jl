
using PProf, Profile 
# include("debugging_main.jl")
include("main_gyre2.jl")

main()
Profile.@profile begin 
debug_main()
end

PProf.pprof()