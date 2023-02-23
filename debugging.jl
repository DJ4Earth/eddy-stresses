
using PProf, Profile 
include("debugging_main.jl")

debug_main()
Profile.@profile begin 
debug_main()
end

PProf.pprof()