"""
Placing the checkpointed integraton, loop, and all of the packages here, so that when running 
include("technical_paper.jl") we never re-include them. This should help with Enzyme compile times

When running the experiments for the technical paper only ever include("technical_paper_integration.jl") once
"""

