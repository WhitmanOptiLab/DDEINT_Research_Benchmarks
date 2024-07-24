println("Running DDE tests")

include("BC_Model.jl")
include("CV_Model.jl")
include("RS_Model.jl")

using DataFrames, CSV

bc_et = BC_Model.Run_BC_Model()
cv_et = CV_Model.Run_CV_Model()
rs_et = RS_Model.Run_Repressilator_Model()

println("DDE tests completed successfully")

data = DataFrame(model=["BC", "CV", "RS"], execution_time=[bc_et, cv_et, rs_et])
CSV.write("dde_elapsed_times_jl.csv", data)
