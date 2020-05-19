"""
    setResidualVariance(LRes::Array{Float64,2})

Set the residual variance used in the simulation.
Dimensions must match the number of QTL traits.
The matrix can be either the dense covariance matrix or
can be the transpose of the covariance matrix cholesky decomposition upper factor.
"""
function setResidualVariance(LRes::Float64)
    setResidualVariance(fill(LRes,1,1))
end
function setResidualVariance(LRes::Array{Float64,2})
    setResidualVariance(cholesky(LRes).U')
end
function setResidualVariance(chol_fact_Ud::Adjoint{Float64,UpperTriangular{Float64,Array{Float64,2}}})
    if get_num_traits(common.G) != size(chol_fact_Ud)[1]
        throw(DimensionMismatch("Number of Residual variance traits ($(size(chol_fact_Ud)[1])) != Number of qtl_traits ($(get_num_traits(common.G)))."))
    end
    common.LRes = chol_fact_Ud
end
