using Printf

function find_non_zeros(A)

    for j = 1:size(A)[1]
        for k = 1:size(A)[2]

            if A[j,k] != 0.0
                @printf "(%i, %i) = %e" j k A[j,k]
            end

        end
    end

end