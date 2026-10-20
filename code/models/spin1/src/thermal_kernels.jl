# Off-diagonal coefficients are +/-1 or +/-2. Exposing these constants avoids
# general Int64 multiplication in the inner SIMD loops, especially on ARM.
# Floor each output row while it is hot in cache. No arithmetic is approximated.
function multiply_real!(y, rows, x, shift)
    lanes = size(x, 1)
    @inbounds for i in eachindex(rows)
        row = rows[i]
        _, diagonal, _ = row[1]
        @simd for lane in 1:lanes
            y[lane, i] = diagonal * x[lane, i]
        end
        for entry in 2:length(row)
            j, a, _ = row[entry]
            if a == -2
                @simd for lane in 1:lanes
                    y[lane, i] -= x[lane, j] << 1
                end
            elseif a == 2
                @simd for lane in 1:lanes
                    y[lane, i] += x[lane, j] << 1
                end
            elseif a == -1
                @simd for lane in 1:lanes
                    y[lane, i] -= x[lane, j]
                end
            elseif a == 1
                @simd for lane in 1:lanes
                    y[lane, i] += x[lane, j]
                end
            else
                error("unexpected real off-diagonal coefficient")
            end
        end
        @simd for lane in 1:lanes
            y[lane, i] >>= shift
        end
    end
end

function multiply_eisenstein!(y, z, rows, u, v, shift)
    lanes = size(u, 1)
    @inbounds for i in eachindex(rows)
        row = rows[i]
        _, diagonal, _ = row[1]
        @simd for lane in 1:lanes
            y[lane, i] = diagonal * u[lane, i]
            z[lane, i] = diagonal * v[lane, i]
        end
        for entry in 2:length(row)
            j, a, b = row[entry]
            if (a, b) == (-1, 0)
                @simd for lane in 1:lanes
                    y[lane, i] -= u[lane, j]
                    z[lane, i] -= v[lane, j]
                end
            elseif (a, b) == (0, -1)
                @simd for lane in 1:lanes
                    y[lane, i] += v[lane, j]
                    z[lane, i] += v[lane, j] - u[lane, j]
                end
            elseif (a, b) == (1, 1)
                @simd for lane in 1:lanes
                    y[lane, i] += u[lane, j] - v[lane, j]
                    z[lane, i] += u[lane, j]
                end
            else
                error("unexpected Eisenstein off-diagonal coefficient")
            end
        end
        @simd for lane in 1:lanes
            y[lane, i] >>= shift
            z[lane, i] >>= shift
        end
    end
end
