using DataStructures

tab = [3.5, 7.2, 1.8, 9.0, 2.1]

# Create min-heap of (value, index) tuples
h = BinaryMinHeap{Tuple{Float64, Int}}()
for (i, val) in enumerate(tab)
    push!(h, (val, i))
end

# Use `top()` instead of `peek()` to avoid method clash
min_val, min_index = top(h)

println("Minimum value: $min_val at index: $min_index")



# Implement a heap by myself
# To allow for insertion, replacement, etc