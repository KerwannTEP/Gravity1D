# Is there a mistake somewhere ?
# Test for low N


mutable struct MinHeap
    data::Vector{BigFloat} # Data
    index_HP::Vector{Int64} # Particle index at for each heap index
    index_PH::Vector{Int64} # Heap index at for each particle index
end

# Create an empty heap
function MinHeap()
    return MinHeap(BigFloat[], Int64[], Int64[])
end

# Helper functions

function parent(i::Int64)
    return div(i, 2)
end

function left(i::Int64)
    return 2*i
end

function right(i::Int64)
    return 2*i+1
end

# Swap two elements
function swap!(v::Vector, i::Int64, j::Int64)
    v[i], v[j] = v[j], v[i]
    return nothing
end

# Bubble up to restore heap after insert
function bubble_up!(heap::MinHeap, i::Int64)
    while((i>1) && (heap.data[i] < heap.data[parent(i)]))
        index_particle = heap.index_HP[i]
        index_particle_parent = heap.index_HP[parent(i)]
        # Swap data between parent and child
        swap!(heap.data, i, parent(i))
        # Swap particle index between parent and child
        swap!(heap.index_HP, i, parent(i))
        # Swap heap index between parent and child
        swap!(heap.index_PH, index_particle, index_particle_parent)

        i = parent(i)
    end
    return nothing
end

# Bubble down to restore heap after removal
function bubble_down!(heap::MinHeap, i::Int)
    n = length(heap.data)
    current = i
    while true
        l = left(current)
        r = right(current)
        smallest = current

        if ((l <= n) && (heap.data[l] < heap.data[smallest]))
            smallest = l
        end
        if ((r <= n) && (heap.data[r] < heap.data[smallest]))
            smallest = r
        end

        if (smallest == current) # Parent is smaller than children
            break
        end

        index_particle = heap.index_HP[current]
        index_particle_parent = heap.index_HP[smallest]
        # Swap data between parent and child
        swap!(heap.data, current, smallest)
        # Swap particle index between parent and child
        swap!(heap.index_HP, current, smallest)
        # Swap heap index between parent and child
        swap!(heap.index_PH, index_particle, index_particle_parent)

        current = smallest
    end
    return nothing
end

# Insert a new value
function push!(heap::MinHeap, val::BigFloat, index::Int64)
    Base.push!(heap.data, val)
    Base.push!(heap.index_HP, index)
    Base.push!(heap.index_PH, index)
    bubble_up!(heap, index)
    return nothing
end

# Peek at the minimum value
function top(heap::MinHeap)
    return (heap.data[1], heap.index_HP[1])
end

# Pop the minimum value
function pop!(heap::MinHeap)
    n = length(heap.data)
    if n == 0
        error("Heap is empty")
    elseif n == 1
        return (Base.pop!(heap.data), Base.pop!(heap.index_HP), Base.pop!(heap.index_PH))
    else
        min_val = heap.data[1]
        min_index_HP = heap.index_HP[1]
        min_index_PH = heap.index_PH[1]
        heap.data[1] = Base.pop!(heap.data)  # Move last to root
        heap.index_HP[1] = Base.pop!(heap.index_HP)  # Move last to root
        heap.index_PH[1] = Base.pop!(heap.index_PH)  # Move last to root
        bubble_down!(heap, 1)
        return min_val, min_index_HP, min_index_PH
    end
end

# Replace value at known heap index
function replace!(heap::MinHeap, i::Int, new_val::BigFloat)
    old_val = heap.data[i]
    heap.data[i] = new_val
    if new_val < old_val
        bubble_up!(heap, i)
       
    elseif new_val > old_val
        bubble_down!(heap, i)
    end
    return nothing
end

# Check if heap is empty
function isempty(heap::MinHeap)
    return Base.isempty(heap.data)
end

