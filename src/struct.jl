mutable struct Stack{T}
    last :: Int
    data :: Array{T,1}
    function Stack{T}(N::Int) where T
        new(0,Array{T,1}(undef,N))
    end 
 end

function push!(s::Stack, sth)
    s.last+=1
    s.data[s.last] = sth
end

function pop!(s::Stack)
    s.last-=1
    s.data[s.last+1]
end

function isempty(s::Stack)
    s.last==0
end