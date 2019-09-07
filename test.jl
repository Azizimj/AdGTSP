#=
test:
- Julia version: 
- Author: MJazizi
- Date: 2019-09-02
=#



# function powerset(x::Vector{T}) where T
#     result = Vector{T}[[]]
#     for elem in x, j in eachindex(result)
#         push!(result, [result[j] ; elem])
#     end
#     result
# end
#
# a = powerset([1,2,3])
# for aa in a:
#     print(aa)

# import Pkg; Pkg.add("Combinatorics")
#
# using Combinatorics
# a = collect(powerset([1,2,3]))
# print(a)

# using Random
# print(rand(1)*10)

if size(ARGS)[1]>0
print(ARGS[1]=="sss")
end

# a =1
# print(string(a, 2)+" -")
# x =


