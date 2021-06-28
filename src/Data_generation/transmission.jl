using JLD: length
using MAT, DelimitedFiles, JLD, LinearAlgebra



# path to your set output folder of GlobalEnergyGIS data
inputdata = "D:\\Eigene Dateien\\Studium\\Master\\RA\\Copernicus\\output"

# Name of the regionset you defined in Inputdata.jl
regionset = "EuropeTest"

TransmissionData = matread(joinpath(inputdata, "distances_$regionset.mat"))

Distances = typeof(TransmissionData["distances"]) == Float64 ? [TransmissionData["distances"]] : TransmissionData["distances"]


# get every unique distance by taking the lower left triangle of the distance matrix
Dist = Distances[tril!(trues(size(Distances)), -1)]

n = 5
# get the combinations of the pair of nodes

Nodes = [1:n;]
Pairs = collect(combinations(Nodes, 2))

g = hcat(Pairs, Dist)

# Create Matrix for remaining parameters of transmission

TransData = [
   # :cost       :converter_cost         :M      :C      :B      :efficiency         :lifetime       :tcap_min       :tcap_max
    400         150000                  0.02    0       1       0.95                40              0               1000000000
]

# Fill the matrix according to the number of pairs in the system
TransDatafull = repeat(TransData; outer=[length(Pairs)])

# combine the tables, add header and create CSV file
TransmissionData = hcat(Pairs, Dist, TransDatafull)
Header = ["line" "dist"  "cost" "converter_cost" "M" "C"  "B" "efficiency" "lifetime" "tcap_min" "tcap_max"]
TransmissionTable = vcat(Header, TransmissionData)
writedlm("transmission.csv", TransmissionTable, ',')
