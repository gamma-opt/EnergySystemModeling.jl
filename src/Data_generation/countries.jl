using MAT, DelimitedFiles, JLD, JSON, CSV, LinearAlgebra, Combinatorics, GlobalEnergyGIS, DataFrames

Dataset_EU27_CH_NO_UK = [
  # "Country"                   "GADM"                      "RoR_data"       "PHS_data"
    "Austria"                   GADM("Austria")             8120             5231
    "Belgium"                   GADM("Belgium")             112              1310
    "Bulgaria"                  GADM("Bulgaria")            2206             1013
    "Croatia"                   GADM("Croatia")             1915             293
    "Cyprus"                    GADM("Cyprus")              0                0
    "Czech_Republic"            GADM("Czech Republic")      1088             1172
    "Denmark"                   GADM("Denmark")             7                0
    "Estonia"                   GADM("Estonia")             6                0
    "Finland"                   GADM("Finland")             3249             0
    "France"                    GADM("France")              18163            7115
    "Germany"                   GADM("Germany")             4577             6822
    "Greece"                    GADM("Greece")              2693             699
    "Hungary"                   GADM("Hungary")             57               0
    "Ireland"                   GADM("Ireland")             237              292
    "Italy"                     GADM("Italy")               14628            7592
    "Latvia"                    GADM("Latvia")              1589             0
    "Lithuania"                 GADM("Lithuania")           117              760
    "Luxembourg"                GADM("Luxembourg")          34               1296
    "Malta"                     GADM("Malta")               0                0
    "Netherlands"               GADM("Netherlands")         37               0
    "Norway"                    GADM("Norway")              29939            1397
    "Poland"                    GADM("Poland")              588              1782
    "Portugal"                  GADM("Portugal")            4379             1789
    "Romania"                   GADM("Romania")             6359             371
    "Slovakia"                  GADM("Slovakia")            1606             916
    "Slovenia"                  GADM("Slovenia")            1115             180
    "Spain"                     GADM("Spain")               14086            5967
    "Sweden"                    GADM("Sweden")              16230            99
    "Switzerland"               GADM("Switzerland")         11850            1839
    "United_Kingdom"            GADM("United Kingdom")      1759             2744
]
writedlm("Countries_EU_CH_NO.CSV", Dataset_EU27_CH_NO_UK)



function receive_countries(Country_names)

    Country_List = Dataset_EU27_CH_NO_UK[:,1:2]
    n = length(Country_List[:,1])
    m = length(Country_names)
    Dataset_Countries = []
    GADM_Name = []

    for j in 1:m
        for i in 1:n
                if occursin.(Country_names[j,1], Country_List[i,1])
                    Dataset_Countries = [Dataset_Countries; permutedims(Country_names[j,:])]
                    GADM_Name = [GADM_Name; Country_List[i,2]]
                end
        end
    end
    Countries = hcat(Dataset_Countries, GADM_Name)
end
