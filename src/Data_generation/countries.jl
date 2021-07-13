using MAT, DelimitedFiles, JLD, JSON, CSV, LinearAlgebra, Combinatorics, GlobalEnergyGIS, DataFrames

Europe = [
    # "Country"                   "GADM"    
      "Albania"                   GADM("Albania") 
      "Armenia"                   GADM("Armenia") 
      "Austria"                   GADM("Austria") 
      "Azerbaijan"                GADM("Azerbaijan") 
      "Belarus"                   GADM("Belarus")                
      "Belgium"                   GADM("Belgium") 
      "Bosnia and Herzegovina"    GADM("Bosnia and Herzegovina")            
      "Bulgaria"                  GADM("Bulgaria")           
      "Croatia"                   GADM("Croatia")        
      "Cyprus"                    GADM("Cyprus")             
      "Czech Republic"            GADM("Czech Republic")     
      "Denmark"                   GADM("Denmark")          
      "Estonia"                   GADM("Estonia")           
      "Finland"                   GADM("Finland")         
      "France"                    GADM("France")              
      "Germany"                   GADM("Germany")         
      "Greece"                    GADM("Greece")           
      "Hungary"                   GADM("Hungary")     
      "Iceland"                   GADM("Iceland")             
      "Ireland"                   GADM("Ireland")           
      "Italy"                     GADM("Italy")  
      "Kosovo"                    GADM("Kosovo")                        
      "Latvia"                    GADM("Latvia")           
      "Lithuania"                 GADM("Lithuania")        
      "Luxembourg"                GADM("Luxembourg")        
      "Malta"                     GADM("Malta")     
      "Moldova"                   GADM("Moldova") 
      "Montenegro"                GADM("Montenegro")                              
      "Netherlands"               GADM("Netherlands")   
      "North_Macedonia"           GADM("Macedonia")    
      "Norway"                    GADM("Norway")          
      "Poland"                    GADM("Poland")             
      "Portugal"                  GADM("Portugal")           
      "Romania"                   GADM("Romania") 
      "Russia"                    GADM("Russia")
      "Serbia"                    GADM("Serbia")              
      "Slovakia"                  GADM("Slovakia")           
      "Slovenia"                  GADM("Slovenia")          
      "Spain"                     GADM("Spain")              
      "Sweden"                    GADM("Sweden")            
      "Switzerland"               GADM("Switzerland")  
      "Turkey"                    GADM("Turkey")
      "Ukraine"                   GADM("Ukraine")         
      "United Kingdom"            GADM("United Kingdom")     
  ]  

Regions_dict = Dict( "Nordics" => (["Finland", "Sweden", "Norway", "Denmark"]),
  "Central" => (["Germany", "Austria",  "Switzerland", "Czech Republic"]),
  "Western" => (["France", "United Kingdom", "Ireland", "Netherlands", "Belgium", "Luxembourg"]),
  "Mediterranian" => (["Spain", "Portugal", "Italy", "Greece", "Croatia", "Malta", "Albania", "Bosnia and Herzegovina"]),
  "Eastern" => (["Poland", "Slowakia", "Hungary", "Lithuania", "Latvia", "Estonia"]))

# function if one wants specific countries
# Check if the regions specified in run_data_generation.jl are in Europe or in Regions_dict
# Get corresponding values
function get_countries(Regions)
    if sum(occursin.(Regions[1,1], Europe[:,1]))>0
        n = length(Europe[:,1])
        m = length(Regions)
        Dataset_Countries = []
        GADM_Name = []

        for j in 1:m
            for i in 1:n
                    if occursin.(Regions[j,1], Europe[i,1])
                        Dataset_Countries = [Dataset_Countries; permutedims(Regions[j,:])]
                        GADM_Name = [GADM_Name; GADM(Europe[i])]
                    end
            end
        end
        Countries = hcat(Dataset_Countries, GADM_Name)

    elseif sum(occursin.(Regions[1,1], collect(keys(Regions_dict))))>0
        Keys = []
        Values = []
        m = length(Regions)
        R = getindex.(Ref(Regions_dict),(Regions))

        for i in 1:m
            RVal = getindex.(Ref(Regions_dict),(Regions))[i]
            Keys = [Keys; Regions[i]]
            Values = [Values; GADM(RVal)]
        end
        GADM_List = hcat(Keys, Values)
    end
end



# function if one wants to create whole regions


# if you want all countries in Europe
function get_countries()
   Europe
end



