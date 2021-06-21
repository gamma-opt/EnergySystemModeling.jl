using GlobalEnergyGIS

create_scenario_datasets("SSP2", 2050)

EuropeSmall = [
    "Mediterranian"     GADM("Spain", "Portugal", "Italy", "Croatia", "Greece", "Bosnia and Herzegovina", "Slowenia", "Serbia", "Kosovo", "Albania", "Macedonia")
    "Nordics"           GADM("Norway", "Sweden", "Finland", "Denmark")
    "Western"           GADM("France", "Belgium", "Netherland", "Luxemburg")
    "Central"           GADM("Germany", "Austria", "Switzerland")    
    "Eastern"           GADM("Poland", "Czechia", "Slovakia", "Hungary")
]

saveregions("EuropeSmall", EuropeSmall)

makedistances("EuropeSmall")

createmaps("EuropeSmall")

