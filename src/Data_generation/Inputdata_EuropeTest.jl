using GlobalEnergyGIS

function create_InputData()
# set scenario and target year

# define regions and countries and name the regionset
Dataset_Nordics = [
    "Norway"           GADM("Norway")  
    "Sweden"           GADM("Sweden")
    "Finland"          GADM("Finland")
    "Denmark"          GADM("Denmark") 
]

Dataset_EU27 = [
    "Austria"                   GADM("Austria") 
    "Belgium"                   GADM("Belgium") 
    "Bulgaria"                  GADM("Bulgaria") 
    "Croatia"                   GADM("Croatia") 
    "Cyprus"                    GADM("Cyprus") 
    "Czech_Republic"            GADM("Czech Republic") 
    "Denmark"                   GADM("Denmark") 
    "Estonia"                   GADM("Estonia") 
    "Finland"                   GADM("Finland") 
    "France"                    GADM("France") 
    "Germany"                   GADM("Germany") 
    "Greece"                    GADM("Greece") 
    "Hungary"                   GADM("Hungary") 
    "Ireland"                   GADM("Ireland") 
    "Italy"                     GADM("Italy") 
    "Latvia"                    GADM("Latvia") 
    "Lithuania"                 GADM("Lithuania") 
    "Luxembourg"                GADM("Luxembourg") 
    "Malta"                     GADM("Malta") 
    "Netherlands"               GADM("Netherlands") 
    "Poland"                    GADM("Poland") 
    "Portugal"                  GADM("Portugal") 
    "Romania"                   GADM("Romania") 
    "Slovakia"                  GADM("Slovakia") 
    "Slovenia"                  GADM("Slovenia") 
    "Spain"                     GADM("Spain") 
    "Sweden"                    GADM("Sweden") 
]

Dataset_EU27_CH_NO_UK = [
    "Austria"                   GADM("Austria") 
    "Belgium"                   GADM("Belgium") 
    "Bulgaria"                  GADM("Bulgaria") 
    "Croatia"                   GADM("Croatia") 
    "Cyprus"                    GADM("Cyprus") 
    "Czech_Republic"            GADM("Czech Republic") 
    "Denmark"                   GADM("Denmark") 
    "Estonia"                   GADM("Estonia") 
    "Finland"                   GADM("Finland") 
    "France"                    GADM("France") 
    "Germany"                   GADM("Germany") 
    "Greece"                    GADM("Greece") 
    "Hungary"                   GADM("Hungary") 
    "Ireland"                   GADM("Ireland") 
    "Italy"                     GADM("Italy") 
    "Latvia"                    GADM("Latvia") 
    "Lithuania"                 GADM("Lithuania") 
    "Luxembourg"                GADM("Luxembourg") 
    "Malta"                     GADM("Malta") 
    "Netherlands"               GADM("Netherlands") 
    "Norway"                    GADM("Norway")
    "Poland"                    GADM("Poland") 
    "Portugal"                  GADM("Portugal") 
    "Romania"                   GADM("Romania") 
    "Slovakia"                  GADM("Slovakia") 
    "Slovenia"                  GADM("Slovenia") 
    "Spain"                     GADM("Spain") 
    "Sweden"                    GADM("Sweden") 
    "Switzerland"               GADM("Switzerland")
    "United_Kingdom"            GADM("United Kingdom")
]

saveregions("EuropeTest", Dataset_EU27_CH_NO_UK)
makedistances("EuropeTest")
createmaps("EuropeTest")

# generate VRE data and demand
GISsolar(gisregion = "EuropeTest")
GISwind(gisregion = "EuropeTest")
GIShydro(gisregion = "EuropeTest")
predictdemand(gisregion = "EuropeTest", sspscenario="ssp2-26", sspyear=2050, era_year=2018)

end


# GIS options

solaroptions() = Dict(
    :gisregion => "EuropeTest",            # "Europe8", "Eurasia38", "Scand3"
    :filenamesuffix => "",              # e.g. "_landx2" to save high land availability data as "GISdata_solar2018_Europe8_landx2.mat" 

    :pv_density => 45,                  # Solar PV land use 45 Wp/m2 = 45 MWp/km2 (includes PV efficiency & module spacing, add latitude dependency later)
    :csp_density => 35,                 # CSP land use 35 W/m2

    :pvroof_area => .05,                # area available for rooftop PV after the masks have been applied
    :plant_area => .05,                 # area available for PV or CSP plants after the masks have been applied

    :distance_elec_access => 300,       # max distance to grid [km] (for solar classes of category B)
    :plant_persons_per_km2 => 150,      # not too crowded, max X persons/km2 (both PV and CSP plants)
    :pvroof_persons_per_km2 => 200,     # only in populated areas, so AT LEAST x persons/km2
                                        # US census bureau requires 1000 ppl/mile^2 = 386 ppl/km2 for "urban" (half in Australia)
                                        # roughly half the people of the world live at density > 300 ppl/km2
    :exclude_landtypes => [0,1,2,3,4,5,8,12],       # exclude water, forests and croplands. See codes in table below.
    :protected_codes => [1,2,3,4,5,8],  # IUCN codes to be excluded as protected areas. See codes in table below.

    :scenarioyear => "ssp2_2050",       # default scenario and year for population and grid access datasets
    :era_year => 2018,                  # which year of the ERA5 time series to use 

    :res => 0.01,                       # resolution of auxiliary datasets [degrees per pixel]
    :erares => 0.28125,                 # resolution of ERA5 datasets [degrees per pixel]

    :pvclasses_min => [0.08,0.14,0.18,0.22,0.26],   # lower bound on annual PV capacity factor for class X    [0:0.01:0.49;]
    :pvclasses_max => [0.14,0.18,0.22,0.26,1.00],   # upper bound on annual PV capacity factor for class X    [0.01:0.01:0.50;]
    :cspclasses_min => [0.10,0.18,0.24,0.28,0.32],  # lower bound on annual CSP capacity factor for class X
    :cspclasses_max => [0.18,0.24,0.28,0.32,1.00]  # upper bound on annual CSP capacity factor for class X
)

windoptions() = Dict(
    :gisregion => "EuropeSmall",            # "Europe8", "Eurasia38", "Scand3"
    :filenamesuffix => "",              # e.g. "_landx2" to save high land availability data as "GISdata_solar2018_Europe8_landx2.mat" 

    :onshore_density => 5,              # about 30% of existing farms have at least 5 W/m2, will become more common
    :offshore_density => 8,             # varies a lot in existing parks (4-18 W/m2)
                                        # For reference: 10D x 5D spacing of 3 MW turbines (with 1D = 100m) is approximately 6 MW/km2 = 6 W/m2
    :area_onshore => .08,               # area available for onshore wind power after the masks have been applied
    :area_offshore => .33,              # area available for offshore wind power after the masks have been applied

    :distance_elec_access => 300,       # max distance to grid [km] (for wind classes of category B and offshore)
    :persons_per_km2 => 150,            # not too crowded, max X persons/km2
                                        # US census bureau requires 1000 ppl/mile^2 = 386 ppl/km2 for "urban" (half in Australia)
                                        # roughly half the people of the world live at density > 300 ppl/km2
    :max_depth => 40,                   # max depth for offshore wind [m]
    :min_shore_distance => 5,           # minimum distance to shore for offshore wind [km]
    :exclude_landtypes => [0,11,13],    # exclude water, wetlands and urban areas. See codes in table below.
    :protected_codes => [1,2,3,4,5,8],  # IUCN codes to be excluded as protected areas. See codes in table below.

    :scenarioyear => "ssp2_2050",       # default scenario and year for population and grid access datasets
    :era_year => 2018,                  # which year of the ERA5 time series to use 
    :rescale_to_wind_atlas => true,     # rescale the ERA5 time series to fit annual wind speed averages from the Global Wind Atlas

    :res => 0.01,                       # resolution of auxiliary datasets [degrees per pixel]
    :erares => 0.28125,                 # resolution of ERA5 datasets [degrees per pixel]

    :onshoreclasses_min => [2,5,6,7,8],     # lower bound on annual onshore wind speeds for class X    [0:0.25:12.25;]
    :onshoreclasses_max => [5,6,7,8,99],    # upper bound on annual onshore wind speeds for class X    [0.25:0.25:12.5;]
    :offshoreclasses_min => [3,6,7,8,9],    # lower bound on annual offshore wind speeds for class X
    :offshoreclasses_max => [6,7,8,9,99]    # upper bound on annual offshore wind speeds for class X
)

hydrooptions() = Dict(
    :gisregion => "Europe8",                    # "Europe8", "Eurasia38", "Scand3"

    :costclasses_min => [ 0,  50, 100],         # US $/MWh
    :costclasses_max => [50, 100, 999],

    :storageclasses_min => [   0, 1e-6,  12],   # weeks (discharge time)
    :storageclasses_max => [1e-6,   12, 9e9]
)

