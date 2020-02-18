# EnergySystemModel.jl
Documentation for EnergySystemModel.jl


## Model
Mathematical reference for the energy system model. The model presented here is based on the model in [^1]. The units are expressed in the square brackets.

### Indices and Sets
*  $g∈G$: Generation technologies
*  $G^r⊆G$: Renewable generation technologies
*  $n∈N$: Nodes
*  $l∈L$: Transmission lines, bidimensional vectors $(i,j)$ where $i,j∈N$
*  $t∈T$: Time steps, depending on the number of clusters per month
*  $s∈S$: Storage technologies

### Parameters
Constant parameters

*  $κ∈[0,1]$: Renewables participation required by the system
*  $C$: Shedding cost [€/MWh]
*  $\bar{C}$: Shedding capacity [€/MWh]

Time clustered parameters

*  $τ$: Number of time slices per month
*  $τ_{t}$: Duration of time period $t$ [h]
*  $Q_{g,n}$: Initial capacity [MW]
*  $A_{g,n,t}∈[0,1]$: Availability of technology $g$ per node $n$ at time step $t$
*  $D_{n,t}$: Clustered demand per node $n$ per time step $t$  [MWh]

Generation technology parameters

*  $I_g^G$: Annualised investment cost for generation per MW of technology $g$ [€/MW]
*  $M_g^G$: Annualised maintenance cost for generation per MW of technology $g$ [€/MW]
*  $C_g^G$: Operational cost per MWh of technology $g$ [€/MWh]
*  $r_g^{-}$: Relative ramp-down limit of technology $g$
*  $r_g^{+}$: Relative ramp-up limit of technology $g$

Transmission parameters

*  $I_l^F$: Annualised investment cost for transmission per line $l$ [€/MW]
*  $M_l^F$: Annualised maintenance cost for transmission per line $l$ [€/MW]
*  $C_l^F$: Transmission cost per line $l$ [€/MWh]
*  $B_l$: Susceptance per line $l$

Storage parameters

*  $I_s^S$: Annualised investment cost of storage technology $s$ per MW [€/MWh]
*  $C_s^S$: Storage operational cost of storage technology $s$ [€/MW]
*  $b_{s,n}^0$: Initial capacity of storage $s$ at node $n$ [MWh]
*  $ξ_s$: Round-trip efficiency of storage technology $s$


### Variables
Generation technology variables

*  $p_{g,n,t}≥0$: Dispatch from technology $g$ at node $n$ in each time step $t$ [MWh]
*  $\bar{p}_{g,n}≥0$: Generation capacity invested in each technology $g$ at node $n$ [MW]

Shedding variables

*  $σ_{n,t}≥0$: Loss of load at node $n$ in each time step $t$ [MWh]

Transmission variables

*  $f_{l,t}$: Transmission flow per line $l$ in each time step $t$ [MWh]
*  $|f_{l,t}|$: Absolute value of transmission flow per line $l$ in each time step $t$ [MWh]
*  $\bar{f}_l$: Transmission capacity per line $l$ [MW]

Storage variables

*  $b_{s,n,t}≥0$: Storage level of storage $s$ at node $n$ in each time step $t$ [MWh]
*  $\bar{b}_{s,n}≥0$: Storage capacity of storage $s$ at node $n$ [MWh]
*  $b_{s,n,t}^{+}≥0$: Charging of storage $s$ at node $n$ in each time step $t$ [MWh]
*  $b_{s,n,t}^{-}≥0$: Discharging of storage $s$ at node $n$ in each time step $t$ [MWh]

Voltage angle variables

*  $θ_{n,t}≥0$: Voltage angle at node $n$ in each time step $t$
*  $θ'_{n,t}≥0$: Voltage angle at node $n$ in each time step $t$

### Objective
The objective is

$$\mathrm{minimize}_{p_{g,t}, \bar{p}_g, σ_{t}, f_{l,t}, \bar{f}_l, b_{s,n,t}^{+}, b_{s,n,t}^{-}} (f_1 + ... + f_7),$$

where

$$f_1=\sum_{g,n} (I_g^G+M_g^G)\bar{p}_{g,n}$$

$$f_2=\sum_{g,n,t} C_g^G p_{g,t,n} τ_{t}$$

$$f_3=\sum_{n,t} C σ_{n,t} τ_{t}$$

$$f_4=\sum_{l} (I_l^F+M_l^F) \bar{f}_l$$

$$f_5=\sum_{l,t} C_l^F ⋅ |f_{l,t}| τ_{t}$$

$$f_6=\sum_{s,n} I_s^S \bar{b}_{s,n}$$

$$f_7=\sum_{s,n,t} C_s^S (b_{s,n,t}^{+}+b_{s,n,t}^{-}) τ_{t}$$


### Constraints
#### Balance
Energy balance $t=1$

$$\sum_{g} p_{g,n,t} + σ_{n,t} + \sum_{(i,j)=l∈L∣j=n} f_{l,t} - \sum_{(i,j)=l∈L∣i=n} f_{l,t} + ξ_s b_{s,n,t} = D_{n,t},\quad ∀s,n,t=1$$

Energy balance $t>1$

$$\sum_{g} p_{g,n,t} + σ_{n,t} + \sum_{(i,j)=l∈L∣j=n} f_{l,t} - \sum_{(i,j)=l∈L∣i=n} f_{l,t} + ξ_s (b_{s,n,t}-b_{s,n,t-1}) = D_{n,t},\quad ∀s,n,t>1$$

#### Generation / Shedding
Generation capacity

$$p_{g,n,t} ≤ A_{g,n,t} (Q_{g,n} + \bar{p}_{g,n}),\quad ∀g,n,t$$

Minimum renewables share

$$\sum_{g∈G^r,n,t} p_{g,n,t} ≥ κ \sum_{g,n,t} p_{g,n,t}$$

Shedding upper bound

$$σ_{n,t} ≤ \bar{C} D_{n,t},\quad ∀n,t$$

#### Transmission
Transmission capacity

$$f_{l,t} ≤ \bar{f}_l,\quad ∀l,t$$

$$f_{l,t} ≥ -\bar{f}_l,\quad ∀l,t$$

The absolute value of the transmission

$$|f_{l,t}|≥f_{l,t},\quad ∀l,t$$

$$|f_{l,t}|≥-f_{l,t},\quad ∀l,t$$

#### Storage
Charge and discharge at $t=1$

$$\begin{aligned}
& b_{s,n,t}^{+}≥b_{s,n,t} - b_{s,n}^0,\quad ∀s,n,t=1 \\
& b_{s,n,t}^{-}≥b_{s,n,t} - b_{s,n}^0,\quad ∀s,n,t=1
\end{aligned}$$

Charge and discharge at $t>1$

$$\begin{aligned}
& b_{s,n,t}^{+}≥b_{s,n,t} - b_{s,n,t-1},\quad ∀s,n,t>1 \\
& b_{s,n,t}^{-}≥b_{s,n,t} - b_{s,n,t-1},\quad ∀s,n,t>1
\end{aligned}$$

Storage capacity

$$b_{s,n,t}≤\bar{b}_{s,n},\quad ∀s,n,t$$

Storage continuity

$$b_{s,n,t=1} = b_{s,n,t=t_{end}},\quad ∀s,n$$


#### Ramping Limits

$$\begin{aligned}
p_{g,n,t} - p_{g,n,t-1} &≥ r_g^{+}, \quad ∀g,n,t>1 \\
p_{g,n,t} - p_{g,n,t-1} &≤ -r_g^{-}, \quad ∀g,n,t>1
\end{aligned}$$

#### Voltage Angles

$$(θ_{n,t} - θ_{n',t}') B_l = p_{g,n,t} - p_{g,n',t}, \quad ∀g,l,n,n',t>1$$


## Input
Users can provide input parameters for different instances as a directory containing CSV and JSON files, and also include a README file, which describes the instance. Users can distribute instances as `.zip` archives.

We have `instance` directory, with files:

- `indices.json` -- Indices.
- `nodes/` -- Time clustered data from the nodes.
  - `1.csv`
  - `2.csv`
  - ...
- `constants.json` -- Constant parameters.
- `transmission.csv` -- Transmission parameters.
- `technology.csv` -- Technology parameters.
- `storage.csv` -- Storage parameters.
- `README.md` -- Description about the instance.

The parameters naming convention is documented in the parameters section.


## API
```@docs
Parameters
Specs
load_parameters
energy_system_model
```

## References

[^1]: Pineda, S., & Morales, J. M. (2018). Chronological time-period clustering for optimal capacity expansion planning with storage. IEEE Transactions on Power Systems, 33(6), 7162–7170. https://doi.org/10.1109/TPWRS.2018.2842093
