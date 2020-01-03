# EnergySystemModel.jl
Documentation for EnergySystemModel.jl


## Model
Mathematical reference of the energy system model. The model presented here is similar to the model in [^1].

### Indexes and Sets
*  $g∈G$: Generation technologies
*  $G^r⊆G$: Renewable generation technologies
*  $n∈N$: Nodes
*  $l∈L$: Transmission lines, bidimensional vectors $(i,j)$ where $i,j∈N$
*  $t∈T$: Time steps, depending on the number of clusters per month

### Parameters
*  $A_g∈\{0,1\}$: Availability of technology $g$
*  $I_g^G$: Annualised investment cost for generation per MW of technology $g$ [€/MW]
*  $M_g^G$: Annualised maintenance cost for generation per MW of technology $g$ [€/MW]
*  $C_g^G$: Operational cost per MWh of technology $g$ [€/MWh]
*  $τ_{t}$: Cluster size of $t$
*  $C$: Shedding cost [€/MWh]
*  $D_{n,t}$: Clustered demand per node $n$ per time step $t$  [MWh]
*  $κ∈[0,1]$: Renewables participation required by the system
*  $I_l^F$: Annualised investment cost for transmission per line $l$ [€/MW]
*  $M_l^F$: Annualised maintenance cost for transmission per line $l$ [€/MW]
*  $C_l^F$: Transmission cost per line $l$
*  $ξ$: Battery's roundtrip efficiency
*  $I^S$: Annualised investment cost for storage per MW [€/MWh]
*  $C^S$: Storage operational cost [€/MW]
*  $b_{n}^0$: Initial storage at node $n$ [MWh]
*  $r_g^{-}$: Relative ramp-down limit of technology $g$
*  $r_g^{+}$: Relative ramp-up limit of technology $g$
*  $B_l$: Susceptance per line $l$

### Variables
*  $p_{g,n,t}≥0$: Dispatch from technology $g$ at node $n$ in each time step $t$ [MWh]
*  $\bar{p}_{g,n}≥0$: Generation capacity invested in each technology $g$ at node $n$ [MW]
*  $σ_{n,t}≥0$: Loss of load at node $n$ in each time step $t$ [MWh]
*  $f_{l,t}≥0$: Transmission flow per line $l$ in each time step $t$ [MWh]
*  $\bar{f}_l≥0$: Transmission capacity per line $l$ [MW]
*  $b_{n,t}≥0$: Storage level at node $n$ in each time step $t$ [MWh]
*  $\bar{b}_{n}≥0$: Storage capacity at node $n$ [MWh]
*  $b_{n,t}^{+}≥0$: Charging at node $n$ in each time step $t$ [MWh]
*  $b_{n,t}^{-}≥0$: Discharging at node $n$ in each time step $t$ [MWh]
*  $θ_{n,t}≥0$: Voltage angle at node $n$ in each time step $t$
*  $θ'_{n,t}≥0$: Voltage angle at node $n$ in each time step $t$

### Objective
Minimize for $p_{g,t}, \bar{p}_g, σ_{t}, f_{l,t}, \bar{f}_l, b_{n,t}^{+}, b_{n,t}^{-}$

$$\begin{aligned}
& \sum_{g,n} (I_g^G+M_g^G)\bar{p}_{g,n} + \\
& \sum_{g,n,t} C_g^G p_{g,t,n} τ_{t} + \\
& \sum_{n,t} C σ_{t} τ_{t} + \\
& \sum_{l} (I_l^F+M_l^F) \bar{f}_l + \\
& \sum_{l,t} C_l^F ⋅ f_{t,l} τ_{t} + \\
& \sum_{n} I^S \bar{b}_n + \\
& \sum_{n,t} C^S (b_{t,n}^{+}+d_{t,n}^{-}) τ_{t}
\end{aligned}$$


### Constraints
#### Balance
Energy balance $t=1$

$$\sum_{g} p_{g,n,t} + σ_{n,t} + \sum_{(i,j)=l∈L∣j=n} f_{l,t} - \sum_{(i,j)=l∈L∣i=n} f_{l,t} + ξ b_{n,t} = D_{n,t},\quad ∀n,t=1$$

Energy balance $t>1$

$$\sum_{g} p_{g,n,t} + σ_{n,t} + \sum_{(i,j)=l∈L∣j=n} f_{l,t} - \sum_{(i,j)=l∈L∣i=n} f_{l,t} + ξ (b_{n,t}-b_{n,t-1}) = D_{n,t},\quad ∀n,t>1$$

#### Generation / Shedding
Generation capacity

$$p_{g,n,t} ≤ A_g \bar{p}_g,\quad ∀g,n,t$$

Minimum renewables share

$$\sum_{g∈G^r,n,t} p_{g,n,t} ≥ κ \sum_{g,n,t} p_{g,n,t}$$

Shedding upper bound

$$σ_{n,t} ≤ C D_{n,t},\quad ∀n,t$$

#### Transmission
Transmission capacity

$$f_{l,t} ≤ \bar{f}_l,\quad ∀l,t$$

#### Storage
Charge and discharge at $t=1$

$$\begin{aligned}
& b_{n,t}^{+}≥b_{n,t} - b_{n}^0,\quad ∀n,t=1 \\
& b_{n,t}^{-}≥b_{n,t} - b_{n}^0,\quad ∀n,t=1
\end{aligned}$$

Charge and discharge at $t>1$

$$\begin{aligned}
& b_{n,t}^{+}≥b_{n,t} - b_{n,t-1},\quad ∀n,t>1 \\
& b_{n,t}^{-}≥b_{n,t} - b_{n,t-1},\quad ∀n,t>1
\end{aligned}$$

Storage capacity

$$b_{n,t}≤\bar{b}_n,\quad ∀n,t$$

Storage

$$b_{n,t=t_{end}} = b_{n,t=0},\quad ∀n$$


#### Ramping Limits

$$\begin{aligned}
p_{g,n,t} - p_{g,n,t-1} &≥ r_g^{+}, \quad ∀g,n,t>1 \\
p_{g,n,t} - p_{g,n,t-1} &≤ -r_g^{-}, \quad ∀g,n,t>1
\end{aligned}$$

#### Voltage Angles

$$(θ_{n,t} - θ_{n',t}') B_l = p_{g,n,t} - p_{g,n',t}, \quad ∀g,l,n,n',t>1$$


## API
```@docs
Specs
```

```@docs
energy_system_model
```

## References

[^1]: Pineda, S., & Morales, J. M. (2018). Chronological time-period clustering for optimal capacity expansion planning with storage. IEEE Transactions on Power Systems, 33(6), 7162–7170. https://doi.org/10.1109/TPWRS.2018.2842093
