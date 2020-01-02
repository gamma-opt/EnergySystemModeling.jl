# EnergySystemModel.jl
Documentation for EnergySystemModel.jl


## Model
Mathematical reference of the energy system model. The model presented here is similar to the model in [^1].

### Indexes and Sets
*  $g∈G$: Generation technologies
*  $G^r⊆G$: Renewable generation technologies
*  $n∈N$: Nodes
*  $l∈L$: Transmission lines, bidimensional vectors $(i,j)$ where $i,j∈V$
*  $t∈T$: Time steps, depending on the numbber of clusters per month

### Parameters
*  $A_g∈\{0,1\}$: Availability of technology $g$
*  $I_g^G$: Annualised investment cost for generation per MW of technology $g$ [€/MW]
*  $M_g^G$: Annualised maintenance cost for generation per MW of technology $g$ [€/MW]
*  $C_g^G$: Operational cost per MWh of technology $g$ [€/MWh]
*  $τ_{t}$: Cluster size of $t$
*  $C$: Shedding cost [€/MWh]
*  $D_{t,n}$: Clustered demand per time step $t$ per node $n$ [MWh]
*  $κ∈[0,1]$: Renewables participation required by the system (computed over the generation)
*  $I_l^F$: Annualised investment cost for transmission per line $l$ [€/MW]
*  $M_l^F$: Annualised maintenance cost for transmission per line $l$ [€/MW]
*  $ξ$: Battery's roundtrip efficiency
*  $I^S$: Annualised investment cost for storage per MW [€/MWh]
*  $C^S$: Storage operational cost [€/MW]
*  $b_{n}^0$: Initial storage at node $n$ [MWh]
*  $r_g^{-}$: Relative ramp-down limit of technology $g$
*  $r_g^{+}$: Relative ramp-up limit of technology $g$
*  $B_l$: Susceptance per line $l$

### Variables
*  $p_{g,t,n}≥0$: Dispatch from technology $g$ in each time step $t$ at node $n$ [MWh]
*  $\bar{p}_{g,n}≥0$: Generation capacity invested in each technology $g$ at node $n$ [MW]
*  $σ_{t,n}≥0$: Loss of load in each time step $t$ at node $n$ [MWh]
*  $f_{t,l}≥0$: Transmission flow in each time step $t$ per line $l$ [MWh]
*  $\bar{f}_l≥0$: Transmission capacity per line $l$ [MW]
*  $b_{t,n}≥0$: Storage level in each time step $t$ at node $n$ [MWh]
*  $\bar{b}_{n}≥0$: Storage capacity at node $n$ [MWh]
*  $b_{t,n}^{+}≥0$: Charging in each time step $t$ at node $n$ [MWh]
*  $b_{t,n}^{-}≥0$: Discharging in each time step $t$ at node $n$ [MWh]
*  $θ^1_{t,n}≥0$: Voltage angle in each time step $t$ at node $n$
*  $θ^2_{t,n}≥0$: Voltage angle in each time step $t$ at node $n$

### Objective
Minimize for $p_{g,t}, \bar{p}_g, σ_{t}, f_{t,l}, \bar{f}_l, b_{t,n}^{+}, b_{t,n}^{-}$

$$\begin{aligned}
& \sum_{n} \sum_{g} (I_g^G+M_g^G)\bar{p}_{g,n} + \\
& \sum_{n} \sum_{t} \sum_{g} C_g^G p_{g,t,n} τ_{t} + \\
& \sum_{n} \sum_{t} C σ_{t} τ_{t} + \\
& \sum_{l} (I_l^F+M_l^F) \bar{f}_l + \\
& \sum_{t} \sum_{l} KT ⋅ f_{t,l} τ_{t} + \\
& \sum_{n} I^S \bar{b}_n + \\
& \sum_{n} \sum_{t} C^S (b_{t,n}^{+}+d_{t,n}^{-}) τ_{t}
\end{aligned}$$

FIXME: what is KT?

### Constraints
#### Balance
Energy balance $t=1$

$$\sum_{g∈G} p_{g,t,n} + σ_{t,n} + \sum_{(i,j)=l∈L∣j=n} f_{t,l} - \sum_{(i,j)=l∈L∣i=n} f_{t,l} + ξ b_{t,n} = D_{t,n},\quad ∀t=1,n$$

Energy balance $t>1$

$$\sum_{g∈G} p_{g,t,n} + σ_{t,n} + \sum_{(i,j)=l∈L∣j=n} f_{t,l} - \sum_{(i,j)=l∈L∣i=n} f_{t,l} + ξ (b_{t,n}-b_{t-1,n}) = D_{t,n},\quad ∀t>1,n$$

#### Generation / Shedding
Generation capacity

$$p_{g,t,n} ≤ A_g \bar{p}_g,\quad ∀g,t,n$$

Min RES

$$\sum_{n} \sum_{t} \sum_{g∈G^r} p_{g,t,n} ≥ κ \sum_{n} \sum_{t} \sum_{g} p_{g,t,n}$$

Shedding upper bound

$$σ_{t,n} ≤ C D_{t,n},\quad ∀t,n$$

#### Transmission
Transmission capacity

$$f_{t,l} ≤ \bar{f}_l,\quad ∀l,t$$

#### Storage
Charge / Discharge $t=1$

$$\begin{aligned}
& b_{t,n}^{+}≥b_{t,n} - b_{n}^0 \\
& b_{t,n}^{-}≥b_{t,n} - b_{n}^0,\quad ∀t=1, n
\end{aligned}$$

Charge / Discharge $t>1$

$$\begin{aligned}
& b_{t,n}^{+}≥b_{t,n} - b_{t-1,n} \\
& b_{t,n}^{-}≥b_{t,n} - b_{t-1,n},\quad ∀t>1, n
\end{aligned}$$

Storage capacity

$$b_{t,n}≤\bar{b}_n,\quad ∀t,n$$

Storage

$$b_{t=T[end], n} = b_{0, n},\quad ∀n$$

FIXME: what is T[end]?

#### Ramping Limits

$$\begin{aligned}
p_{g,t,n} - p_{g,t-1,n} &≥ r_g^{+}, \quad ∀t>1, n, g \\
p_{g,t,n} - p_{g,t-1,n} &≤ -r_g^{-}, \quad ∀t>1, n, g
\end{aligned}$$

#### Voltage Angles

$$(θ_{t,n}^1 - θ_{t,n'}^2) B_l = p_{g,t,n} - p_{g,t,n'}, \quad ∀t>1,g,l,n,n'$$


## API
```@docs
Specs
```

```@docs
energy_system_model
```

## References

[^1]: Pineda, S., & Morales, J. M. (2018). Chronological time-period clustering for optimal capacity expansion planning with storage. IEEE Transactions on Power Systems, 33(6), 7162–7170. https://doi.org/10.1109/TPWRS.2018.2842093
