# EnergySystemModel.jl
Documentation for EnergySystemModel.jl


## Model
Mathematical reference of the energy system model. The model presented here is similar to the model in [^1].

### Indexes and Sets
*  $n∈N$: Nodes
*  $l∈L$: Transmission lines, bidimensional vectors $(i,j)$ where $i,j∈V$
*  $t∈T$: Time steps, depending on the numbber of clusters per month
*  $g∈G$: Generation technologies
*  $G^r⊆G$: Renewable generation technologies

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

### Variables
*  $p_{g,t,n}$: Dispatch from technology $g$ in each time step $t$ at node $n$ [MWh]
*  $\bar{p}_{g,n}$: Generation capacity invested in each technology $g$ at node $n$ [MW]
*  $σ_{t,n}$: Loss of load in each time step $t$ at node $n$ [MWh]
*  $f_{t,l}$: Transmission flow in each time step $t$ per line $l$ [MWh]
*  $\bar{f}_l$: Transmission capacity per line $l$ [MW]
*  $b_{t,n}$: Storage level in each time step $t$ at node $n$ [MWh]
*  $\bar{b}_{n}$: Storage capacity at node $n$ [MWh]
*  $b_{t,n}^{+}$: Charging in each time step $t$ at node $n$ [MWh]
*  $b_{t,n}^{-}$: Discharging in each time step $t$ at node $n$ [MWh]

### Objective
Minimize for $p_{g,t}, \bar{p}_g, σ_{t}, f_{t,l}, \bar{f}_l, b_{t,n}^{+}, b_{t,n}^{-}$

$$\begin{aligned}
& ∑_{n∈N} ∑_{g∈G} (I_g^G+M_g^G)\bar{p}_{g,n} + \\
& ∑_{n∈N} ∑_{t∈T} ∑_{g∈G} C_g^G p_{g,t,n} τ_{t} + \\
& ∑_{n∈N} ∑_{t∈T} C σ_{t} τ_{t} + \\
& ∑_{l∈L} (I_l^F+M_l^F) \bar{f}_l + \\
& ∑_{t∈T} ∑_{l∈L} KT ⋅ f_{t,l} τ_{t} + \\
& ∑_{n∈N} I^S \bar{b}_n + \\
& ∑_{n∈N} ∑_{t∈T} C^S (b_{t,n}^{+}+d_{t,n}^{-}) τ_{t}
\end{aligned}$$

FIXME: what is KT?

### Constraints
#### Balance
Energy balance $t=1$: $∀t∈\{1\}, n∈N$

$$∑_{g∈G} p_{g,t,n} + σ_{t,n} + ∑_{(i,j)=l∈L∣j=n} f_{t,l} - ∑_{(i,j)=l∈L∣i=n} f_{t,l} + ξ b_{t,n} = D_{t,n}$$

Energy balance $t>1$: $∀t∈t∖\{1\}, n∈N$

$$∑_{g∈G} p_{g,t,n} + σ_{t,n} + ∑_{(i,j)=l∈L∣j=n} f_{t,l} - ∑_{(i,j)=l∈L∣i=n} f_{t,l} + ξ (b_{t,n}-b_{t-1,n}) = D_{t,n}$$

#### Generation / Shedding
Generation capacity: $∀g∈G, t∈T, n∈N$

$$p_{g,t,n} ≤ A_g \bar{p}_g$$

Min RES

$$∑_{n∈N} ∑_{t∈T} ∑_{g∈G^r} p_{g,t,n} ≥ κ ∑_{n∈N} ∑_{t∈T} ∑_{g∈G^r} p_{g,t,n}$$

Shedding upper bound: $∀t∈T, n∈N$

$$σ_{t,n} ≤ C D_{t,n}$$

#### Transmission
Transmission capacity: $∀l∈L, t∈T$

$$f_{t,l} ≤ \bar{f}_l$$

#### Storage
Charge / Discharge $t=1$: $∀t∈\{1\}, n∈N$

$$\begin{aligned}
b_{t,n}^{+}≥b_{t,n} - b_{n}^0 \\
b_{t,n}^{-}≥b_{t,n} - b_{n}^0
\end{aligned}$$

Charge / Discharge $t>1$: $∀t∈T∖\{1\}, n∈N$

$$\begin{aligned}
b_{t,n}^{+}≥b_{t,n} - b_{t-1,n} \\
b_{t,n}^{-}≥b_{t,n} - b_{t-1,n}
\end{aligned}$$

Storage capacity: $∀t∈T, n∈N$

$$b_{t,n}≤\bar{b}_n$$

Storage: $∀n∈N$

$$b_{t=T[end], n} = b_{0, n}$$


## API
```@docs
Specs
```

```@docs
energy_system_model
```

## References

[^1]: Pineda, S., & Morales, J. M. (2018). Chronological time-period clustering for optimal capacity expansion planning with storage. IEEE Transactions on Power Systems, 33(6), 7162–7170. https://doi.org/10.1109/TPWRS.2018.2842093
