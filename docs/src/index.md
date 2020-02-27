# Energy System Model
Mathematical reference for the energy system model. The model presented here is based on the model in [^1]. We express units using square brackets.

## Utility
We calculate annualized costs using [*equivalent annual cost (EAC)*](https://en.wikipedia.org/wiki/Equivalent_annual_cost) formula

$$EAC(c,r,n) = \frac{c}{a_{n,r}},\quad a_{n,r} = \frac{1-(1+r)^{-n}}{r},\quad a_{n,0}=n,$$

where $c$ is the net present cost of the project, $n$ is the number of payments, and $r$ is the interest rate.

## Indices and Sets
Indices and sets define the different objects and dimensions in the model.

*  $g∈G$: Generation technologies
*  $G^r⊆G$: Renewable generation technologies
*  $n∈N$: Nodes
*  $l∈L$: Transmission lines, bidimensional vectors $(i,j)$ where $i,j∈N$
*  $t∈T$: Time steps, depending on the number of clusters per month
*  $s∈S$: Storage technologies

In the code, we store both indices and parameters in the [`Parameters`](@ref) struct.

## Parameters
Constant parameters

*  $κ∈[0,1]$: Renewables participation required by the system
*  $C$: Shedding cost [€/MWh]
*  $\bar{C}$: Shedding capacity [€/MWh]
*  $r≥0$: Interest rate

Time clustered parameters

*  $τ$: Number of time slices per month
*  $τ_{t}$: Duration of time period $t$ [h]
*  $Q_{g,n}$: Initial capacity [MW]
*  $A_{g,n,t}∈[0,1]$: Availability of technology $g$ per node $n$ at time step $t$
*  $D_{n,t}$: Clustered demand per node $n$ per time step $t$  [MWh]

Generation technology parameters

*  $I_g^G$: Annualised investment cost for generation per MW of technology $g$ [€/MW]. Calculated as $I_g^G=EAC(c_g, t_g, r)$ where $c_g$ is the cost and $t_g$ is lifetime of technology $g$.
*  $M_g^G$: Annualised maintenance cost for generation per MW of technology $g$ [€/MW]
*  $C_g^G$: Operational cost per MWh of technology $g$ [€/MWh]. Calculated as $c_g/c'_g/1000$ where $c_g$ is fuel cost 1 and $c'_g$ fuel cost 2 of technology $g.$
*  $r_g^{-}$: Relative ramp-down limit of technology $g$
*  $r_g^{+}$: Relative ramp-up limit of technology $g$

Transmission parameters

*  $I_l^F$: Annualised investment cost for transmission per line $l$ [€/MW]. Calculated as $I_l^F=EAC(c_l⋅d_l + M_l^F, t_l, r)$ where $c_l$ is cost per kilometer and $d_l$ distance in kilometers, $t_l$ the lifetime of transmission line $l.$
*  $M_l^F$: Annualised maintenance cost for transmission per line $l$ [€/MW]
*  $C_l^F$: Transmission cost per line $l$ [€/MWh]
*  $B_l$: Susceptance per line $l$

Storage parameters

*  $I_s^S$: Annualised investment cost of storage technology $s$ per MW [€/MWh]. Calculated as $I_s^S=EAC(c_s, t_s, r)$ where $c_s$ is the cost and $t_s$ the lifetime of storage $s.$
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
*  $\bar{b}_{s,n}≥0$: Storage capacity of storage $s$ at node $n$ [MW]
*  $b_{s,n,t}^{+}≥0$: Charging of storage $s$ at node $n$ in each time step $t$ [MW]
*  $b_{s,n,t}^{-}≥0$: Discharging of storage $s$ at node $n$ in each time step $t$ [MW]

Voltage angle variables

*  $θ_{n,t}≥0$: Voltage angle at node $n$ in each time step $t$
*  $θ'_{n,t}≥0$: Voltage angle at node $n$ in each time step $t$

## Objective
We define the objective as cost minimization

$$\mathrm{minimize}_{p_{g,t}, \bar{p}_g, σ_{t}, f_{l,t}, \bar{f}_l, b_{s,n,t}^{+}, b_{s,n,t}^{-}} (f_1 + ... + f_7).$$

The individual objectives are defined as follows.

Investment and maintenance cost of generation capacity

$$f_1=\sum_{g,n} (I_g^G+M_g^G)\bar{p}_{g,n}$$

The operational cost of generation dispatch

$$f_2=\sum_{g,n,t} C_g^G p_{g,t,n} τ_{t}$$

Shedding cost

$$f_3=\sum_{n,t} C σ_{n,t} τ_{t}$$

Investment and maintenance cost of transmission capacity

$$f_4=\sum_{l} (I_l^F+M_l^F) \bar{f}_l$$

The operational cost of transmission flow

$$f_5=\sum_{l,t} C_l^F ⋅ |f_{l,t}| τ_{t}$$

Investment cost of storage capacity

$$f_6=\sum_{s,n} I_s^S \bar{b}_{s,n}$$

The operational cost of storage

$$f_7=\sum_{s,n,t} C_s^S (b_{s,n,t}^{+}+b_{s,n,t}^{-}) τ_{t}$$


## Constraints
Within the code, we use [`Specs`](@ref) struct to control whether we include certain constraints in the model.

### Balance
Energy balance $t=1$

$$\sum_{g} p_{g,n,t} + σ_{n,t} + \sum_{(i,j)=l∈L∣j=n} f_{l,t} - \sum_{(i,j)=l∈L∣i=n} f_{l,t} + ξ_s b_{s,n,t} = D_{n,t},\quad ∀s,n,t=1$$

Energy balance $t>1$

$$\sum_{g} p_{g,n,t} + σ_{n,t} + \sum_{(i,j)=l∈L∣j=n} f_{l,t} - \sum_{(i,j)=l∈L∣i=n} f_{l,t} + ξ_s (b_{s,n,t}-b_{s,n,t-1}) = D_{n,t},\quad ∀s,n,t>1$$

### Generation / Shedding
Generation capacity

$$p_{g,n,t} ≤ A_{g,n,t} (Q_{g,n} + \bar{p}_{g,n}),\quad ∀g,n,t$$

Minimum renewables share

$$\sum_{g∈G^r,n,t} p_{g,n,t} ≥ κ \sum_{g,n,t} p_{g,n,t}$$

Shedding upper bound

$$σ_{n,t} ≤ \bar{C} D_{n,t},\quad ∀n,t$$

### Transmission
Transmission capacity

$$f_{l,t} ≤ \bar{f}_l,\quad ∀l,t$$

$$f_{l,t} ≥ -\bar{f}_l,\quad ∀l,t$$

The absolute value of the transmission

$$|f_{l,t}|≥f_{l,t},\quad ∀l,t$$

$$|f_{l,t}|≥-f_{l,t},\quad ∀l,t$$

### Storage
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


### Ramping Limits
Ramping limit up and down

$$\begin{aligned}
p_{g,n,t} - p_{g,n,t-1} &≥ r_g^{+}, \quad ∀g,n,t>1 \\
p_{g,n,t} - p_{g,n,t-1} &≤ -r_g^{-}, \quad ∀g,n,t>1
\end{aligned}$$

### Voltage Angles
Faraday law for accounting voltage angles

$$(θ_{n,t} - θ_{n',t}') B_l = p_{g,n,t} - p_{g,n',t}, \quad ∀g,l,n,n',t>1$$


## Instance
Users can provide input parameters for different instances as a directory containing CSV and JSON files, and also include a README file, which describes the instance. Users can distribute instances as `.zip` archives. We recommend checking out the example instance in `examples/instance` for reference. We describe the input format in [`load_parameters`](@ref) and the resulting output format in [`save_results`](@ref).

## References
[^1]: Pineda, S., & Morales, J. M. (2018). Chronological time-period clustering for optimal capacity expansion planning with storage. IEEE Transactions on Power Systems, 33(6), 7162–7170. https://doi.org/10.1109/TPWRS.2018.2842093
