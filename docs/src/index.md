# Energy System Modeling
Mathematical reference for the energy system model. The model presented here is based on the model in [^1]. We express units of parameters and variables using square brackets. In the code, we implement the model as [`EnergySystemModel(::Params, ::Specs)`](@ref) method, which constructs an [`EnergySystemModel`](@ref) instance.

The model considers that each time period might encompass varaible number of hours, as a result of using input data that has been clustered. Otherwise, each time period refers to one hour.

## Indices and Sets
Indices and sets define the different objects and dimensions in the model.

*  $g∈G$: Generation technologies
*  $G^r⊆G$: Renewable generation technologies
*  $n∈N$: Nodes
*  $l∈L$: Transmission lines, refers to a pair of nodes $(i,j)$ where $i,j∈N$
*  $t∈T$: Time periods
*  $s∈S$: Storage technologies

## Parameters
General parameters

*  $κ∈[0,1]$: Fraction of the demand that must be met by renewable generation
*  $C$: Shedding cost [€/MWh]
*  $\bar{C}$: Shedding capacity [MWh]
*  $r≥0$: Interest rate
*  $τ_{t}$: Number of hours clustered in time period $t$ [h]
*  $D_{n,t}$: Demand per node $n$ per time period $t$ [MWh]

Generation technology parameters

*  $I_g^G$: Annualised investment cost for generation per MW of technology $g$ [€/MW]. Calculated as $I_g^G=EAC(c_g, t_g, r)$ where $c_g$ is the cost and $t_g$ is lifetime of technology $g$.
*  $M_g^G$: Annualised maintenance cost for generation per MW of technology $g$ [€/MW]
*  $C_g^G$: Operational cost per MWh of technology $g$ [€/MWh]. Calculated as $c_g/c'_g/1000$ where $c_g$ is fuel cost 1 and $c'_g$ fuel cost 2 of technology $g.$
*  $r_g^{-}$: Relative ramp-down limit of technology $g$
*  $r_g^{+}$: Relative ramp-up limit of technology $g$
*  $Q_{g,n}$: Initial capacity [MW]
*  $A_{g,n,t}∈[0,1]$: Availability of technology $g$ per node $n$ at time step $t$

Transmission parameters

*  $I_l^F$: Annualised investment cost for transmission per line $l$ [€/MW]. Calculated as $I_l^F=EAC(c_l⋅d_l + M_l^F, t_l, r)$ where $c_l$ is cost per kilometer and $d_l$ distance in kilometers, $t_l$ the lifetime of transmission line $l.$
*  $M_l^F$: Annualised maintenance cost for transmission per line $l$ [€/MW]
*  $C_l^F$: Transmission cost per line $l$ [€/MWh]
*  $B_l$: Susceptance per line $l$

Storage parameters

*  $I_s^S$: Annualised investment cost of storage technology $s$ per MW [€/MW]. Calculated as $I_s^S=EAC(c_s, t_s, r)$ where $c_s$ is the (upfront) investment cost and $t_s$ the lifetime of storage $s.$
*  $C_s^S$: Storage operational cost of storage technology $s$ [€/MWh]
*  $b_{s,n}^0$: Initial capacity of storage $s$ at node $n$ [MWh]
*  $ξ_s$: Round-trip efficiency of storage technology $s$


Remarks:

*  In the code, we store both indices and parameters in the [`Params`](@ref) struct. 
*  All annualized costs are calculated using [*equivalent annual cost (EAC)*](https://en.wikipedia.org/wiki/Equivalent_annual_cost) given by

$$EAC(c,r,n) = \frac{c}{a_{n,r}},\quad a_{n,r} = \frac{1-(1+r)^{-n}}{r},\quad a_{n,0}=n,$$

where $c$ is the net present cost of the project, $n$ is the number of periods, and $r$ is the interest rate.

## Variables
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

We use [`Variables`](@ref) struct to store the variable values after optimization. We can query the values from the model using [`Variables(::EnergySystemModel)`](@ref) method.

## Objective
We define the objective as cost minimization

$$\mathrm{minimize}_{p_{g,t}, \bar{p}_g, σ_{t}, f_{l,t}, \bar{f}_l, b_{s,n,t}^{+}, b_{s,n,t}^{-}} (f_1 + ... + f_7).$$

The individual objectives are defined as follows.

Investment and maintenance cost of generation capacity

$$f_1=\sum_{g,n} (I_g^G+M_g^G)\bar{p}_{g,n} \tag{f1}$$

The operational cost of generation dispatch

$$f_2=\sum_{g,n,t} C_g^G p_{g,t,n} τ_{t} \tag{f2}$$

Shedding cost

$$f_3=\sum_{n,t} C σ_{n,t} τ_{t} \tag{f3}$$

Investment and maintenance cost of transmission capacity

$$f_4=\sum_{l} (I_l^F+M_l^F) \bar{f}_l \tag{f4}$$

The operational cost of transmission flow

$$f_5=\sum_{l,t} C_l^F ⋅ |f_{l,t}| τ_{t} \tag{f5}$$

Investment cost of storage capacity

$$f_6=\sum_{s,n} I_s^S \bar{b}_{s,n} \tag{f6}$$

The operational cost of storage

$$f_7=\sum_{s,n,t} C_s^S (b_{s,n,t}^{+}+b_{s,n,t}^{-}) τ_{t} \tag{f7}$$

We use [`Objectives`](@ref) struct to store the objetive values after optimization. We can query the values from the model using [`Objectives(::EnergySystemModel)`](@ref) method.

## Constraints
In this section, we list all the constraints in the energy system model and explain their function. Each constraint is named in the code using the equation labels. We can then access the individual constraints using the standard JuMP syntax, for example, `model[:b1]`.

We use the [`Specs`](@ref) struct to control whether we include certain constraints in the model.

### Balance
Transmission lines to node $n$

$$L_n^-=\{l∈L∣i∈N,(i,n)=l\}$$

Transmission lines from node $n$

$$L_n^+=\{l∈L∣j∈N,(n,j)=l\}$$

Energy balance $t=1$

$$\sum_{g} p_{g,n,t} + σ_{n,t} + \sum_{l∈L_n^-} f_{l,t} - \sum_{l∈L_n^+} f_{l,t} + ξ_s b_{s,n,t} = D_{n,t},\quad ∀s,n,t=1 \tag{b1}$$

Energy balance $t>1$

$$\sum_{g} p_{g,n,t} + σ_{n,t} + \sum_{l∈L_n^-} f_{l,t} - \sum_{l∈L_n^+} f_{l,t} + ξ_s (b_{s,n,t}-b_{s,n,t-1}) = D_{n,t},\quad ∀s,n,t>1 \tag{b2}$$

### Generation
Generation capacity

$$p_{g,n,t} ≤ A_{g,n,t} (Q_{g,n} + \bar{p}_{g,n}),\quad ∀g,n,t \tag{g1}$$

Minimum renewables share

$$\sum_{g∈G^r,n,t} p_{g,n,t} ≥ κ \sum_{g,n,t} p_{g,n,t} \tag{g2}$$

### Shedding
Shedding upper bound

$$σ_{n,t} ≤ \bar{C} D_{n,t},\quad ∀n,t \tag{g3}$$

### Transmission
Transmission capacity

$$f_{l,t} ≤ \bar{f}_l,\quad ∀l,t \tag{t1}$$

$$f_{l,t} ≥ -\bar{f}_l,\quad ∀l,t \tag{t2}$$

The absolute value of the transmission

$$|f_{l,t}|≥f_{l,t},\quad ∀l,t \tag{t3}$$

$$|f_{l,t}|≥-f_{l,t},\quad ∀l,t \tag{t4}$$

### Storage
Charge and discharge at $t=1$

$$b_{s,n,t}^{+}≥b_{s,n,t} - b_{s,n}^0,\quad ∀s,n,t=1 \tag{s1}$$

$$b_{s,n,t}^{-}≥b_{s,n,t} - b_{s,n}^0,\quad ∀s,n,t=1 \tag{s2}$$

Charge and discharge at $t>1$

$$b_{s,n,t}^{+}≥b_{s,n,t} - b_{s,n,t-1},\quad ∀s,n,t>1 \tag{s3}$$

$$b_{s,n,t}^{-}≥b_{s,n,t} - b_{s,n,t-1},\quad ∀s,n,t>1 \tag{s4}$$

Storage capacity

$$b_{s,n,t}≤\bar{b}_{s,n},\quad ∀s,n,t \tag{s5}$$

Storage continuity

$$b_{s,n,t=1} = b_{s,n,t=t_{end}},\quad ∀s,n \tag{s6}$$


### Ramping Limits
Ramping limit up and down

$$p_{g,n,t} - p_{g,n,t-1} ≥ r_g^{+}, \quad ∀g,n,t>1 \tag{r1}$$

$$p_{g,n,t} - p_{g,n,t-1} ≤ -r_g^{-}, \quad ∀g,n,t>1 \tag{r2}$$

### Voltage Angles
Faraday law for accounting voltage angles

$$(θ_{n,t} - θ_{n',t}') B_l = p_{g,n,t} - p_{g,n',t}, \quad ∀g,l,n,n',t>1 \tag{v1}$$


## Instances
Users can provide input parameters for different instances as a directory containing CSV and JSON files, and also include a README file, which describes the instance. Users can distribute instances as `.zip` archives. We provide an example instance in `examples/instance` as a reference. We describe the input format in [`Params(::AbstractString)`](@ref).

We can write [`Specs`](@ref), [`Params`](@ref), [`Variables`](@ref), and [`Objectives`](@ref) structs into JSON files using [`save_json`](@ref) and read them from JSON files using [`load_json`](@ref).


## References
[^1]: Pineda, S., & Morales, J. M. (2018). Chronological time-period clustering for optimal capacity expansion planning with storage. IEEE Transactions on Power Systems, 33(6), 7162–7170. https://doi.org/10.1109/TPWRS.2018.2842093
