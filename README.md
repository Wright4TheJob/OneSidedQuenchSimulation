# One Sided Quench Simulation
This small program simulates the transient temperatures through a thick material while quenching one side. This is a 1D transient simulation, so plate edge boundary conditions are ignored. This simulation uses rectangular, not radial, coordinates, so is appropriate for flat material or large-radius curved objects like pipe, but will lose accuracy for geometries where the plate thickness is a significant fraction of the outer diameter.

The program currently assumes quenching of carbon steel but more materials may be added in the future.

The boundary conditions are laid out according to this convention, where $dx = t/n$, and $t$ is the thickness of the material:

| convective boundary | cell | cell | ... | cell | cell | Insulated layer (no heat flow)|
|--|--|--|--|--|--|--|
|x dimension | $dx/2$ | $dx/2 +dx$ ||||$t-dx/2$|

The convective boundary heat flow is driven by the equation:

$$q = h*(T_fluid - T_cell)$$

where negative heat flow means heat leaving the cell.
