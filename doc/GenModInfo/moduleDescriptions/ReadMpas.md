
# ReadMpas
Read MPAS files

<svg width="60.0em" height="8.6em" >
<style>.text { font: normal 1.0em sans-serif;}tspan{ font: italic 1.0em sans-serif;}.moduleName{ font: bold 1.0em sans-serif;}</style>
<rect x="0em" y="0.8em" width="6.0em" height="3.0em" rx="0.1em" ry="0.1em" style="fill:#64c8c8ff;" />
<text x="0.2em" y="2.6500000000000004em" class="moduleName" >ReadMpas</text><rect x="0.2em" y="2.8em" width="1.0em" height="1.0em" rx="0em" ry="0em" style="fill:#c8c81eff;" >
<title>grid_out</title></rect>
<rect x="0.7em" y="3.8em" width="0.03333333333333333em" height="4.0em" rx="0em" ry="0em" style="fill:#000000;" />
<rect x="0.7em" y="7.8em" width="1.0em" height="0.03333333333333333em" rx="0em" ry="0em" style="fill:#000000;" />
<text x="1.9em" y="7.8999999999999995em" class="text" >grid<tspan> (grid_out)</tspan></text>
<rect x="1.4em" y="2.8em" width="1.0em" height="1.0em" rx="0em" ry="0em" style="fill:#c8c81eff;" >
<title>data_out_0</title></rect>
<rect x="1.9em" y="3.8em" width="0.03333333333333333em" height="3.0em" rx="0em" ry="0em" style="fill:#000000;" />
<rect x="1.9em" y="6.8em" width="1.0em" height="0.03333333333333333em" rx="0em" ry="0em" style="fill:#000000;" />
<text x="3.0999999999999996em" y="6.8999999999999995em" class="text" >scalar data<tspan> (data_out_0)</tspan></text>
<rect x="2.5999999999999996em" y="2.8em" width="1.0em" height="1.0em" rx="0em" ry="0em" style="fill:#c8c81eff;" >
<title>data_out_1</title></rect>
<rect x="3.0999999999999996em" y="3.8em" width="0.03333333333333333em" height="2.0em" rx="0em" ry="0em" style="fill:#000000;" />
<rect x="3.0999999999999996em" y="5.8em" width="1.0em" height="0.03333333333333333em" rx="0em" ry="0em" style="fill:#000000;" />
<text x="4.3em" y="5.8999999999999995em" class="text" >scalar data<tspan> (data_out_1)</tspan></text>
<rect x="3.8em" y="2.8em" width="1.0em" height="1.0em" rx="0em" ry="0em" style="fill:#c8c81eff;" >
<title>data_out_2</title></rect>
<rect x="4.3em" y="3.8em" width="0.03333333333333333em" height="1.0em" rx="0em" ry="0em" style="fill:#000000;" />
<rect x="4.3em" y="4.8em" width="1.0em" height="0.03333333333333333em" rx="0em" ry="0em" style="fill:#000000;" />
<text x="5.5em" y="4.8999999999999995em" class="text" >scalar data<tspan> (data_out_2)</tspan></text>
</svg>

## Parameters
|name|description|type|
|-|-|-|
|first_step|first timestep to read|Int|
|last_step|last timestep to read (-1: last)|Int|
|step_increment|number of steps to increment|Int|
|first_rank|rank for first partition of first timestep|Int|
|check_convexity|whether to check convexity of grid cells|Int|
|grid_file|File containing the grid|String|
|data_file|File containing data|String|
|zGrid_file_dir|File containing the vertical coordinates (elevation of cell from mean sea level|String|
|part_file|File containing the grid partitions|String|
|numParts|Number of Partitions|Int|
|numLevels|Number of vertical levels to read|Int|
|altitudeScale|value to scale the grid altitude (zGrid)|Float|
|cellsOnCell|List of neighboring cells (for ghosts))|String|
|var_dim|Dimension of variables (2D, 3D, other)|String|
|Variable_0|Variable_0 ((NONE))|String|
|Variable_1|Variable_1 ((NONE))|String|
|Variable_2|Variable_2 ((NONE))|String|