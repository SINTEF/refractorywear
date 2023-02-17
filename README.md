# Sidenor refractory wear model
Python code to model the refractory wear of a steel ladle

## Running model
pip install -r requirements.txt
python example.py
## Inputs
The inputs are all located in various input files and read using pandas. The data files contains more data than is currently used, but these data points might be useful in the future

### Time dependent inputs

Variable                   location             Description
----------                 ----------           -----------
Temperatura                \<heat number\>.gz   Measurement of temperature in steel. 
                                                First value is used as a backup if no EAF 
                                                temperature is given, all values are used 
                                                to compare the model to data
time                       \<heat number\>.gz   Preprocessed data
power                      \<heat number\>.gz   Preprocessed data
Consumo_electrico          \<heat number\>.gz   Consumed electricity, energy will go into 
                                                slag (and then steel). This will give a 
                                                temperature increase 
Caudal\_gas                \<heat number\>.gz   The gas rate is used to calculate utau, 
                                                which is used to calculate heat transfer 
                                                coefficients in steel and slag
Fechahora                  \<heat number\>.gz   Time, used to determine what happens when
Heat number                Alloys addition      Need the heat number to connect the rest of data
                           time and weights.csv
Addition time              Alloys addition      At what time additions are added, used to 
                           time and weights.csv determine when to add additions
Mass added per heat (ton)  Alloys addition      Mass of additions
                           time and weights.csv

### Time independent inputs
Variables listed in no particular order

Variable       location          Description
----------     ----------        -----------
heat number    tt.csv            
tapping T (ºC) tt.csv            Temperature from the EAF used as an input temperature
ncol           2019\_Process     Heat number, needed to find other data matching
               Parameters.csv
ncuch          2019\_Process     Ladle number
               Parameters.csv
nusos          2019\_Process     Number of heats since last total rebuild
               Parameters.csv
acero_liquido  2019\_Process     Amount of liquid steel, used as the amount of 
               Parameters.csv    steel in the ladle for simulation
Heatnumber     Ladlewaittime.csv Heat number
Inactive (min) Ladlewaittime.csv

## Time dependent Outputs
Variable       Description
-------------  -----------
TsteelV        Temperature (C) in steel 
TslagV         Temperature (C) in slag 
TInnerW        Temperature (C) in inner most wall layer at brick 5
TOuterW        Temperature (C) in outer most wall layer at brick 5
massstelV      Mass (kg) of steel in ladle
massslagV      Mass (kg) of slag in ladle 
TimeV          Time (s)
TmeasV         Measured temperature (also given in data, this is for plotting convenience)
PowerV         Power added to system
ElectV         Electric power added 
qslagV         Heat flux in slag 
qsteelV        Heat flux in steel
qwallV         Heat flux in wall
qbottomV       Heat flux in bottom
qtot           Total heat flux
TsideV         Temperature in wall
TbottomV       Temperature in bottom
yerode         Fraction of erosion for each brick layer


