**This is a trivially small compilation of MATLAB functions to calculate 'long-period' teleseismic surface-wave observables using classical ray theory, as described in Woodhouse and Wong, 1986: "Amplitude, phase and path anomalies of mantle waves", GJI. But I thought I'd upload these in case they were of use, since I haven't found nice implementations of these calculations anywhere else.**

Anant Hariharan, June 2026

### Amplitudes. 

To generate focusing-predicted amplitudes for any phase velocity map, call the function Calculate_FocusingAmps_WW86_AllinOne_ArbitraryModel.m
This has the following syntax: [Amplist] = Calculate_FocusingAmps_WW86_AllinOne_ArbitraryModel(lon_src,lat_src,lon_stas,lat_stas,lon_c,lat_c,c)

The 7 inputs are as follows: 
-lon_src = The longitude of the Earthquake source
- lat_src: The latitude of the Earthquake source
-lon_stas: The longitudes of the stations, or coordinates at which you
 want to predict the amplitudes
- lat_stas:The latitudes of the stations, or coordinates at which you
% want to predict the amplitudes
-lon_c: Longitudes at which the phase velocity map is defined.
-lat_c: Latitudes at whcih the phase velocity map is defined.
-c: Phase velocity map (km/s or m/s; since this is eventually non-dimensionalized, the units don't matter) at each of the coordinates in lon_c and lat_c

The output, Amplist, is a vector of focusing-predicted amplitudes corresponding to the station locations in  lon_stas and lat_stas.

The code works, following the workflow described in studies such as Laske and Masters and Larson and Ekstrom, by (for every source-station path), first rotating the coordinate system so that the source-receiver path is along the equator. Then, the second derivative of the phase velocity maps transverse to the source-receiver path is calculated using a central difference approach. This, along with a few other trigonometric terms, is used as the integrand over the integral from source to receiver- exactly as in the equation below:

$$ \ln A_{ij}^{F} =
+ \frac{1}{2\sin\Delta} \int_{0}^{\Delta} \left[ \sin(\Delta-\phi)\sin\phi\, \frac{\partial^2}{\partial\theta^2} \left(\frac{\delta c}{c_0}\right) - \cos(\Delta-2\phi) \left(\frac{\delta c}{c_0}\right) \right] \,d\phi $$


Note: If you don't want to specify your own phase velocity map, just use the function Calculate_FocusingAmps_WW86_AllinOne_PreloadModel.m 
In this case, instead of specifying the phase velocity map (lon_c,lat_c,c), just specify the period you're interested in and the earth model. 'GDM52' or 'Ma2014'. The code will then attempt to load in the GDM52 (Ekström, 2011) or the Ma et al., 2014 global phase velocity maps, which are provided for some periods in files with the naming convention '*.pix' or 'map_rayleigh*'. 


### Arrival Angles. 

To generate focusing-predicted amplitudes for any phase velocity map, call the function Calculate_Arrival_Angle_WW86_AllinOne_ArbitraryModel
This has the following syntax: [AAlist] = Calculate_FocusingAmps_WW86_AllinOne_ArbitraryModel(lon_src,lat_src,lon_stas,lat_stas,lon_c,lat_c,c)

The 7 inputs are as follows: 
- lon_src = The longitude of the Earthquake source
- lat_src: The latitude of the Earthquake source
- lon_stas: The longitudes of the stations, or coordinates at which you
 want to predict the amplitudes
- lat_stas:The latitudes of the stations, or coordinates at which you want to predict the amplitudes
- lon_c: Longitudes at which the phase velocity map is defined.
- lat_c: Latitudes at whcih the phase velocity map is defined.
- c: Phase velocity map (km/s or m/s; since this is eventually non-dimensionalized, the units don't matter) at the coordinates in lon_c and lat_c

The output, AAlist, is a vector of Arrival angles (that is, deviations from the direction corresponding to the great-circle path) corresponding to the station locations in  lon_stas and lat_stas.

The code works as described previously, but we use the first derivative of the phase velocity map instead of the second derivative, and some of the other terms in the integrand are different. The convention for the arrival angles is positive clockwise.  

Note: If you don't want to specify your own phase velocity map and are happy just using the global phase velocity maps GDM52 or Ma et al., 2014, which are great, just use the function Calculate_Arrival_Angle_WW86_AllinOne_ArbitraryModel.m 
This works as explained previously for amplitudes. 


## Getting Started

There are three examples you could run to get familiar with the outputs and their visualization!

Example1_FocusingAmplitudes_LowVelocityAnomaly.m predicts surface-wave amplitudes due to a single gaussian anomaly in an otherwise homogeneous phase velocity map, and visualizes these amplitudes. See below for what it should look like, and note the streak of high amplitudes forming downstream of the low-velocity anomaly. 

![alt text](https://github.com/AnantHariharan1996/Woodhouse_and_Wong_RayTheory/blob/main/Example1.png "Example of amplitude predctions")

Example2_FocusingAmplitudes_GDM52.m predicts amplitudes using the GDM52 map at 100s instead of a single anomaly. 

Finally, Example3_WW86_ArrivalAngles_DiffModels.m generates arrival angles for a single event recorded at the TA, using both sets of phase velocity maps at 100 s, and compares the outputs. 


### In progress: Implementation of exact ray-tracing using the formalism in Woodhouse & Wong, 1986.
