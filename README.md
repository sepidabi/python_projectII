# Project II codes
THIS REPOSITORY INCLUDES (for now) THE PYTHON SCRIPTS THAT I AM WRITING FOR MY SECOND PHD PROJECT: **The time evolution and oscillations of the Ca II K bright fibrils.**

>Note: We assume here that the cuts across the fibrils have already been defined via crispex. So we are picking up from there. Also the inversions on the space-time intensity maps of the cuts have been done and stored

## fibril_big_oscillation.py: A user interactive code to trace the priliminary path of the fibrillar oscillations
- The time vs. intensity image along the each cut (defined by user) get extracted and enhaced to bring out the oscillations.
- User clicks along the path of of each oscillation
- The coords of the traced path and smoothed path is saved in .txt format
- Once done with tracing, the enhanced image of time-intensity with overploted smoothed paths and their code is saved.

## fibril_oscillation_properties.py: An automated object-oriented code to extract/calculate the fibrillar oscillation properties by using the coords of the traced paths and the inversion results.
- Extacts the inversion results along the **traced smoothed** path of each oscillation
- Calculates the pos velocity
- Calculates the self-correlation of the oscillations
- Detects the extremum points in the velocity oscillations of both los and pos
- Saves the corresponding information in object class in cut*.obj format
> Note: no oscillation at this point is excluded and the characteristics that are saved are purely calculated

The class under which each object assigned to each oscillation follows as:
>class oscillation:
    def init(self, time, vpos_f, vlos_f, ncc_pos, ncc_los, pos_xmax, pos_ymax, pos_xmin, pos_ymin, pos_ncc_xmax, pos_ncc_ymax, los_xmax, los_ymax, los_xmin, los_ymin, los_ncc_xmax, los_ncc_ymax, osc_mid_coord, phase, osc_cat):

## fibril_object_corrected.py: checks the oscillations automatically and excludes the ones that are not acceptable and saves the results in new object format
- Extract the objects save by "object_oscillation_properties.py" and check if
   - the oscillation is valid
   - if there is oscillation in los direction
   - if the oscillation in los and pos in each case are correlated
- the oscillations info such as Period, phase and amplitude in new object class as following:
>class oscillation_result:
    def init(self,
                 per_los,
                 per_los_m,
                 per_los_dom,
                 per_pos,
                 per_pos_m,
                 per_pos_dom,
                 A_los_m,
                 A_pos_m,
                 dphi,
                 ncc_min,
                 ncc_max,
                 ncc,
                 mid_coord,
                 cor_rate,
                 theta,
                 per_m,
                 per_dom,
                 ratio,
                 dist,
                 ):

## fibril_oscillation_final_check.py: Provides the chance to visually inspect the traced fibrils to make sure that the automatic judjemt in oscillations validity is actually correct
# So the oscillations in general iclude the following categories:
- the ones that oscillate long enough to define a dphi for them in [-1,1]
- The ones that oscillate in both direcitions long enough but are not correlated (dphi = -5)
- The ones where the oscillations are so short (less than one full period) that it is not possible to define dphi (dphi=-10)

The purpose of this script is to:
- make sure that the oscillations where dphi is automatically set to -10 are valid oscillations.
- ensure that the oscillations have defined dphi [-1,1] or -5, also have convincing oscillations path 
> The oscillations that were not subjectively acceptable, to the "removed_oscillations" folder

## fibril_object_results.py: plots the statistical anlysis results

## fibril_single_plot_new.py: plotting the proto-oscillations from each categories in 2&3D

## fibril_hor_oscilation.py: tracking oscillation in the one roi example

## fibril_single_oscillation.py: plot the one specific oscillation inversion and tracking results