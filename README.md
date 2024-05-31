# II
Project II codes
THIS REPOSITORY INCLUDES (for now) THE PYTHON SCRIPTS THAT I AM WRITING FOR MY SECOND PHD PROJECT:The time evolution and oscillations of the Ca II K bright fibrils.
# python_projectI

Note: We assume here that the cuts across the fibrils have already been defined via crispex. So we are picking up from there. Also the inversions on the space-time intensity maps of the cuts have been done and stored

###fibril_big_oscillation.py
## A user interactive code to trace the priliminary path of the fibrillar oscillations
1- The time vs. intensity image along the each cut (defined by user) get extracted and enhaced to bring out the oscillations.
2- User clicks along the path of of each oscillation
3- The coords of the traced path and smoothed path is saved in .txt format
4- Once done with tracing, the enhanced image of time-intensity with overploted smoothed paths and their code is saved.

###fibril_oscillation_properties.py
## An automated object-oriented code to extract/calculate the fibrillar oscillation properties by using the coords of the traced paths and the inversion results.
1- Extacts the inversion results along the **traced smoothed** path of each oscillation
2- Calculates the pos velocity
3- Calculates the self-correlation of the oscillations
4- Detects the extremum points in the velocity oscillations of both los and pos
5- Saves the corresponding information in object class in cut*.obj format
> Note: no oscillation at this point is excluded and the characteristics that are saved are purely calculated

# The class under which each object assigned to each oscillation follows as:
>class oscillation:
    def __init__(self, time, vpos_f, vlos_f, ncc_pos, ncc_los, pos_xmax, pos_ymax, pos_xmin, pos_ymin, pos_ncc_xmax, pos_ncc_ymax, los_xmax, los_ymax, los_xmin, los_ymin, los_ncc_xmax, los_ncc_ymax, osc_mid_coord, phase, osc_cat):
        self.time = time
        self.vpos_f = vpos_f
        self.vlos_f = vlos_f
        self.ncc_pos = ncc_pos
        self.ncc_los = ncc_los
        self.pos_xmax = pos_xmax
        self.pos_ymax = pos_ymax
        self.pos_xmin = pos_xmin
        self.pos_ymin = pos_ymin
        self.pos_ncc_xmax = pos_ncc_xmax
        self.pos_ncc_ymax = pos_ncc_ymax
        self.los_xmax = los_xmax
        self.los_ymax = los_ymax
        self.los_xmin = los_xmin
        self.los_ymin = los_ymin
        self.los_ncc_xmax = los_ncc_xmax
        self.los_ncc_ymax = los_ncc_ymax
        self.osc_mid_coord = osc_mid_coord
        self.phase = phase
        self.osc_cat = osc_cat

###fibril_object_corrected.py
#checks the oscillations automatically and excludes the ones that are not acceptable and saves the results in new object format
1- Extract the objects save by "object_oscillation_properties.py" and check if
   - the oscillation is valid
   - if there is oscillation in los direction
   - if the oscillation in los and pos in each case are correlated
2- the oscillations info such as Period, phase and amplitude in new object class as following:
>class oscillation_result:
    def __init__(self,
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

###fibril_oscillation_final_check.py
#Provides the chance to visually inspect the traced fibrils to make sure that the automatic judjemt in oscillations validity is actually correct
1- Move the files corresponding to the oscillations that were not subjectively acceptable, to the "removed_oscillations" folder

###fibril_object_results.py
#plot the results

