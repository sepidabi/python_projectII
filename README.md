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

class oscillation:
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
