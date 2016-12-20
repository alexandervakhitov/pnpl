Synthetic experiments list

In the folder 'main' you will find

A) Fig. 3 in the paper (Points+lines. Accuracy w.r.t. image noise level)

General case:

main_lines_pts_noise.m accuracy w.r.t. noise 
main_lines_pts_eq_constr_noise.m accuracy w.r.t. noise for equal number of constraints

Planar case:

main_planar_lines_pts_noise.m accuracy w.r.t. noise 
main_planar_eq_constr_lines_pts_noise.m accuracy w.r.t. noise for equal number of constraints

B) Fig. 4 in the paper (Points+lines. Accuracy w.r.t. increasing number
of point or line correspondences)

General case:

main_lines_pts_small_number.m accuracy w.r.t. number of points and lines
main_lines_pts_eq_constr_small_num.m accuracy w.r.t. number of points and lines for equal number of constraints

Planar case:

main_planar_lines_pts_small_num.m accuracy w.r.t. number of points and lines
main_planar_lines_pts_eq_constr_small_num.m accuracy w.r.t. number of points and lines for equal number of constraints

C) Fig. 5 in the paper (Lines only. Accuracy w.r.t. image noise level)

General case:

main_lines_noise.m accuracy w.r.t. noise 

Planar case:

main_planar_lines_noise.m accuracy w.r.t. noise 


D) Fig. 6 in the paper (Timing of the algorithms)
main_lines_pts_number.m - timing of all the evaluated algorithms
main_detailed_timing.m - detailed timing of the EPnPL and OPnPL

E) Fig. 7 in the paper (Robustness to the amount of shift between the 3D reference lines and the
projected ones.)

main_lines_shiftvar.m 

F) Fig. 1 in the Supplemental (Pose from only line correspondences. Accuracy (mean and median error) w.r.t. increasing number of point or line correspondences. The image noise std is set to 1 pixel for all experiments.)

main_lines_small_num.m general case
main_planar_lines_small_num.m planar case



