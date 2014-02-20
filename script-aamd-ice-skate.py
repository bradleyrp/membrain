# Pseudocode.
# Needs: surface pickle + ion trajectory
# 1. For each ion, from start to end by time (triple loop)
# 		Record the distance from start point by interval length
#		Record if in the zone in separate array (1 or 0)
# 2. Filter array (1 or 0) by some % time in the zone
# 		Only include delta t distances if the mean in the 1/0 array is above some level
