set font ",12"

#============================================
#
# Set these variables before plotting!
#
#============================================

# dx = (x_max - x_min)/n_points
dx = 0.01

# dt = out_every_time
dt = 0.01

# tf = final_time
tf = 0.5





imax = int(tf/dt)
print imax

do for [it=0:imax] {
    t=it*dt
    plot \
	 'pres.dat'    i it  u 2:3 w l lw 3 t sprintf( "presion, t=%.3f", t ), \
	 'vel.dat'     i it  u 2:3 w l lw 3 t sprintf( "velocidad, t=%.3f", t ), \
	 'density.dat' i it  u 2:3 w l lw 3 t sprintf( "densidad, t=%.3f", t )	 
    print t
    pause -1
}
