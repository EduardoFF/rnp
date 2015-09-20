

# set terminal png transparent nocrop enhanced font arial 8 size 420,320 
# set output 'scatter.1.png'
set term postscript
set output "data.ps"
set nokey
#set grid

plot "nodes.dat" with points 1 1, "base.dat" with points 3 21, "relay.dat" with points 2 0
