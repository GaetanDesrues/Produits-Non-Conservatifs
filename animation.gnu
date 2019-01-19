it = it+1
file = sprintf('Output/Sol/Sol_it=%d.txt', it)
set yrange [0:1]
plot file using 1:2 w l
pause 0.05
if (it<it_max) reread
