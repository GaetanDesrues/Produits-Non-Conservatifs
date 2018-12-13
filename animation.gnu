it = it+1
file = sprintf('Output/Sol_it=%d.txt', it)
plot file using 1:4 w l
pause 0.05
if (it<it_max) reread
