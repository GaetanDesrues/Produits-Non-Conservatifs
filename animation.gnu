it = it+1
file = sprintf('Output/Sol_it=%d.txt', it)
plot file using 1:3 w l
pause 0.2
if (it<it_max) reread
