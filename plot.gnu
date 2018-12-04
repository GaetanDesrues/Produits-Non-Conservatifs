#set yrange [0:2]

plot "Output/Sol_it=0.txt" using 1:2 with lines title 't=0.00'
replot "Output/Sol_it=500.txt" using 1:2 with lines title 't=1'
