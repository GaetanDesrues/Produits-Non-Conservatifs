#set yrange [0:2]

plot "Output/Sol_it=0.txt" using 1:3 with lines title 't=0.00'
replot "Output/Sol_it=20.txt" using 1:3 with lines title 't=0.1'
