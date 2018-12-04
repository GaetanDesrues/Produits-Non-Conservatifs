#set yrange [0:2]

plot "Output/Sol_it=0.txt" using 1:2 with lines title 't=0.00'
# replot "Output/Sol_it=1.txt" using 1:2 with lines title 't=0.40'
# replot "Output/Sol_it=2.txt" using 1:2 with lines title 't=0.40'
replot "Output/Sol_it=20.txt" using 1:2 with lines title 't=0.40'
# replot "Output/Sol_it=30.txt" using 1:2 with lines title 't=0.40'
#replot "Output/Sol_it=30.txt" using 1:2 with lines title 't=0.30'
#replot "Output/Sol_it=4.txt" using 1:2 with lines title 't=0.04'
#replot "Output/Sol_it=14.txt" using 1:2 with lines title 't=0.05'
#replot "Output/Sol_t=0.20.txt" with lines title 't=0.02'
#replot "Output/Sol_t=0.10.txt" with lines title 't=0.10'
#replot "Output/Sol_t=0.20.txt" with lines title 't=0.20'


#replot "Output/Sol_t=1.10.txt" with lines title 'Sol-exacte-t=0.10'
#replot "Output/Sol_t=1.20.txt" with lines title 'Sol-exacte-t=0.20'

#replot "Output/Sol_t=0.50.txt" with lines title 't=0.50'
#replot "Output/Sol_t=1.50.txt" with lines title 'Sol-exacte-t=0.50'
#replot "Output/Sol_t=1.00.txt" with lines title 't=1.00'
