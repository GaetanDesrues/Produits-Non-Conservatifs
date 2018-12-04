#set yrange [0:2]

plot "Output/Sol_t=0.00.txt" using 1:2 with lines title 't=0.00'
replot "Output/Sol_t=0.01.txt" using 1:2 with lines title 't=0.01'
replot "Output/Sol_t=0.02.txt" using 1:2 with lines title 't=0.02'
replot "Output/Sol_t=0.03.txt" using 1:2 with lines title 't=0.03'
replot "Output/Sol_t=0.04.txt" using 1:2 with lines title 't=0.04'
replot "Output/Sol_t=0.05.txt" using 1:2 with lines title 't=0.05'
#replot "Output/Sol_t=0.20.txt" with lines title 't=0.02'
#replot "Output/Sol_t=0.10.txt" with lines title 't=0.10'
#replot "Output/Sol_t=0.20.txt" with lines title 't=0.20'


#replot "Output/Sol_t=1.10.txt" with lines title 'Sol-exacte-t=0.10'
#replot "Output/Sol_t=1.20.txt" with lines title 'Sol-exacte-t=0.20'

#replot "Output/Sol_t=0.50.txt" with lines title 't=0.50'
#replot "Output/Sol_t=1.50.txt" with lines title 'Sol-exacte-t=0.50'
#replot "Output/Sol_t=1.00.txt" with lines title 't=1.00'
