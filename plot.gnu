set yrange [0:1.2]

plot "Output/Sol_t=0.00.txt" with lines title 't=0.00'
#replot "Output/Sol_t=0.01.txt" with lines title 't=0.01'
#replot "Output/Sol_t=0.02.txt" with lines title 't=0.02'
replot "Output/Sol_t=0.10.txt" with lines title 't=0.10'
replot "Output/Sol_t=0.20.txt" with lines title 't=0.20'


#replot "Output/Sol_t=1.10.txt" with lines title 'Sol-exacte-t=0.10'
#replot "Output/Sol_t=1.20.txt" with lines title 'Sol-exacte-t=0.20'

#replot "Output/Sol_t=0.50.txt" with lines title 't=0.50'
#replot "Output/Sol_t=1.50.txt" with lines title 'Sol-exacte-t=0.50'
#replot "Output/Sol_t=1.00.txt" with lines title 't=1.00'
