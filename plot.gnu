
syl = 'Densité Vitesse Pression'
col = '6 7 10'


set multiplot layout 2, 2 title "Equations d'Euler en coordonnées Lagrangiennes sous forme conservative" font ",14"
set yrange[0:1.2]
set title "n = 0"
plot for [c=2:4] "Output/Sol_it=0.txt" using 1:c with lines lw 1.5 lt word(col,c-1) title word(syl,c-1)
set title "n = 10"
plot for [c=2:4] "Output/Sol_it=10.txt" using 1:c with lines lw 1.5 lt word(col,c-1) title word(syl,c-1)
set title "n = 20"
plot for [c=2:4] "Output/Sol_it=20.txt" using 1:c with lines lw 1.5 lt word(col,c-1) title word(syl,c-1)
set title "n = 30"
plot for [c=2:4] "Output/Sol_it=30.txt" using 1:c with lines lw 1.5 lt word(col,c-1) title word(syl,c-1)
