# Fichier gnuplot pour l'affichage des courbes


syl = 'Densité Vitesse Pression'
col = '6 7 10'
a4 = `find ./Output/Sol -type f | wc -l` - 1
a3 = a4/3
a2 = 2*a3
a1 = 0

set multiplot layout 2, 2 title "Equations d'Euler en coordonnées Lagrangiennes sous forme conservative" font ",14"
# set yrange[0:3]
# set xrange[0.4:0.6]
# set yrange[-10:10]
set title "n = ".a1
plot for [c=2:4] "Output/Sol/Sol_it=".a1.".txt" using 1:c with lines lw 1.5 lt word(col,c-1) title word(syl,c-1)
set title "n = ".a3
plot for [c=2:4] "Output/Sol/Sol_it=".a3.".txt" using 1:c with lines lw 1.5 lt word(col,c-1) title word(syl,c-1)
set title "n = ".a2
plot for [c=2:4] "Output/Sol/Sol_it=".a2.".txt" using 1:c with lines lw 1.5 lt word(col,c-1) title word(syl,c-1)
set title "n = ".a4
plot for [c=2:4] "Output/Sol/Sol_it=".a4.".txt" using 1:c with lines lw 1.5 lt word(col,c-1) title word(syl,c-1)
