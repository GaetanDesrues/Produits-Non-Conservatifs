# Fichier gnuplot pour l'affichage des courbes

set terminal pdf
set output 'Output/ParesRus2_VonLeer.pdf'

syl = 'Densité Vitesse Pression'
col = '6 7 10'
fileName = 'PCR2_VonLeer'
titre = "Forme non-conservative - Schéma de Pares-Rusanov d'ordre 2 (Von Leer)"

set multiplot layout 2, 2 title titre font ",14"
set yrange[0:2.5]

do for [k=0:3]{
    a = k*33
    t = 2*a # ms
    set title 't = '.t.' ms'
    plot for [c=2:4] "Output/Non-Conservative/".fileName."_it=".a.".txt" using 1:c with lines lw 1.5 lt word(col,c-1) title word(syl,c-1)
}








###
# Fichier gnuplot pour l'affichage des courbes


# syl = 'Densité Vitesse Pression'
# col = '6 7 10'
# a4 = 100#`find ./Output/Sol -type f | wc -l` - 1
# a3 = a4/3
# a2 = 2*a3
# a1 = 0
# fileName = 'PCR2_VonLeer'
#
# set multiplot layout 2, 2 title "Forme non-conservative - Schéma de Pares-Rusanov Von Leer" font ",14"
# # set yrange[0:2.5]
# # set xrange[0.4:0.6]
# # set yrange[-10:10]
# set title "n = ".a1
# plot for [c=2:4] "Output/Sol/".fileName."_it=".a1.".txt" using 1:c with lines lw 1.5 lt word(col,c-1) title word(syl,c-1)
# set title "n = ".a3
# plot for [c=2:4] "Output/Sol/".fileName."_it=".a3.".txt" using 1:c with lines lw 1.5 lt word(col,c-1) title word(syl,c-1)
# set title "n = ".a2
# plot for [c=2:4] "Output/Sol/".fileName."_it=".a2.".txt" using 1:c with lines lw 1.5 lt word(col,c-1) title word(syl,c-1)
# set title "n = ".a4
# plot for [c=2:4] "Output/Sol/".fileName."_it=".a4.".txt" using 1:c with lines lw 1.5 lt word(col,c-1) title word(syl,c-1)
