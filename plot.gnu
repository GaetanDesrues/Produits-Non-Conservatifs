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
    t = 2*a
    set title 't = '.t.' ms'
    plot for [c=2:4] "Output/Non-Conservative/".fileName."_it=".a.".txt" using 1:c with lines lw 1.5 lt word(col,c-1) title word(syl,c-1)
}
