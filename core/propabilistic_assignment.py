
"""

For each single peak, set a probability distribution over all the possibile assignment.
Given a function phi(p,p') = {1/dist(p, p')} / {\sum_s 1/dist(p, s)}


definiamo le probabilità in base allo scarto con il valor medio

anzichè assegnare il particolare picco tale per cui lo scarto con il valor medio sia minimo,
ovvero che minimizzi la seguente funzione di costo: cost(p,p') = ||csp(i) - avg_csp(i)||,

noi definiamo una distribuzione di probabilità  P = {1/cost(p, p')} / {\sum_s 1/cost(p, s)}

Poi realizziamo gli assegnamenti usando beam-search
Quindi aggiorniamo tutte le distribuzioni
Ripetiamo la cosa fino a convergenza.



"""