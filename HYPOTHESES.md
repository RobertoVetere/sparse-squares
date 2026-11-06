# Registro de hipótesis

| # | Hipótesis / experimento | Resultado | Decisión |
|---|---|---|---|
| H0 | El exponente sale de una torre p(x)p(x^b)… con semilla densa; el soporte de P_K² son sumas de dígitos con acarreo de supp(p²). | Autómata de acarreos = fuerza bruta (300+400 casos aleatorios). Coppersmith–Davenport grado 12 → λ = 9,87, exp 0,893. | Base de todo lo demás. |
| H0' | Corrección de mi primer resumen: dije que las semillas de Verdenius no podían ser densas. | Falso: su cota es para polinomios completos; la semilla es densa y el mecanismo es el de H0. | Corregido. |
| H1 | Recocido sobre patrones de ceros Z con m = n+1. | |Z| = n reproduce exactamente log 6/log 9 (n = 8) y log 8/log 13 (n = 12). | Mecanismo de Verdenius identificado. |
| H2 | Ansatz simétrico x^n p(1/x) = C p(ρx); enumeración exhaustiva de sistemas cuadrados, n ≤ 20, m = n+1. | Óptimo = Verdenius n = 12; familia λ = n/2+2. Semilla entera n = 8 reproducida exactamente. | Familia agotada. |
| H3 | Mezclar semillas distintas por nivel (mixed radix). | ρ(M₁M₂) nunca mejora al mejor individual (mejor 0,81071). | Descartada. |
| H4 | Cancelaciones de coeficientes en P_K² más allá del soporte. | Verdenius n = 8: conteo exacto = cota hasta K = 5. CD12 y A6 tampoco. | Descartada para estas semillas. |
| H5 | Sistemas simétricos sobredeterminados en 1 ecuación ("milagros"); exponentes potenciales 0,747–0,827. | Levenberg–Marquardt complejo, 150–300 arranques: 0 soluciones no degeneradas. | Descartada (no exhaustivo, pero sin indicios). |
| H6 | **Base solapada m < n+1.** | Exhaustivo n ≤ 18: n = 14, m = 12 → 0,78309. Familia de tres bloques S = A + m{0,1,2}; óptimo n = 20, m = 17 → **0,77552**. | **Seguida.** |
| H6a | ¿Existen semillas reales para esos patrones? | Homotopía de grado total: n = 14 → 28 sol. complejas, 8 reales densas; n = 20 → 259/262 complejas, 42 reales densas. Krawczyk certifica una de cada. | Confirmado. |
| H6b | Densidad de la torre solapada. | Vía de cocientes consecutivos + Möbius: no concluyente (los cocientes se acercan a los valores prohibidos). **Vía ℓ-ádica + lema de suma de dígitos mínima única: demostrada para todo K.** | Resuelto. |
| H6c | Super-semilla P_K apilada sin solapamiento. | Converge a 0,7755, pero necesitaría K ≈ 8–9 explícito. | Innecesaria tras H6b. |
| H7 | Recocido más amplio (n ≤ 32, todas las bases). | Siempre cae en la familia de tres bloques; mejor 0,77552. | Óptimo del marco (heurístico). |
| H8 | Semillas **racionales** con exponente < 0,8155 (récord entero). | Todas las soluciones aisladas, n ≤ 16: solo Verdenius n = 8 es racional. | Negativo. |
| H9 | Curva n = 18, m = 15 (exp 0,81137): puntos racionales vía barrido en ρ ∈ ℚ. | ~20 valores de ρ: 70–78 soluciones complejas, ~20 reales cada uno, ninguna racional (los "candidatos" eran falsos positivos de tolerancia). | Negativo (búsqueda limitada). |
| H10 | ¿Puede este marco llegar a 1/2 (objetivo de Epoch)? | Argumento de conteo heurístico: no; óptimo ≈ 0,77. | Haría falta otra idea (README §5). |
