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
| H11 | Semillas racionales conocidas (V8, A6, CD12 con α ∈ {−110, −253, 15625}) con base solapada. | Mejor V8 a m=9 (0,8155); con solapamiento ninguna mejora: A6 0,847, CD12 0,850. | Negativo. |
| H12 | Racionalidad de todas las soluciones aisladas para n = 19 y n = 20 (patrones de 0,7925 y 0,7755). | 894 + 521 soluciones complejas, 0 racionales. | Negativo. |
| H13 | Haz racional: P = A(x)·(v + u x^k) con u/v raíz racional de una condición de cancelación, 4 rondas, ancho 40, grado ≤ 40. | ~16 000 candidatos; mejor 0,850 (es CD12); empeora con la profundidad. | Negativo. |
| H14 | Muestreo de puntos racionales en las 5 curvas (n = 18, 22, 23, 24, 26; exponentes potenciales 0,800–0,815): fijar ρ racional, resolver el sistema cuadrado por homotopía y comprobar racionalidad a 60 dígitos (`src/rational_point_sampler.py`; el control recupera la semilla de Verdenius). | Tanda local de prueba: ~100 s por muestra y núcleo; 5 muestras (n = 22, 23): 2247 soluciones, 252 reales densas, 0 racionales. | Abandonada: ritmo de ~1 muestra/min con 2 núcleos y probabilidad a priori baja (no hay indicios de que las curvas sean racionales). Sería ejecutable en una máquina propia. |
| H15 | Búsqueda agnóstica a la estructura: minimizar una aproximación suave de \|supp p²\| (continuación en δ más barrera) sin fijar el patrón de ceros (`src/sparsify.py`). | Grado 12, ~9 min: el mínimo hallado es \|p²\| = 14 = n+2 (conteo genérico). **No** redescubre los ejemplos de Coppersmith–Davenport (12): los "milagros" ocupan regiones minúsculas del espacio y el descenso no las alcanza. | Tal cual, peor que la teoría. Haría falta sembrarla con estructura (productos o composiciones) para que sirva. |
