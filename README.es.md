# sparse-squares

Ataque computacional al problema abierto de Epoch AI *"Dense Polynomials Whose Squares Are Sparse"*
(FrontierMath: Open Problems): hallar un polinomio **denso** p (todos los coeficientes enteros no nulos,
grado n) con |p²| < √n términos.

**Resumen honesto**

| Pregunta | Estado |
|---|---|
| Problema de Epoch (entero, denso, \|p²\| < √n) | **No resuelto.** No hay candidato. Argumento heurístico de por qué este método no llega a 1/2 (§5). |
| Récord **entero/racional** (Verdenius 1949: exponente log 6/log 9 ≈ 0,8155) | **No mejorado.** Reproducido exactamente; búsqueda de semillas racionales negativa (§4). |
| Récord **real** para polinomios completos (Verdenius 1949: log 8/log 13 ≈ 0,81071) | **Mejorado a log 9/log 17 ≈ 0,77552**, con semilla real certificada por aritmética de intervalos y prueba de densidad (§2–3). *Pendiente de revisión por expertos.* |

---

## 1. El mecanismo (Verdenius, reconstruido)

Sea p denso de grado n y S = supp(p²) ⊂ [0, 2n]. Para una base m, la "torre"

  P_K(x) = p(x) · p(x^m) · p(x^{m²}) ··· p(x^{m^{K−1}})

cumple P_K² = ∏ p²(x^{m^j}), así que **supp(P_K²) ⊆ { Σ s_j m^j : s_j ∈ S }**. Con acarreos, ese conjunto
crece como λ^K con λ ≪ |S|. λ se calcula exactamente con un autómata de acarreos (`src/digits.py`,
`src/digits2.py`, contrastado con fuerza bruta). Con m = n+1 (sin solapamiento) P_K es denso trivialmente y
el exponente es log λ / log(n+1).

Para que p² tenga muchos ceros con pocos grados de libertad se usa un **ansatz simétrico**:
x^n p(1/x) = C·p(ρx), con p ≡ √(parte baja de p²) mod x^{…}. Un patrón de soporte es "resoluble
genéricamente" cuando el sistema simétrico es cuadrado.

* Exhaustivo para n ≤ 20 con m = n+1 (`src/enum_sym.py`): el óptimo es **exactamente Verdenius** — n = 8
  (λ = 6, racional) y n = 12 (λ = 8, real). En esa familia λ = n/2 + 2, minimizado en n = 12. Esto explica
  por qué el récord no se movió desde 1949 con este método.
* Semilla entera de Verdenius reproducida: `p = [1,−2,−2,−4,−10,−28,−98,−686,2401]`, |p²| = 9,
  soporte {0,1,6,…,10,15,16}, λ = 6. Torre densa, conteo exacto = autómata hasta K = 5.

## 2. Resultado nuevo: base solapada m < n+1 (exponente real 0,77552)

**Familia de tres bloques.** Para r par, n = 3r + 2, m = n − r/2, el patrón

  S = [0, r] ∪ [m, m+r] ∪ [2m, 2m+r]

da un sistema simétrico cuadrado. Los dígitos son (bloque bajo + acarreo ∈ {0,1,2}) < m, sin acarreos
posteriores, así que **|supp P_K²| ≤ (r+1)·(r+3)^{K−1}·3**, de forma exacta y elemental. El exponente es
log(r+3)/log(m). El óptimo está en r = 6: **n = 20, m = 17, λ = 9**. Un recocido sobre todos los patrones
simétricos con n ≤ 32 y todas las bases converge siempre a esta familia (`src/anneal_overlap.py`).

**Semilla** (`seeds/seed_n20_m17.json`): p de grado 20 con coeficientes reales, palindrómico
(a_{20−i} = a_i), p² con exactamente 21 términos, soporte [0,6] ∪ [17,23] ∪ [34,40]. Es una de las 42
soluciones reales densas que da la homotopía de grado total (`src/homotopy.py`).

### Teorema (demostrado salvo error; ver §3 para lo verificado por ordenador)

Existe p ∈ ℝ[x] denso de grado 20 con supp(p²) = [0,6] ∪ [17,23] ∪ [34,40]. Para todo K ≥ 1 existe un
polinomio real **completo** (denso) de grado D_K = 20(17^K − 1)/16 cuyo cuadrado tiene **como mucho
21·9^{K−1}** términos. En consecuencia, Q̂(k) ≪ k^{log 9/log 17} ≈ k^{0,77552} a lo largo de esa sucesión,
frente a k^{0,81071} de Verdenius (1949). La ventaja numérica empieza en K = 7 (≈ 5·10⁸ términos).

**Prueba.**
1. *Existencia de p:* certificado de Krawczyk (abajo). Los ceros de p² fuera de S son exactos por construcción.
2. *Cota superior:* supp(P_K²) ⊆ sumas de dígitos; contarlas da 7·9^{K−1}·3.
3. *Densidad para todo K (lema nuevo):* los a_i son algebraicos no nulos, así que existe un primo ℓ tal que
   todos son unidades en un primo 𝔭 | ℓ de su cuerpo de números. Reescalando p(x) → p(ℓx) (no cambia el
   soporte de p²), el coeficiente de x^N en P_K es Σ_reps (∏ a_{i_j}) ℓ^{Σ i_j}, sumando sobre las
   escrituras N = Σ i_j m^j con 0 ≤ i_j ≤ n. Cada término tiene valoración e·(suma de dígitos).
   **Lema:** si m − 1 ≤ n < 2m, la escritura de suma de dígitos mínima es única. Prueba: si r ≠ r' son
   mínimas y j es el primer dígito en que difieren, entonces i_j = i'_j + m. Si algún dígito de r por
   encima de j es < n, una cadena de acarreos reduce la suma de dígitos en un múltiplo de (m − 1),
   contradicción. Si todos son n, el valor de r' por encima de j excedería el máximo representable,
   contradicción. Por tanto hay un único término de valoración mínima y el coeficiente es ≠ 0. ∎

(La misma construcción con n = 14, m = 12, r = 4 da log 7/log 12 ≈ 0,78309; también certificada.)

**Advertencias:** (i) es una mejora para coeficientes **reales**; con coeficientes racionales o enteros
(el marco de Epoch y de Erdős #485 en ℚ[x]) **no hay mejora**. (ii) Es una cota a lo largo de una
sucesión, como la de Verdenius; la extensión a todo k usaría el argumento de interpolación de Erdős
(1949), que no he formalizado. (iii) No es revisión por pares: conviene que un especialista revise el
lema de densidad y el certificado.

## 3. Verificación (programas independientes del buscador)

* `certify/krawczyk.py`: aritmética de intervalos (mpmath.iv, 110 dígitos). Test de Krawczyk sobre el
  sistema simétrico: **solución exacta única** en una caja de radio 10⁻⁶⁰. Además comprueba que todos los
  a_i son ≠ 0 y que [x^j]p² ≠ 0 en S. Salidas: `results/krawczyk_*.json` (todo `true` para n = 14 y n = 20).
* `certify/tower_check.py`: no importa `src/`. Re-deriva la semilla, cuenta las sumas de dígitos por
  enumeración directa (189, 1701, …, 1 240 029 = 21·9^{K−1} hasta K = 6) y calcula P_K² en alta precisión:
  K = 2 denso, |P_2²| = 189 = cota. Ver `results/tower_check.txt`.
* `src/minds.py`: comprobación por fuerza bruta del lema de unicidad (K ≤ 4); el lema está demostrado arriba.
* `verify.py`: verificador exacto con enteros del criterio de Epoch (para entregas en formato JSON).

## 4. Intentos negativos (ver `HYPOTHESES.es.md`)

* Mezclar semillas distintas por nivel: solo promedia los exponentes.
* Cancelaciones extra de coeficientes en P_K²: no aparecen (el conteo exacto coincide con la cota).
* Sistemas sobredeterminados ("milagros" tipo Coppersmith–Davenport): ninguna solución no degenerada en
  los patrones que habrían dado entre 0,747 y 0,796.
* Semillas racionales: todas las soluciones aisladas (n ≤ 16, exponente < 0,816) calculadas por homotopía.
  La única racional es la de Verdenius n = 8. En la curva n = 18 (exponente 0,8114, que sí batiría el
  récord entero) un barrido en ρ racional no encuentra puntos racionales.

## 5. ¿Por qué no se llega a 1/2? (heurístico, no demostrado)

En esta familia, un patrón resoluble necesita |S| ≈ n+1 (conteo de dimensiones). Con S = A + m·E,
λ ≈ |A + E| ≥ |A| + |E| − 1 y m ≈ 2n/max E. Equilibrar esas restricciones deja el exponente cerca de 0,77,
lejos de 1/2. Bajar de 1/2 exigiría patrones con |S| ≪ n que sean resolubles: sistemas muy
sobredeterminados con soluciones, de los que no he encontrado ninguno. O bien una construcción distinta
de las torres de una sola semilla.

## Estructura

```
verify.py                 verificador entero independiente (formato Epoch)
certify/                  certificación rigurosa (Krawczyk) y comprobación independiente de torres
seeds/                    semillas (JSON); verdenius_integer_n8.json en formato Epoch
results/                  salidas de las verificaciones
src/                      código de búsqueda (autómatas, enumeraciones, homotopía, recocido)
HYPOTHESES.es.md             registro de hipótesis y resultados
```

Requisitos: Python 3, numpy, scipy, mpmath, python-flint.

Referencias: Rényi (1947); Erdős, Nieuw Arch. Wisk. 23 (1949); Verdenius, Indag. Math. 11 (1949);
Coppersmith–Davenport, Acta Arith. 58 (1991); Abbott, Math. Comp. 71 (2002); Schinzel–Zannier (2009);
Erdős Problem #485 (erdosproblems.com).
