# EDA de resistencias antimicrobianas (AMR) de la historia clínica de Stanford Health Care

Análisis exploratorio y estadístico de qué factores sociodemográficos, clínicos y hospitalarios se asocian a la resistencia antimicrobiana en cultivos microbiológicos positivos, y de qué bacterias y antibióticos concentran más resistencia, utilizando una población adulta de Stanford. 

La resistencia a los antimicrobianos (RAM) constituye uno de los principales retos actuales de salud pública, debido a su creciente impacto sobre la eficacia de los tratamientos y el manejo de las enfermedades infecciosas. Su aparición y propagación están relacionadas con múltiples factores, entre ellos el uso de antimicrobianos y las características clínicas y epidemiológicas de los pacientes. Por ello, resulta fundamental identificar los grupos de población y los contextos asistenciales en los que se observa una mayor prevalencia de microorganismos no sensibles. Este análisis pretende estudiar la relación entre los factores demográficos, socioeconómicos y clínicos y la presencia de resistencia antimicrobiana. Los resultados pueden contribuir a comprender mejor su distribución epidemiológica e identificar posibles grupos prioritarios para las estrategias de vigilancia y prevención.

**Autoría:** _Jorge Sánchez Hernández_ · **Notebook principal:** `notebooks/memoria.ipynb`

## Preguntas de investigación

1. ¿Qué factores clínicos, demográficos, hospitalarios y socioeconómicos se asocian a la AMR en los cultivos positivos?
2. ¿Qué bacterias son las más resistentes?
3. ¿A qué antibióticos?

## Datos

Se usa **ARMD (Antibiotic Resistance Microbiology Dataset)**, datos desidentificados de la historia clínica
electrónica de Stanford Health Care (pacientes adultos, 1999-2024). Son públicos:

- Descarga: <https://datadryad.org/dataset/doi:10.5061/dryad.jq2bvq8kp>
- Artículo: Nateghi Haredasht F. et al., *Antibiotic Resistance Microbiology Dataset (ARMD): A Resource for Antimicrobial
  Resistance from EHRs*, Scientific Data. Preprint: <https://arxiv.org/abs/2503.07664>

**Los datos no están en este repositorio** (debido a su gran tamaño). Descarga las cinco tablas siguientes del enlace de descarga y colócalas en `data/`
**con estos nombres** (en Dryad llevan el prefijo `microbiology_cultures_`; quítalo o cambia las rutas del notebook):

| Fichero esperado | Contenido que se usa |
|---|---|
| `data/cohort.csv` | Cultivo, ámbito, tipo de cultivo, positividad, organismo, antibiótico, sensibilidad |
| `data/demographics.csv` | Edad (grupos) y género (codificado 0/1) |
| `data/microbial_resistance.csv` | Resistencias previas y días hasta el cultivo |
| `data/ward_info.csv` | Indicadores de ámbito: IP, OP, ER, ICU |
| `data/adi_scores.csv` | Area Deprivation Index (score y rango estatal de exclusión social) |

## Estructura del repositorio

```
.
├── README.md
├── requirements.txt
├── notebooks/
│   └── memoria.ipynb          # análisis completo (secciones 1-5)
├── utils/
│   └── amr_def.py             # funciones y estilo de gráficos (el notebook hace `from utils.amr_def import ...`)
├── data/         # CSV de ARMD + ficheros intermedios (no se versionan)
├── graphs/       #Gráficos
│    └── grafico 1
│    └── grafico 2...
├── venv/ 
└── img/
    └── diagrama_flujo_cohorte_amr.svg
```

El notebook añade la carpeta padre a `sys.path` y lee `../data` y `../img`, así que **hay que ejecutarlo desde `notebooks/`**.

## Cómo ejecutarlo

```bash
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt
jupyter lab notebooks/memoria.ipynb
```

Ejecuta las celdas en orden (Kernel → Restart & Run All). Las secciones 2 a 4 generan ficheros intermedios en `data/`
que usan las secciones siguientes:

| Fichero generado | Qué contiene |
|---|---|
| `cohort_amr_summary.csv` | Una fila por cultivo con AMR, nº de antibióticos, de resistencias y de organismos |
| `cohort_preanalysis.csv` | Cohorte unida a edad, género, ámbito y ADI (con y sin categorizar) |
| `clinical_sociodemographic_resistance.csv` | Tabla de análisis por cultivo (sección 5.1-5.2) |
| `microorganism_antibiotic_resistance.csv` | Antibiogramas por cultivo, organismo y antibiótico (sección 5.3) |
| `cohorte_combinaciones_resistentes.csv` | Nº de combinaciones resistentes por paciente |

## Introducción y carga de datos.

Las secciones 1 y 2 introducen el por qué de nuestro análisis y se procede a la carga de los datos.

## Estructura del análisis

| Sección | Contenido | Unidad de análisis |
|---|---|---|
| 3 | Calidad de datos: nulos, duplicados, IDs huérfanos, consistencia | Tablas originales |
| 4 | Cohorte de estudio, cruce de tablas, tratamiento de ADI ausente y quintiles | Cultivo |
| 5.1 | ¿Cuántas resistencias acumula cada paciente? Tests no paramétricos por edad, género y ADI | **Paciente** |
| 5.2 | Factores asociados a AMR: % con IC95 % de Wilson, V de Cramér, correlaciones | **Cultivo** |
| 5.3 | Bacterias y antibióticos más resistentes, heatmap organismo × antibiótico | Cultivo / paciente × organismo × antibiótico |
| 5.4 | Regresión logística multivariable con errores estándar agrupados por paciente | Cultivo |
| 5.5 | Conclusiones y limitaciones | — |

### Definiciones clave

- **Cohorte:** cultivos con `was_positive == 1` y al menos un resultado Susceptible, Intermediate o Resistant.
- **AMR (sección 5.2 y 5.4):** el cultivo tiene al menos un antibiótico con resultado *Resistant*. *Intermediate* cuenta como no resistente.
- **Combinaciones resistentes (sección 5.1):** combinaciones únicas paciente-bacteria-antibiótico con resultado
  *Resistant* . 
- **Ámbitos (IP, OP, ER, ICU):** no son excluyentes; se mantienen como indicadores 0/1.
- **ADI:** quintiles calculados sobre la muestra de estudio; los valores ausentes (~20 %) se conservan como categoría `Missing`.

## Resultados principales


- **Cohorte final:** 113.991 cultivos de 64.648 pacientes (de 751.075 cultivos y 283.715 pacientes en ARMD).
- **Prevalencia de AMR:** 67,6 % de los cultivos tiene al menos una resistencia.
- **Acumulación por paciente:** mediana de 2 combinaciones resistentes (media 3,24, máximo 121). El 27,7 % de los pacientes
  no tiene ninguna y el 7,1 % acumula 10 o más. Las diferencias por edad, género y ADI son significativas pero pequeñas (ε² ≤ 0,013).
- **Factores asociados a AMR:** todas las asociaciones son débiles (V de Cramér ≤ 0,09). Las mayores diferencias son de 8-9 puntos
  porcentuales: género (74,1 % frente a 65,2 % entre los grupos codificados), tipo de cultivo (respiratorio 74,6 % frente a orina 66,2 %),
  ADI (quintil más desfavorecido 73,1 % frente a 65,1 %) y edad (62,1 % a los 18-24 años frente a ≈ 69-70 % desde los 55).
- **Modelo multivariable:** AUC 0,58 y pseudo-R² de McFadden 0,012. Las variables disponibles explican muy poco de la AMR. El modelo de regresión logística presentó un pseudo-R² de McFadden de 0,0123, lo que indica una mejora limitada del ajuste respecto al modelo nulo. Este resultado sugiere que las variables sociodemográficas y clínicas incluidas aportan una capacidad explicativa conjunta reducida para la presencia de resistencia antimicrobiana.
- **Antibióticos de uso habitual (resistencia por paciente-organismo-antibiótico):** Conviene dejar claro que son porcentajes de resultados resistentes entre las combinaciones evaluadas, no prevalencias generales de resistencia en la población:  trimetoprim/sulfametoxazol 23,4 %,
  levofloxacino 17,3 %, ciprofloxacino 17,2 %, ceftriaxona 11,0 %, meropenem 1,4 % y ertapenem 0,5 %.

## Limitaciones

Detalladas en la sección 5.5 del notebook. Las más importantes:

- Diseño observacional: se describen asociaciones, no causas.
- Un paciente aporta varios cultivos: los chi-cuadrado y los IC de Wilson asumen independencia y son optimistas;
  la regresión corrige esto con errores agrupados por paciente.
- La AMR se define como "al menos una resistencia entre los antibióticos probados", que depende de cuántos y cuáles se prueben
  y que incluye resistencia intrínseca (p. ej., ampicilina en muchas enterobacterias).
- No se incluyen comorbilidades ni exposición antibiótica previa.
- El género viene codificado 0/1 sin indicar cuál es cuál; el ADI falta en ~20 % y no al azar.
- Un solo centro (Stanford): puede no generalizarse a otros sistemas sanitarios. Distintas zonas básicas de salud tienen diferentes poblaciones y por tanto diferentes perfiles de resistencias a antimicrobianos 