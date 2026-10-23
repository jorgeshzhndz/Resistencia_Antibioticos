# EDA de Resistencia Antimicrobiana (Antimicrobial Resistance - AMR)

Análisis Exploratorio de Datos (EDA) eficiente sobre los patrones de resistencia bacteriana, combinando datos clínicos, demográficos, hospitalarios y socioeconómicos.

# Pregunta principal 
¿Qué factores clínicos, demográficos, hospitalarios y socioeconómicos están más fuertemente asociados a la aparición de resistencia antimicrobiana (AMR) en los pacientes?

---

## 📂 Dataset y ficheros utilizados 

El proyecto se divide en un **Análisis general** para un análisis rápido y completo de la base de datos, y **Análisis extra** para profundizar en determinantes socioeconómicos y exposición farmacológica previa. En análisis futuros se podría analizar la comorbilidad de los paciente y su relación con la aparición de resistencias a antibióticos. 

### 1. Análisis general _ 

*   **`cohort.csv`**: Define la cohorte base de pacientes/episodios de estudio.
*   **`microbial_resistance.csv`**: Contiene los aislamientos microbiológicos y sus perfiles de sensibilidad/resistencia a distintos antibióticos.
*   **`demographics.csv`**: Variables demográficas clave (edad, género, etc.).
*   **`ward_info.csv`**: Localización hospitalaria del paciente (unidades de hospitalización: IP, OP, ER, ICU).

### 2. Extras
*   **`adi_scores.csv` (48 MB)**: Añade el Índices de Deprivación del Área (**ADI**). Permite correlacionar la resistencia con factores socioeconómicos del paciente.
*   **`antibiotic_class_exposure.csv` (540 MB)**: Historial de exposición previa a clases de antibióticos. *Nota técnica: Requiere agregación previa por cultivo/paciente antes de hacer el `merge` con la tabla principal.*

##  Uso y ética

El dataset es público en Dryad con fines de investigación y está desidentificado (identificadores anonimizados, edad agrupada, fechas alteradas). El estudio original fue aprobado por el IRB de Stanford (eProtocol #70466). No se debe intentar reidentificar a ningún paciente.

##  Referencia

Nateghi Haredasht, F., Amrollahi, F., Maddali, M., Marshall, N., Ma, S., Chang, A., Banaei, N., Deresinski, S., Asch, S., Goldstein, M., & Chen, J. (2025). *Antibiotic Resistance Microbiology Dataset (ARMD): A resource for antimicrobial resistance from EHRs* [Dataset]. Dryad. https://doi.org/10.5061/dryad.jq2bvq8kp
