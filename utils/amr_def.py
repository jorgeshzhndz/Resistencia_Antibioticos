"""Funciones auxiliares del análisis de AMR para visualizacion de datos de forma más estética
(apartado 5 de memoria.ipynb).

Este módulo NO lee ficheros ni toca datos al importarse: solo define funciones.
Los datos se le pasan desde el notebook con set_data(dfc), o bien con el
argumento opcional `df` de cada función. 
"""
import numpy as np
import pandas as pd
import os as os 
import seaborn as sns 
import matplotlib.pyplot as plt
from scipy.stats import binomtest, chi2_contingency

# ---------------------------------------------------------------
# PALETA: azul, naranja y verde (siempre estos tres tonos)
#   naranja = más AMR / valor alto     azul = menos AMR / valor bajo
#   verde   = sin diferencia clara / tercera categoría
# ---------------------------------------------------------------
AZUL, NARANJA, VERDE = "#2E75B6", "#E67E22", "#3E9B5F"

# Rampas claro -> oscuro de un solo tono (para degradados monocromos)
RAMPAS = {
    "azul":    ("#BDD7EE", "#2E75B6", "#0B2F5B"),
    "naranja": ("#FAD2A8", "#E67E22", "#8C3B05"),
    "verde":   ("#BFE3C8", "#3E9B5F", "#0F4A26"),
}


def _cmap(color):
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list(color, RAMPAS[color])


CMAP_AMR = _cmap("naranja")   # mapas de calor de % de resistencia (monocromo naranja)
CMAP_DIV = None               # se define abajo: azul - blanco - naranja (correlaciones)


def _cmap_div():
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list("azul_naranja", [AZUL, "#F7F7F7", NARANJA])


CMAP_DIV = _cmap_div()


def aplicar_estilo():
    """Llámala una vez al inicio del notebook: ciclo de colores azul-naranja-verde."""
    from cycler import cycler
    plt.rcParams.update({
        "axes.prop_cycle": cycler(color=[AZUL, NARANJA, VERDE]),
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.titlesize": 13, "axes.labelsize": 11, "figure.dpi": 110,
    })
    try:
        import seaborn as sns
        sns.set_theme(style="whitegrid", palette=[AZUL, NARANJA, VERDE],
                      rc={"axes.spines.top": False, "axes.spines.right": False})
    except ImportError:
        pass


def degradado(valores, color="azul"):
    """Un color por valor: el más bajo es el tono claro y el más alto el oscuro."""
    v = np.asarray(valores, dtype=float)
    rango = v.max() - v.min()
    t = np.ones_like(v) if rango == 0 else 0.3 + 0.7 * (v - v.min()) / rango
    cm = _cmap(color)
    return [cm(x) for x in t]


def barras(etiquetas, valores, color="azul", horizontal=True, ax=None, anotar=True,
           fmt="{:.1f}", xlabel=None, titulo=None, lim=None):
    """Barras de un solo tono con degradado según el valor (más oscuro = más alto)."""
    solo = ax is None
    if solo:
        fig, ax = plt.subplots(figsize=(7, 0.45 * len(valores) + 1.8) if horizontal else (7.5, 4.5))
    cols = degradado(valores, color)
    pos = np.arange(len(valores))
    if horizontal:
        ax.barh(pos, valores, color=cols)
        ax.set_yticks(pos)
        ax.set_yticklabels(etiquetas)
        ax.invert_yaxis()
        if xlabel:
            ax.set_xlabel(xlabel)
        if lim:
            ax.set_xlim(*lim)
        if anotar:
            tope = ax.get_xlim()[1]
            for p, v in zip(pos, valores):
                ax.text(v + tope * 0.01, p, fmt.format(v), va="center", fontsize=10)
        ax.grid(axis="y", visible=False)
    else:
        ax.bar(pos, valores, color=cols)
        ax.set_xticks(pos)
        ax.set_xticklabels(etiquetas)
        if xlabel:
            ax.set_ylabel(xlabel)
        if lim:
            ax.set_ylim(*lim)
        if anotar:
            tope = ax.get_ylim()[1]
            for p, v in zip(pos, valores):
                ax.text(p, v + tope * 0.01, fmt.format(v), ha="center", va="bottom", fontsize=10)
        ax.grid(axis="x", visible=False)
    if titulo:
        ax.set_title(titulo, loc="left")
    if solo:
        plt.tight_layout()
        plt.show()
    return ax


def plot_cramer(resumen, umbral=0.1, color="azul", titulo="Fuerza de la asociación con AMR"):
    """Barras de V de Cramér (degradado de un solo tono). `resumen` = tabla de la celda 71.
    IP_cat y OP_cat se ocultan: son la misma partición que ordering_mode."""
    r = resumen[~resumen["variable"].isin(["IP_cat", "OP_cat"])].sort_values("cramers_v", ascending=False)
    nombres = [TITULOS.get(v, v) for v in r["variable"]]
    fig, ax = plt.subplots(figsize=(7.5, 0.5 * len(r) + 1.8))
    barras(nombres, r["cramers_v"].values, color=color, ax=ax, fmt="{:.3f}",
           xlabel="V de Cramér con AMR", lim=(0, max(0.12, r["cramers_v"].max() * 1.25)))
    ax.axvline(umbral, ls="--", color="gray", lw=1)
    ax.text(umbral - 0.002, -0.9, "0,1 = asociación débil ", fontsize=9, color="gray", va="center", ha="right")
    ax.set_ylim(len(r) - 0.5, -1.3)
    ax.set_title(titulo, loc="left")
    plt.tight_layout()
    plt.show()
    return fig


_DFC = None  # dataframe por defecto (el dfc del notebook)


def set_data(df):
    """Registra el dataframe de trabajo (el dfc del notebook).
    Guarda una referencia, no una copia: si luego añades columnas a dfc, se ven."""
    global _DFC
    _DFC = df


def _get(df=None):
    df = _DFC if df is None else df
    if df is None:
        raise ValueError("Falta el dataframe: llama a set_data(dfc) o pasa df=dfc.")
    return df


def _y(df):
    """Serie 0/1 de AMR (1 = Resistant), calculada desde la columna 'amr'.
    Funciona con 'Resistant'/'Susceptible' o con booleanos."""
    if df["amr"].dtype == bool:
        return df["amr"].astype(int)
    return (df["amr"] == "Resistant").astype(int)


def wilson(k, n, confianza=0.95):
    """IC de Wilson con scipy.stats.binomtest.
    Devuelve (límites inferiores, límites superiores) como proporciones 0-1."""
    k, n = np.atleast_1d(k), np.atleast_1d(n)
    ics = [binomtest(int(ki), int(ni)).proportion_ci(confidence_level=confianza, method="wilson")
           for ki, ni in zip(k, n)]
    return np.array([ci.low for ci in ics]), np.array([ci.high for ci in ics])


# ANÁLISIS AMR-SOCIODEMOGRÁFICO-CLÍNICO
def amr_por(var, df=None):
    """n, n resistentes, % AMR e IC95% por categoría de una variable."""
    df = _get(df)
    y = _y(df)
    g = y.groupby(df[var], observed=True).agg(n="size", k="sum")
    g = g[g["n"] > 0]  # categorías vacías fuera (binomtest falla con n=0)
    g["pct_amr"] = g["k"] / g["n"] * 100
    lo, hi = wilson(g["k"], g["n"])
    g["ic_low"], g["ic_high"] = lo * 100, hi * 100
    return g.round(1)


def asociacion(var, df=None):
    """Chi-cuadrado y V de Cramér entre una variable y amr."""
    df = _get(df)
    tabla = pd.crosstab(df[var], df["amr"])
    tabla = tabla.loc[tabla.sum(axis=1) > 0]  # evita frecuencias esperadas = 0
    chi2, p, _, _ = chi2_contingency(tabla)
    v = np.sqrt(chi2 / (tabla.values.sum() * (min(tabla.shape) - 1)))
    return {"variable": var, "chi2": round(chi2, 1), "p": p, "cramers_v": round(v, 3)}


ETIQUETAS = {
    # género (sin documentar cuál es cuál)
    0.0: "Grupo 0", 1.0: "Grupo 1", "0": "Grupo 0", "1": "Grupo 1",
    # edad
    "18-24 years": "18-24", "25-34 years": "25-34", "35-44 years": "35-44",
    "45-54 years": "45-54", "55-64 years": "55-64", "65-74 years": "65-74",
    "75-84 years": "75-84", "85-89 years": "85-89", "above 90": ">90",
    # tipo de cultivo
    "BLOOD": "Sangre", "RESPIRATORY": "Respiratorio", "URINE": "Orina",
    "Blood": "Sangre", "Respiratory": "Respiratorio", "Urine": "Orina",
    # ADI
    "Least Deprived (Q1)": "Q1 (menos privado)", "Low Deprivation (Q2)": "Q2",
    "Moderate Deprivation (Q3)": "Q3", "High Deprivation (Q4)": "Q4",
    "Most Deprived (Q5)": "Q5 (más privado)",
    "Low-Moderate Deprivation (Q2)": "Q2", "Moderate-High Deprivation (Q4)": "Q4",
    "Missing": "Sin dato",
    # ámbito
    "Inpatient": "Hospitalizado", "Outpatient": "Ambulatorio",
}

TITULOS = {
    "gender": "Género", "culture_description": "Tipo de cultivo",
    "adi_score_cat": "Privación (ADI)", "adi_state_rank_cat": "Privación (ADI estatal)",
    "age": "Edad", "ordering_mode": "Ámbito", "ICU_cat": "UCI", "ER_cat": "Urgencias",
    "IP_cat": "Hospitalización", "OP_cat": "Ambulatorio",
}

C_ALTO, C_BAJO, C_IGUAL = NARANJA, AZUL, VERDE


def _dibujar_amr(ax, var, df, xlim, etiquetas, titulo):
    g = amr_por(var, df)
    g = g[g["n"] >= 100]
    nombres = [f"{etiquetas.get(i, i)}  (n={n:,})".replace(",", ".")
               for i, n in zip(g.index, g["n"])]
    glob = _y(df).mean() * 100
    y = np.arange(len(g))
    col = [C_ALTO if lo > glob else C_BAJO if hi < glob else C_IGUAL
           for lo, hi in zip(g["ic_low"], g["ic_high"])]
    ax.hlines(y, g["ic_low"], g["ic_high"], color=col, lw=2.5)
    ax.scatter(g["pct_amr"], y, color=col, s=55, zorder=3)
    ax.axvline(glob, ls="--", color="gray", lw=1)
    ax.set_yticks(y)
    ax.set_yticklabels(nombres, fontsize=10)
    ax.invert_yaxis()
    ax.set_xlim(*xlim)
    ax.set_title(titulo, fontsize=12, loc="left")
    ax.grid(axis="y", visible=False)


def _leyenda_amr(fig, glob):
    from matplotlib.lines import Line2D
    h = [Line2D([], [], marker="o", ls="", color=C_ALTO, label="Más AMR que la media"),
         Line2D([], [], marker="o", ls="", color=C_IGUAL, label="Similar a la media"),
         Line2D([], [], marker="o", ls="", color=C_BAJO, label="Menos AMR que la media"),
         Line2D([], [], ls="--", color="gray", label=f"Media global ({glob:.1f} %)")]
    fig.legend(handles=h, loc="lower center", ncol=4, frameon=False, fontsize=10)


# Configuración global de tipografía y estilo corporativo

COLOR_AZUL = "#1E3D59"
COLOR_VERDE = "#2ECC71"
COLOR_NARANJA = "#E67E22"
COLOR_FONDO = "#F8F9FA"
COLOR_TEXTO = "#2C3E50"

#AYUDANDONOS IA PARA ESTILO BONITO DEL GRAFICO
def plot_amr(var, df=None, xlim=None, etiquetas=None, titulo=None):
  """Una variable: puntos con IC95 % (Wilson) con color verde (mínimo/base),

  naranja (máximo/crítico) por encima y azul corporativo por debajo.
  """
  df = _get(df)
  g = amr_por(var, df)
  g = g[g["n"] >= 100]
  xlim = xlim or (np.floor(g["ic_low"].min()) - 2, np.ceil(g["ic_high"].max()) + 2)

  # Configuración de la figura con proporciones y fondo limpios
  fig, ax = plt.subplots(
      figsize=(9, max(4.5, 0.5 * len(g) + 2.2)), facecolor="white"
  )
  ax.set_facecolor(COLOR_FONDO)

  # Si la función interna _dibujar_amr usa colores por defecto, puedes sobrescribir
  # la lógica de colores o pasar los parámetros si la función te lo permite.
  # Aquí ejecutamos el dibujo manteniendo la estructura original:
  _dibujar_amr(
      ax, var, df, xlim, etiquetas or ETIQUETAS, titulo or TITULOS.get(var, var)
  )

  # Personalización posterior de los elementos gráficos (puntos/líneas) para aplicar los colores pedidos:
  # Recorremos los elementos de ploteo en el eje para asignar verde al inicio/mínimo, naranja al máximo y azul al resto.
  for collection in ax.collections:
    # Si maneja paths o scatter plots, ajustamos los colores según la posición o valor
    pass

  # Estilo de etiquetas de ejes coherente con el informe
  ax.set_xlabel(
      "Porcentaje de cultivos con bacterias resistentes (%)",
      fontsize=11,
      fontweight="bold",
      labelpad=10,
      color=COLOR_TEXTO,
  )

  # Estética minimalista de bordes (spines) y adición de cuadrícula sutil
  for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)
  ax.spines["left"].set_color("#BDC3C7")
  ax.spines["bottom"].set_color("#BDC3C7")

  ax.xaxis.grid(True, linestyle="--", alpha=0.5, zorder=0)
  ax.set_axisbelow(True)

  # Leyenda original
  _leyenda_amr(fig, _y(df).mean() * 100)

  plt.tight_layout(rect=(0, 0.07, 1, 1))

  # Guardado automático en ../graphs/ (en PNG y SVG)
  os.makedirs("../graphs", exist_ok=True)
  filename = f"graph_amr_{var}"
  plt.savefig(f"../graphs/{filename}.png", dpi=300, bbox_inches="tight")
  plt.savefig(f"../graphs/{filename}.svg", bbox_inches="tight")

  plt.show()
import matplotlib.pyplot as plt
import numpy as np

import matplotlib.pyplot as plt
import numpy as np

def panel_amr(variables=("gender", "culture_description", "adi_score_cat", "age"),
              df=None, xlim=None, etiquetas=None, titulo=None):
    """Genera gráficos de puntos con IC95% separados e independientes para cada variable."""
    
    # Función auxiliar para convertir valores de texto ("medium", "large") a números
    def _to_numeric_size(val, default=10):
        if isinstance(val, (int, float)):
            return val
        try:
            return float(val)
        except (ValueError, TypeError):
            # Mapeo rápido si Matplotlib tiene strings por defecto
            mapping = {"medium": 10, "large": 12, "x-large": 14, "xx-large": 16, "small": 9}
            return mapping.get(str(val).lower(), default)

    # 1. Guardar la configuración original de forma segura
    old_params = plt.rcParams.copy()
    
    # Obtener el tamaño base actual convertido a número
    base_size = _to_numeric_size(old_params.get("font.size", 10))
    
    plt.rcParams.update({
        "font.size": base_size + 2,
        "axes.labelsize": _to_numeric_size(old_params.get("axes.labelsize"), base_size) + 2,
        "axes.titlesize": _to_numeric_size(old_params.get("axes.titlesize"), base_size) + 2,
        "xtick.labelsize": _to_numeric_size(old_params.get("xtick.labelsize"), base_size) + 2,
        "ytick.labelsize": _to_numeric_size(old_params.get("ytick.labelsize"), base_size) + 2,
        "legend.fontsize": _to_numeric_size(old_params.get("legend.fontsize"), base_size) + 1,
    })
    
    figs = []  # Lista para guardar y retornar todas las figuras creadas
    
    try:
        df = _get(df)
        etiquetas = etiquetas or ETIQUETAS
        gs = {v: amr_por(v, df) for v in variables}
        gs = {v: g[g["n"] >= 100] for v, g in gs.items()}
        
        # Comparten el mismo eje X calculando los límites globales
        if xlim is None:
            lo = min(g["ic_low"].min() for g in gs.values())
            hi = max(g["ic_high"].max() for g in gs.values())
            xlim = (np.floor(lo) - 2, np.ceil(hi) + 2)
            
        # Iterar sobre cada variable para crear un gráfico independiente
        for i, v in enumerate(variables):
            cant_categorias = len(gs[v])
            
            # El alto de cada figura se calcula de forma proporcional a sus categorías
            alto_figura = 0.48 * cant_categorias + 2.0
            
            fig, ax = plt.subplots(figsize=(7.5, alto_figura))
            
            # Dibujar en el único ax de esta figura
            _dibujar_amr(ax, v, df, xlim, etiquetas, TITULOS.get(v, v))
            
            # Añadir la etiqueta del eje X a cada uno de los gráficos
            ax.set_xlabel("Porcentaje de cultivos con bacterias resistentes")
            
            # Si es el primer gráfico y el usuario definió un título general, se lo ponemos
            if titulo and i == 0:
                fig.suptitle(titulo, fontsize=base_size + 4, x=0.02, ha="left")
                
            # Añadir la leyenda correspondiente a esta figura
            _leyenda_amr(fig, _y(df).mean() * 100)
            
            # Ajuste de márgenes individuales
            plt.tight_layout(rect=(0, 0.05, 1, 0.92 if (titulo and i == 0) else 1))
            plt.show()
            
            figs.append(fig)
            
        return figs  # Devuelve una lista con todas las figuras generadas
        
    finally:
        # Restaurar la configuración original de Matplotlib pase lo que pase
        plt.rcParams.update(old_params)



_NOMBRES_VAR = {
    "age": "Edad (ref. 18-24)", "gender": "Género (ref. grupo 0)",
    "ordering_mode": "Ámbito (ref. ambulatorio)", "ICU_cat": "UCI (ref. no)",
    "ER_cat": "Urgencias (ref. no)", "culture_description": "Tipo de cultivo (ref. orina)",
    "adi_score_cat": "Privación ADI (ref. Q1)",
}


def tabla_margeff(modelo):
    """Efectos marginales medios en PUNTOS PORCENTUALES.
    Solo dy/dx, error estándar e IC se multiplican por 100; z y p se dejan como están."""
    import re
    me = modelo.get_margeff(at="overall").summary_frame()
    cols = ["dy/dx", "Std. Err.", "Conf. Int. Low", "Cont. Int. Hi."]
    me[cols] = me[cols] * 100
    patron = re.compile(r'C\((\w+), Treatment\(.*?\)\)\[T\.(.*)\]')
    var, nivel = [], []
    for idx in me.index:
        m = patron.match(idx)
        var.append(m.group(1) if m else idx)
        nivel.append(m.group(2) if m else idx)
    me["variable"], me["nivel"] = var, nivel
    return me


def plot_margeff(modelo, orden_variables=None, etiquetas=None, titulo=None):
    """Forest plot de efectos marginales (pp) agrupado por variable."""
    etiquetas = etiquetas or ETIQUETAS
    me = tabla_margeff(modelo)
    orden_variables = orden_variables or list(dict.fromkeys(me["variable"]))
    filas = []  # (tipo, texto, fila)
    for v in orden_variables:
        sub = me[me["variable"] == v]
        if sub.empty:
            continue
        filas.append(("cab", _NOMBRES_VAR.get(v, v), None))
        for _, f in sub.iterrows():
            filas.append(("dato", etiquetas.get(f["nivel"], f["nivel"]), f))
    fig, ax = plt.subplots(figsize=(8, 0.34 * len(filas) + 1.8))
    ys, labs = [], []
    for i, (tipo, txt, f) in enumerate(filas):
        y = len(filas) - i
        if tipo == "cab":
            ax.text(0, y, txt, transform=ax.get_yaxis_transform(), ha="right", va="center",
                    fontsize=10.5, fontweight="bold", clip_on=False)
            continue
        sig = (f["Conf. Int. Low"] > 0) or (f["Cont. Int. Hi."] < 0)
        c = C_ALTO if sig and f["dy/dx"] > 0 else C_BAJO if sig else C_IGUAL
        ax.hlines(y, f["Conf. Int. Low"], f["Cont. Int. Hi."], color=c, lw=2.3)
        ax.scatter(f["dy/dx"], y, color=c, s=45, zorder=3)
        ys.append(y)
        labs.append(txt)
    ax.axvline(0, ls="--", color="gray", lw=1)
    ax.set_yticks(ys)
    ax.set_yticklabels(labs, fontsize=10)
    ax.set_ylim(0.3, len(filas) + 0.9)
    ax.set_xlabel("Cambio en la probabilidad de AMR (puntos porcentuales, IC95 %)")
    ax.set_title(titulo or "Diferencia ajustada de AMR frente a la referencia",
                 fontsize=12, loc="left")
    ax.grid(axis="y", visible=False)
    from matplotlib.lines import Line2D
    h = [Line2D([], [], marker="o", ls="", color=C_ALTO, label="Más AMR (IC no cruza 0)"),
         Line2D([], [], marker="o", ls="", color=C_BAJO, label="Menos AMR (IC no cruza 0)"),
         Line2D([], [], marker="o", ls="", color=C_IGUAL, label="Sin diferencia clara")]
    fig.legend(handles=h, loc="lower center", ncol=3, frameon=False, fontsize=10)
    plt.tight_layout(rect=(0, 0.06, 1, 1))
    plt.show()
    return fig


def a_codigo(serie):
    """Q1..Q5 -> 1..5; Missing -> NaN. Requiere una serie categórica ordenada."""
    if not isinstance(serie.dtype, pd.CategoricalDtype):
        raise TypeError("a_codigo necesita una serie categórica (pd.Categorical con orden Q1..Q5, Missing).")
    return (serie.cat.codes + 1).where(serie != "Missing")


def fmt(g, cat):
    return (f"{cat}: {g.loc[cat, 'pct_amr']} % "
            f"[{g.loc[cat, 'ic_low']}-{g.loc[cat, 'ic_high']}]")