"""Funciones auxiliares del análisis de AMR (apartado 5 de memoria.ipynb).
 
Este módulo NO lee ficheros ni toca datos al importarse: solo define funciones.
Los datos se le pasan desde el notebook con set_data(dfc), o bien con el
argumento opcional `df` de cada función.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import binomtest, chi2_contingency
 
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
 
 
def plot_amr(var, df=None):
    df = _get(df)
    g = amr_por(var, df)
    fig, ax = plt.subplots(figsize=(7, 0.45 * len(g) + 1.5))
    err = [(g["pct_amr"] - g["ic_low"]).clip(lower=0),
           (g["ic_high"] - g["pct_amr"]).clip(lower=0)]
    ax.barh(g.index.astype(str), g["pct_amr"], xerr=err, capsize=3)
    ax.axvline(_y(df).mean() * 100, ls="--", color="gray")
    ax.set_xlabel("% de cultivos con AMR")
    ax.set_title(f"AMR según {var}")
    ax.invert_yaxis()
    plt.tight_layout()
    plt.show()
 
 
def a_codigo(serie):
    """Q1..Q5 -> 1..5; Missing -> NaN. Requiere una serie categórica ordenada."""
    if not isinstance(serie.dtype, pd.CategoricalDtype):
        raise TypeError("a_codigo necesita una serie categórica (pd.Categorical con orden Q1..Q5, Missing).")
    return (serie.cat.codes + 1).where(serie != "Missing")
 
 
def fmt(g, cat):
    return (f"{cat}: {g.loc[cat, 'pct_amr']} % "
            f"[{g.loc[cat, 'ic_low']}-{g.loc[cat, 'ic_high']}]")
