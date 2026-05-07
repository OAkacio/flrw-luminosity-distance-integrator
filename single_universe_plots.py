#
# * =============================================================================
# * DEPENDÊNCIAS
# * =============================================================================

# ? -----------------------------------------------------------------------------
# ?         BIBLIOTECAS
# ? -----------------------------------------------------------------------------

import numpy as np

#
# ? -----------------------------------------------------------------------------
# ?         MÓDULOS LOCAIS
# ? -----------------------------------------------------------------------------

from src.constants import c, H0
from src.parameters import *
from pytools import graphs as gp
from pytools import system as sy
from pytools import saveload as sl

# * =============================================================================
# * ROTINA PRINCIPAL
# * =============================================================================

sy.header("Análise de Universo Simulado", Folder="/data/...")

# ? -----------------------------------------------------------------------------
# ?         CARREGANDO DADOS DO UNIVERSO SIMULADO
# ? -----------------------------------------------------------------------------

try:
    sy.status("Iniciando carregamento de dados...")
    DLvectorX = sl.loadtable(f"data/DLdados.txt")[0]
    DLAPvectorX = sl.loadtable(f"data/DLAPdados.txt")[0]
    DIFvectorX = sl.loadtable(f"data/DIFdados.txt")[0]
    MUvectorX = sl.loadtable(f"data/MUdados.txt")[0]
    ANALMUX = sl.loadtable(f"data/ANALMU_Mdados.txt")[0]
    MagApX = sl.loadtable(f"data/MagAp_list.txt")[0]
    astropyX = sl.loadtable(f"data/astropy_MU.txt")[0]
    DLvectorY = sl.loadtable(f"data/DLdados.txt")[1]
    DLAPvectorY = sl.loadtable(f"data/DLAPdados.txt")[1]
    DIFvectorY = sl.loadtable(f"data/DIFdados.txt")[1]
    MUvectorY = sl.loadtable(f"data/MUdados.txt")[1]
    ANALMU_MY = sl.loadtable(f"data/ANALMU_Mdados.txt")[1]
    ANALMU_EEY = sl.loadtable(f"data/ANALMU_EEdados.txt")[1]
    ANALMU_VY = sl.loadtable(f"data/ANALMU_Vdados.txt")[1]
    MagApY = sl.loadtable(f"data/MagAp_list.txt")[1]
    astropy_DL = sl.loadtable(f"data/astropy_DL.txt")[1]
    astropy_MU = sl.loadtable(f"data/astropy_MU.txt")[1]
    infos = sl.loadtable(f"data/infos.txt")[0]
    sy.ok(
        (
            "infos",
            "DLdados",
            "MUdados",
            "DLAPdados",
            "DIFdados",
            "ANALMU_Mdados",
            "ANALMU_EEdados",
            "ANALMU_Vdados",
            "astropy_DL",
            "astropy_MU",
        )
    )
    sy.status("Dados carregados com sucesso!")
    sy.param(
        ("Omega_M", infos[0]), ("Omega_EE", infos[1]), ("w", infos[2]), ("z", infos[3])
    )

    # ? -----------------------------------------------------------------------------
    # ?         GERANDO GRÁFICOS DE PARÂMETROS
    # ? -----------------------------------------------------------------------------

    sy.status("Iniciando criação do gráficos...")
    gp.basic(
        x_data=DLvectorX,
        y_data=DLvectorY / (c / H0),
        title="",
        x_label=r"$z$",
        y_label=r"$d_L \; /\; \left(\frac{c}{H_0}\right)$",
        save_fig=True,
        filename="DLdistribuicao",
        show_plot=False,
    )
    gp.multi(
        x_list=(
            DLvectorX,
            DLAPvectorX,
        ),
        y_list=(DLvectorY / (c / H0), DLAPvectorY / (c / H0)),
        title="",
        x_label=r"$z$",
        y_label=r"$d_L \; /\; \left(\frac{c}{H_0}\right)$",
        curve_names=(
            rf"Curva Exata ( $\Omega_m={Omega_M}, \Omega_\Lambda={Omega_EE}$ )",
            "Curva Aproximada",
        ),
        save_fig=True,
        filename="ExataAproximadadistribuicao",
        show_plot=False,
    )
    gp.multi(
        x_list=(DIFvectorX, DIFvectorX),
        y_list=(DIFvectorY, np.full(len(DIFvectorY), 0)),
        title="",
        x_label=r"$z$",
        y_label=r"$d_L \; /\; \left(\frac{c}{H_0}\right)$",
        curve_names=(
            rf"Diferença Exato vs. Aproximado ( $\Omega_m={Omega_M}, \Omega_\Lambda={Omega_EE}$ )",
            "Zero",
        ),
        save_fig=True,
        filename="DIFdistribuicao",
        show_plot=False,
    )
    gp.basic(
        x_data=MUvectorX,
        y_data=MUvectorY,
        title="",
        x_label=r"$z$",
        y_label=r"$\mu$ (mag)",
        save_fig=True,
        filename="MUdistribuicao",
        show_plot=False,
    )
    gp.basic(
        x_data=MagApX,
        y_data=MagApY,
        title="",
        x_label=r"$z$",
        y_label=r"$m$ (mag)",
        save_fig=True,
        filename="MagApdistribuicao",
        show_plot=False,
    )
    gp.multi(
        x_list=[ANALMUX, MUvectorX],
        y_list=[ANALMU_MY, MUvectorY],
        title="",
        x_label=r"$z$",
        y_label=r"$\mu$ (mag)",
        curve_names=(
            "Curva Analítica - Universo de Somente Matéria",
            rf"Curva Numérica ( $\Omega_m={Omega_M}, \Omega_\Lambda={Omega_EE}$ )",
        ),
        save_fig=True,
        filename="AnalitVSNumer_M",
        show_plot=False,
    )
    gp.multi(
        x_list=[ANALMUX, MUvectorX],
        y_list=[ANALMU_EEY, MUvectorY],
        title="",
        x_label=r"$z$",
        y_label=r"$\mu$ (mag)",
        curve_names=(
            "Curva Analítica - Universo de Somente Energia",
            rf"Curva Numérica ( $\Omega_m={Omega_M}, \Omega_\Lambda={Omega_EE}$ )",
        ),
        save_fig=True,
        filename="AnalitVSNumer_EE",
        show_plot=False,
    )
    gp.multi(
        x_list=[ANALMUX, MUvectorX],
        y_list=[ANALMU_VY, MUvectorY],
        title="",
        x_label=r"$z$",
        y_label=r"$\mu$ (mag)",
        curve_names=(
            "Curva Analítica - Universo Vazio",
            rf"Curva Numérica ( $\Omega_m={Omega_M}, \Omega_\Lambda={Omega_EE}$ )",
        ),
        save_fig=True,
        filename="AnalitVSNumer_V",
        show_plot=False,
    )
    gp.multi(
        x_list=(
            DLvectorX,
            astropyX,
        ),
        y_list=(DLvectorY / (c / H0), astropy_DL / (c / H0)),
        title="",
        x_label=r"$z$",
        y_label=r"$d_L \; /\; \left(\frac{c}{H_0}\right)$",
        curve_names=(
            rf"Curva Exata ( $\Omega_m={Omega_M}, \Omega_\Lambda={Omega_EE}$ )",
            "Curva de Validação Astropy",
        ),
        save_fig=True,
        filename="validacaoDLastropy",
        show_plot=False,
    )
    gp.multi(
        x_list=(
            MUvectorX,
            astropyX,
        ),
        y_list=(MUvectorY, astropy_MU),
        title="",
        x_label=r"$z$",
        y_label=r"$\mu$",
        curve_names=(
            rf"Curva Exata ( $\Omega_m={Omega_M}, \Omega_\Lambda={Omega_EE}$ )",
            "Curva de Validação Astropy",
        ),
        save_fig=True,
        filename="validacaoMUastropy",
        show_plot=False,
    )
    sy.ok(
        (
            "Gráfico de Distância de Luminosidade",
            "Gráfico de Distância de Luminosidade Aproximada",
            "Gráfico de Anlálise de Erro Aproximação vs. Exato",
            "Gráfico de Módulo de Distância",
            "Gráfico de comparação de Universo Only-Matter vs. Simulado",
            "Gráfico de comparação de Universo Only-Energy vs. Simulado",
            "Gráfico de comparação de Universo Vazio vs. Simulado",
        )
    )
    sy.fim()
except Exception as e:
    sy.ok(f"Falha no processo de plotagem dos gráficos! Erro: {e}", False)
