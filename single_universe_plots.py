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
    DLvectorY = sl.loadtable(f"data/DLdados.txt")[1]
    DLAPvectorY = sl.loadtable(f"data/DLAPdados.txt")[1]
    DIFvectorY = sl.loadtable(f"data/DIFdados.txt")[1]
    MUvectorY = sl.loadtable(f"data/MUdados.txt")[1]
    infos = sl.loadtable(f"data/infos.txt")[0]
    sy.ok(("infos", "DLdados", "MUdados", "DLAPdados", "DIFdados"))
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
        curve_names=("Curva Exata", "Curva Aproximada"),
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
        curve_names=("Diferença Exato vs. Aproximado", "Zero"),
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
    sy.ok(
        (
            "Gráfico de Distância de Luminosidade",
            "Gráfico de Distância de Luminosidade Aproximada",
            "Gráfico de Anlálise de Erro Aproximação vs. Exato",
            "Gráfico de Módulo de Distância",
        )
    )
    sy.fim()
except Exception as e:
    sy.ok(f"Falha no processo de plotagem dos gráficos! Erro: {e}", False)
