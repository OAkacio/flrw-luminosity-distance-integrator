#
# * =============================================================================
# * DEPENDÊNCIAS
# * =============================================================================

# ? -----------------------------------------------------------------------------
# ?         BIBLIOTECAS
# ? -----------------------------------------------------------------------------

import numpy as np

# ? -----------------------------------------------------------------------------
# ?         MÓDULOS LOCAIS
# ? -----------------------------------------------------------------------------

from src.parameters import *
from src.core import *
from src.constants import *
from pytools import system as sy
from pytools import saveload as sl

# * =============================================================================
# * ROTINA PRINCIPAL
# * =============================================================================


def main(Omega_M, Omega_EE, w, z, type="return"):
    sy.header(
        "iniciando HoggCosmoMeasures...", Omega_M=Omega_M, Omega_EE=Omega_EE, w=w, z=z
    )

    # ? -----------------------------------------------------------------------------
    # ?         CÁLCULO DE PARÂMETROS PONTUAIS
    # ? -----------------------------------------------------------------------------

    try:
        sy.status("Iniciando processo de integração numérica para parâmetros pontuais")
        resintlist = integracao(integral, Omega_M, Omega_EE, z, w)
        sy.param(
            ("Integração Numérica", resintlist[0], "Mpc"),
            ("Erro Estimado", resintlist[1], "Mpc"),
            (
                "Distância de Luminosidade para (dL(z))",
                dL(Omega_M, Omega_EE, resintlist[0], z),
                "Mpc",
            ),
            (
                "Módulo de Distância (mu(z))",
                mu(Omega_M, Omega_EE, resintlist[0], z),
                "mag",
            ),
        )
    except Exception as e:
        sy.ok(f"Processo de integração numérica falhou! Erro: {e}...", False)
    try:
        sy.status("Iniciando cálculo de parâmetros do universo")
        sy.param(
            ("Tipo de universo", UniType(Omega_K(Omega_M, Omega_EE)), "adm."),
            ("Constante de curvatura espacial (k)", k(Omega_M, Omega_EE), "adm."),
            (
                "Parâmetro derivado de curvatura (Omega_K)",
                Omega_K(Omega_M, Omega_EE),
                "adm.",
            ),
            ("Distância comóvel radial (dC)", dC(resintlist[0]), "Mpc"),
            ("Parâmetro de desaceleração (q0)", q0(Omega_M, Omega_EE, w), "adm."),
        )
    except Exception as e:
        sy.ok(
            f"Processo de cálculo de parâmetros do universo falhou! Erro: {e}...", False
        )

    # ? -----------------------------------------------------------------------------
    # ?         CÁLCULO DE PARÂMETROS PARA TODO O INTERVALO
    # ? -----------------------------------------------------------------------------

    try:
        sy.status("Iniciando integração por todo o intervalo de redshift...")
        sollist = solution(Omega_M, Omega_EE, z, z_step, w)
        DLvectorX = sollist[0]
        DLvectorY = sollist[1]
        MUvectorX = sollist[2]
        MUvectorY = sollist[3]
        DLAPvectorX = sollist[4]
        DLAPvectorY = sollist[5]
        DIFvectorX = sollist[6]
        DIFvectorY = sollist[7]

        # ? -----------------------------------------------------------------------------
        # ?         EXPORTAÇÃO DE DADOS
        # ? -----------------------------------------------------------------------------

        if type == "custom":
            sy.status("Iniciando exportação de dados...")
            sl.savetable("infos", ((Omega_M, Omega_EE, w, z), ("", "", "", "")))
            sl.savetable("DLdados", (DLvectorX, DLvectorY))
            sl.savetable("MUdados", (MUvectorX, MUvectorY))
            sl.savetable("DLAPdados", (DLAPvectorX, DLAPvectorY))
            sl.savetable("DIFdados", (DIFvectorX, DIFvectorY))
            sy.ok(("infos", "DLdados", "MUdados", "DLAPdados", "DIFdados"))
        elif type == "M":
            sy.status("Iniciando exportação de dados")
            sy.status("Iniciando exportação de dados...")
            sl.savetable("infosM", ((Omega_M, Omega_EE, w, z), ("", "", "", "")))
            sl.savetable("DLdadosM", (DLvectorX, DLvectorY))
            sl.savetable("MUdadosM", (MUvectorX, MUvectorY))
            sl.savetable("DLAPdadosM", (DLAPvectorX, DLAPvectorY))
            sl.savetable("DIFdadosM", (DIFvectorX, DIFvectorY))
            sy.ok(("infosM", "DLdadosM", "MUdadosM", "DLAPdadosM", "DIFdadosM"))
        elif type == "EE":
            sy.status("Iniciando exportação de dados")
            sy.status("Iniciando exportação de dados...")
            sl.savetable("infosEE", ((Omega_M, Omega_EE, w, z), ("", "", "", "")))
            sl.savetable("DLdadosEE", (DLvectorX, DLvectorY))
            sl.savetable("MUdadosEE", (MUvectorX, MUvectorY))
            sl.savetable("DLAPdadosEE", (DLAPvectorX, DLAPvectorY))
            sl.savetable("DIFdadosEE", (DIFvectorX, DIFvectorY))
            sy.ok(("infosEE", "DLdadosEE", "MUdadosEE", "DLAPdadosEE", "DIFdadosEE"))
        elif type == "return":
            sy.status("EXECUÇÃO FINALIZADA!")
            return [
                dL(Omega_M, Omega_EE, resintlist[0], z),
                mu(Omega_M, Omega_EE, resintlist[0], z),
            ]
        sy.fim()
    except Exception as e:
        sy.ok(f"Falha no processo de salvamento! Erro: {e}...", False)


# ? -----------------------------------------------------------------------------
# ?         MAIN GUARD
# ? -----------------------------------------------------------------------------


if __name__ == "__main__":
    main(Omega_M, Omega_EE, w, z, type="custom")
