#
# * =============================================================================
# * DEPENDÊNCIAS
# * =============================================================================

# ? -----------------------------------------------------------------------------
# ?         MÓDULOS LOCAIS
# ? -----------------------------------------------------------------------------

from src.core import *
from main import *
from src.parameters import *
from pytools import graphs as gp
from pytools import system as sy
from pytools import saveload as sl

# ? -----------------------------------------------------------------------------
# ?         BIBLIOTECAS
# ? -----------------------------------------------------------------------------

import numpy as np

# * =============================================================================
# * ROTINA PRINCIPAL
# * =============================================================================

sy.header("Inferência de dados", Dados="obs_data.txt")

sy.status("Iniciando carregamento de dados...")
try:
    obsdatalist = sl.loadtable(
        "obs_data.txt", separator=sl.identificadora("obs_data.txt")
    )
    z_list = np.array(obsdatalist[0])
    mu_obs_list = np.array(obsdatalist[1])
    ERROmu_obs_list = np.array(obsdatalist[2])
    sy.ok(("z_list", "mu_obs_list", "ERROmu_obs_list"))
except Exception as e:
    sy.ok(f"Um erro foi encontrado ao tentar ler os dados observados. Erro: {e}", False)

# ? -----------------------------------------------------------------------------
# ?         VARREDURA UNIDIMENSIONAL
# ? -----------------------------------------------------------------------------

try:
    sy.status("Iniciando a varredura 1D...")
    omega_list = []
    omega_var = mesh_inter_omega[0]
    while omega_var < mesh_inter_omega[1]:
        omega_list.append(omega_var)
        omega_var = omega_var + meshgrid_step
    var1d = varredura_1D(omega_list, mu_obs_list, ERROmu_obs_list, z_list, w)
    chi_min = np.min(var1d[0])
    INDchi_min = np.argmin(var1d[0])
    sy.status("Varredura 1D concluida!")
except Exception as e:
    sy.ok(f"Um erro foi encontrado ao tentar realizar a varredura 1D. Erro: {e}", False)
try:
    sy.param(
        ("Qui-quadrado mínimo (chi2_min)", chi_min, "adm."),
        ("Densidade de materia Bestfit (BF_Omega_M)", var1d[1][INDchi_min], "adm."),
        ("Densidade de energia Bestfit (BF_Omega_EE)", var1d[2][INDchi_min], "adm."),
    )
    sy.status("Iniciando cálculo de incertezas...")
    chi_sigma1 = chi_min + 1
    chi_sigma2 = chi_min + 4
    chi_sigma3 = chi_min + 9
    R_INDchi_sigma1 = np.abs(var1d[0][INDchi_min:] - chi_sigma1).argmin() + INDchi_min
    L_INDchi_sigma1 = np.abs(var1d[0][:INDchi_min] - chi_sigma1).argmin()
    R_INDchi_sigma2 = np.abs(var1d[0][INDchi_min:] - chi_sigma2).argmin() + INDchi_min
    L_INDchi_sigma2 = np.abs(var1d[0][:INDchi_min] - chi_sigma2).argmin()
    R_INDchi_sigma3 = np.abs(var1d[0][INDchi_min:] - chi_sigma3).argmin() + INDchi_min
    L_INDchi_sigma3 = np.abs(var1d[0][:INDchi_min] - chi_sigma3).argmin()
    sy.param(
        ("Qui-quadrado 1_sigma (chi2_1_sigma)", chi_sigma1, "adm."),
        ("Qui-quadrado 2_sigma (chi2_2_sigma)", chi_sigma2, "adm."),
        ("Qui-quadrado 3_sigma (chi2_3_sigma)", chi_sigma3, "adm."),
    )
except Exception as e:
    sy.ok(f"Um erro foi encontrado ao tentar calcular as incertezas. Erro: {e}", False)
try:
    sy.status("Iniciando criação de gráfico de distribuição qui-quadrado...")
    deltachi1=var1d[1][INDchi_min]-var1d[1][L_INDchi_sigma1]
    deltachi2=var1d[1][INDchi_min]-var1d[1][L_INDchi_sigma2]
    deltachi3=var1d[1][INDchi_min]-var1d[1][L_INDchi_sigma3]
    gp.basicstyle(
        x_data=var1d[1],
        y_data=var1d[0],
        highlight_point=(var1d[1][INDchi_min], chi_min),
        title="",
        x_label=r"$\Omega_M$",
        y_label=r"$\chi^2$",
        save_fig=True,
        linestyle="--",
        highlight_label="Densidade de Matéria Bestfit",
        highlight_marker="|",
        highlight_size=200,
        curve_label=r"Distribuição de $\chi^2$",
        filename="CHIdistribuicao",
        show_plot=False,
        sigma_intervals=((var1d[1][INDchi_min]-deltachi1,var1d[1][INDchi_min]+deltachi1),(var1d[1][INDchi_min]-deltachi2,var1d[1][INDchi_min]+deltachi2),(var1d[1][INDchi_min]-deltachi3,var1d[1][INDchi_min]+deltachi3)),
        show_sigma_lines=True,
    )
    sy.ok(("Gráfico de distribuição qui-quadrado"))
except Exception as e:
    sy.ok(
        f"Um erro foi encontrado ao tentar plotar o gráfico de distribuição qui-quadrado. Erro: {e}",
        False,
    )

# ? -----------------------------------------------------------------------------
# ?         PROBABILIDADE DE ENERGIA ESCURA
# ? -----------------------------------------------------------------------------

try:
    sy.status("Iniciando cálculo da probabilidade de (Omega_EE > 0.5)...")
    prob_ee = veross1d(var1d[0], chi_min, var1d[2])
    if prob_ee > 0.5:
        ee1d = "confirmada!"
    else:
        ee1d = "não confirmada!"
    sy.param(
        ("Probabilidade de Energia Escura (P_ee)", prob_ee, "adm."),
        ("Energia Escura do universo superior a 50%", ee1d),
    )
except Exception as e:
    sy.ok(
        f"Um erro foi encontrado ao tentar calcular a probabilidade de energia escura. Erro: {e}",
        False,
    )

# ? -----------------------------------------------------------------------------
# ?         VARREDURA BIDIMENSIONAL
# ? -----------------------------------------------------------------------------

try:
    sy.status("Iniciando a varredura 2D...")
    omegaM_list = []
    omegaM_var = mesh_inter_omega[0]
    while omegaM_var < mesh_inter_omega[1]:
        omegaM_list.append(omegaM_var)
        omegaM_var = omegaM_var + meshgrid_step
    omegaEE_list = omegaM_list
    var2d = varredura_2D(
        omegaM_list, omegaEE_list, mu_obs_list, ERROmu_obs_list, z_list, w
    )
    chi2d_min = np.min(var2d)
    MINDchi2d_min, EEINDchi2d_min = np.unravel_index(np.argmin(var2d), var2d.shape)
    sy.status("Varredura 2D concluida...")
except Exception as e:
    sy.ok(f"Um erro foi encontrado ao tentar realizar a varredura 2D. Erro: {e}", False)
try:
    sy.param(
        ("Qui-quadrado mínimo (2dchi2_min)", chi2d_min, "adm."),
        (
            "Densidade de materia Bestfit (2dBF_Omega_M)",
            omegaM_list[MINDchi2d_min],
            "adm.",
        ),
        (
            "Densidade de energia Bestfit (2dBF_Omega_EE)",
            omegaEE_list[EEINDchi2d_min],
            "adm.",
        ),
    )
    sy.status("Iniciando cálculo de incertezas...")
    chi2d_sigma1 = chi2d_min + 2.30
    chi2d_sigma2 = chi2d_min + 6.18
    chi2d_sigma3 = chi2d_min + 11.83
    sy.param(
        ("Qui-quadrado 1_sigma (chi2_1_sigma)", chi2d_sigma1),
        ("Qui-quadrado 2_sigma (chi2_2_sigma)", chi2d_sigma2),
        ("Qui-quadrado 3_sigma (chi2_3_sigma)", chi2d_sigma3),
    )
    niveis = [chi2d_sigma1, chi2d_sigma2, chi2d_sigma3]
except Exception as e:
    sy.ok(f"Um erro foi encontrado ao tentar calcular as incertezas. Erro: {e}", False)
try:
    sy.status("Iniciando criação de gráfico de distribuição qui-quadrado...")
    gp.elipse(
        x_data=omegaM_list,
        y_data=omegaEE_list,
        z_data=np.transpose(var2d),
        highlight_point=(
            round(omegaM_list[MINDchi2d_min], 2),
            round(omegaEE_list[EEINDchi2d_min], 2),
        ),
        ellipse_levels=niveis,
        sigma_names=(r"$\sigma_1$", r"$\sigma_2$", r"$\sigma_3$"),
        x_label=r"$\Omega_m$",
        y_label=r"$\Omega_\Lambda$",
        highlight_size=125,
        ellipse_styles=["-", "-", "-"],
        save_fig=True,
        show_grid=False,
        extra_line_x=omegaM_list,
        extra_line_y=1 - np.array(omegaM_list),
        extra_line_label="Universo Plano",
        legend_frame=True,
        colorbar_format="max",
        legend_alpha=0.15,
        legend_fontsize=10,
        colorbar_ticks=3,
        highlight_label="Melhor Ajuste",
        filename="mapadecalorMEE",
        show_plot=False,
    )
    sy.ok(("Gráfico de Distribuição de Qui-Quadrado"))
except Exception as e:
    sy.ok(
        f"Um erro foi encontrado ao tentar criar o gráfico de distribuição qui-quadrado. Erro: {e}",
        False,
    )

# ? -----------------------------------------------------------------------------
# ?         QUEBRA DE DEGENERECÊNCIA (PRIOR DA CMB)
# ? -----------------------------------------------------------------------------
try:

    sy.status(
        "Iniciando quebra de degenerescência: Aplicando Prior da CMB na malha 2D..."
    )
    var2dPRIOR = quebra_degenerecencia(
        omegaM_list, omegaEE_list, var2d, Omega_K_obs, ERROOmega_K_obs
    )
    chi2d_minPRIOR = np.min(var2dPRIOR)
    MINDchi2d_minPRIOR, EEINDchi2d_minPRIOR = np.unravel_index(
        np.argmin(var2dPRIOR), var2dPRIOR.shape
    )
    sy.param(
        (
            "Qui-quadrado mínimo com quebra de degenerecência (2dchi2_minPRIOR)",
            chi2d_minPRIOR,
        ),
        (
            "Densidade de materia Bestfit com quebra de degenerecência (2dBF_Omega_MPRIOR)",
            omegaM_list[MINDchi2d_minPRIOR],
        ),
        (
            "Densidade de energia Bestfit com quebra de degenerecência (2dBF_Omega_EEPRIOR)",
            omegaEE_list[EEINDchi2d_minPRIOR],
        ),
    )
except Exception as e:
    sy.ok(
        f"Um erro foi encontrado ao tentar realizar a quebra de degenerecência. Erro: {e}",
        False,
    )
try:
    sy.status("Iniciando cálculo de incertezas")
    chi2d_sigma1PRIOR = chi2d_minPRIOR + 2.30
    chi2d_sigma2PRIOR = chi2d_minPRIOR + 6.18
    chi2d_sigma3PRIOR = chi2d_minPRIOR + 11.83
    sy.param(
        ("Qui-quadrado 1_sigma (chi2_1_sigma)", chi2d_sigma1PRIOR),
        ("Qui-quadrado 2_sigma (chi2_2_sigma)", chi2d_sigma2PRIOR),
        ("Qui-quadrado 3_sigma (chi2_3_sigma)", chi2d_sigma3PRIOR),
    )
    niveisPRIOR = [chi2d_sigma1PRIOR, chi2d_sigma2PRIOR, chi2d_sigma3PRIOR]
except Exception as e:
    sy.ok(f"Um erro foi encontrado ao tentar calcular as incertezas. Erro: {e}", False)
try:
    sy.status("Iniciando criação de gráfico de distribuição qui-quadrado corrigido...")

    gp.elipse(
        x_data=omegaM_list,
        y_data=omegaEE_list,
        z_data=np.transpose(var2dPRIOR),
        highlight_point=(
            round(omegaM_list[MINDchi2d_minPRIOR], 2),
            round(omegaEE_list[EEINDchi2d_minPRIOR], 2),
        ),
        ellipse_levels=niveisPRIOR,
        sigma_names=(r"$\sigma_1$", r"$\sigma_2$", r"$\sigma_3$"),
        x_label=r"$\Omega_m$",
        y_label=r"$\Omega_\Lambda$",
        highlight_size=125,
        ellipse_styles=["-", "-", "-"],
        save_fig=True,
        show_grid=False,
        extra_line_x=omegaM_list,
        extra_line_y=1 - np.array(omegaM_list),
        extra_line_label="Universo Plano",
        legend_frame=True,
        colorbar_format="max",
        legend_alpha=0.15,
        legend_fontsize=10,
        colorbar_ticks=3,
        highlight_label="Melhor Ajuste",
        filename="mapadecalormEEPRIOR",
        show_plot=False,
    )

    sy.ok(("Gráfico de Distribuição de Qui-Quadrado Corrigido"))
except Exception as e:
    sy.ok(
        f"Um erro foi encontrado ao tentar plotar o gráfico de distribuição qui-quadrado corrigido. Erro: {e}",
        False,
    )

# ? -----------------------------------------------------------------------------
# ?         PROBABILIDADE DE ACELERAÇÃO 2D
# ? -----------------------------------------------------------------------------

try:
    sy.status("Iniciando cálculo da probabilidade de aceleração")
    prob_aceleracao = veross2d(var2dPRIOR, chi2d_minPRIOR, omegaM_list, omegaEE_list, w)
    if prob_aceleracao > 0.997:
        acel = "confirmada!"
    else:
        acel = "não confirmada!"
    sy.param(
        ("Probabilidade de Aceleração (P_acel)", prob_aceleracao),
        ("Aceleração do universo superior a 3-sigma", acel),
    )
except Exception as e:
    sy.ok(
        f"Um erro foi encontrado ao tentar calcular a probabilidade de aceleração. Erro: {e}",
        False,
    )

# ? -----------------------------------------------------------------------------
# ?         VARREDURA BIDIMENSIONAL - OMEGA_M VS W
# ? -----------------------------------------------------------------------------

try:
    sy.status("Iniciando varredura bidimensional (OMEGA_M vs W)...")
    w_list = []
    w_var = mesh_inter_w[0]
    while w_var < mesh_inter_w[1]:
        w_list.append(w_var)
        w_var = w_var + meshgrid_step
    var2dOW = varreduraOW(w_list, omegaM_list, mu_obs_list, ERROmu_obs_list, z_list)
    chi2dOW_min = np.min(var2dOW)
    WINDchi2dOW_min, MINDchi2dOW_min = np.unravel_index(
        np.argmin(var2dOW), var2dOW.shape
    )
except Exception as e:
    sy.ok(
        f"Um erro foi encontrado ao tentar realizar a varredura bidimensional (OMEGA_M vs W). Erro: {e}",
        False,
    )
try:
    sy.status("Iniciando cálculo de parâmetros da varredura 2D (w vs Omega_M)...")
    sy.param(
        ("Qui-quadrado mínimo (2dchi2_min)", chi2dOW_min),
        ("Parâmetro da equação de estado Bestfit (2dBF_w)", w_list[WINDchi2dOW_min]),
        ("Densidade de matéria Bestfit (2dBF_Omega_M)", omegaM_list[MINDchi2dOW_min]),
    )
    sy.status("Iniciando cálculo de incertezas")
    chi2dOW_sigma1 = chi2dOW_min + 2.30
    chi2dOW_sigma2 = chi2dOW_min + 6.18
    chi2dOW_sigma3 = chi2dOW_min + 11.83
    sy.param(
        ("Qui-quadrado 1_sigma (chi2_1_sigma)", chi2dOW_sigma1),
        ("Qui-quadrado 2_sigma (chi2_2_sigma)", chi2dOW_sigma2),
        ("Qui-quadrado 3_sigma (chi2_3_sigma)", chi2dOW_sigma3),
    )
    niveis = [chi2dOW_sigma1, chi2dOW_sigma2, chi2dOW_sigma3]
except Exception as e:
    sy.ok(
        f"Um erro foi encontrado ao tentar calcular as incertezas. Erro: {e}",
        False,
    )
try:
    sy.status("Iniciando criação de gráfico de distribuição qui-quadrado")

    gp.elipse(
        x_data=w_list,
        y_data=omegaM_list,
        z_data=np.transpose(var2dOW),
        highlight_point=(
            round(w_list[WINDchi2dOW_min], 2),
            round(omegaM_list[MINDchi2dOW_min], 2),
        ),
        ellipse_levels=niveis,
        sigma_names=(r"$\sigma_1$", r"$\sigma_2$", r"$\sigma_3$"),
        y_label=r"$\Omega_m$",
        x_label=r"$w$",
        highlight_size=125,
        ellipse_styles=["-", "-", "-"],
        save_fig=True,
        show_grid=False,
        legend_frame=True,
        colorbar_format="max",
        legend_alpha=0.15,
        legend_fontsize=10,
        colorbar_ticks=3,
        highlight_label="Melhor Ajuste",
        filename="mapadecalorWM",
        show_plot=False,
    )
except Exception as e:
    sy.ok(
        f"Um erro foi encontrado ao tentar criar o gráfico de distribuição qui-quadrado. Erro: {e}",
        False,
    )

sy.fim()
