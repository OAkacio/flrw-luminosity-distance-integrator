#
# * =============================================================================
# * DEPENDÊNCIAS
# * =============================================================================

# ? -----------------------------------------------------------------------------
# ?         MÓDULOS LOCAIS
# ? -----------------------------------------------------------------------------

from src.constants import *

# * =============================================================================
# * PARÂMETROS INICIAIS
# * =============================================================================

Omega_M = 0.315  # Densidade de matéria [adm.]
Omega_EE = 0.685  # Densidade de energia escura [adm.]
w = -1  # Parâmetro de equação de estado da energia escura [adm.]
z_max = 10  # Redshift máximo [adm.]
MagABS = -19.05  # Magnitude absoluta de vela padrão [mag]

# * =============================================================================
# * PARÂMETROS DE ANÁLISE
# * =============================================================================

# ? -----------------------------------------------------------------------------
# ?         MAIN
# ? -----------------------------------------------------------------------------

z_step = 1e-5  # Intervalo entre um ponto e outro nos dados

# ? -----------------------------------------------------------------------------
# ?         INFERENCE
# ? -----------------------------------------------------------------------------

mesh_inter_omega = [0, 1]  # Intervalo de análise [omega_ini, omega_fin]
mesh_inter_w = [-1.25, 0.25]  # Intervalo de análise [w_ini, w_fin]
meshgrid_step = 5e-3  # Intervalo entre um ponto e outro nos dados
Omega_K_obs = -0.06  # Parâmetro de Densidade da Curvatura observado
ERROOmega_K_obs = 0.05  # Erro do Parâmetro de Densidade da Curvatura observado
