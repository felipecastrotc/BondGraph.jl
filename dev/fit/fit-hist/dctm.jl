using Statistics

function dctm(f, range; ϵ=1e-14, l=NaN, k_max=2000)
  #dctm Algoritmo de Busca Dicotômica
  #   Essse algoritmo de Busca Dicotômica foi desenvolvido como parte da disciplina
  #   'Métodos de Otimização de Sistemas' da Faculdade de Engenharia Mecânica
  #   da Universidade Estadual de Campinas.

  #   Autor: Felipe de Castro Teixeira Carvalho
  #   RA: 190823

  # Valores padrão das variáveis
  isnan(l) && (l = 2.1*ϵ)
  # Determina-se um valor inicial do intervalo atual (l_0)
  l_0 = Inf
  k = 0
  while k <= k_max && l_0 > l
      # Calcula-se o ponto médio
      x_M = (range[1] + range[2])/2
      # Calcula-se os valores de λ e ϵ
      λ = x_M - ϵ
      μ = x_M + ϵ
      # Cria-se um novo intervalo baseado em f(λ) e f(μ)
      if f(λ)[1] > f(μ)[1]
          range = [λ range[2]]
      else
          range = [range[1] μ]
      end
      # Adiciona-se 1 na variável contadora
      k += 1
      l_0 = abs(range[2] - range[1])
  end
  erro = l_0
  # Calcula-se a média do intervalo
  out = mean(range)
  return out, k-1, erro
end
